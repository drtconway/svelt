use std::{
    collections::HashMap,
    io::{BufReader, BufWriter},
};

use autocompress::{CompressionLevel, autodetect_create, autodetect_open};
use datafusion::{
    arrow::{
        array::{Array, Float64Array, GenericStringArray, PrimitiveArray, RecordBatch},
        datatypes::Float64Type,
    },
    common::JoinType,
    functions_aggregate::{count::count, expr_fn::first_value, min_max::max, sum::sum},
    prelude::{DataFrame, SessionContext, col, lit, sqrt},
    scalar::ScalarValue,
};
use itertools::Itertools;
use noodles::fasta::{
    self,
    record::{Definition, Sequence},
};

use crate::{
    disjoint_set::DisjointSet,
    expressions::prefix_cols,
    features::FeatureIndex,
    kmers::Kmer,
    kmers_table::kmer_frequencies_to_table,
    options::{CommonOptions, QueryOptions, make_session_context},
    sequence::{SequenceIterator, fasta::FastaSequenceIterator, vcf::VcfSequenceIterator},
};

pub async fn find_similar(
    features: &str,
    query: &QueryOptions,
    k: usize,
    common: &CommonOptions,
) -> std::io::Result<()> {
    let ctx = make_session_context(common);

    if let Some(query) = &query.query {
        return find_similar_single(features, query, k, common).await;
    }

    let idx = FeatureIndex::load(features, &ctx).await?;

    if let Some(query_file) = &query.query_file {
        let itr = FastaSequenceIterator::new(&query_file)?;
        find_similar_inner(itr, k, &idx, &ctx).await?;
    }

    if let Some(vcf_file) = &query.vcf {
        let itr = VcfSequenceIterator::new(&vcf_file)?;
        find_similar_inner(itr, k, &idx, &ctx).await?;
    }

    Ok(())
}

async fn find_similar_inner<Itr: SequenceIterator>(
    itr: Itr,
    k: usize,
    idx: &FeatureIndex,
    ctx: &SessionContext,
) -> std::io::Result<()> {
    let mut curr: Option<DataFrame> = None;
    let mut curr_count: usize = 0;
    for rec in itr {
        let (name, sequence) = rec?;
        log::info!("scanning {} ({}bp)", name, sequence.len());

        let mut fwd = Vec::new();
        let mut rev = Vec::new();
        Kmer::with_many_both(k, &sequence, |x, y| {
            fwd.push(x.0);
            rev.push(y.0);
        });

        fwd.sort();
        let fwd: Vec<(u64, usize)> = fwd
            .into_iter()
            .chunk_by(|x| *x)
            .into_iter()
            .map(|(x, xs)| (x, xs.count()))
            .collect();
        curr_count += fwd.len();
        let fwd = kmer_frequencies_to_table(&fwd, &ctx)
            .await?
            .with_column("strand", lit("+"))?;

        rev.sort();
        let rev: Vec<(u64, usize)> = rev
            .into_iter()
            .chunk_by(|x| *x)
            .into_iter()
            .map(|(x, xs)| (x, xs.count()))
            .collect();
        curr_count += rev.len();
        let rev = kmer_frequencies_to_table(&rev, &ctx)
            .await?
            .with_column("strand", lit("-"))?;

        let block = fwd.union(rev)?.with_column("name", lit(name))?;
        curr = if let Some(block0) = curr {
            Some(block0.union(block)?)
        } else {
            Some(block)
        };

        if curr_count > 500000 {
            log::info!("ranking over {} k-mers.", curr_count);
            let tbl = curr.take().unwrap();
            let res = idx.rank_many(tbl).await?;
            res.aggregate(
                vec![col("name")],
                vec![
                    max(col("dot")).alias("dot"),
                    first_value(col("idx_name"), Some(vec![col("dot").sort(false, false)]))
                        .alias("feature"),
                    first_value(col("class"), Some(vec![col("dot").sort(false, false)]))
                        .alias("class"),
                ],
            )?
            .show()
            .await?;
            curr_count = 0;
        }
    }

    if curr.is_none() {
        return Ok(());
    }

    log::info!("ranking over {} k-mers.", curr_count);
    let all = curr.unwrap();
    let res = idx.rank_many(all).await?;

    res.aggregate(
        vec![col("name")],
        vec![
            first_value(col("idx_name"), Some(vec![col("dot").sort(false, false)]))
                .alias("feature"),
            first_value(col("class"), Some(vec![col("dot").sort(false, false)])).alias("class"),
        ],
    )?
    .show()
    .await?;

    Ok(())
}

async fn find_similar_single(
    features: &str,
    query: &str,
    k: usize,
    common: &CommonOptions,
) -> std::io::Result<()> {
    let ctx = make_session_context(common);

    let idx = FeatureIndex::load(features, &ctx).await?;

    let mut fwd = Vec::new();
    let mut rev = Vec::new();
    Kmer::with_many_both(k, &query, |x, y| {
        fwd.push(x.0);
        rev.push(y.0);
    });

    fwd.sort();
    let fwd: Vec<(u64, usize)> = fwd
        .into_iter()
        .chunk_by(|x| *x)
        .into_iter()
        .map(|(x, xs)| (x, xs.count()))
        .collect();
    let fwd = kmer_frequencies_to_table(&fwd, &ctx).await?;
    let fwd = idx.rank(fwd).await?.with_column("strand", lit("+"))?;

    rev.sort();
    let rev: Vec<(u64, usize)> = rev
        .into_iter()
        .chunk_by(|x| *x)
        .into_iter()
        .map(|(x, xs)| (x, xs.count()))
        .collect();
    let rev = kmer_frequencies_to_table(&rev, &ctx).await?;
    let rev = idx.rank(rev).await?.with_column("strand", lit("-"))?;

    fwd.union(rev)?
        .sort(vec![col("distance").sort(true, false)])?
        .show()
        .await?;

    Ok(())
}

pub async fn cluster_sequences(
    sequences: &str,
    out: &str,
    k: usize,
    cutoff: f64,
    common: &CommonOptions,
) -> std::io::Result<()> {
    let ctx = make_session_context(common);

    let reader = autodetect_open(sequences)?;
    let reader = BufReader::new(reader);
    let mut reader = fasta::io::reader::Builder::default().build_from_reader(reader)?;

    let mut sequence_index = HashMap::new();
    let mut curr: Option<DataFrame> = None;
    let mut curr_count: usize = 0;
    for rec in reader.records() {
        let rec = rec?;

        let name = String::from_utf8(rec.name().to_vec()).unwrap();
        log::info!("scanning {}", name);

        let seq = rec.sequence();

        sequence_index.insert(
            name.clone(),
            String::from_utf8(seq.as_ref().to_vec()).unwrap(),
        );

        let mut fwd = Vec::new();
        Kmer::with_many_both(k, seq, |x, _y| {
            fwd.push(x.0);
        });

        fwd.sort();
        let fwd: Vec<(u64, usize)> = fwd
            .into_iter()
            .chunk_by(|x| *x)
            .into_iter()
            .map(|(x, xs)| (x, xs.count()))
            .collect();
        curr_count += fwd.len();
        let fwd = kmer_frequencies_to_table(&fwd, &ctx).await?;

        let block = fwd.with_column("name", lit(name))?;
        curr = if let Some(block0) = curr {
            Some(block0.union(block)?)
        } else {
            Some(block)
        };
    }

    if curr.is_none() {
        log::warn!("no sequences to cluster!");
        return Ok(());
    }

    log::info!("ranking over {} k-mers.", curr_count);
    let all = curr.unwrap();

    let _mag = all
        .clone()
        .clone()
        .aggregate(
            vec![col("name")],
            vec![
                sum(col("count") * col("count")).alias("mag"),
                sum(col("count")).alias("mag_count"),
            ],
        )?
        .with_column("mag", sqrt(col("mag")))?
        .with_column_renamed("name", "mag_name")?;

    let lhs = prefix_cols(all.clone(), "lhs")?;
    let rhs = prefix_cols(all.clone(), "rhs")?;

    let tbl = lhs.join(
        rhs,
        JoinType::Inner,
        &["lhs_kmer"],
        &["rhs_kmer"],
        Some(lit(true).and(col("lhs_name").lt(col("rhs_name")))),
    )?;

    let lhs_missing = all
        .clone()
        .join(
            tbl.clone(),
            JoinType::LeftAnti,
            &["name", "kmer"],
            &["lhs_name", "lhs_kmer"],
            None,
        )?
        .aggregate(
            vec![col("name")],
            vec![sum(col("count")).alias("lhs_missing")],
        )?;

    let rhs_missing = all
        .clone()
        .join(
            tbl.clone(),
            JoinType::LeftAnti,
            &["name", "kmer"],
            &["rhs_name", "rhs_kmer"],
            None,
        )?
        .aggregate(
            vec![col("name")],
            vec![sum(col("count")).alias("rhs_missing")],
        )?;

    let tbl = tbl
        .with_column("count_difference", col("lhs_count") - col("rhs_count"))?
        .with_column("count_sum", col("lhs_count") + col("rhs_count"))?;

    if false {
        tbl.clone().show().await?;
    }

    let tbl = tbl
        .aggregate(
            vec![col("lhs_name"), col("rhs_name")],
            vec![
            sum(col("lhs_count") * col("rhs_count")).alias("raw_dot"),
            sum(lit(1.0) * col("count_difference") * col("count_difference") / col("count_sum"))
                .alias("chi_sq"),
            count(col("count_sum")).alias("df")
        ],
        )?
        .join(lhs_missing, JoinType::Left, &["lhs_name"], &["name"], None)?
        .drop_columns(&["name"])?
        .join(rhs_missing, JoinType::Left, &["rhs_name"], &["name"], None)?
        .drop_columns(&["name"])?
        .fill_null(
            ScalarValue::from(0i32),
            vec![String::from("lhs_missing"), String::from("rhs_missing")],
        )?
        .with_column(
            "chi_sq",
            lit(0.5) * (col("chi_sq") + col("lhs_missing") + col("rhs_missing")),
        )?
        .drop_columns(&["lhs_missing", "rhs_missing"])?;

    let tbl = tbl.clone().sort(vec![col("score").sort(true, false)])?;

    if true {
        tbl.clone().show().await?;
    }

    /*
        if use_cosine {
            let tbl = tbl
                .join(
                    mag.clone(),
                    JoinType::Left,
                    &["lhs_name"],
                    &["mag_name"],
                    None,
                )?
                .with_column_renamed("mag", "lhs_mag")?
                .with_column_renamed("mag_count", "lhs_count")?
                .drop_columns(&["mag_name"])?
                .join(
                    mag.clone(),
                    JoinType::Left,
                    &["rhs_name"],
                    &["mag_name"],
                    None,
                )?
                .with_column_renamed("mag", "rhs_mag")?
                .with_column_renamed("mag_count", "rhs_count")?
                .drop_columns(&["mag_name"])?
                .with_column(
                    "dot",
                    col("raw_dot") * lit(1.0) / (col("lhs_mag") * col("rhs_mag")),
                )?
                .with_column(
                    "ratio",
                    lit(1.0) * least(vec![col("lhs_count"), col("rhs_count")])
                        / greatest(vec![col("lhs_count"), col("rhs_count")]),
                )?
                .with_column("score", col("dot") * col("ratio"))?;

            let tbl = tbl.clone().sort(vec![
                col("score").sort(false, false),
                col("dot").sort(false, false),
                col("ratio").sort(false, false),
            ])?;

            if true {
                tbl.clone()
                    .sort(vec![
                        col("score").sort(false, false),
                        col("lhs_name").sort(true, false),
                        col("rhs_name").sort(true, false),
                    ])?
                    .show()
                    .await?;
            }
        }
    */

    let batch = tbl.clone().collect().await?;
    let itr = batch.iter().flat_map(|recs| MergeIterator::new(recs));

    let mut uf: DisjointSet<usize> = DisjointSet::new();
    let mut names = Vec::new();
    let mut name_index: HashMap<String, usize> = HashMap::new();
    let mut edge_count: HashMap<usize, usize> = HashMap::new();
    for item in itr {
        if !name_index.contains_key(item.0) {
            let n = name_index.len();
            names.push(String::from(item.0));
            name_index.insert(String::from(item.0), n);
        }
        if !name_index.contains_key(item.1) {
            let n = name_index.len();
            names.push(String::from(item.1));
            name_index.insert(String::from(item.1), n);
        }

        let d = item.2;
        if d < cutoff {
            let l = *name_index.get(item.0).unwrap();
            let r = *name_index.get(item.1).unwrap();

            *edge_count.entry(l).or_default() += 1;
            *edge_count.entry(r).or_default() += 1;

            let a = uf.find(l);
            let b = uf.find(r);
            if a != b {
                uf.union(a, b);
            }
        }
    }

    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut group_index: Vec<usize> = (0..name_index.len()).map(|_| 0).collect();
    for i in 0..name_index.len() {
        let j = uf.find(i);
        group_index[i] = j;
        groups.entry(j).or_default().push(i);
    }

    log::info!("at cutoff {} there were {} groups", cutoff, groups.len());

    let out = autodetect_create(out, CompressionLevel::default())?;
    let out = BufWriter::new(out);
    let mut out = fasta::io::Writer::new(out);

    for group_item in groups.iter() {
        let mut max_member = *group_item.0;
        let mut max_link = 0;
        for j in group_item.1.iter() {
            let c = if let Some(c) = edge_count.get(j) {
                *c
            } else {
                0
            };
            if c > max_link {
                max_member = *j;
                max_link = c;
            }
        }
        let name = &names[max_member];
        let seq = sequence_index.get(name).unwrap();
        let seq = Sequence::from(seq.as_bytes().to_vec());

        let rec = fasta::Record::new(Definition::new(name as &str, None), seq);
        out.write_record(&rec)?;
    }

    Ok(())
}

struct MergeIterator<'a> {
    lhs_name: &'a GenericStringArray<i32>,
    rhs_name: &'a GenericStringArray<i32>,
    score: &'a PrimitiveArray<Float64Type>,
    i: usize,
}

impl<'a> MergeIterator<'a> {
    pub fn new(recs: &'a RecordBatch) -> MergeIterator<'a> {
        let lhs_name = Self::get_array::<GenericStringArray<i32>>(recs, "lhs_name");
        let rhs_name = Self::get_array::<GenericStringArray<i32>>(recs, "rhs_name");
        let score = Self::get_array::<Float64Array>(recs, "score");
        MergeIterator {
            lhs_name,
            rhs_name,
            score,
            i: 0,
        }
    }

    fn get_array<Type: 'static>(recs: &'a RecordBatch, name: &str) -> &'a Type {
        recs.column_by_name(name)
            .unwrap()
            .as_any()
            .downcast_ref::<Type>()
            .unwrap()
    }
}

impl<'a> Iterator for MergeIterator<'a> {
    type Item = (&'a str, &'a str, f64);

    fn next(&mut self) -> Option<Self::Item> {
        if self.i < self.lhs_name.len() {
            let i = self.i;
            self.i += 1;

            let lhs_name = self.lhs_name.value(i);
            let rhs_name = self.rhs_name.value(i);
            let dot = self.score.value(i);
            Some((lhs_name, rhs_name, dot))
        } else {
            None
        }
    }
}
