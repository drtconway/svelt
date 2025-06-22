use std::{
    collections::HashMap,
    io::{BufReader, BufWriter, Error, ErrorKind},
    sync::Arc,
};

use autocompress::{CompressionLevel, autodetect_create, autodetect_open};
use datafusion::{
    arrow::{
        array::{
            Array, Float64Array, GenericStringArray, PrimitiveArray, PrimitiveBuilder, RecordBatch,
            StringDictionaryBuilder,
        },
        datatypes::{DataType, Field, Float64Type, Schema, UInt32Type, UInt64Type},
    },
    common::JoinType,
    config::{ParquetColumnOptions, TableParquetOptions},
    dataframe::DataFrameWriteOptions,
    functions_aggregate::{count::count, expr_fn::first_value, min_max::max, sum::sum},
    prelude::{DataFrame, ParquetReadOptions, SessionContext, cast, col, lit, sqrt},
    scalar::ScalarValue,
};
use itertools::Itertools;
use noodles::fasta::{
    self,
    record::{Definition, Sequence},
};
use regex::Regex;

use crate::{
    disjoint_set::DisjointSet,
    distance::{DistanceMetric, distance},
    either::Either::{self, Left, Right},
    errors::{SveltError, as_io_error},
    expressions::prefix_cols,
    kmers::Kmer,
    kmers_table::kmer_frequencies_to_table,
    options::{CommonOptions, IndexingOptions, QueryOptions, make_session_context},
    sequence::{FastaSequenceIterator, SequenceIterator, VcfSequenceIterator},
};

pub async fn index_features(
    features_name: &str,
    out: &str,
    options: &IndexingOptions,
    common: &CommonOptions,
) -> std::io::Result<()> {
    let reader = autodetect_open(features_name)?;
    let reader = BufReader::new(reader);
    let mut reader = fasta::io::reader::Builder::default().build_from_reader(reader)?;

    let mut names = Vec::new();
    let mut kmer_index: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();

    for rec in reader.records() {
        let rec = rec?;
        let name = rec.definition().to_string().split_off(1);

        let n = names.len();
        names.push(name.clone());

        let mut fwd = Vec::new();
        Kmer::with_many_both(options.k, rec.sequence(), |x, _y| {
            fwd.push(x.0);
        });
        fwd.sort();
        for (x, c) in fwd
            .into_iter()
            .chunk_by(|x| *x)
            .into_iter()
            .map(|(x, xs)| (x, xs.count()))
        {
            kmer_index.entry(x).or_default().push((n, c));
        }
    }

    let ctx = make_session_context(common);

    let rx = Regex::new(&options.pattern).map_err(|e| Error::new(ErrorKind::Other, e))?;
    let nm: Either<usize, String> = if let Ok(x) = options.name.parse() {
        Left(x)
    } else {
        Right(options.name.clone())
    };
    let cls: Either<usize, String> = if let Ok(x) = options.class.parse() {
        Left(x)
    } else {
        Right(options.class.clone())
    };

    // Write out the names

    let mut name_builder = StringDictionaryBuilder::<UInt32Type>::new();
    let mut class_builder = StringDictionaryBuilder::<UInt32Type>::new();
    let mut nix_builder = PrimitiveBuilder::<UInt32Type>::new();

    for (nix, orig) in names.into_iter().enumerate() {
        if let Some(m) = rx.captures(&orig) {
            let name = match &nm {
                Left(n) => String::from(&m[*n]),
                Right(v) => String::from(&m[v as &str]),
            };
            let class = match &cls {
                Left(n) => String::from(&m[*n]),
                Right(v) => String::from(&m[v as &str]),
            };
            name_builder.append_value(name);
            class_builder.append_value(class);
            nix_builder.append_value(nix as u32);
        } else {
            log::warn!("Couldn't parse feature sequence name: '{}'", orig);
            name_builder.append_value(orig);
            class_builder.append_value(String::from("<unknown>"));
            nix_builder.append_value(nix as u32);
        }
    }

    let name_array = name_builder.finish();
    let class_array = class_builder.finish();
    let nix_array = nix_builder.finish();

    let name_schema = Arc::new(Schema::new(vec![
        Field::new_dictionary("name", DataType::UInt32, DataType::Utf8, false),
        Field::new_dictionary("class", DataType::UInt32, DataType::Utf8, false),
        Field::new("nix", DataType::UInt32, false),
    ]));

    let recs = RecordBatch::try_new(
        name_schema,
        vec![
            Arc::new(name_array),
            Arc::new(class_array),
            Arc::new(nix_array),
        ],
    )
    .map_err(|e| Error::new(ErrorKind::Other, e))?;

    let df = ctx.read_batch(recs)?;

    let opts = DataFrameWriteOptions::default();

    df.sort_by(vec![col("class"), col("name")])?
        .write_parquet(&format!("{}-nidx.parquet", out), opts, None)
        .await?;

    // Write out the kmer index

    let mut kmer_builder = PrimitiveBuilder::<UInt64Type>::new();
    let mut nix_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut count_builder = PrimitiveBuilder::<UInt32Type>::new();

    for (x, postings) in kmer_index.into_iter() {
        for (n, c) in postings.into_iter() {
            kmer_builder.append_value(x);
            nix_builder.append_value(n as u32);
            count_builder.append_value(c as u32);
        }
    }

    let kmer_array = kmer_builder.finish();
    let nix_array = nix_builder.finish();
    let count_array = count_builder.finish();

    let kmer_schema = Arc::new(Schema::new(vec![
        Field::new("kmer", DataType::UInt64, false),
        Field::new("nix", DataType::UInt32, false),
        Field::new("count", DataType::UInt32, false),
    ]));

    let recs = RecordBatch::try_new(
        kmer_schema,
        vec![
            Arc::new(kmer_array),
            Arc::new(nix_array),
            Arc::new(count_array),
        ],
    )
    .map_err(|e| Error::new(ErrorKind::Other, e))?;

    let df = ctx.read_batch(recs)?;

    let opts = DataFrameWriteOptions::default();

    let mut kidx_opts = TableParquetOptions::default();
    let mut kmer_opts = ParquetColumnOptions::default();
    kmer_opts.compression = Some(String::from("zstd(19)"));
    kidx_opts
        .column_specific_options
        .insert(String::from("kmer"), kmer_opts);
    kidx_opts
        .key_value_metadata
        .insert(String::from("k"), Some(options.k.to_string()));

    df.sort_by(vec![col("kmer"), col("nix")])?
        .write_parquet(&format!("{}-kidx.parquet", out), opts, Some(kidx_opts))
        .await?;

    Ok(())
}

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

    let idx = FeatureIndex::new(features, &ctx).await?;

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

    let idx = FeatureIndex::new(features, &ctx).await?;

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

    let chi_squared_cdf = datafusion_statrs::distribution::chi_squared::cdf();

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
        .drop_columns(&["lhs_missing", "rhs_missing"])?
        .with_column(
            "score",
            chi_squared_cdf.call(vec![col("chi_sq"), cast(col("df"), DataType::Float64)]),
        )?;

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

pub struct FeatureIndex {
    k: usize,
    kmers: DataFrame,
    names: DataFrame,
}

impl FeatureIndex {
    pub async fn new(features: &str, ctx: &SessionContext) -> std::io::Result<FeatureIndex> {
        log::info!("loading index '{}'", features);
        let opts = ParquetReadOptions::default().skip_metadata(false);
        let kmers = ctx
            .read_parquet(&format!("{}-kidx.parquet", features), opts.clone())
            .await?;
        let meta = kmers.schema().metadata();
        let k: usize = if let Some(k_str) = meta.get("k") {
            if let Ok(k) = k_str.parse() {
                k
            } else {
                return Err(as_io_error(SveltError::MissingK(String::from(features))));
            }
        } else {
            return Err(as_io_error(SveltError::MissingK(String::from(features))));
        };

        let names = ctx
            .read_parquet(&format!("{}-nidx.parquet", features), opts)
            .await?;
        let names = names
            .with_column_renamed("name", "idx_name")?
            .with_column_renamed("nix", "idx_nix")?;

        log::info!("loading index done.");

        Ok(FeatureIndex { k, kmers, names })
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub async fn rank(&self, query: DataFrame) -> std::io::Result<DataFrame> {
        let query = query.with_column("name", lit("query"))?;

        let subject = self.kmers.clone().with_column_renamed("nix", "name")?;

        let res = distance(query, subject, DistanceMetric::Cosine).await?;
        let res = res
            .join(
                self.names.clone(),
                JoinType::Left,
                &["subject_name"],
                &["idx_nix"],
                None,
            )?
            .drop_columns(&["subject_name", "idx_nix"])?
            .with_column_renamed("idx_name", "subject_name")?;

        Ok(res)
    }

    pub async fn rank_many(&self, queries: DataFrame) -> std::io::Result<DataFrame> {
        let subject = self.kmers.clone().with_column_renamed("nix", "name")?;

        let res = distance(queries, subject, DistanceMetric::Cosine).await?;
        let res = res
            .join(
                self.names.clone(),
                JoinType::Left,
                &["subject_name"],
                &["idx_nix"],
                None,
            )?
            .drop_columns(&["subject_name", "idx_nix"])?
            .with_column_renamed("idx_name", "subject_name")?;
        Ok(res)
    }
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
