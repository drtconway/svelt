use std::{
    collections::HashMap,
    io::{BufReader, Error, ErrorKind},
    sync::Arc,
};

use autocompress::autodetect_open;
use datafusion::{
    arrow::{
        array::{Float64Array, PrimitiveBuilder, RecordBatch, StringDictionaryBuilder},
        datatypes::{DataType, Field, Schema, UInt32Type, UInt64Type},
    },
    common::JoinType,
    config::{ParquetColumnOptions, TableParquetOptions},
    dataframe::DataFrameWriteOptions,
    functions_aggregate::sum::sum,
    prelude::{DataFrame, ParquetReadOptions, SessionContext, col, lit, sqrt},
};
use itertools::Itertools;
use noodles::fasta;
use regex::Regex;

use crate::{
    either::Either::{self, Left, Right},
    kmers::Kmer,
    options::{CommonOptions, IndexingOptions, QueryOptions, make_session_context},
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
    let mut name_index = HashMap::new();

    let mut kmer_index: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();

    for rec in reader.records() {
        let rec = rec?;
        let name = rec.definition().to_string().split_off(1);

        let n = names.len();
        names.push(name.clone());
        name_index.insert(name, n);

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

    for (orig, nix) in name_index.into_iter() {
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

    let idx = FeatureIndex::new(features, &ctx).await?;

    if let Some(query) = &query.query {
        return find_similar_single(features, query, k, common).await;
    }

    let query_file = query.query_file.as_ref().unwrap();

    let reader = autodetect_open(query_file)?;
    let reader = BufReader::new(reader);
    let mut reader = fasta::io::reader::Builder::default().build_from_reader(reader)?;

    let mut curr: Option<DataFrame> = None;
    let mut curr_count: usize = 0;
    for rec in reader.records() {
        let rec = rec?;

        let name = String::from_utf8(rec.name().to_vec()).unwrap();
        log::info!("scanning {}", name);

        let seq = rec.sequence();
        if seq.len() > 10000 {
            log::info!("skipping long sequence {} ({} bp)", name, seq.len());
            continue;
        }
        let mut fwd = Vec::new();
        let mut rev = Vec::new();
        Kmer::with_many_both(k, seq, |x, y| {
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

        if curr_count > 10000 {
            log::info!("ranking over {} k-mers.", curr_count);
            let tbl = curr.take().unwrap();
            idx.rank_many(tbl).await?;
            curr_count = 0;
        }
    }

    if curr.is_none() {
        return Ok(());
    }

    log::info!("ranking over {} k-mers.", curr_count);
    let all = curr.unwrap();
    idx.rank_many(all).await?;

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
        .sort(vec![col("dot").sort(false, false)])?
        .show()
        .await?;

    Ok(())
}

struct FeatureIndex {
    kmers: DataFrame,
    freqs: DataFrame,
    names: DataFrame,
}

impl FeatureIndex {
    pub async fn new(features: &str, ctx: &SessionContext) -> std::io::Result<FeatureIndex> {
        log::info!("loading index '{}'", features);
        let opts = ParquetReadOptions::default();
        let kmers = ctx
            .read_parquet(&format!("{}-kidx.parquet", features), opts.clone())
            .await?;
        let kmers = kmers
            .with_column_renamed("kmer", "idx_kmer")?
            .with_column_renamed("count", "idx_count")?;
        let names = ctx
            .read_parquet(&format!("{}-nidx.parquet", features), opts)
            .await?;
        let names = names
            .with_column_renamed("name", "idx_name")?
            .with_column_renamed("nix", "idx_nix")?;

        let freqs = kmers
            .clone()
            .aggregate(
                vec![col("nix")],
                vec![sum(col("idx_count") * col("idx_count")).alias("idx_mag")],
            )?
            .with_column("idx_mag", sqrt(col("idx_mag")))?
            .with_column_renamed("nix", "idx_nix")?;

        log::info!("loading index done.");

        Ok(FeatureIndex {
            kmers,
            freqs,
            names,
        })
    }

    pub async fn rank(&self, other: DataFrame) -> std::io::Result<DataFrame> {
        let mag = other
            .clone()
            .aggregate(vec![], vec![sum(col("count") * col("count")).alias("mag")])?
            .with_column("mag", sqrt(col("mag")))?;
        let mag = mag.collect().await?[0]
            .column(0)
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap()
            .value(0);

        let tbl = self.kmers.clone().join(
            other.clone(),
            JoinType::Right,
            &["idx_kmer"],
            &["kmer"],
            None,
        )?;

        if false {
            tbl.clone().show().await?;
        }

        let tbl = tbl.aggregate(
            vec![col("nix")],
            vec![sum(col("idx_count") * col("count")).alias("raw_dot")],
        )?;

        if false {
            tbl.clone()
                .sort(vec![col("raw_dot").sort(false, false)])?
                .show()
                .await?;
        }

        let tbl = tbl
            .join(
                self.freqs.clone(),
                JoinType::Left,
                &["nix"],
                &["idx_nix"],
                None,
            )?
            .drop_columns(&["idx_nix"])?
            .with_column("mag", lit(mag))?
            .with_column(
                "dot",
                col("raw_dot") * lit(1.0) / (col("idx_mag") * col("mag")),
            )?
            .join(
                self.names.clone(),
                JoinType::Left,
                &["nix"],
                &["idx_nix"],
                None,
            )?
            .drop_columns(&["idx_nix"])?
            .filter(col("nix").is_not_null())?
            .with_column_renamed("idx_name", "name")?;

        if false {
            tbl.clone()
                .sort(vec![col("dot").sort(false, false)])?
                .show()
                .await?;
        }

        Ok(tbl)
    }

    pub async fn rank_many(&self, queries: DataFrame) -> std::io::Result<DataFrame> {
        if false {
            queries.clone().show().await?;
        }

        let mag = queries
            .clone()
            .aggregate(
                vec![col("name")],
                vec![sum(col("count") * col("count")).alias("mag")],
            )?
            .with_column("mag", sqrt(col("mag")))?
            .with_column_renamed("name", "mag_name")?;

        let tbl = self.kmers.clone().join(
            queries.clone(),
            JoinType::Right,
            &["idx_kmer"],
            &["kmer"],
            None,
        )?;

        if false {
            tbl.clone().show().await?;
        }

        let tbl = tbl.aggregate(
            vec![col("name"), col("nix"), col("strand")],
            vec![sum(col("idx_count") * col("count")).alias("raw_dot")],
        )?;

        if false {
            tbl.clone()
                .sort(vec![col("raw_dot").sort(false, false)])?
                .show()
                .await?;
        }

        let tbl = tbl
            .join(
                self.freqs.clone(),
                JoinType::Left,
                &["nix"],
                &["idx_nix"],
                None,
            )?
            .drop_columns(&["idx_nix"])?
            .join(mag, JoinType::Left, &["name"], &["mag_name"], None)?
            .drop_columns(&["mag_name"])?
            .with_column(
                "dot",
                col("raw_dot") * lit(1.0) / (col("idx_mag") * col("mag")),
            )?
            .join(
                self.names.clone(),
                JoinType::Left,
                &["nix"],
                &["idx_nix"],
                None,
            )?
            .drop_columns(&["idx_nix"])?
            .filter(col("nix").is_not_null())?
            .filter(col("dot").gt_eq(lit(0.25)))?;

        if true {
            tbl.clone()
                .sort(vec![col("dot").sort(false, false)])?
                .show()
                .await?;
        }

        Ok(tbl)
    }
}

async fn kmer_frequencies_to_table(
    items: &Vec<(u64, usize)>,
    ctx: &SessionContext,
) -> std::io::Result<DataFrame> {
    let mut kmer_builder = PrimitiveBuilder::<UInt64Type>::new();
    let mut count_builder = PrimitiveBuilder::<UInt32Type>::new();

    for (x, c) in items.iter() {
        kmer_builder.append_value(*x);
        count_builder.append_value(*c as u32);
    }

    let kmer_array = kmer_builder.finish();
    let count_array = count_builder.finish();

    let kmer_schema = Arc::new(Schema::new(vec![
        Field::new("kmer", DataType::UInt64, false),
        Field::new("count", DataType::UInt32, false),
    ]));

    let recs = RecordBatch::try_new(
        kmer_schema,
        vec![Arc::new(kmer_array), Arc::new(count_array)],
    )
    .map_err(|e| Error::new(ErrorKind::Other, e))?;

    let df = ctx.read_batch(recs)?;

    Ok(df)
}

fn _dot(lhs: &Vec<(u64, usize)>, rhs: &Vec<(u64, usize)>) -> f64 {
    let mut i = 0;
    let mut j = 0;
    let mut d = 0.0;
    let mut lmag = 0.0;
    let mut rmag = 0.0;
    while i < lhs.len() && j < rhs.len() {
        if lhs[i].0 < rhs[j].0 {
            lmag += (lhs[i].1 * lhs[i].1) as f64;
            i += 1;
            continue;
        }
        if lhs[i].0 > rhs[j].0 {
            rmag += (rhs[j].1 * rhs[j].1) as f64;
            j += 1;
            continue;
        }
        d += (lhs[i].1 * rhs[j].1) as f64;
        lmag += (lhs[i].1 * lhs[i].1) as f64;
        rmag += (rhs[j].1 * rhs[j].1) as f64;
        i += 1;
        j += 1;
    }
    d / (lmag.sqrt() * rmag.sqrt())
}
