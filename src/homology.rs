use std::{
    collections::HashMap,
    io::{BufReader, Error, ErrorKind},
    sync::Arc,
};

use autocompress::autodetect_open;
use datafusion::{
    arrow::{
        array::{PrimitiveBuilder, RecordBatch, StringDictionaryBuilder},
        datatypes::{DataType, Field, Schema, UInt32Type, UInt64Type},
    },
    config::{ParquetColumnOptions, TableParquetOptions},
    dataframe::DataFrameWriteOptions,
    prelude::col,
};
use itertools::Itertools;
use noodles::fasta;
use regex::Regex;

use crate::{
    either::Either::{self, Left, Right},
    kmers::Kmer,
    options::{CommonOptions, IndexingOptions, make_session_context},
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
        let name = String::from_utf8(rec.name().to_vec()).unwrap();

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
    .write_parquet(&format!("{}.nidx", out), opts, None)
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
        Field::new("seq_num", DataType::UInt32, false),
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

    df.sort_by(vec![col("kmer"), col("seq_num")])?
        .write_parquet(&format!("{}.kidx", out), opts, Some(kidx_opts))
        .await?;

    Ok(())
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
