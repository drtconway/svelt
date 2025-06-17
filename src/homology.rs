use std::{
    collections::HashMap,
    io::{BufReader, Error, ErrorKind},
    sync::Arc,
};

use autocompress::autodetect_open;
use datafusion::{
    arrow::{
        array::{PrimitiveBuilder, RecordBatch},
        datatypes::{DataType, Field, Schema, UInt32Type, UInt64Type},
    },
    config::CsvOptions,
    dataframe::DataFrameWriteOptions,
    prelude::{SessionContext, col},
};
use itertools::Itertools;
use noodles::fasta;

use crate::kmers::Kmer;

const K: usize = 11;

pub async fn index_features(features_name: &str) -> std::io::Result<()> {
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
        Kmer::with_many_both(K, rec.sequence(), |x, _y| {
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

    let mut kmer_builder = PrimitiveBuilder::<UInt64Type>::new();
    let mut seq_num_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut count_builder = PrimitiveBuilder::<UInt32Type>::new();

    for (x, postings) in kmer_index.into_iter() {
        for (n, c) in postings.into_iter() {
            kmer_builder.append_value(x);
            seq_num_builder.append_value(n as u32);
            count_builder.append_value(c as u32);
        }
    }

    let kmer_array = kmer_builder.finish();
    let seq_num_array = seq_num_builder.finish();
    let count_array = count_builder.finish();

    let schema = Arc::new(Schema::new(vec![
        Field::new("kmer", DataType::UInt64, false),
        Field::new("seq_num", DataType::UInt32, false),
        Field::new("count", DataType::UInt32, false),
    ]));

    let recs = RecordBatch::try_new(
        schema,
        vec![
            Arc::new(kmer_array),
            Arc::new(seq_num_array),
            Arc::new(count_array),
        ],
    )
    .map_err(|e| Error::new(ErrorKind::Other, e))?;

    let ctx = SessionContext::new();
    let df = ctx.read_batch(recs)?;

    let opts = DataFrameWriteOptions::default();
    let csv_opts = CsvOptions::default().with_delimiter(b'\t');
    let csv_opts = Some(csv_opts);

    df.sort_by(vec![col("kmer"), col("seq_num")])?
        .write_csv("feature-index.tsv", opts, csv_opts)
        .await?;

    Ok(())
}

fn dot(lhs: &Vec<(u64, usize)>, rhs: &Vec<(u64, usize)>) -> f64 {
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
