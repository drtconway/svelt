use std::{
    io::{Error, ErrorKind},
    sync::Arc,
};

use datafusion::{
    arrow::{
        array::{PrimitiveBuilder, RecordBatch},
        datatypes::{DataType, Field, Int32Type, Schema, UInt64Type},
    },
    prelude::{lit, DataFrame, SessionContext},
};
use itertools::Itertools as _;

use crate::kmers::Kmer;

pub async fn kmers_table(seq: &str, k: usize, ctx: &SessionContext) -> std::io::Result<(DataFrame, usize)> {
    let mut fwd = Vec::new();
    let mut rev = Vec::new();
    Kmer::with_many_both(k, &seq, |x, y| {
        fwd.push(x.0);
        rev.push(y.0);
    });

    let mut kmer_count = 0;

    fwd.sort();
    let fwd: Vec<(u64, usize)> = fwd
        .into_iter()
        .chunk_by(|x| *x)
        .into_iter()
        .map(|(x, xs)| (x, xs.count()))
        .collect();
    kmer_count += fwd.len();
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
    kmer_count += rev.len();
    let rev = kmer_frequencies_to_table(&rev, &ctx)
        .await?
        .with_column("strand", lit("-"))?;

    let res = fwd.union(rev)?;

    Ok((res, kmer_count))
}

pub async fn kmer_frequencies_to_table(
    items: &Vec<(u64, usize)>,
    ctx: &SessionContext,
) -> std::io::Result<DataFrame> {
    let mut kmer_builder = PrimitiveBuilder::<UInt64Type>::new();
    let mut count_builder = PrimitiveBuilder::<Int32Type>::new();

    for (x, c) in items.iter() {
        kmer_builder.append_value(*x);
        count_builder.append_value(*c as i32);
    }

    let kmer_array = kmer_builder.finish();
    let count_array = count_builder.finish();

    let kmer_schema = Arc::new(Schema::new(vec![
        Field::new("kmer", DataType::UInt64, false),
        Field::new("count", DataType::Int32, false),
    ]));

    let recs = RecordBatch::try_new(
        kmer_schema,
        vec![Arc::new(kmer_array), Arc::new(count_array)],
    )
    .map_err(|e| Error::new(ErrorKind::Other, e))?;

    let df = ctx.read_batch(recs)?;

    Ok(df)
}
