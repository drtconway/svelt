use std::{
    io::{BufReader, Error, ErrorKind},
    sync::Arc,
};

use autocompress::{Processor, io::ProcessorReader};
use datafusion::{
    arrow::{
        array::{GenericStringBuilder, RecordBatch},
        datatypes::{DataType, Field, Schema},
    },
    functions_aggregate::count::count,
    prelude::{DataFrame, SessionContext, col, lit},
};

use crate::kmers::kmerize::kmers_fwd;

pub trait SequenceIterator: Iterator<Item = std::io::Result<(String, String)>> {}

type InnerReader = BufReader<
    ProcessorReader<Box<dyn Processor + Send + Unpin + 'static>, BufReader<std::fs::File>>,
>;

pub mod fasta;
pub mod vcf;

pub async fn make_kmer_table<Itr: SequenceIterator>(
    k: usize,
    itr: Itr,
    ctx: &SessionContext,
) -> std::io::Result<DataFrame> {
    let kmers = kmers_fwd();

    let mut name_builder = GenericStringBuilder::<i32>::new();
    let mut seq_builder = GenericStringBuilder::<i32>::new();
    for item in itr {
        let (name, seq) = item?;
        name_builder.append_value(name);
        seq_builder.append_value(seq);
    }

    let name_array = name_builder.finish();
    let seq_array = seq_builder.finish();

    let schema = Arc::new(Schema::new(vec![
        Field::new("name", DataType::Utf8, true),
        Field::new("seq", DataType::Utf8, true),
    ]));

    let recs = RecordBatch::try_new(schema, vec![Arc::new(name_array), Arc::new(seq_array)])
        .map_err(|e| Error::new(ErrorKind::Other, e))?;

    let df = ctx.read_batch(recs)?.distinct()?;

    let df = df
        .with_column("kmer", kmers.call(vec![lit(k as u32), col("seq")]))?
        .drop_columns(&["seq"])?
        .unnest_columns(&["kmer"])?
        .aggregate(
            vec![col("name"), col("kmer")],
            vec![count(lit(1)).alias("count")],
        )?;

    df.clone().show().await?;

    Ok(df)
}
