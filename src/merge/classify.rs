use crate::features::FeatureIndex;
use datafusion::{
    arrow::{
        array::{
            Array as _, GenericStringArray, GenericStringBuilder, PrimitiveBuilder, RecordBatch,
        },
        datatypes::{DataType, Field, Float64Type, Schema},
    },
    prelude::{DataFrame, SessionContext},
};
use std::{
    io::{Error, ErrorKind}, sync::Arc, time::Instant, u32
};

pub(crate) async fn find_classifications(
    batch: Vec<RecordBatch>,
    features: &str,
    ctx: &SessionContext,
) -> std::io::Result<DataFrame> {
    let n: usize = batch.iter().map(|recs| recs.num_rows()).sum();
    log::info!("number of sequences to classify: {}", n);

    let idx = FeatureIndex::load(features, ctx).await?;

    log::info!("classifying insertion sequences with '{}'", features);

    let now = Instant::now();

    let itr = batch.iter().flat_map(|recs| MergeIterator::new(recs));

    let mut seq_hash_builder = GenericStringBuilder::<i32>::new();
    let mut class_builder = GenericStringBuilder::<i32>::new();
    let mut strand_builder = GenericStringBuilder::<i32>::new();
    let mut distance_builder = PrimitiveBuilder::<Float64Type>::new();

    for (seq_hash, sequence) in itr {
        let (fwd, rev) = idx.rank(sequence)?;
        let mut best_nix = u32::MAX;
        let mut best_score = -1.0;
        let mut best_strand = true;
        for (nix, score) in fwd {
            if score > best_score {
                best_nix = nix;
                best_score = score;
                best_strand = true;
            }
        }
        for (nix, score) in rev {
            if score > best_score {
                best_nix = nix;
                best_score = score;
                best_strand = false;
            }
        }

        if best_score > 0.5 {
            let class = &idx.names[best_nix as usize];
            seq_hash_builder.append_value(seq_hash);
            class_builder.append_value(class);
            strand_builder.append_value(if best_strand { "+" } else { "-" });
            distance_builder.append_value(1.0 - best_score);
        }
    }

    let dt = now.elapsed().as_secs_f64();
    let sps = (n as f64) / dt;
    log::info!("insertion sequence classification took {:2.2}s ({:2.2} sequences/s)", dt, sps);

    let seq_hash_array = seq_hash_builder.finish();
    let class_array = class_builder.finish();
    let strand_array = strand_builder.finish();
    let distance_array = distance_builder.finish();

    let schema = Arc::new(Schema::new(vec![
        Field::new("query_name", DataType::Utf8, false),
        Field::new("class", DataType::Utf8, false),
        Field::new("strand", DataType::Utf8, false),
        Field::new("distance", DataType::Float64, false),
    ]));

    let recs = RecordBatch::try_new(
        schema,
        vec![
            Arc::new(seq_hash_array),
            Arc::new(class_array),
            Arc::new(strand_array),
            Arc::new(distance_array),
        ],
    )
    .map_err(|e| Error::new(ErrorKind::Other, e))?;

    let result = ctx.read_batch(recs)?;

    Ok(result)
}

pub(crate) async fn _materialise(df: DataFrame, ctx: &SessionContext) -> std::io::Result<DataFrame> {
    let batches = df.collect().await?;
    let df = ctx.read_batches(batches)?;
    Ok(df)
}

pub(crate) struct MergeIterator<'a> {
    pub(crate) seq_hash: &'a GenericStringArray<i32>,
    pub(crate) alt_seq: &'a GenericStringArray<i32>,
    pub(crate) i: usize,
}

impl<'a> MergeIterator<'a> {
    pub fn new(recs: &'a RecordBatch) -> MergeIterator<'a> {
        for field in recs.schema_ref().fields().iter() {
            log::debug!("field: {:?}", field);
        }
        let seq_hash = Self::get_array::<GenericStringArray<i32>>(recs, "seq_hash");
        let alt_seq = Self::get_array::<GenericStringArray<i32>>(recs, "alt_seq");
        MergeIterator {
            seq_hash,
            alt_seq,
            i: 0,
        }
    }

    pub(crate) fn get_array<Type: 'static>(recs: &'a RecordBatch, name: &str) -> &'a Type {
        recs.column_by_name(name)
            .unwrap()
            .as_any()
            .downcast_ref::<Type>()
            .unwrap()
    }
}

impl<'a> Iterator for MergeIterator<'a> {
    type Item = (&'a str, &'a str);

    fn next(&mut self) -> Option<Self::Item> {
        if self.i < self.seq_hash.len() {
            let i = self.i;
            self.i += 1;

            let seq_hash = self.seq_hash.value(i);
            let alt_seq = self.alt_seq.value(i);
            Some((seq_hash, alt_seq))
        } else {
            None
        }
    }
}
