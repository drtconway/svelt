use crate::features::FeatureIndex;
use crate::kmers_table::kmers_table;
use datafusion::{
    arrow::array::{Array as _, GenericStringArray, RecordBatch},
    functions_aggregate::{expr_fn::first_value, min_max::min},
    prelude::{DataFrame, SessionContext, col, concat, lit, substr, substring},
};
use std::time::Instant;

pub(crate) async fn find_classifications(
    batch: Vec<RecordBatch>,
    features: &str,
    ctx: &SessionContext,
) -> std::io::Result<Option<DataFrame>> {
    let n: usize = batch.iter().map(|recs| recs.num_rows()).sum();
    log::info!("number of sequences to classify: {}", n);
    if n == 0 {
        return Ok(None);
    }

    let idx = FeatureIndex::load(features, &ctx).await?;

    let k = idx.k();

    log::info!("classifying insertion sequences with '{}'", features);
    let now = Instant::now();
    let mut last = Instant::now();

    let itr = batch.iter().flat_map(|recs| MergeIterator::new(recs));

    let mut result: Option<DataFrame> = None;
    let mut total_ins_sequences: usize = 0;

    let mut curr: Option<DataFrame> = None;
    let mut seq_count: usize = 0;
    let mut curr_count: usize = 0;
    for rec in itr {
        if total_ins_sequences & 0xF == 0xF {
            if last.elapsed().as_secs_f64() > 10.0 {
                log::info!(
                    "classified {} sequences in {}s ({} seq/sec).",
                    total_ins_sequences,
                    now.elapsed().as_secs_f32(),
                    total_ins_sequences as f32 / now.elapsed().as_secs_f32(),
                );
                last = Instant::now();
            }
        }
        total_ins_sequences += 1;
        seq_count += 1;
        let (name, sequence) = rec;
        if sequence.len() > 16384 {
            log::info!("long insertion: {} ({}bp)", name, sequence.len());
        }

        let (block, count) = kmers_table(sequence, k, ctx).await?;
        let block = block
            .with_column("name", concat(vec![col("strand"), lit(name)]))?
            .drop_columns(&["strand"])?;

        curr = if let Some(block0) = curr {
            Some(block0.union(block)?)
        } else {
            Some(block)
        };
        curr_count += count;

        if curr_count > 1000000 {
            log::debug!(
                "ranking over {} sequences with {} k-mers.",
                seq_count,
                curr_count
            );
            let tbl = curr.take().unwrap();
            let res = idx.rank_many(tbl).await?;
            let res = materialise(res, &ctx).await?;
            let res = res
                .with_column("strand", substring(col("query_name"), lit(1), lit(1)))?
                .with_column("query_name", substr(col("query_name"), lit(2)))?;
            let res = res.aggregate(
                vec![col("query_name")],
                vec![
                    min(col("distance")).alias("distance"),
                    first_value(
                        col("subject_name"),
                        Some(vec![col("distance").sort(true, false)]),
                    )
                    .alias("feature"),
                    first_value(col("class"), Some(vec![col("distance").sort(true, false)]))
                        .alias("class"),
                    first_value(col("strand"), Some(vec![col("distance").sort(true, false)]))
                        .alias("strand"),
                ],
            )?;
            if let Some(result0) = result {
                result = Some(result0.union(res)?);
            } else {
                result = Some(res);
            }
            seq_count = 0;
            curr_count = 0;
        }
    }

    if let Some(tbl) = curr {
        log::debug!(
            "ranking over {} sequences with {} k-mers.",
            seq_count,
            curr_count
        );

        let res = idx.rank_many(tbl).await?;
        let res = materialise(res, &ctx).await?;
        let res = res
            .with_column("strand", substring(col("query_name"), lit(1), lit(1)))?
            .with_column("query_name", substr(col("query_name"), lit(2)))?;
        let res = res.aggregate(
            vec![col("query_name")],
            vec![
                min(col("distance")).alias("distance"),
                first_value(
                    col("subject_name"),
                    Some(vec![col("distance").sort(true, false)]),
                )
                .alias("feature"),
                first_value(col("class"), Some(vec![col("distance").sort(true, false)]))
                    .alias("class"),
                first_value(col("strand"), Some(vec![col("distance").sort(true, false)]))
                    .alias("strand"),
            ],
        )?;
        if let Some(result0) = result {
            result = Some(result0.union(res)?);
        } else {
            result = Some(res);
        }
    }

    if false {
        if let Some(tbl) = &result {
            tbl.clone().show().await?;
        } else {
            log::info!("no results!");
        }
    }

    log::info!(
        "classified {} sequences in {}s.",
        total_ins_sequences,
        now.elapsed().as_secs_f32()
    );

    Ok(result)
}

pub(crate) async fn materialise(df: DataFrame, ctx: &SessionContext) -> std::io::Result<DataFrame> {
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
