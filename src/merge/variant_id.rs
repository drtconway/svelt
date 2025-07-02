use std::{collections::HashMap, sync::Arc};

use datafusion::{
    arrow::{
        array::{
            ArrayRef, GenericStringArray, GenericStringBuilder, PrimitiveBuilder, RecordBatch,
            UInt32Array,
        },
        datatypes::{DataType, Field, Schema, UInt32Type},
    },
    common::{JoinType, cast::as_binary_array},
    error::DataFusionError,
    functions_aggregate::expr_fn::first_value,
    logical_expr::{ColumnarValue, ScalarUDF, Volatility},
    prelude::{DataFrame, SessionContext, SimpleScalarUDF, col, concat_ws, left, lit, sha256},
};

/// Generate values to populate the ID column.
///
/// There are competing priorities in constructing IDs:
///
/// 1. The IDs within a VCF must be unique.
///
/// 2. The IDs should not be too long, just for practical reasons.
///
/// 3. It is very desirable for the IDs to be such that two different variants in
///    different VCFs do not lead to the same ID.
///
/// The first constraint is the only hard constraint, but the other two factors
/// have substantial practical consequences. While it would be ideal for a given
/// variant to generate the same ID in any VCF in which it occurs (i.e. make them
/// globally unique and deterministic) this is impossible because some tools
/// produce duplicate IDs, so a globally unique and deterministic naming system
/// would lead to violations of the VCF specification (i.e. #1 above).
///
/// The compromise solution we employ is to generate an ID with an integer component
/// that allows multiple possible IDs to be generated for a given variant, while
/// making it vanishingly unlikely that two different variants will have the same ID.
///
pub async fn construct_variant_ids(
    orig: DataFrame,
    ctx: &SessionContext,
) -> std::io::Result<DataFrame> {
    let enc = base62();

    // First, select the variant that will be the representative variant reported
    // in the output, which is the "leftmost" variant for each row_key.
    let rhs = orig
        .clone()
        .aggregate(
            vec![col("row_key")],
            vec![
                first_value(col("row_id"), Some(vec![col("vix").sort(true, false)]))
                    .alias("left_row_id"),
            ],
        )?
        .drop_columns(&["row_key"])?;

    let primary = orig
        .clone()
        .join(rhs, JoinType::Inner, &["row_id"], &["left_row_id"], None)?
        .drop_columns(&["left_row_id"])?;

    // Now generate the basic ID values corresponding to the variant characteristics.
    let ids = primary
        .clone()
        .with_column(
            "vid",
            concat_ws(
                lit("_"),
                vec![
                    col("kind"),
                    col("chrom"),
                    col("start"),
                    col("end"),
                    col("length"),
                    col("chrom2"),
                    col("end2"),
                    col("seq_hash"),
                ],
            ),
        )?
        .with_column("vid_hash", left(enc.call(vec![sha256(col("vid"))]), lit(7)))?
        .select_columns(&["row_key", "vid", "vid_hash"])?
        .sort_by(vec![col("row_key")])?;

    // Now we're going to drop out of DataFusion to track the occurence number
    // of each ID, so we can number-apart any collisions.
    let id_batches = ids.collect().await?;

    let mut seen: HashMap<String, u32> = HashMap::new();

    let mut row_key_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut vid_rep_builder = PrimitiveBuilder::<UInt32Type>::new();

    for batch in id_batches.into_iter() {
        let row_keys = get_array::<UInt32Array>(&batch, "row_key");
        let vid_hashes = get_array::<GenericStringArray<i32>>(&batch, "vid_hash");

        for i in 0..row_keys.len() {
            let row_key = row_keys.value(i);
            let vid_hash = vid_hashes.value(i).to_string();
            let count = seen.entry(vid_hash).or_default();
            let so_far = *count;
            *count += 1;

            row_key_builder.append_value(row_key);
            vid_rep_builder.append_value(so_far);
        }
    }

    let row_key_array = row_key_builder.finish();
    let vid_rep_array = vid_rep_builder.finish();

    let schema = Arc::new(Schema::new(vec![
        Field::new("vid_row_key", DataType::UInt32, false),
        Field::new("vid_rep", DataType::UInt32, false),
    ]));

    let recs = RecordBatch::try_new(
        schema,
        vec![Arc::new(row_key_array), Arc::new(vid_rep_array)],
    )
    .unwrap();

    let vid_reps = ctx.read_batch(recs)?;

    // Now we're back in DataFusion, join the occurrence number back on to
    // the ID generation table, and make the final IDs.
    let ids = primary
        .join(
            vid_reps,
            JoinType::Inner,
            &["row_key"],
            &["vid_row_key"],
            None,
        )?
        .with_column(
            "vid",
            concat_ws(
                lit("_"),
                vec![
                    col("kind"),
                    col("chrom"),
                    col("start"),
                    col("end"),
                    col("length"),
                    col("chrom2"),
                    col("end2"),
                    col("seq_hash"),
                    col("vid_rep"),
                ],
            ),
        )?
        .with_column("vid_hash", left(enc.call(vec![sha256(col("vid"))]), lit(7)))?
        .with_column(
            "variant_id",
            concat_ws(lit("_"), vec![lit("SVELT"), col("kind"), col("vid_hash")]),
        )?
        .select(vec![col("row_key").alias("rhs_row_key"), col("variant_id")])?;

    if false {
        ids.clone()
            .sort_by(vec![col("rhs_row_key")])?
            .show()
            .await?;
    }

    let res = orig
        .join(ids, JoinType::Left, &["row_key"], &["rhs_row_key"], None)?
        .drop_columns(&["rhs_row_key"])?;

    Ok(res)
}

fn get_array<'a, Type: 'static>(recs: &'a RecordBatch, name: &str) -> &'a Type {
    if false {
        log::info!("getting {}", name);
    }
    recs.column_by_name(name)
        .unwrap()
        .as_any()
        .downcast_ref::<Type>()
        .unwrap()
}

fn base62() -> ScalarUDF {
    ScalarUDF::from(SimpleScalarUDF::new(
        "base62",
        vec![DataType::Binary],
        DataType::Utf8,
        Volatility::Immutable,
        Arc::new(|args| do_base62(args)),
    ))
}

fn do_base62(args: &[ColumnarValue]) -> std::result::Result<ColumnarValue, DataFusionError> {
    let args = ColumnarValue::values_to_arrays(args)?;

    let blobs = as_binary_array(&args[0])?;

    let mut result_builder = GenericStringBuilder::<i32>::new();

    for blob in blobs.iter() {
        if let Some(blob) = blob {
            let mut s = Vec::new();
            let mut n: usize = 0;
            let mut x: u128 = 0;
            for b in blob {
                x = (x << 8) | (*b as u128);
                n += 1;
                if n == 16 {
                    let v = base62::encode(x);
                    s.push(v);
                    n = 0;
                    x = 0;
                }
            }
            if n > 0 {
                let v = base62::encode(x);
                s.push(v);
            }
            let s = s.join("");
            result_builder.append_value(s);
        } else {
            result_builder.append_null();
        }
    }

    let result_array = result_builder.finish();

    Ok(ColumnarValue::from(Arc::new(result_array) as ArrayRef))
}
