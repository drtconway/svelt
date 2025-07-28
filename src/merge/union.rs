use std::{
    collections::{HashMap, HashSet},
    io::{Error, ErrorKind},
    sync::Arc,
};

use datafusion::{
    arrow::{
        array::{PrimitiveBuilder, RecordBatch, UInt32Array, UInt64Array},
        datatypes::{DataType, Field, Schema, UInt32Type, UInt64Type},
    },
    common::JoinType,
    prelude::{DataFrame, SessionContext, coalesce, col, concat_ws, lit, nullif},
};

use crate::disjoint_set::DisjointSet;

pub async fn merge_with(
    tbl: DataFrame,
    union: DataFrame,
    ctx: &SessionContext,
    criterion: &str,
) -> std::io::Result<DataFrame> {
    let updates = make_merge_table(union, ctx)
        .await?
        .with_column("new_criterion", lit(criterion))?;

    if false {
        for field in tbl.schema().fields().iter() {
            log::info!("tbl field: {:?}", field);
        }

        for field in updates.schema().fields().iter() {
            log::info!("update field: {:?}", field);
        }
    }

    let tbl = tbl
        .join(
            updates,
            JoinType::Left,
            &["row_key"],
            &["orig_row_key"],
            None,
        )?
        .with_column(
            "row_key",
            coalesce(vec![col("new_row_key"), col("row_key")]),
        )?
        .with_column(
            "vix_set",
            coalesce(vec![col("new_vix_set"), col("vix_set")]),
        )?
        .with_column(
            "vix_count",
            coalesce(vec![col("new_vix_count"), col("vix_count")]),
        )?
        .with_column(
            "criteria",
            concat_ws(
                lit(","),
                vec![nullif(col("criteria"), lit("")), col("new_criterion")],
            ),
        )?
        .drop_columns(&[
            "orig_row_key",
            "new_row_key",
            "new_vix_set",
            "new_vix_count",
            "new_criterion",
        ])?;

    if false {
        tbl.clone()
            .drop_columns(&["alt_seq", "seq_hash"])?
            .sort_by(vec![col("chrom_id"), col("start"), col("row_id")])?
            .show()
            .await?;
    }

    Ok(tbl)
}

async fn make_merge_table(union: DataFrame, ctx: &SessionContext) -> std::io::Result<DataFrame> {
    let union =
        union.select_columns(&["lhs_row_key", "lhs_vix_set", "rhs_row_key", "rhs_vix_set"])?;

    if false {
        union.clone().show().await?;
    }

    let union = union.collect().await?;

    let itr = union.iter().flat_map(|recs| MergeIterator::new(recs));

    let mut updated_row_keys = HashSet::new();
    let mut vix_set_index = HashMap::new();
    let mut sets = DisjointSet::new();
    for (lhs_row_key, lhs_vix_set, rhs_row_key, rhs_vix_set) in itr {
        let x = sets.find(lhs_row_key);
        let y = sets.find(rhs_row_key);
        if x != y {
            // First, make sure we're not about to merge variants from the same VCF
            //
            let x_vix_set = if let Some(vix_set) = vix_set_index.get(&x) {
                *vix_set
            } else {
                lhs_vix_set
            };
            let y_vix_set = if let Some(vix_set) = vix_set_index.get(&y) {
                *vix_set
            } else {
                rhs_vix_set
            };
            if x_vix_set & y_vix_set != 0 {
                continue;
            }

            updated_row_keys.insert(lhs_row_key);
            updated_row_keys.insert(rhs_row_key);

            let z = sets.union(x, y);
            let vix_set = x_vix_set | y_vix_set;

            vix_set_index.remove(&x);
            vix_set_index.remove(&y);
            vix_set_index.insert(z, vix_set);
        }
    }

    log::info!(
        "updating merge information for {} entries",
        updated_row_keys.len()
    );

    let mut orig_row_key_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut new_row_key_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut new_vix_set_builder = PrimitiveBuilder::<UInt64Type>::new();
    let mut new_vix_count_builder = PrimitiveBuilder::<UInt32Type>::new();
    for orig_row_key in updated_row_keys.into_iter() {
        let new_row_key = sets.find(orig_row_key);
        let new_vix_set = *vix_set_index.get(&new_row_key).unwrap();
        orig_row_key_builder.append_value(orig_row_key);
        new_row_key_builder.append_value(new_row_key);
        new_vix_set_builder.append_value(new_vix_set);
        new_vix_count_builder.append_value(new_vix_set.count_ones());
    }
    let orig_row_key_array = orig_row_key_builder.finish();
    let new_row_key_array = new_row_key_builder.finish();
    let new_vix_set_array = new_vix_set_builder.finish();
    let new_vix_count_array = new_vix_count_builder.finish();

    let schema = Arc::new(Schema::new(vec![
        Field::new("orig_row_key", DataType::UInt32, false),
        Field::new("new_row_key", DataType::UInt32, false),
        Field::new("new_vix_set", DataType::UInt64, false),
        Field::new("new_vix_count", DataType::UInt32, false),
    ]));

    let recs = RecordBatch::try_new(
        schema,
        vec![
            Arc::new(orig_row_key_array),
            Arc::new(new_row_key_array),
            Arc::new(new_vix_set_array),
            Arc::new(new_vix_count_array),
        ],
    )
    .unwrap();
    //.map_err(|e| Error::new(ErrorKind::Other, e))?;

    ctx.read_batch(recs)
        .map_err(|e| Error::new(ErrorKind::Other, e))
}

pub(crate) struct MergeIterator<'a> {
    pub(crate) lhs_row_key: &'a UInt32Array,
    pub(crate) lhs_vix_set: &'a UInt64Array,
    pub(crate) rhs_row_key: &'a UInt32Array,
    pub(crate) rhs_vix_set: &'a UInt64Array,
    pub(crate) i: usize,
}

impl<'a> MergeIterator<'a> {
    pub fn new(recs: &'a RecordBatch) -> MergeIterator<'a> {
        for field in recs.schema_ref().fields().iter() {
            log::debug!("field: {:?}", field);
        }
        let lhs_row_key = Self::get_array::<UInt32Array>(recs, "lhs_row_key");
        let lhs_vix_set = Self::get_array::<UInt64Array>(recs, "lhs_vix_set");
        let rhs_row_key = Self::get_array::<UInt32Array>(recs, "rhs_row_key");
        let rhs_vix_set = Self::get_array::<UInt64Array>(recs, "rhs_vix_set");
        MergeIterator {
            lhs_row_key,
            lhs_vix_set,
            rhs_row_key,
            rhs_vix_set,
            i: 0,
        }
    }

    pub(crate) fn get_array<Type: 'static>(recs: &'a RecordBatch, name: &str) -> &'a Type {
        if false {
            log::info!("getting {}", name);
        }
        recs.column_by_name(name)
            .unwrap()
            .as_any()
            .downcast_ref::<Type>()
            .unwrap()
    }
}

impl<'a> Iterator for MergeIterator<'a> {
    type Item = (u32, u64, u32, u64);

    fn next(&mut self) -> Option<Self::Item> {
        if self.i < self.lhs_vix_set.len() {
            let i = self.i;
            self.i += 1;

            let lhs_row_key = self.lhs_row_key.value(i);
            let lhs_vix_set = self.lhs_vix_set.value(i);
            let rhs_row_key = self.rhs_row_key.value(i);
            let rhs_vix_set = self.rhs_vix_set.value(i);
            Some((lhs_row_key, lhs_vix_set, rhs_row_key, rhs_vix_set))
        } else {
            None
        }
    }
}
