use std::{
    collections::{HashMap, HashSet},
    io::{Error, ErrorKind},
    sync::Arc,
};

use datafusion::{
    arrow::{
        array::{PrimitiveBuilder, RecordBatch, UInt32Array},
        datatypes::{DataType, Field, Int64Type, Schema, UInt32Type},
    },
    common::JoinType,
    prelude::{DataFrame, SessionContext, lit},
};

use crate::disjoint_set::DisjointSet;

pub async fn merge_with(
    tbl: DataFrame,
    union: DataFrame,
    ctx: &SessionContext,
) -> std::io::Result<DataFrame> {
    let updates = make_merge_table(union, ctx).await?;
    let tbl = tbl
        .join(
            updates,
            JoinType::Left,
            &["row_key"],
            &["orig_row_key"],
            None,
        )?
        .drop_columns(&["orig_row_key"])?
        .with_column_renamed("new_row_key", "row_key")?
        .with_column_renamed("new_vix_set", "vix_set")?
        .with_column_renamed("new_vix_count", "vix_count")?;

    Ok(tbl)
}

async fn make_merge_table(union: DataFrame, ctx: &SessionContext) -> std::io::Result<DataFrame> {
    let union =
        union.select_columns(&["lhs_row_key", "lhs_vix_set", "rhs_row_key", "rhs_vix_set"])?;
    let union = union.collect().await?;

    let itr = union.iter().flat_map(|recs| MergeIterator::new(recs));

    let mut all = HashMap::new();
    let mut sets = DisjointSet::new();
    for (lhs_row_key, lhs_vix_set, rhs_row_key, rhs_vix_set) in itr {
        all.insert(lhs_row_key, lhs_vix_set);
        all.insert(rhs_row_key, rhs_vix_set);
        let x = sets.find(lhs_row_key);
        let y = sets.find(rhs_row_key);
        if x != y {
            sets.union(x, y);
        }
    }

    let mut final_sets: HashMap<u32, u32> = HashMap::new();
    for item in all.iter() {
        let row_key = *item.0;
        let vix_set = *item.1;
        let x = sets.find(row_key);
        *final_sets.entry(x).or_default() |= vix_set;
    }

    let mut orig_row_key_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut new_row_key_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut new_vix_set_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut new_vix_count_builder = PrimitiveBuilder::<UInt32Type>::new();
    for (orig_row_key, _vix_set) in all.into_iter() {
        let new_row_key = sets.find(orig_row_key);
        let new_vix_set = *final_sets.get(&new_row_key).unwrap();
        if orig_row_key != new_row_key {
            orig_row_key_builder.append_value(orig_row_key);
            new_row_key_builder.append_value(new_row_key);
            new_vix_set_builder.append_value(new_vix_set);
            new_vix_count_builder.append_value(new_vix_set.count_ones());
        }
    }
    let orig_row_key_array = orig_row_key_builder.finish();
    let new_row_key_array = new_row_key_builder.finish();
    let new_vix_set_array = new_vix_set_builder.finish();
    let new_vix_count_array = new_vix_set_builder.finish();

    let schema = Arc::new(Schema::new(vec![
        Field::new("orig_row_key", DataType::UInt32, false),
        Field::new("new_row_key", DataType::UInt32, false),
        Field::new("new_vix_set", DataType::UInt32, false),
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
    .map_err(|e| Error::new(ErrorKind::Other, e))?;

    ctx.read_batch(recs)
        .map_err(|e| Error::new(ErrorKind::Other, e))
}

pub(crate) struct MergeIterator<'a> {
    pub(crate) lhs_row_key: &'a UInt32Array,
    pub(crate) lhs_vix_set: &'a UInt32Array,
    pub(crate) rhs_row_key: &'a UInt32Array,
    pub(crate) rhs_vix_set: &'a UInt32Array,
    pub(crate) i: usize,
}

impl<'a> MergeIterator<'a> {
    pub fn new(recs: &'a RecordBatch) -> MergeIterator<'a> {
        for field in recs.schema_ref().fields().iter() {
            log::debug!("field: {:?}", field);
        }
        let lhs_row_key = Self::get_array::<UInt32Array>(recs, "lhs_row_key");
        let lhs_vix_set = Self::get_array::<UInt32Array>(recs, "lhs_vix_set");
        let rhs_row_key = Self::get_array::<UInt32Array>(recs, "rhs_row_key");
        let rhs_vix_set = Self::get_array::<UInt32Array>(recs, "rhs_vix_set");
        MergeIterator {
            lhs_row_key,
            lhs_vix_set,
            rhs_row_key,
            rhs_vix_set,
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
    type Item = (u32, u32, u32, u32);

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
