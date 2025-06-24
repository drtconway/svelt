use std::{
    collections::HashMap,
    io::{Error, ErrorKind},
    sync::Arc,
};

use datafusion::{
    arrow::{
        array::{PrimitiveBuilder, RecordBatch, UInt32Array},
        datatypes::{DataType, Field, Schema, UInt32Type},
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
    let flip = criterion == "flip";

    let mut updates = make_merge_table(union, ctx)
        .await?
        .with_column("new_criterion", lit(criterion))?;
    if flip {
        updates = updates.with_column("new_flip", lit(true))?;
    }

    if false {
        for field in tbl.schema().fields().iter() {
            log::info!("tbl field: {:?}", field);
        }

        for field in updates.schema().fields().iter() {
            log::info!("update field: {:?}", field);
        }
    }

    let mut tbl = tbl
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

    // If we're flipping, then we need to swap over the chrom/chrom2 end/end2
    //
    if flip {
        tbl = tbl
            .with_column("flip", coalesce(vec![col("new_flip"), col("flip")]))?
            .drop_columns(&["new_flip"])?;

        let rhs = tbl.clone().filter(col("flip").eq(lit(true)))?.select(vec![
            col("row_id").alias("flip_row_id"),
            col("chrom").alias("flip_chrom"),
            col("chrom_id").alias("flip_chrom_id"),
            col("end").alias("flip_end"),
            col("chrom2").alias("flip_chrom2"),
            col("chrom2_id").alias("flip_chrom2_id"),
            col("end2").alias("flip_end2"),
        ])?;

        tbl = tbl
            .join(rhs, JoinType::Left, &["row_id"], &["flip_row_id"], None)?
            .with_column("chrom", coalesce(vec![col("flip_chrom2"), col("chrom")]))?
            .with_column(
                "chrom_id",
                coalesce(vec![col("flip_chrom2_id"), col("chrom_id")]),
            )?
            .with_column("start", coalesce(vec![col("flip_end2"), col("start")]))?
            .with_column("end", coalesce(vec![col("flip_end2"), col("end")]))?
            .with_column("chrom2", coalesce(vec![col("flip_chrom"), col("chrom2")]))?
            .with_column(
                "chrom2_id",
                coalesce(vec![col("flip_chrom_id"), col("chrom2_id")]),
            )?
            .with_column("end2", coalesce(vec![col("flip_end"), col("end2")]))?
            .drop_columns(&[
                "flip_row_id",
                "flip_chrom",
                "flip_chrom_id",
                "flip_end",
                "flip_chrom2",
                "flip_chrom2_id",
                "flip_end2",
            ])?;
    }

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

    log::info!("updating merge information for {} entries", all.len());

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
    .unwrap();
    //.map_err(|e| Error::new(ErrorKind::Other, e))?;

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
