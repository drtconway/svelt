use std::{
    collections::HashMap,
    io::{Error, ErrorKind},
    sync::Arc,
};

use datafusion::{
    arrow::{
        array::{Int64Array, PrimitiveBuilder, RecordBatch},
        datatypes::{DataType, Field, Schema, UInt32Type},
    }, common::JoinType, prelude::{coalesce, col, concat_ws, lit, nullif, DataFrame}
};

use crate::disjoint_set::DisjointSet;

pub async fn resolve_groups(tbl: DataFrame) -> std::io::Result<RecordBatch> {
    let batch = tbl.collect().await?;

    let mut uf: DisjointSet<u32> = DisjointSet::new();
    let mut vix_sets: HashMap<u32, u32> = HashMap::new();
    let mut row_ids = HashMap::new();

    for recs in batch.into_iter() {
        for field in recs.schema_ref().fields().iter() {
            log::debug!("field: {:?}", field);
        }
        let lhs_row_id = recs
            .column_by_name("lhs_row_id")
            .unwrap()
            .as_any()
            .downcast_ref::<Int64Array>()
            .unwrap();
        let lhs_row_key = recs
            .column_by_name("lhs_row_key")
            .unwrap()
            .as_any()
            .downcast_ref::<Int64Array>()
            .unwrap();
        let lhs_vix_set = recs
            .column_by_name("lhs_vix_set")
            .unwrap()
            .as_any()
            .downcast_ref::<Int64Array>()
            .expect("downcast failed");
        let rhs_row_id = recs
            .column_by_name("rhs_row_id")
            .unwrap()
            .as_any()
            .downcast_ref::<Int64Array>()
            .unwrap();
        let rhs_row_key = recs
            .column_by_name("rhs_row_key")
            .unwrap()
            .as_any()
            .downcast_ref::<Int64Array>()
            .unwrap();
        let rhs_vix_set = recs
            .column_by_name("rhs_vix_set")
            .unwrap()
            .as_any()
            .downcast_ref::<Int64Array>()
            .unwrap();
        let z = lhs_row_key.len();
        for i in 0..z {
            let ln = lhs_row_id.value(i) as u32;
            let lk = lhs_row_key.value(i) as u32;
            let ls = lhs_vix_set.value(i) as u32;
            let rn = rhs_row_id.value(i) as u32;
            let rk = rhs_row_key.value(i) as u32;
            let rs = rhs_vix_set.value(i) as u32;

            let a = uf.find(lk);
            let b = uf.find(rk);
            if a == b {
                // They already got joined.
                continue;
            }

            let av = *vix_sets.entry(a).or_insert(ls);
            let bv = *vix_sets.entry(b).or_insert(rs);

            if av & bv != 0 {
                continue;
            }

            row_ids.insert(ln, lk);
            row_ids.insert(rn, rk);

            let c = uf.union(a, b);
            let cv = av | bv;
            vix_sets.insert(c, cv);
            if c == a {
                vix_sets.remove(&b);
            } else {
                vix_sets.remove(&a);
            }
        }
    }

    let mut orig_row_id_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut new_row_key_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut vix_count_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut vix_set_builder = PrimitiveBuilder::<UInt32Type>::new();

    for x in row_ids.into_iter() {
        let a = uf.find(x.1);
        let s = *vix_sets.get(&a).unwrap();
        orig_row_id_builder.append_value(x.0);
        new_row_key_builder.append_value(a);
        vix_count_builder.append_value(s.count_ones());
        vix_set_builder.append_value(s);
    }

    let orig_row_id_array = orig_row_id_builder.finish();
    let new_row_key_array = new_row_key_builder.finish();
    let vix_count_array = vix_count_builder.finish();
    let vix_set_array = vix_set_builder.finish();

    let schema = Arc::new(Schema::new(vec![
        Field::new("orig_row_id", DataType::UInt32, false),
        Field::new("new_row_key", DataType::UInt32, false),
        Field::new("new_vix_count", DataType::UInt32, false),
        Field::new("new_vix_set", DataType::UInt32, false),
    ]));

    let recs = RecordBatch::try_new(
        schema,
        vec![
            Arc::new(orig_row_id_array),
            Arc::new(new_row_key_array),
            Arc::new(vix_count_array),
            Arc::new(vix_set_array),
        ],
    )
    .map_err(|e| Error::new(ErrorKind::Other, e))?;

    Ok(recs)
}

pub async fn update_tables(
    tbl: DataFrame,
    resolution: DataFrame,
    criterion: &str,
) -> std::io::Result<DataFrame> {
    let resolution = resolution.with_column("new_criterion", lit(criterion))?;
    let tbl = tbl
        .join(
            resolution,
            JoinType::Left,
            &["row_id"],
            &["orig_row_id"],
            None,
        )?
        .with_column(
            "row_key",
            coalesce(vec![col("new_row_key"), col("row_key")]),
        )?
        .with_column(
            "vix_count",
            coalesce(vec![col("new_vix_count"), col("vix_count")]),
        )?
        .with_column(
            "vix_set",
            coalesce(vec![col("new_vix_set"), col("vix_set")]),
        )?
        .with_column(
            "criteria",
            nullif(concat_ws(lit(","), vec![col("criteria"), col("new_criterion")]), lit("")),
        )?
        .drop_columns(&["orig_row_id", "new_row_key", "new_vix_count", "new_vix_set", "new_criterion"])?;

    if false {
        tbl.clone()
            .sort_by(vec![
                col("chrom_id"),
                col("start"),
                col("end"),
                col("row_id"),
            ])?
            .show()
            .await?;
    }

    Ok(tbl)
}
