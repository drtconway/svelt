use std::{
    collections::HashMap,
    io::{Error, ErrorKind},
    sync::Arc,
};

use datafusion::{
    arrow::{
        array::{Int64Array, PrimitiveBuilder, RecordBatch},
        datatypes::{DataType, Field, Schema, UInt32Type},
    },
    common::HashSet,
    functions_aggregate::count::count,
    prelude::{DataFrame, Expr, SessionContext, abs, case, cast, coalesce, col, lit},
};

use crate::disjoint_set::DisjointSet;

pub async fn find_almost_exact(
    ctx: &SessionContext,
    tbl: DataFrame,
    n: u32,
) -> std::io::Result<DataFrame> {
    let tbl = find_almost_exact_non_bnd(ctx, tbl, n).await?;
    find_almost_exact_bnd(ctx, tbl, n).await
}

pub async fn find_almost_exact_non_bnd(
    ctx: &SessionContext,
    tbl: DataFrame,
    n: u32,
) -> std::io::Result<DataFrame> {
    let lhs = tbl
        .clone()
        .filter(col("vix_count").lt(lit(n)))?
        .select(vec![
            col("row_id").alias("lhs_row_id"),
            col("chrom_id").alias("lhs_chrom_id"),
            col("start").alias("lhs_start"),
            col("end").alias("lhs_end"),
            col("kind").alias("lhs_kind"),
            col("length").alias("lhs_length"),
            abs(col("length")).alias("lhs_abs_length"),
            col("row_key").alias("lhs_row_key"),
            col("vix_count").alias("lhs_vix_count"),
            col("vix_set").alias("lhs_vix_set"),
        ])?;
    let rhs = tbl
        .clone()
        .filter(col("vix_count").lt(lit(n)))?
        .select(vec![
            col("row_id").alias("rhs_row_id"),
            col("chrom_id").alias("rhs_chrom_id"),
            col("start").alias("rhs_start"),
            col("end").alias("rhs_end"),
            col("kind").alias("rhs_kind"),
            col("length").alias("rhs_length"),
            abs(col("length")).alias("rhs_abs_length"),
            col("row_key").alias("rhs_row_key"),
            col("vix_count").alias("rhs_vix_count"),
            col("vix_set").alias("rhs_vix_set"),
        ])?;

    //lhs.clone().show().await?;

    let almost_exact = lhs
        .join(
            rhs,
            datafusion::common::JoinType::Inner,
            &["lhs_chrom_id", "lhs_kind"],
            &["rhs_chrom_id", "rhs_kind"],
            Some(
                lit(true)
                    .and(col("lhs_row_id").lt(col("rhs_row_id")))
                    .and(abs(col("lhs_start") - col("rhs_start")).lt(lit(25)))
                    .and(abs(col("lhs_end") - col("rhs_end")).lt(lit(25)))
                    .and((col("lhs_vix_set") & col("rhs_vix_set")).eq(lit(0))),
            ),
        )?
        //.with_column("min_start", pmin(col("lhs_start"), col("rhs_start"))?)?
        //.with_column("max_start", pmax(col("lhs_start"), col("rhs_start"))?)?
        //.with_column("min_end", pmin(col("lhs_end"), col("rhs_end"))?)?
        //.with_column("max_end", pmax(col("lhs_end"), col("rhs_end"))?)?
        .with_column(
            "max_abs_length",
            pmax(col("lhs_abs_length"), col("rhs_abs_length"))?,
        )?
        .with_column(
            "min_abs_length",
            pmin(col("lhs_abs_length"), col("rhs_abs_length"))?,
        )?
        .with_column(
            "length_ratio",
            cast(col("min_abs_length"), DataType::Float64)
                / cast(col("max_abs_length"), DataType::Float64),
        )?
        .drop_columns(&[
            //"min_start",
            //"max_start",
            //"min_end",
            //"max_end",
            "max_abs_length",
            "min_abs_length",
        ])?
        .filter(col("length_ratio").gt(lit(0.9)))?
        .sort(vec![
            (col("lhs_vix_count") + col("rhs_vix_count")).sort(false, false),
            col("length_ratio").sort(false, false),
            col("lhs_row_key").sort(true, false),
            col("rhs_row_key").sort(true, false),
        ])?;

    if false {
        almost_exact.clone().show().await?;
    }

    let batch: Vec<RecordBatch> = almost_exact.collect().await?;

    let resolution = resolve_groups(batch).await?;
    let resolution = ctx
        .read_batch(resolution)
        .map_err(|e| Error::new(ErrorKind::Other, e))?;

    if false {
        resolution.clone().show().await?;
    }

    if false {
        tbl.clone()
            .aggregate(
                vec![col("row_id")],
                vec![count(col("row_key")).alias("count")],
            )?
            .filter(col("count").gt(lit(1)))?
            .sort_by(vec![col("row_id")])?
            .show()
            .await?;
    }

    let tbl = tbl
        .join(
            resolution,
            datafusion::common::JoinType::Left,
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
        .drop_columns(&["orig_row_id", "new_row_key", "new_vix_count", "new_vix_set"])?;

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

pub async fn find_almost_exact_bnd(
    ctx: &SessionContext,
    tbl: DataFrame,
    n: u32,
) -> std::io::Result<DataFrame> {
    let lhs = tbl
        .clone()
        .filter(col("vix_count").lt(lit(n)).and(col("kind").eq(lit("BND"))))?
        .select(vec![
            col("row_id").alias("lhs_row_id"),
            col("chrom_id").alias("lhs_chrom_id"),
            col("start").alias("lhs_start"),
            col("end").alias("lhs_end"),
            col("chrom2_id").alias("lhs_chrom2_id"),
            col("end2").alias("lhs_end2"),
            col("row_key").alias("lhs_row_key"),
            col("vix_count").alias("lhs_vix_count"),
            col("vix_set").alias("lhs_vix_set"),
        ])?;
    let rhs = tbl
        .clone()
        .filter(col("vix_count").lt(lit(n)).and(col("kind").eq(lit("BND"))))?
        .select(vec![
            col("row_id").alias("rhs_row_id"),
            col("chrom_id").alias("rhs_chrom_id"),
            col("start").alias("rhs_start"),
            col("end").alias("rhs_end"),
            col("chrom2_id").alias("rhs_chrom2_id"),
            col("end2").alias("rhs_end2"),
            col("row_key").alias("rhs_row_key"),
            col("vix_count").alias("rhs_vix_count"),
            col("vix_set").alias("rhs_vix_set"),
        ])?;

    //lhs.clone().show().await?;

    let almost_exact = lhs
        .join(
            rhs,
            datafusion::common::JoinType::Inner,
            &["lhs_chrom_id", "lhs_chrom2_id"],
            &["rhs_chrom_id", "rhs_chrom2_id"],
            Some(
                lit(true)
                    .and(col("lhs_row_id").lt(col("rhs_row_id")))
                    .and(abs(col("lhs_start") - col("rhs_start")).lt(lit(25)))
                    .and(abs(col("lhs_end") - col("rhs_end")).lt(lit(25)))
                    .and((col("lhs_vix_set") & col("rhs_vix_set")).eq(lit(0))),
            ),
        )?
        .sort(vec![
            (col("lhs_vix_count") + col("rhs_vix_count")).sort(false, false),
            col("lhs_row_key").sort(true, false),
            col("rhs_row_key").sort(true, false),
        ])?;

    if false {
        almost_exact.clone().show().await?;
    }

    let batch: Vec<RecordBatch> = almost_exact.collect().await?;

    let resolution = resolve_groups(batch).await?;
    let resolution = ctx
        .read_batch(resolution)
        .map_err(|e| Error::new(ErrorKind::Other, e))?;

    if true {
        resolution.clone().sort_by(vec![col("new_row_key")])?.show().await?;
    }

    if false {
        tbl.clone()
            .aggregate(
                vec![col("row_id")],
                vec![count(col("row_key")).alias("count")],
            )?
            .filter(col("count").gt(lit(1)))?
            .sort_by(vec![col("row_id")])?
            .show()
            .await?;
    }

    let tbl = tbl
        .join(
            resolution,
            datafusion::common::JoinType::Left,
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
        .drop_columns(&["orig_row_id", "new_row_key", "new_vix_count", "new_vix_set"])?;

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

fn ifelse(cond: Expr, then: Expr, otherwise: Expr) -> datafusion::common::Result<Expr> {
    case(cond)
        .when(lit(true), then)
        .when(lit(false), otherwise)
        .end()
}

fn pmin(lhs: Expr, rhs: Expr) -> datafusion::common::Result<Expr> {
    ifelse(lhs.clone().lt_eq(rhs.clone()), lhs, rhs)
}

fn pmax(lhs: Expr, rhs: Expr) -> datafusion::common::Result<Expr> {
    ifelse(lhs.clone().gt_eq(rhs.clone()), lhs, rhs)
}

async fn resolve_groups(batch: Vec<RecordBatch>) -> std::io::Result<RecordBatch> {
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
            .unwrap();
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
