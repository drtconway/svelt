use std::io::{Error, ErrorKind};

use datafusion::{
    arrow::datatypes::DataType,
    functions_aggregate::count::count,
    prelude::{DataFrame, SessionContext, abs, cast, col, lit},
};

use crate::resolve::resolve_groups;
use crate::{
    expressions::{pmax, pmin},
    resolve::update_tables,
};

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
    log::info!("resolving non-BND inexact matches");

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

    let resolution = resolve_groups(almost_exact).await?;
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

    let tbl = update_tables(tbl, resolution).await?;

    Ok(tbl)
}

pub async fn find_almost_exact_bnd(
    ctx: &SessionContext,
    tbl: DataFrame,
    n: u32,
) -> std::io::Result<DataFrame> {
    log::info!("resolving BND inexact matches");

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

    let resolution = resolve_groups(almost_exact).await?;
    let resolution = ctx
        .read_batch(resolution)
        .map_err(|e| Error::new(ErrorKind::Other, e))?;

    if false {
        resolution
            .clone()
            .sort_by(vec![col("new_row_key")])?
            .show()
            .await?;
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

    let tbl = update_tables(tbl, resolution).await?;

    Ok(tbl)
}
