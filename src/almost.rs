use std::{
    io::{Error, ErrorKind},
    time::Instant,
};

use datafusion::{
    common::JoinType,
    prelude::{DataFrame, SessionContext, abs, col, greatest, least, lit},
};

use crate::{
    expressions::prefix_cols,
    near_merge::near_merge_table,
    options::MergeOptions,
    resolve::{resolve_groups, update_tables},
};

pub async fn find_almost_exact(
    ctx: &SessionContext,
    tbl: DataFrame,
    n: u32,
    options: &MergeOptions,
) -> std::io::Result<DataFrame> {
    let tbl = find_almost_exact_non_bnd(ctx, tbl, n, options).await?;
    find_almost_exact_bnd(ctx, tbl, n, options).await
}

pub async fn find_almost_exact_non_bnd(
    ctx: &SessionContext,
    tbl: DataFrame,
    n: u32,
    options: &MergeOptions,
) -> std::io::Result<DataFrame> {
    log::info!("resolving non-BND inexact matches");
    let now = Instant::now();
    let candidates = tbl
        .clone()
        .filter(
            lit(true)
                .and(col("kind").not_eq(lit("BND")))
                .and(col("vix_count").lt(lit(n))),
        )?
        .with_column("abs_length", abs(col("length")))?;

    // WIP code for ~linear range join.
    let middle = near_merge_table(candidates.clone()).await?;
    let middle = ctx
        .read_batch(middle)
        .map_err(|e| Error::new(ErrorKind::Other, e))?;

    let lhs = prefix_cols(candidates.clone(), "lhs")?;
    let rhs = prefix_cols(candidates.clone(), "rhs")?;

    let almost_exact = lhs
        .join(
            middle,
            JoinType::Left,
            &["lhs_row_id"],
            &["lhs_merge_row_id"],
            None,
        )?
        .join(
            rhs,
            JoinType::Left,
            &["rhs_merge_row_id"],
            &["rhs_row_id"],
            None,
        )?
        .with_column(
            "length_ratio",
            (lit(1.0) * least(vec![col("lhs_abs_length"), col("rhs_abs_length")]))
                / greatest(vec![col("lhs_abs_length"), col("rhs_abs_length")]),
        )?;

    let almost_exact = almost_exact
        .filter(
            lit(true)
                .and(col("lhs_row_id").lt(col("rhs_row_id")))
                .and(col("lhs_row_key").not_eq(col("rhs_row_key")))
                .and((col("lhs_vix_set") & col("rhs_vix_set")).eq(lit(0)))
                .and(abs(col("lhs_start") - col("rhs_start")).lt(lit(options.position_window)))
                .and(abs(col("lhs_end") - col("rhs_end")).lt(lit(options.position_window)))
                .and(col("length_ratio").gt_eq(lit(options.length_ratio))),
        )?
        .sort(vec![
            (col("lhs_vix_count") + col("rhs_vix_count")).sort(false, false),
            col("length_ratio").sort(false, false),
            col("lhs_row_key").sort(true, false),
            col("rhs_row_key").sort(true, false),
        ])?;

    let resolution = resolve_groups(almost_exact).await?;
    let resolution = ctx
        .read_batch(resolution)
        .map_err(|e| Error::new(ErrorKind::Other, e))?;

    let tbl = update_tables(tbl, resolution, "near").await?;

    log::info!("near_merge too {}s", now.elapsed().as_secs_f32());

    Ok(tbl)
}

pub async fn find_almost_exact_bnd(
    ctx: &SessionContext,
    tbl: DataFrame,
    n: u32,
    options: &MergeOptions,
) -> std::io::Result<DataFrame> {
    log::info!("resolving BND inexact matches");

    let candidates = tbl.clone().filter(col("vix_count").lt(lit(n)))?;
    let lhs = prefix_cols(candidates.clone(), "lhs")?;
    let rhs = prefix_cols(candidates.clone(), "rhs")?;
    //lhs.clone().show().await?;

    let almost_exact = lhs
        .join(
            rhs,
            JoinType::Inner,
            &["lhs_chrom_id", "lhs_chrom2_id"],
            &["rhs_chrom_id", "rhs_chrom2_id"],
            Some(
                lit(true)
                    .and(col("lhs_row_id").lt(col("rhs_row_id")))
                    .and(abs(col("lhs_start") - col("rhs_start")).lt(lit(options.position_window)))
                    .and(abs(col("lhs_end") - col("rhs_end")).lt(lit(options.position_window)))
                    .and(abs(col("lhs_end2") - col("rhs_end2")).lt(lit(options.end2_window)))
                    .and((col("lhs_vix_set") & col("rhs_vix_set")).eq(lit(0))),
            ),
        )?
        .sort(vec![
            (col("lhs_vix_count") + col("rhs_vix_count")).sort(false, false),
            col("lhs_row_key").sort(true, false),
            col("rhs_row_key").sort(true, false),
        ])?;

    let resolution = resolve_groups(almost_exact).await?;
    let resolution = ctx
        .read_batch(resolution)
        .map_err(|e| Error::new(ErrorKind::Other, e))?;

    let tbl = update_tables(tbl, resolution, "near-BND").await?;

    Ok(tbl)
}
