use std::io::{Error, ErrorKind};

use datafusion::{
    common::JoinType,
    prelude::{DataFrame, SessionContext, col, lit},
};

use crate::{
    expressions::prefix_cols,
    options::MergeOptions,
    resolve::{resolve_groups, update_tables},
};

pub async fn find_exact(
    ctx: &SessionContext,
    tbl: DataFrame,
    n: u32,
    options: &MergeOptions,
) -> std::io::Result<DataFrame> {
    let tbl = find_exact_non_bnd(ctx, tbl, n, options).await?;
    find_exact_bnd(ctx, tbl, n, options).await
}

pub async fn find_exact_non_bnd(
    ctx: &SessionContext,
    tbl: DataFrame,
    n: u32,
    _options: &MergeOptions,
) -> std::io::Result<DataFrame> {
    log::info!("resolving exact non-BND matches");

    let candidates = tbl.clone().filter(col("vix_count").lt(lit(n)))?;
    let lhs = prefix_cols(candidates.clone(), "lhs")?;
    let rhs = prefix_cols(candidates.clone(), "rhs")?;

    let exact = lhs
        .join(
            rhs,
            JoinType::Inner,
            &[
                "lhs_chrom_id",
                "lhs_start",
                "lhs_end",
                "lhs_kind",
                "lhs_length",
            ],
            &[
                "rhs_chrom_id",
                "rhs_start",
                "rhs_end",
                "rhs_kind",
                "rhs_length",
            ],
            Some(
                lit(true)
                    .and(col("lhs_row_id").lt(col("rhs_row_id")))
                    .and((col("lhs_vix_set") & col("rhs_vix_set")).eq(lit(0)))
                    .and(col("lhs_kind").not_eq(lit("INS")).or(col("lhs_seq_hash").eq(col("rhs_seq_hash"))))
            ),
        )?
        .sort(vec![
            (col("lhs_vix_count") + col("rhs_vix_count")).sort(false, false),
            col("lhs_row_key").sort(true, false),
            col("rhs_row_key").sort(true, false),
        ])?;

    let resolution = resolve_groups(exact).await?;
    let resolution = ctx
        .read_batch(resolution)
        .map_err(|e| Error::new(ErrorKind::Other, e))?;

    let tbl = update_tables(tbl, resolution, "exact").await?;

    Ok(tbl)
}

pub async fn find_exact_bnd(
    ctx: &SessionContext,
    tbl: DataFrame,
    n: u32,
    _options: &MergeOptions,
) -> std::io::Result<DataFrame> {
    log::info!("resolving exact BND matches");

    let candidates = tbl.clone().filter(col("vix_count").lt(lit(n)))?;
    let lhs = prefix_cols(candidates.clone(), "lhs")?;
    let rhs = prefix_cols(candidates.clone(), "rhs")?;
    //lhs.clone().show().await?;

    let exact = lhs
        .join(
            rhs,
            JoinType::Inner,
            &[
                "lhs_chrom_id",
                "lhs_start",
                "lhs_end",
                "lhs_kind",
                "lhs_chrom2_id",
                "lhs_end2",
            ],
            &[
                "rhs_chrom_id",
                "rhs_start",
                "rhs_end",
                "rhs_kind",
                "rhs_chrom2_id",
                "rhs_end2",
            ],
            Some(
                lit(true)
                    .and(col("lhs_row_id").lt(col("rhs_row_id")))
                    .and((col("lhs_vix_set") & col("rhs_vix_set")).eq(lit(0))),
            ),
        )?
        .sort(vec![
            (col("lhs_vix_count") + col("rhs_vix_count")).sort(false, false),
            col("lhs_row_key").sort(true, false),
            col("rhs_row_key").sort(true, false),
        ])?;

    let resolution = resolve_groups(exact).await?;
    let resolution = ctx
        .read_batch(resolution)
        .map_err(|e| Error::new(ErrorKind::Other, e))?;

    let tbl = update_tables(tbl, resolution, "exact-BND").await?;

    Ok(tbl)
}
