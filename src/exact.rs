use std::io::{Error, ErrorKind};

use datafusion::{
    common::JoinType,
    functions_aggregate::count::count,
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

    if false {
        lhs.clone()
            .drop_columns(&[
                "lhs_chrom",
                "lhs_chrom2_id",
                "lhs_chrom2",
                "lhs_end2",
                "lhs_alt_seq",
            ])?
            .sort_by(vec![col("lhs_chrom_id"), col("lhs_start"), col("lhs_end")])?
            .show()
            .await?;
    }

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

    if false {
        exact
            .clone()
            .drop_columns(&[
                "lhs_chrom",
                "lhs_chrom2_id",
                "lhs_chrom2",
                "lhs_end2",
                "lhs_alt_seq",
                "rhs_chrom",
                "rhs_chrom2_id",
                "rhs_chrom2",
                "rhs_end2",
                "rhs_alt_seq",
            ])?
            .show()
            .await?;
    }

    let resolution = resolve_groups(exact).await?;
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

    if false {
        exact.clone().show().await?;
    }

    let resolution = resolve_groups(exact).await?;
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

    let tbl = update_tables(tbl, resolution, "exact-BND").await?;

    Ok(tbl)
}
