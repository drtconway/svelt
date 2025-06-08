use std::io::{Error, ErrorKind};

use datafusion::{
    common::JoinType,
    functions_aggregate::count::count,
    prelude::{DataFrame, SessionContext, col, lit},
};

use crate::
    resolve::{resolve_groups, update_tables}
;

pub async fn find_exact(
    ctx: &SessionContext,
    tbl: DataFrame,
    n: u32,
) -> std::io::Result<DataFrame> {
    let tbl = find_exact_non_bnd(ctx, tbl, n).await?;
    find_exact_bnd(ctx, tbl, n).await
} 

pub async fn find_exact_non_bnd(
    ctx: &SessionContext,
    tbl: DataFrame,
    n: u32,
) -> std::io::Result<DataFrame> {
    log::info!("resolving exact non-BND matches");

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
            col("chrom2_id").alias("lhs_chrom2_id"),
            col("end2").alias("lhs_end2"),
            col("seq_hash").alias("lhs_seq_hash"),
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
            col("chrom2_id").alias("rhs_chrom2_id"),
            col("end2").alias("rhs_end2"),
            col("seq_hash").alias("rhs_seq_hash"),
            col("row_key").alias("rhs_row_key"),
            col("vix_count").alias("rhs_vix_count"),
            col("vix_set").alias("rhs_vix_set"),
        ])?;

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
                "lhs_length",
                "lhs_seq_hash",
            ],
            &[
                "rhs_chrom_id",
                "rhs_start",
                "rhs_end",
                "rhs_kind",
                "rhs_length",
                "rhs_seq_hash",
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

    let tbl = update_tables(tbl, resolution).await?;

    Ok(tbl)
}

pub async fn find_exact_bnd(
    ctx: &SessionContext,
    tbl: DataFrame,
    n: u32,
) -> std::io::Result<DataFrame> {
    log::info!("resolving exact BND matches");

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
            col("chrom2_id").alias("lhs_chrom2_id"),
            col("end2").alias("lhs_end2"),
            col("seq_hash").alias("lhs_seq_hash"),
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
            col("chrom2_id").alias("rhs_chrom2_id"),
            col("end2").alias("rhs_end2"),
            col("seq_hash").alias("rhs_seq_hash"),
            col("row_key").alias("rhs_row_key"),
            col("vix_count").alias("rhs_vix_count"),
            col("vix_set").alias("rhs_vix_set"),
        ])?;

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

    let tbl = update_tables(tbl, resolution).await?;

    Ok(tbl)
}
