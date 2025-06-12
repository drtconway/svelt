use std::io::{Error, ErrorKind};

use datafusion::{
    common::JoinType,
    prelude::{DataFrame, SessionContext, abs, coalesce, col, lit},
};

use crate::{
    expressions::prefix_cols,
    options::Options,
    resolve::{resolve_groups, update_tables},
};

pub async fn find_backwards_bnds(
    ctx: &SessionContext,
    tbl: DataFrame,
    n: u32,
    options: &Options,
) -> std::io::Result<DataFrame> {
    log::info!("resolving backward BND inexact matches");

    let candidates = tbl
        .clone()
        .filter(col("vix_count").lt(lit(n)).and(col("kind").eq(lit("BND"))))?;
    let lhs = prefix_cols(candidates.clone(), "lhs")?;
    let rhs = prefix_cols(candidates.clone(), "rhs")?;
    //lhs.clone().show().await?;

    let flip = lhs
        .join(
            rhs,
            JoinType::Inner,
            &["lhs_chrom_id", "lhs_chrom2_id"],
            &["rhs_chrom_id", "rhs_chrom2_id"],
            Some(
                lit(true)
                    .and(abs(col("lhs_end") - col("rhs_end")).lt(lit(options.end2_window)))
                    .and(abs(col("lhs_end") - col("rhs_end")).gt_eq(lit(options.position_window)))
                    .and(abs(col("lhs_end2") - col("rhs_end2")).lt(lit(options.position_window)))
                    .and(col("lhs_row_key").not_eq(col("rhs_row_key")))
                    .and((col("lhs_vix_set") & col("rhs_vix_set")).eq(lit(0))),
            ),
        )?
        .sort(vec![
            (col("lhs_vix_count") + col("rhs_vix_count")).sort(false, false),
            col("lhs_row_key").sort(true, false),
            col("rhs_row_key").sort(true, false),
        ])?;

    if false {
        flip.clone()
            .drop_columns(&[
                "lhs_length",
                "lhs_alt_seq",
                "lhs_seq_hash",
                "rhs_length",
                "rhs_alt_seq",
                "rhs_seq_hash",
            ])?
            .show()
            .await?;
    }

    let resolution = resolve_groups(flip).await?;
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

    let tbl = update_tables(tbl, resolution.clone(), "flip-BND").await?;

    let tbl = tbl
        .join(
            resolution.select(vec![col("orig_row_id"), lit(true).alias("new_flip")])?,
            JoinType::Left,
            &["row_id"],
            &["orig_row_id"],
            None,
        )?
        .with_column("flip", col("flip").or(col("new_flip")))?
        .drop_columns(&["orig_row_id", "new_flip"])?;

    let rhs = tbl.clone().filter(col("flip").eq(lit(true)))?.select(vec![
        col("row_id").alias("flip_row_id"),
        col("chrom").alias("flip_chrom"),
        col("chrom_id").alias("flip_chrom_id"),
        col("end").alias("flip_end"),
        col("chrom2").alias("flip_chrom2"),
        col("chrom2_id").alias("flip_chrom2_id"),
        col("end2").alias("flip_end2"),
    ])?;

    let tbl = tbl
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

    if false {
        tbl.clone()
            .drop_columns(&["alt_seq"])?
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
