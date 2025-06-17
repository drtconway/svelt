use std::{
    io::{Error, ErrorKind},
    time::Instant,
};

use datafusion::{
    arrow::datatypes::DataType,
    common::JoinType,
    functions_aggregate::count::count,
    prelude::{DataFrame, SessionContext, abs, cast, col, greatest, least, lit},
};

use crate::{
    expressions::{ifelse, prefix_cols},
    near_merge::near_merge_table,
    options::Options,
    range_joins::overlap_join,
    resolve::{resolve_groups, update_tables},
};

pub async fn find_almost_exact(
    ctx: &SessionContext,
    tbl: DataFrame,
    n: u32,
    options: &Options,
) -> std::io::Result<DataFrame> {
    let tbl = find_almost_exact_non_bnd(ctx, tbl, n, options).await?;
    find_almost_exact_bnd(ctx, tbl, n, options).await
}

pub async fn find_almost_exact_non_bnd(
    ctx: &SessionContext,
    tbl: DataFrame,
    n: u32,
    options: &Options,
) -> std::io::Result<DataFrame> {
    log::info!("resolving non-BND inexact matches");
    if true {
        let now = Instant::now();
        let candidates = tbl
            .clone()
            .filter(col("kind").not_eq(lit("BND")))?
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
            .filter(
                lit(true)
                    .and(col("lhs_row_id").lt(col("rhs_row_id")))
                    .and(col("lhs_row_key").not_eq(col("rhs_row_key")))
                    .and((col("lhs_vix_set") & col("rhs_vix_set")).eq(lit(0)))
                    .and(abs(col("lhs_start") - col("rhs_start")).lt(lit(options.position_window)))
                    .and(abs(col("lhs_end") - col("rhs_end")).lt(lit(options.position_window)))
                    .and(
                        ((lit(1.0) * least(vec![col("lhs_abs_length"), col("rhs_abs_length")]))
                            / greatest(vec![col("lhs_abs_length"), col("rhs_abs_length")]))
                        .gt_eq(lit(options.length_ratio)),
                    ),
            )?
            .with_column(
                "length_ratio",
                (lit(1.0) * least(vec![col("lhs_abs_length"), col("rhs_abs_length")]))
                    / greatest(vec![col("lhs_abs_length"), col("rhs_abs_length")]),
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

        return Ok(tbl);
    }
    if false {
        let now = Instant::now();

        let candidates = tbl
            .clone()
            .filter(col("kind").not_eq(lit("BND")))?
            .with_column("abs_length", abs(col("length")))?
            .with_column(
                "expanded_start",
                col("start") - lit(options.position_window),
            )?
            .with_column("expanded_end", col("end") + lit(options.position_window))?;
        let lhs = prefix_cols(candidates.clone(), "lhs")?;
        let rhs = prefix_cols(candidates.clone(), "rhs")?;

        overlap_join(
            lhs,
            rhs,
            &["lhs_kind", "lhs_chrom_id"],
            &["rhs_kind", "rhs_chrom_id"],
            (col("lhs_expanded_start"), col("lhs_expanded_end")),
            (col("rhs_expanded_start"), col("rhs_expanded_end")),
            Some(
                lit(true)
                    .and(col("lhs_row_id").lt(col("rhs_row_id")))
                    .and(col("lhs_row_key").not_eq(col("rhs_row_key")))
                    .and((col("lhs_vix_set") & col("rhs_vix_set")).eq(lit(0)))
                    .and(abs(col("lhs_start") - col("rhs_start")).lt(lit(options.position_window)))
                    .and(abs(col("lhs_end") - col("rhs_end")).lt(lit(options.position_window)))
                    .and(
                        ((lit(1.0) * least(vec![col("lhs_abs_length"), col("rhs_abs_length")]))
                            / greatest(vec![col("lhs_abs_length"), col("rhs_abs_length")]))
                        .gt_eq(lit(options.length_ratio)),
                    ),
            ),
        )
        .await?
        .drop_columns(&[
            "lhs_alt_seq",
            "lhs_chrom2_id",
            "lhs_chrom2",
            "lhs_end2",
            "lhs_flip",
            "lhs_seq_hash",
            "rhs_alt_seq",
            "rhs_chrom2_id",
            "rhs_chrom2",
            "rhs_end2",
            "rhs_flip",
            "rhs_seq_hash",
        ])?
        .distinct()?
        .sort_by(vec![col("lhs_chrom_id"), col("lhs_start")])?
        .show()
        .await?;
        log::info!("overlap_join took {}s", now.elapsed().as_secs_f32());
    }

    let candidates = tbl
        .clone()
        .filter(col("vix_count").lt(lit(n)))?
        .with_column("abs_length", abs(col("length")))?;
    let lhs = prefix_cols(candidates.clone(), "lhs")?;
    let rhs = prefix_cols(candidates.clone(), "rhs")?;
    //lhs.clone().show().await?;

    let almost_exact = lhs
        .join(
            rhs,
            JoinType::Inner,
            &["lhs_chrom_id", "lhs_kind"],
            &["rhs_chrom_id", "rhs_kind"],
            Some(
                lit(true)
                    .and(col("lhs_row_id").lt(col("rhs_row_id")))
                    .and(col("lhs_row_key").not_eq(col("rhs_row_key")))
                    .and(abs(col("lhs_start") - col("rhs_start")).lt(lit(options.position_window)))
                    .and(abs(col("lhs_end") - col("rhs_end")).lt(lit(options.position_window)))
                    .and((col("lhs_vix_set") & col("rhs_vix_set")).eq(lit(0))),
            ),
        )?
        .with_column("min_start", least(vec![col("lhs_start"), col("rhs_start")]))?
        .with_column(
            "max_start",
            greatest(vec![col("lhs_start"), col("rhs_start")]),
        )?
        .with_column("min_end", least(vec![col("lhs_end"), col("rhs_end")]))?
        .with_column("max_end", greatest(vec![col("lhs_end"), col("rhs_end")]))?
        .with_column("start_displacement", col("max_start") - col("min_start"))?
        .with_column("end_displacement", col("max_end") - col("min_end"))?
        .with_column(
            "displacement",
            ifelse(
                col("lhs_kind").eq(lit("INS")),
                col("start_displacement"),
                col("start_displacement") + col("end_displacement"),
            )?,
        )?
        .with_column(
            "max_abs_length",
            greatest(vec![col("lhs_abs_length"), col("rhs_abs_length")]),
        )?
        .with_column(
            "min_abs_length",
            least(vec![col("lhs_abs_length"), col("rhs_abs_length")]),
        )?
        .with_column(
            "length_ratio",
            cast(col("min_abs_length"), DataType::Float64)
                / cast(col("max_abs_length"), DataType::Float64),
        )?
        .drop_columns(&[
            "min_start",
            "max_start",
            "min_end",
            "max_end",
            "max_abs_length",
            "min_abs_length",
            "start_displacement",
            "end_displacement",
        ])?
        .filter(col("length_ratio").gt(lit(options.length_ratio)))?
        .sort(vec![
            (col("lhs_vix_count") + col("rhs_vix_count")).sort(false, false),
            col("length_ratio").sort(false, false),
            col("lhs_row_key").sort(true, false),
            col("rhs_row_key").sort(true, false),
        ])?;

    if false {
        almost_exact
            .clone()
            .drop_columns(&[
                "lhs_chrom",
                "lhs_chrom2_id",
                "lhs_chrom2",
                "lhs_end2",
                "rhs_chrom",
                "rhs_chrom2_id",
                "rhs_chrom2",
                "rhs_end2",
            ])?
            .show()
            .await?;
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

    let tbl = update_tables(tbl, resolution, "near").await?;

    Ok(tbl)
}

pub async fn find_almost_exact_bnd(
    ctx: &SessionContext,
    tbl: DataFrame,
    n: u32,
    options: &Options,
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

    let tbl = update_tables(tbl, resolution, "near-BND").await?;

    Ok(tbl)
}
