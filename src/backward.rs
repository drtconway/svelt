use std::io::{Error, ErrorKind};

use datafusion::{common::JoinType, prelude::{abs, col, lit, DataFrame, SessionContext}};

use crate::{expressions::prefix_cols, options::Options, resolve::{resolve_groups, update_tables}};

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
                    .and(abs(col("lhs_end2") - col("rhs_end2")).lt(lit(options.position_window)))
                    .and((col("lhs_vix_set") & col("rhs_vix_set")).eq(lit(0))),
            ),
        )?
        .sort(vec![
            (col("lhs_vix_count") + col("rhs_vix_count")).sort(false, false),
            col("lhs_row_key").sort(true, false),
            col("rhs_row_key").sort(true, false),
        ])?;

    if true {
        flip.clone().show().await?;
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

    //let tbl = update_tables(tbl, resolution).await?;

    Ok(tbl)
}
