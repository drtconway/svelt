use datafusion::{common::JoinType, prelude::{abs, col, lit, DataFrame}};

use crate::{expressions::prefix_cols, options::MergeOptions};

pub(super) fn approx_flipped_bnd_join(
    orig: DataFrame,
    n: usize,
    options: &MergeOptions,
) -> std::io::Result<DataFrame> {
    let candidates = orig.clone().filter(
        lit(true)
            .and(col("kind").eq(lit("BND")))
            .and(col("vix_count").lt(lit(n as u32))),
    )?;

    let lhs = prefix_cols(candidates.clone(), "lhs")?;
    let rhs = prefix_cols(candidates.clone(), "rhs")?;

    let exact = lhs
        .join(
            rhs,
            JoinType::Inner,
            &["lhs_chrom_id", "lhs_chrom2_id"],
            &["rhs_chrom_id", "rhs_chrom2_id"],
            Some(
                lit(true)
                    .and(col("lhs_row_id").lt(col("rhs_row_id")))
                    .and(abs(col("lhs_start") - col("rhs_start")).lt(lit(options.end2_window)))
                    .and(abs(col("lhs_end") - col("rhs_end")).lt(lit(options.end2_window)))
                    .and(abs(col("lhs_end2") - col("rhs_end2")).lt(lit(options.position_window)))
                    .and((col("lhs_vix_set") & col("rhs_vix_set")).eq(lit(0))),
            ),
        )?
        .sort(vec![
            (col("lhs_vix_count") + col("rhs_vix_count")).sort(false, false),
            col("lhs_row_key").sort(true, false),
            col("rhs_row_key").sort(true, false),
        ])?
        .select_columns(&["lhs_row_key", "lhs_vix_set", "rhs_row_key", "rhs_vix_set"])?
        .distinct()?;

    Ok(exact)
}