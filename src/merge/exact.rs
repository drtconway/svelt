use datafusion::{
    common::JoinType,
    prelude::{DataFrame, abs, col, greatest, least, lit},
};

use crate::{expressions::prefix_cols, options::MergeOptions};

pub(super) fn full_exact_indel_join(orig: DataFrame, n: usize) -> std::io::Result<DataFrame> {
    let candidates = orig.clone().filter(
        lit(true)
            .and(col("kind").not_eq(lit("BND")))
            .and(col("vix_count").lt(lit(n as u32))),
    )?;

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
                    .and(
                        col("lhs_kind")
                            .not_eq(lit("INS"))
                            .or(col("lhs_seq_hash").eq(col("rhs_seq_hash"))),
                    ),
            ),
        )?
        .select_columns(&["lhs_row_key", "lhs_vix_set", "rhs_row_key", "rhs_vix_set"])?
        .distinct()?
        .sort(vec![
            (col("lhs_vix_count") + col("rhs_vix_count")).sort(false, false),
            col("lhs_row_key").sort(true, false),
            col("rhs_row_key").sort(true, false),
        ])?;

    Ok(exact)
}

pub(super) fn full_exact_locus_ins_join(
    orig: DataFrame,
    n: usize,
    options: &MergeOptions,
) -> std::io::Result<DataFrame> {
    let candidates = orig.clone().filter(
        lit(true)
            .and(col("kind").eq(lit("INS")))
            .and(col("vix_count").lt(lit(n as u32))),
    )?;

    let lhs = prefix_cols(candidates.clone(), "lhs")?;
    let rhs = prefix_cols(candidates.clone(), "rhs")?;

    let exact = lhs
        .join(
            rhs,
            JoinType::Inner,
            &["lhs_chrom_id", "lhs_start", "lhs_end", "lhs_kind"],
            &["rhs_chrom_id", "rhs_start", "rhs_end", "rhs_kind"],
            Some(
                lit(true)
                    .and(col("lhs_row_id").lt(col("rhs_row_id")))
                    .and((col("lhs_vix_set") & col("rhs_vix_set")).eq(lit(0)))
                    .and(
                        (least(vec![abs(col("lhs_length")), abs(col("rhs_length"))]) * lit(1.0)
                            / greatest(vec![abs(col("lhs_length")), abs(col("rhs_length"))]))
                        .gt_eq(lit(options.length_ratio)),
                    ),
            ),
        )?
        .select_columns(&["lhs_row_key", "lhs_vix_set", "rhs_row_key", "rhs_vix_set"])?
        .distinct()?
        .sort(vec![
            (col("lhs_vix_count") + col("rhs_vix_count")).sort(false, false),
            col("lhs_row_key").sort(true, false),
            col("rhs_row_key").sort(true, false),
        ])?;

    Ok(exact)
}

pub(super) fn full_exact_bnd(orig: DataFrame, n: usize) -> std::io::Result<DataFrame> {
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
        .select_columns(&["lhs_row_key", "lhs_vix_set", "rhs_row_key", "rhs_vix_set"])?
        .distinct()?
        .sort(vec![
            (col("lhs_vix_count") + col("rhs_vix_count")).sort(false, false),
            col("lhs_row_key").sort(true, false),
            col("rhs_row_key").sort(true, false),
        ])?;

    Ok(exact)
}
