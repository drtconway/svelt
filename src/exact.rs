use std::io::{Error, ErrorKind};

use datafusion::{
    common::JoinType,
    functions_aggregate::{count::count, expr_fn::bit_or, min_max::min},
    prelude::{DataFrame, SessionContext, coalesce, col, lit, nullif},
};

use crate::
    resolve::{resolve_groups, update_tables}
;

pub async fn find_exact(
    ctx: &SessionContext,
    tbl: DataFrame,
    n: u32,
) -> std::io::Result<DataFrame> {
    find_exact_non_bnd(ctx, tbl, n).await
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
            datafusion::common::JoinType::Inner,
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
            datafusion::common::JoinType::Inner,
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

pub async fn find_exact0(tbl: DataFrame, _n: u32) -> std::io::Result<DataFrame> {
    let exact = tbl
        .clone()
        .aggregate(
            vec![
                col("chrom_id"),
                col("start"),
                col("end"),
                col("kind"),
                col("length"),
                col("chrom2_id"),
                col("end2"),
                col("seq_hash"),
            ],
            vec![
                min(col("row_id")).alias("row_key"),
                count(col("vix")).alias("vix_count"),
                bit_or(col("vix")).alias("vix_set"),
            ],
        )?
        .select(vec![
            col("chrom_id").alias("exact_chrom_id"),
            col("start").alias("exact_start"),
            col("end").alias("exact_end"),
            col("kind").alias("exact_kind"),
            col("length").alias("exact_length"),
            col("chrom2_id").alias("exact_chrom2_id"),
            col("end2").alias("exact_end2"),
            col("seq_hash").alias("exact_seq_hash"),
            col("row_key").alias("exact_row_key"),
            col("vix_count").alias("exact_vix_count"),
            col("vix_set").alias("exact_vix_set"),
        ])?;
    //.filter(col("exact_vix_count").eq(lit(n)))?;

    let tbl = tbl
        .with_column("row_key", nullif(lit(1), lit(1)))?
        .with_column("vix_count", nullif(lit(1), lit(1)))?
        .with_column("vix_set", nullif(lit(1), lit(1)))?;

    let exact_ins = exact.clone().filter(col("exact_kind").eq(lit("INS")))?;
    let tbl = find_exact_inner(
        tbl,
        exact_ins,
        &["chrom_id", "start", "end", "kind", "length", "seq_hash"],
        &[
            "exact_chrom_id",
            "exact_start",
            "exact_end",
            "exact_kind",
            "exact_length",
            "exact_seq_hash",
        ],
    )
    .await?;

    let exact_del = exact.clone().filter(col("exact_kind").eq(lit("DEL")))?;
    let tbl = find_exact_inner(
        tbl,
        exact_del,
        &["chrom_id", "start", "end", "kind", "length"],
        &[
            "exact_chrom_id",
            "exact_start",
            "exact_end",
            "exact_kind",
            "exact_length",
        ],
    )
    .await?;

    let exact_dup = exact.clone().filter(col("exact_kind").eq(lit("DUP")))?;
    let tbl = find_exact_inner(
        tbl,
        exact_dup,
        &["chrom_id", "start", "end", "kind", "length"],
        &[
            "exact_chrom_id",
            "exact_start",
            "exact_end",
            "exact_kind",
            "exact_length",
        ],
    )
    .await?;

    let exact_inv = exact.clone().filter(col("exact_kind").eq(lit("INV")))?;
    let tbl = find_exact_inner(
        tbl,
        exact_inv,
        &["chrom_id", "start", "end", "kind", "length"],
        &[
            "exact_chrom_id",
            "exact_start",
            "exact_end",
            "exact_kind",
            "exact_length",
        ],
    )
    .await?;

    let exact_bnd = exact.clone().filter(col("exact_kind").eq(lit("BND")))?;
    let tbl = find_exact_inner(
        tbl,
        exact_bnd,
        &["chrom_id", "start", "end", "kind", "chrom2_id", "end2"],
        &[
            "exact_chrom_id",
            "exact_start",
            "exact_end",
            "exact_kind",
            "exact_chrom2_id",
            "exact_end2",
        ],
    )
    .await?;

    Ok(tbl)
}

async fn find_exact_inner(
    lhs: DataFrame,
    rhs: DataFrame,
    lhs_cols: &[&str],
    rhs_cols: &[&str],
) -> std::io::Result<DataFrame> {
    let tbl = lhs
        .join(rhs, JoinType::Left, lhs_cols, rhs_cols, None)
        .unwrap()
        .drop_columns(&[
            "exact_chrom_id",
            "exact_start",
            "exact_end",
            "exact_kind",
            "exact_length",
            "exact_chrom2_id",
            "exact_end2",
            "exact_seq_hash",
        ])
        .unwrap()
        .with_column(
            "row_key",
            coalesce(vec![col("row_key"), col("exact_row_key")]),
        )
        .unwrap()
        .with_column(
            "vix_count",
            coalesce(vec![col("vix_count"), col("exact_vix_count")]),
        )
        .unwrap()
        .with_column(
            "vix_set",
            coalesce(vec![col("vix_set"), col("exact_vix_set")]),
        )
        .unwrap()
        .drop_columns(&["exact_row_key", "exact_vix_count", "exact_vix_set"])?;

    Ok(tbl)
}
