use datafusion::{
    common::JoinType,
    functions_aggregate::{count::count, expr_fn::bit_or, min_max::min},
    prelude::{DataFrame, coalesce, col, lit, nullif},
};

pub async fn find_exact(tbl: DataFrame, _n: u32) -> std::io::Result<DataFrame> {
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
