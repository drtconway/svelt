use datafusion::{
    common::JoinType,
    functions_aggregate::expr_fn::first_value,
    prelude::{DataFrame, abs, case, col, greatest, least, lit, round},
};

pub async fn make_reporting_table(tbl: DataFrame) -> std::io::Result<DataFrame> {
    let rhs = tbl
        .clone()
        .aggregate(
            vec![col("row_key")],
            vec![
                first_value(col("row_id"), Some(vec![col("vix").sort(true, false)]))
                    .alias("left_row_id"),
            ],
        )?
        .drop_columns(&["row_key"])?;

    let primary = tbl
        .clone()
        .join(rhs, JoinType::Inner, &["row_id"], &["left_row_id"], None)?
        .select(vec![
            col("row_key").alias("primary_row_key"),
            col("start").alias("primary_start"),
            col("end").alias("primary_end"),
            col("end2").alias("primary_end2"),
            abs(col("length")).alias("primary_length"),
        ])?;

    let tbl = tbl
        .clone()
        .join(
            primary,
            JoinType::Left,
            &["row_key"],
            &["primary_row_key"],
            None,
        )?
        .with_column("start_offset", abs(col("start") - col("primary_start")))?
        .with_column("end_offset", abs(col("end") - col("primary_end")))?
        .with_column("end2_offset", abs(col("end2") - col("primary_end2")))?
        .with_column(
            "total_offset",
            case(col("kind"))
                .when(lit("INS"), col("start_offset"))
                .when(lit("BND"), col("start_offset") + col("end2_offset"))
                .otherwise(col("start_offset") + col("end_offset"))?,
        )?
        .with_column(
            "length_ratio",
            least(vec![abs(col("length")), col("primary_length")]) * lit(1.0)
                / greatest(vec![abs(col("length")), col("primary_length")]),
        )?
        .with_column(
            "length_ratio",
            round(vec![col("length_ratio") * lit(100.0)]) / lit(100.0),
        )?;

    Ok(tbl)
}
