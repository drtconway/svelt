use datafusion::{
    config::CsvOptions,
    dataframe::DataFrameWriteOptions,
    prelude::{DataFrame, abs, case, col, greatest, least, lit, round},
};

pub async fn produce_reporting_table(tbl: DataFrame, out: &str) -> std::io::Result<()> {
    let report = tbl
        .clone()
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
            least(vec![abs(col("length")), abs(col("primary_length"))]) * lit(1.0)
                / greatest(vec![abs(col("length")), abs(col("primary_length"))]),
        )?
        .with_column(
            "length_ratio",
            round(vec![col("length_ratio") * lit(100.0)]) / lit(100.0),
        )?;
    let opts = DataFrameWriteOptions::default();
    let csv_opts = CsvOptions::default().with_delimiter(b'\t');
    let csv_opts = Some(csv_opts);
    report
        .clone()
        .sort_by(vec![
            col("chrom_id"),
            col("start"),
            col("end"),
            col("row_key"),
            col("row_id"),
        ])?
        .select_columns(&[
            "variant_id",
            "chrom",
            "start",
            "end",
            "kind",
            "length",
            "start_offset",
            "end_offset",
            "end2_offset",
            "total_offset",
            "length_ratio",
            "chrom2",
            "end2",
            "seq_hash",
            "vix",
            "row_id",
            "row_key",
            "flip",
            "vix_set",
            "vix_count",
            "criteria",
            "alt_seq",
        ])?
        .write_csv(out, opts, csv_opts)
        .await?;

    Ok(())
}
