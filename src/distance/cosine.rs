use datafusion::{
    common::JoinType,
    functions_aggregate::sum::sum,
    prelude::{DataFrame, col, lit, sqrt},
};

use crate::expressions::prefix_cols;

/// Compute the cosine distance metric on a set of queries and subjects.
///
/// Each dataframe is expected to have columns: `name`, `kmer`, `count`.
/// The output is a dataframe with columns `query_name` `subject_name` `distance`.
///
pub async fn cosine(query: DataFrame, subject: DataFrame) -> std::io::Result<DataFrame> {
    log::info!("evaluating cosine measure");

    let query = prefix_cols(query, "query")?;
    let subject = prefix_cols(subject, "subject")?;

    let query_mag = query
        .clone()
        .aggregate(
            vec![col("query_name")],
            vec![sum(col("query_count") * col("query_count")).alias("query_mag")],
        )?
        .with_column("query_mag", sqrt(col("query_mag")))?
        .with_column_renamed("query_name", "query_mag_name")?;

    let subject_mag = subject
        .clone()
        .aggregate(
            vec![col("subject_name")],
            vec![sum(col("subject_count") * col("subject_count")).alias("subject_mag")],
        )?
        .with_column("subject_mag", sqrt(col("subject_mag")))?
        .with_column_renamed("subject_name", "subject_mag_name")?;

    let tbl = query
        .clone()
        .join(
            subject.clone(),
            JoinType::Inner,
            &["query_kmer"],
            &["subject_kmer"],
            None,
        )?
        .aggregate(
            vec![col("query_name"), col("subject_name")],
            vec![sum(col("query_count") * col("subject_count")).alias("distance_numerator")],
        )?
        .join(
            query_mag,
            JoinType::Inner,
            &["query_name"],
            &["query_mag_name"],
            None,
        )?
        .join(
            subject_mag,
            JoinType::Inner,
            &["subject_name"],
            &["subject_mag_name"],
            None,
        )?
        .with_column(
            "distance",
            lit(1.0)
                - col("distance_numerator") * lit(1.0) / (col("query_mag") * col("subject_mag")),
        )?;

    if true {
        tbl.clone()
            .sort_by(vec![col("distance")])?
            .show()
            .await?;
    }

    let tbl = tbl.drop_columns(&[
        "query_mag_name",
        "query_mag",
        "subject_mag_name",
        "subject_mag",
    ])?;

    if false {
        tbl.clone()
            .sort_by(vec![col("query_name"), col("subject_name")])?
            .show()
            .await?;
    }

    Ok(tbl)
}
