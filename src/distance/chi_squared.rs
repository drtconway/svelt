use datafusion::{
    common::JoinType,
    functions_aggregate::{count::count, sum::sum},
    prelude::{DataFrame, coalesce, col, lit},
};

use crate::expressions::prefix_cols;

/// Compute the chi-squared distance metric on a set of queries and subjects.
///
/// Each dataframe is expected to have columns: `name`, `kmer`, `count`.
/// The output is a dataframe with columns `query_name` `subject_name` `distance`.
///
pub async fn chi_squared(query: DataFrame, subject: DataFrame) -> std::io::Result<DataFrame> {
    let query = prefix_cols(query, "query")?;
    let subject = prefix_cols(subject, "subject")?;

    let query_count = query
        .clone()
        .aggregate(
            vec![col("query_name")],
            vec![count(col("query_count")).alias("query_kmer_count")],
        )?
        .with_column_renamed("query_name", "query_count_name")?;

    let subject_count = subject
        .clone()
        .aggregate(
            vec![col("subject_name")],
            vec![count(col("subject_count")).alias("subject_kmer_count")],
        )?
        .with_column_renamed("subject_name", "subject_count_name")?;

    let tbl = query.clone().join(
        subject.clone(),
        JoinType::Inner,
        &["query_kmer"],
        &["subject_kmer"],
        None,
    )?;

    let query_missing = query
        .clone()
        .join(
            tbl.clone()
                .with_column_renamed("query_name", "rhs_query_name")?
                .with_column_renamed("query_kmer", "rhs_query_kmer")?,
            JoinType::LeftAnti,
            &["query_name", "query_kmer"],
            &["rhs_query_name", "rhs_query_kmer"],
            None,
        )?
        .aggregate(
            vec![col("query_name")],
            vec![sum(col("query_count")).alias("query_missing_count")],
        )?
        .with_column_renamed("query_name", "query_missing_name")?;

    let subject_missing = subject
        .clone()
        .join(
            tbl.clone()
                .with_column_renamed("subject_name", "rhs_subject_name")?
                .with_column_renamed("subject_kmer", "rhs_subject_kmer")?,
            JoinType::LeftAnti,
            &["subject_name", "subject_kmer"],
            &["rhs_subject_name", "rhs_subject_kmer"],
            None,
        )?
        .aggregate(
            vec![col("subject_name")],
            vec![sum(col("subject_count")).alias("subject_missing_count")],
        )?
        .with_column_renamed("subject_name", "subject_missing_name")?;

    let tbl = tbl
        .with_column("count_sum", col("query_count") + col("subject_count"))?
        .with_column(
            "count_difference",
            col("query_count") - col("subject_count"),
        )?;

    let tbl = tbl
        .aggregate(
            vec![col("query_name"), col("subject_name")],
            vec![
                sum(
                    (col("count_difference") * col("count_difference")) * lit(1.0)
                        / col("count_sum"),
                )
                .alias("chi_squared"),
                count(col("count_sum")).alias("df"),
            ],
        )?
        .join(
            query_missing,
            JoinType::Left,
            &["query_name"],
            &["query_missing_name"],
            None,
        )?
        .drop_columns(&["query_missing_name"])?
        .with_column(
            "query_missing_count",
            coalesce(vec![col("query_missing_count"), lit(0)]),
        )?
        .join(
            subject_missing,
            JoinType::Left,
            &["subject_name"],
            &["subject_missing_name"],
            None,
        )?
        .drop_columns(&["subject_missing_name"])?
        .with_column(
            "subject_missing_count",
            coalesce(vec![col("subject_missing_count"), lit(0)]),
        )?
        .with_column(
            "chi_squared",
            lit(0.5)
                * (col("chi_squared") + col("query_missing_count") + col("subject_missing_count")),
        )?
        .drop_columns(&["query_missing_count", "subject_missing_count"])?;

    let tbl = tbl
        .join(
            query_count,
            JoinType::Inner,
            &["query_name"],
            &["query_count_name"],
            None,
        )?
        .join(
            subject_count,
            JoinType::Inner,
            &["subject_name"],
            &["subject_count_name"],
            None,
        )?
        .with_column(
            "df",
            lit(1.0) * (col("query_kmer_count") + col("subject_kmer_count") - col("df")),
        )?
        .drop_columns(&[
            "query_count_name",
            //"query_count",
            "subject_count_name",
            //"subject_count",
        ])?
        .with_column("distance", col("chi_squared"))?;

    if true {
        tbl.clone()
            .sort_by(vec![col("chi_squared")])?
            .show()
            .await?;
    }

    Ok(tbl)
}
