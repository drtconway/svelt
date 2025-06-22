use datafusion::prelude::DataFrame;

pub mod chi_squared;
pub mod cosine;

pub enum DistanceMetric {
    ChiSquared,
    Cosine,
}

pub async fn distance(query: DataFrame, subject: DataFrame, metric: DistanceMetric) -> std::io::Result<DataFrame> {
    match metric {
        DistanceMetric::ChiSquared => {
            chi_squared::chi_squared(query, subject).await
        },
        DistanceMetric::Cosine => {
            cosine::cosine(query, subject).await
        },
    }
}