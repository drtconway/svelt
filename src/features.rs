use datafusion::{
    common::JoinType,
    prelude::{DataFrame, ParquetReadOptions, SessionContext, lit},
};

use crate::{
    distance::{DistanceMetric, distance},
    errors::{SveltError, as_io_error},
};

pub(crate) struct FeatureIndex {
    pub(crate) k: usize,
    pub(crate) kmers: DataFrame,
    pub(crate) names: DataFrame,
}

impl FeatureIndex {
    pub async fn new(features: &str, ctx: &SessionContext) -> std::io::Result<FeatureIndex> {
        log::info!("loading index '{}'", features);
        let opts = ParquetReadOptions::default().skip_metadata(false);
        let kmers = ctx
            .read_parquet(&format!("{}-kidx.parquet", features), opts.clone())
            .await?;
        let meta = kmers.schema().metadata();
        let k: usize = if let Some(k_str) = meta.get("k") {
            if let Ok(k) = k_str.parse() {
                k
            } else {
                return Err(as_io_error(SveltError::MissingK(String::from(features))));
            }
        } else {
            return Err(as_io_error(SveltError::MissingK(String::from(features))));
        };

        let names = ctx
            .read_parquet(&format!("{}-nidx.parquet", features), opts)
            .await?;
        let names = names
            .with_column_renamed("name", "idx_name")?
            .with_column_renamed("nix", "idx_nix")?;

        log::info!("loading index done.");

        Ok(FeatureIndex { k, kmers, names })
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub async fn rank(&self, query: DataFrame) -> std::io::Result<DataFrame> {
        let query = query.with_column("name", lit("query"))?;

        let subject = self.kmers.clone().with_column_renamed("nix", "name")?;

        let res = distance(query, subject, DistanceMetric::Cosine).await?;
        let res = res
            .join(
                self.names.clone(),
                JoinType::Left,
                &["subject_name"],
                &["idx_nix"],
                None,
            )?
            .drop_columns(&["subject_name", "idx_nix"])?
            .with_column_renamed("idx_name", "subject_name")?;

        Ok(res)
    }

    pub async fn rank_many(&self, queries: DataFrame) -> std::io::Result<DataFrame> {
        let subject = self.kmers.clone().with_column_renamed("nix", "name")?;

        let res = distance(queries, subject, DistanceMetric::Cosine).await?;
        let res = res
            .join(
                self.names.clone(),
                JoinType::Left,
                &["subject_name"],
                &["idx_nix"],
                None,
            )?
            .drop_columns(&["subject_name", "idx_nix"])?
            .with_column_renamed("idx_name", "subject_name")?;
        Ok(res)
    }
}
