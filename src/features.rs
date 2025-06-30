use std::{
    io::{BufReader, Error, ErrorKind},
    sync::Arc,
};

use autocompress::autodetect_open;
use datafusion::{
    arrow::{
        array::{GenericStringBuilder, PrimitiveBuilder, RecordBatch},
        datatypes::{DataType, Field, Schema, UInt32Type},
    },
    common::JoinType,
    config::{ParquetColumnOptions, TableParquetOptions},
    dataframe::DataFrameWriteOptions,
    functions_aggregate::count::count,
    prelude::{DataFrame, ParquetReadOptions, SessionContext, col, lit, regexp_replace},
};
use noodles::fasta;

use crate::{
    distance::{DistanceMetric, distance, needleman_wunsch::align},
    errors::{SveltError, as_io_error, wrap_file_error},
    kmers::kmerize::kmers_fwd,
    options::IndexingOptions,
};

pub struct FeatureIndex {
    pub(crate) k: usize,
    pub(crate) kmers: DataFrame,
    pub(crate) names: DataFrame,
    pub(crate) sequences: DataFrame,
}

impl FeatureIndex {
    pub async fn build(
        source: &str,
        options: &IndexingOptions,
        ctx: &SessionContext,
    ) -> std::io::Result<FeatureIndex> {
        let reader = autodetect_open(source).map_err(|e| wrap_file_error(e, source))?;
        let reader = BufReader::new(reader);
        let mut reader = fasta::io::reader::Builder::default().build_from_reader(reader)?;

        let mut nix_builder = PrimitiveBuilder::<UInt32Type>::new();
        let mut name_builder = GenericStringBuilder::<i32>::new();
        let mut sequence_builder = GenericStringBuilder::<i32>::new();

        let mut sequence_number: u32 = 0;

        log::info!("reading sequences from '{}'", source);

        for rec in reader.records() {
            let rec = rec?;
            let nix = sequence_number;
            let name = rec.definition().to_string().split_off(1);
            let sequence = String::from_utf8(rec.sequence().as_ref().to_vec()).unwrap();

            nix_builder.append_value(nix);
            name_builder.append_value(name);
            sequence_builder.append_value(sequence);

            sequence_number += 1;
        }

        log::info!("constructing sequences table.");

        let nix_array = nix_builder.finish();
        let name_array = name_builder.finish();
        let sequence_array = sequence_builder.finish();

        let name_schema = Arc::new(Schema::new(vec![
            Field::new("nix", DataType::UInt32, false),
            Field::new("name", DataType::Utf8, true),
            Field::new("sequence", DataType::Utf8, true),
        ]));

        let recs = RecordBatch::try_new(
            name_schema,
            vec![
                Arc::new(nix_array),
                Arc::new(name_array),
                Arc::new(sequence_array),
            ],
        )
        .map_err(|e| Error::new(ErrorKind::Other, e))?;

        let sequences = ctx.read_batch(recs)?;

        log::info!("constructing name/class index");

        let names = sequences.clone().select(vec![
            col("nix"),
            regexp_replace(col("name"), lit(&options.pattern), lit(&options.name), None)
                .alias("name"),
            regexp_replace(
                col("name"),
                lit(&options.pattern),
                lit(&options.class),
                None,
            )
            .alias("class"),
        ])?;

        log::info!("constructing k-mer index");

        let kmerize = kmers_fwd();
        let kmers = sequences
            .clone()
            .select(vec![
                col("nix"),
                kmerize
                    .call(vec![lit(options.k as u32), col("sequence")])
                    .alias("kmer"),
            ])?
            .unnest_columns(&["kmer"])?
            .aggregate(
                vec![col("nix"), col("kmer")],
                vec![count(lit(1)).alias("count")],
            )?;

        log::info!("index construction complete");

        Ok(FeatureIndex {
            k: options.k,
            kmers,
            names,
            sequences,
        })
    }

    pub async fn save(&self, out: &str) -> std::io::Result<()> {
        log::info!("saving sequences table.");
        let opts = DataFrameWriteOptions::default();
        self.sequences
            .clone()
            .sort_by(vec![col("nix")])?
            .write_parquet(&format!("{}-sidx.parquet", out), opts, None)
            .await?;

        log::info!("saving names table.");
        let opts = DataFrameWriteOptions::default();
        self.names
            .clone()
            .sort_by(vec![col("nix")])?
            .write_parquet(&format!("{}-nidx.parquet", out), opts, None)
            .await?;

        log::info!("saving kmers table.");
        let opts = DataFrameWriteOptions::default();
        let mut kidx_opts = TableParquetOptions::default();
        let mut kmer_opts = ParquetColumnOptions::default();
        kmer_opts.compression = Some(String::from("zstd(19)"));
        kidx_opts
            .column_specific_options
            .insert(String::from("kmer"), kmer_opts);
        kidx_opts
            .key_value_metadata
            .insert(String::from("k"), Some(self.k.to_string()));
        self.kmers
            .clone()
            .sort_by(vec![col("kmer"), col("nix")])?
            .write_parquet(&format!("{}-kidx.parquet", out), opts, Some(kidx_opts))
            .await?;

        Ok(())
    }

    pub async fn load(features: &str, ctx: &SessionContext) -> std::io::Result<FeatureIndex> {
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
            .read_parquet(&format!("{}-nidx.parquet", features), opts.clone())
            .await?;
        let names = names
            .with_column_renamed("name", "idx_name")?
            .with_column_renamed("nix", "idx_nix")?;

        let sequences = ctx
            .read_parquet(&format!("{}-sidx.parquet", features), opts)
            .await?;

        log::info!("loading index done.");

        Ok(FeatureIndex {
            k,
            kmers,
            names,
            sequences,
        })
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

    pub async fn score_alignments(&self, queries: DataFrame) -> std::io::Result<DataFrame> {
        let align = align();

        let rhs = self.sequences.clone().select(vec![
            col("nix").alias("rhs_nix"),
            col("sequence").alias("rhs_sequence"),
        ])?;

        let tbl = queries
            .join(rhs, JoinType::Inner, &["nix"], &["rhs_nix"], None)?
            .with_column(
                "score",
                align.call(vec![col("sequence"), col("rhs_sequence")]),
            )?
            .drop_columns(&["rhs_nix", "rhs_sequence"])?;

        Ok(tbl)
    }
}
