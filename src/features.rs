use std::{
    collections::HashMap,
    fmt::format,
    fs::File,
    io::{BufReader, Error, ErrorKind},
    iter::zip,
    sync::Arc,
};

use autocompress::{CompressionLevel, autodetect_open};
use datafusion::{
    arrow::{
        array::{
            Float64Array, GenericStringBuilder, PrimitiveBuilder, RecordBatch, StringArray, UInt32Array, UInt64Array
        },
        datatypes::{DataType, Field, Float64Type, Schema, UInt32Type, UInt64Type},
    },
    common::JoinType,
    parquet::arrow::{arrow_reader::ParquetRecordBatchReaderBuilder, ArrowWriter},
    prelude::{col, lit, DataFrame, ParquetReadOptions, SessionContext},
};
use noodles::fasta;

use crate::{
    distance::{DistanceMetric, distance, needleman_wunsch::align},
    errors::{SveltError, as_io_error, wrap_file_error},
    kmers::{KmerIterator, kmerize::kmers_fwd},
    options::IndexingOptions,
};

pub struct FeatureIndex {
    pub(crate) k: usize,
    pub(crate) kmers: HashMap<u64, Vec<(u32, u32)>>,
    pub(crate) names: Vec<String>,
    pub(crate) sequences: Vec<String>,
    pub(crate) mags: Vec<f64>,
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

        let k = options.k;

        let mut sequence_number: u32 = 0;

        log::info!("reading sequences from '{}'", source);

        let mut names = Vec::new();
        let mut sequences = Vec::new();
        let mut mags = Vec::new();
        let mut kmers: HashMap<u64, Vec<(u32, u32)>> = HashMap::new();

        for rec in reader.records() {
            let rec = rec?;
            let nix = sequence_number;
            let name = rec.definition().to_string().split_off(1);
            let sequence = String::from_utf8(rec.sequence().as_ref().to_vec()).unwrap();

            let mut tmp: HashMap<u64, u32> = HashMap::new();
            for (x, _) in KmerIterator::new(k, sequence.as_bytes().iter()) {
                *tmp.entry(x.0).or_default() += 1;
            }

            let mut mag = 0;
            for (x, count) in tmp.into_iter() {
                mag += count * count;
                kmers.entry(x).or_default().push((nix, count));
            }

            names.push(name);
            sequences.push(sequence);
            mags.push((mag as f64).sqrt());

            sequence_number += 1;
        }

        log::info!("index construction complete");

        Ok(FeatureIndex {
            k: options.k,
            kmers,
            names,
            sequences,
            mags,
        })
    }

    pub async fn save(&self, out: &str) -> std::io::Result<()> {
        log::info!("saving sequences table.");

        let mut name_builder = GenericStringBuilder::<i32>::new();
        for name in self.names.iter() {
            name_builder.append_value(name);
        }
        let name_array = name_builder.finish();

        let mut sequence_builder = GenericStringBuilder::<i32>::new();
        for sequence in self.sequences.iter() {
            sequence_builder.append_value(sequence);
        }
        let sequence_array = sequence_builder.finish();

        let mut mags_builder = PrimitiveBuilder::<Float64Type>::new();
        for mag in self.mags.iter() {
            mags_builder.append_value(*mag);
        }
        let mags_array = mags_builder.finish();

        let name_schema = Arc::new(Schema::new(vec![
            Field::new("name", DataType::Utf8, false),
            Field::new("sequence", DataType::Utf8, false),
            Field::new("mags", DataType::Float64, false),
        ]));

        let recs = RecordBatch::try_new(
            name_schema,
            vec![
                Arc::new(name_array),
                Arc::new(sequence_array),
                Arc::new(mags_array),
            ],
        )
        .map_err(|e| Error::new(ErrorKind::Other, e))?;

        save_record_batch(&recs, &format!("{}-names.parquet", out))?;

        log::info!("saving kmers table.");

        let mut kmers_builder = PrimitiveBuilder::<UInt64Type>::new();
        let mut nixs_builder = PrimitiveBuilder::<UInt32Type>::new();
        let mut counts_builder = PrimitiveBuilder::<UInt32Type>::new();

        for (x, hits) in self.kmers.iter() {
            for (nix, count) in hits.iter() {
                kmers_builder.append_value(*x);
                nixs_builder.append_value(*nix);
                counts_builder.append_value(*count);
            }
        }

        let kmers_array = kmers_builder.finish();
        let nixs_array = nixs_builder.finish();
        let counts_array = counts_builder.finish();

        let kmers_meta: HashMap<String, String> = vec![(String::from("k"), format!("{}", self.k))]
            .into_iter()
            .collect();
        let kmers_schema = Arc::new(Schema::new_with_metadata(
            vec![
                Field::new("kmer", DataType::UInt64, false),
                Field::new("nix", DataType::UInt32, false),
                Field::new("count", DataType::UInt32, false),
            ],
            kmers_meta,
        ));

        let recs = RecordBatch::try_new(
            kmers_schema,
            vec![
                Arc::new(kmers_array),
                Arc::new(nixs_array),
                Arc::new(counts_array),
            ],
        )
        .map_err(|e| Error::new(ErrorKind::Other, e))?;

        save_record_batch(&recs, &format!("{}-kmers.parquet", out))?;

        Ok(())
    }

    pub async fn load(features: &str, ctx: &SessionContext) -> std::io::Result<FeatureIndex> {
        log::info!("loading index '{}'", features);

        let mut names: Vec<String> = Vec::new();
        let mut sequences: Vec<String> = Vec::new();
        let mut mags: Vec<f64> = Vec::new();
        let mut kmers: HashMap<u64, Vec<(u32, u32)>> = HashMap::new();
        let mut k: usize = 0;

        let reader = File::open(&format!("{}-names.parquet", features))?;
        let mut reader = ParquetRecordBatchReaderBuilder::try_new(reader)?.build()?;
        for recs in reader {
            let recs = recs.map_err(|e| wrap_file_error(e, features))?;
            let name_col = recs
                .column(0)
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap();
            let sequence_col = recs
                .column(1)
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap();
                let mags_col = recs
                .column(1)
                .as_any()
                .downcast_ref::<Float64Array>()
                .unwrap();
            for (name, (sequence, mag)) in zip(name_col, zip(sequence_col, mags_col)) {
                let name = name.unwrap();
                let sequence = sequence.unwrap();
                let mag = mag.unwrap();
                names.push(name.to_string());
                sequences.push(sequence.to_string());
                mags.push(mag);
            }
        }

        let reader = File::open(&format!("{}-kmers.parquet", features))?;
        let mut reader = ParquetRecordBatchReaderBuilder::try_new(reader)?.build()?;
        for recs in reader {
            let recs = recs.map_err(|e| wrap_file_error(e, features))?;
            if k == 0 {
                let k_str = recs
                    .schema()
                    .metadata()
                    .get("k")
                    .map(|s| s.to_string())
                    .unwrap();
                k = k_str.parse().unwrap();
            }
            let kmer_col = recs
                .column(0)
                .as_any()
                .downcast_ref::<UInt64Array>()
                .unwrap();
            let nix_col = recs
                .column(1)
                .as_any()
                .downcast_ref::<UInt32Array>()
                .unwrap();
            let count_col = recs
                .column(2)
                .as_any()
                .downcast_ref::<UInt32Array>()
                .unwrap();

            for (kmer, (nix, count)) in zip(kmer_col, zip(nix_col, count_col)) {
                let kmer = kmer.unwrap();
                let nix = nix.unwrap();
                let count = count.unwrap();
                kmers.entry(kmer).or_default().push((nix, count));
            }
        }

        log::info!("loading index done.");

        Ok(FeatureIndex {
            k,
            kmers,
            names,
            sequences,
            mags
        })
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub async fn rank(&self, query: &str) -> std::io::Result<()> {
        let mut fwd: HashMap<u64, u32> = HashMap::new();
        let mut rev: HashMap<u64, u32> = HashMap::new();

        for (x, y) in KmerIterator::new(self.k, query.as_bytes().iter()) {
            *fwd.entry(x.0).or_default() += 1;
            *rev.entry(y.0).or_default() += 1;
        }

        let mut q_mag = 0;
        let mut d_fwd: HashMap<u32, u32> = HashMap::new();
        for (x, q) in fwd.iter() {
            q_mag += *q * *q;
            if let Some(items) = self.kmers.get(x) {
                for (nix, s) in items.iter() {
                    *d_fwd.entry(*nix).or_default() += *q * *s;
                }
            }
        }

        Ok(())
    }

    pub async fn rank_many(&self, queries: DataFrame) -> std::io::Result<DataFrame> {
        todo!()
    }
}

fn save_record_batch(recs: &RecordBatch, path: &str) -> std::io::Result<()> {
    let writer = autocompress::autodetect_create_or_stdout(Some(path), CompressionLevel::best())?;
    let mut writer = ArrowWriter::try_new(writer, recs.schema(), None)?;
    writer.write(recs)?;
    writer.close();
    Ok(())
}
