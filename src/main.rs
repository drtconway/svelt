use std::{
    collections::HashMap,
    io::{BufRead, Error, ErrorKind},
    rc::Rc,
    u32,
};

use autocompress::autodetect_create;
use clap::{Parser, Subcommand};
use datafusion::{
    arrow::array::Int64Array,
    prelude::{DataFrame, SessionContext, col, lit},
};
use noodles::vcf::{self, Header, Record, header::SampleNames, variant::io::Write};
use svelt::{
    almost::find_almost_exact, chroms::ChromSet, construct::construct_record, exact::find_exact,
    record_seeker::RecordSeeker, row_key::RowKey, tables::load_vcf_core, vcf_reader::VcfReader,
};

/// Structuaral Variant (SV) VCF merging
#[derive(Debug, Parser)]
#[command(name = "svelt")]
#[command(about = "Merge structural variants, aligning similar ones.", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Apply annotations from one VCF to another
    #[command(arg_required_else_help = true)]
    Merge {
        /// Force ALTs to be symbolic
        #[arg(short, long)]
        force_alt_tags: bool,

        /// The output filename
        #[arg(short, long)]
        out: String,

        /// SV VCF files to merge
        #[arg(num_args(1..))]
        vcf: Vec<String>,
    },
}

fn load_chroms(path: &str) -> std::io::Result<ChromSet> {
    let reader = autocompress::autodetect_open(path)?;
    let mut reader: vcf::io::Reader<Box<dyn BufRead>> =
        vcf::io::reader::Builder::default().build_from_reader(reader)?;
    let header: Header = reader.read_header()?;

    let mut names: Vec<&str> = Vec::new();
    for contig in header.contigs() {
        let chrom = contig.0;
        names.push(chrom);
    }
    let chroms = ChromSet::from(names.as_ref());
    Ok(chroms)
}

#[tokio::main]
async fn main() -> std::io::Result<()> {
    env_logger::init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Merge {
            out,
            vcf,
            force_alt_tags,
        } => {
            let chroms = load_chroms(&vcf[0])?;
            let chroms = Rc::new(chroms);
            let mut readers = Vec::new();
            for vcf in vcf.iter() {
                let reader = VcfReader::new(vcf, chroms.clone())?;
                readers.push(reader);
            }
            let n = readers.len();

            let ctx = SessionContext::new();

            let mut acc: Option<DataFrame> = None;
            for vix in 0..readers.len() {
                log::info!("reading {}", readers[vix].path);
                let reader: &mut VcfReader = &mut readers[vix];
                let records = load_vcf_core(reader)?;
                let df = ctx
                    .read_batch(records)
                    .map_err(|e| Error::new(ErrorKind::Other, e))?;
                let df = df
                    .with_column("vix", lit(1u32 << vix))?
                    .with_column("row_id", lit(vix as u32) + (col("row_num") * lit(100)))?;

                if let Some(df0) = acc {
                    let df = df0.union(df)?;
                    acc = Some(df);
                } else {
                    acc = Some(df)
                }
            }
            let orig = acc.unwrap();

            let results = find_exact(orig, n as u32).await?;

            let results = find_almost_exact(&ctx, results, n as u32).await?;

            if false {
                results
                    .clone()
                    .sort_by(vec![
                        col("chrom_id"),
                        col("start"),
                        col("end"),
                        col("row_id"),
                    ])?
                    .show()
                    .await?;
            }

            let table = results
                .sort_by(vec![
                    col("chrom_id"),
                    col("start"),
                    col("end"),
                    col("row_key"),
                ])?
                .collect()
                .await?;

            let mut sample_names: Vec<String> = Vec::new();
            let mut vix_samples = Vec::new();

            for h in readers.iter().map(|r| &r.header) {
                for s in h.sample_names() {
                    sample_names.push(s.clone());
                }
                vix_samples.push(h.sample_names().len());
            }

            let mut header = readers[0].header.clone();
            *header.sample_names_mut() =
                SampleNames::from_iter(sample_names.iter().map(|s| s.clone()));

            let mut seekers = Vec::new();
            for path in vcf.iter() {
                let seeker = RecordSeeker::new(path, chroms.clone())?;
                seekers.push(seeker);
            }

            let writer = autodetect_create(out, autocompress::CompressionLevel::Default)?;
            let mut writer = vcf::io::Writer::new(writer);
            writer.write_header(&header)?;

            let mut current_chrom = String::new();
            let mut current_row_key = u32::MAX;
            let mut current_row: Vec<Option<u32>> = (0..n).into_iter().map(|_| None).collect();

            for recs in table.into_iter() {
                for field in recs.schema_ref().fields().iter() {
                    log::debug!("field: {:?}", field);
                }
                let row_ids = recs
                    .column_by_name("row_id")
                    .unwrap()
                    .as_any()
                    .downcast_ref::<Int64Array>()
                    .unwrap();
                let row_keys = recs
                    .column_by_name("row_key")
                    .unwrap()
                    .as_any()
                    .downcast_ref::<Int64Array>()
                    .unwrap();

                for i in 0..row_ids.len() {
                    let row_id = row_ids.value(i) as u32;
                    let row_key = row_keys.value(i) as u32;

                    if row_key != current_row_key {
                        log::debug!("flushing group: {} {:?}", current_row_key, current_row);
                        let mut is_empty = true;
                        let mut recs: Vec<Option<(Rc<Header>, Record)>> =
                            (0..n).into_iter().map(|_| None).collect();
                        let mut m = 0;
                        for vix in 0..n {
                            if let Some(rn) = current_row[vix] {
                                m += 1;
                                is_empty = false;
                                let hnr = seekers[vix].take(rn)?.unwrap();
                                recs[vix] = Some(hnr)
                            }
                        }
                        if !is_empty {
                            let rec =
                                construct_record(&header, recs, &vix_samples, force_alt_tags)?;
                            writer.write_variant_record(&header, &rec)?;
                            if rec.reference_sequence_name() != &current_chrom {
                                current_chrom = String::from(rec.reference_sequence_name());
                                log::info!("scanning {}", current_chrom);
                            }
                        }

                        current_row = (0..n).into_iter().map(|_| None).collect();
                        current_row_key = row_key;
                    }
                    let (vix, rn) = RowKey::decode(row_id);
                    current_row[vix as usize] = Some(rn);
                }
            }
            log::debug!("flushing group: {} {:?}", current_row_key, current_row);
            let mut is_empty = true;
            let mut recs: Vec<Option<(Rc<Header>, Record)>> =
                (0..n).into_iter().map(|_| None).collect();
            let mut m = 0;
            for vix in 0..n {
                if let Some(rn) = current_row[vix] {
                    m += 1;
                    is_empty = false;
                    let hnr = seekers[vix].take(rn)?.unwrap();
                    recs[vix] = Some(hnr)
                }
            }
            if !is_empty {
                let rec = construct_record(&header, recs, &vix_samples, force_alt_tags)?;
                writer.write_variant_record(&header, &rec)?;
            }
        }
    }
    Ok(())
}
