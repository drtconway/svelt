use std::{
    io::{BufRead, Error, ErrorKind},
    rc::Rc,
};

use clap::{Parser, Subcommand};
use datafusion::
    prelude::{
        DataFrame, SessionContext, col,
        lit,
    }
;
use noodles::vcf::{self, Header};
use svelt::{almost::find_almost_exact, chroms::ChromSet, exact::find_exact, tables::load_vcf_core, vcf_reader::VcfReader};

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
        Commands::Merge { out, vcf } => {
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
        }
    }
    Ok(())
}
