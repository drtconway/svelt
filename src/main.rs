use std::{
    io::{BufRead, Error, ErrorKind},
    rc::Rc,
};

use clap::{Parser, Subcommand};
use datafusion::{
    functions_aggregate::{count::count, expr_fn::{array_agg, bit_or}, string_agg::string_agg}, prelude::{col, concat, lit, make_array, DataFrame, SessionContext}
};
use noodles::vcf::{self, Header};
use svelt::{chroms::ChromSet, tables::load_vcf_core, vcf_reader::VcfReader};

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
                let df = df.with_column("vix", lit(1u32 << vix))?
                    .with_column("row_id", make_array(vec![lit(vix as u32), col("row_num")]))?;
                if let Some(df0) = acc {
                    let df = df0.union(df)?;
                    acc = Some(df);
                } else {
                    acc = Some(df)
                }
            }
            let df = acc.unwrap();
            let df = df.aggregate(
                vec![
                    col("chrom_id"),
                    col("start"),
                    col("end"),
                    col("kind"),
                    col("length"),
                    col("chrom2_id"),
                    col("end2"),
                    col("seq_hash"),
                ],
                vec![
                    array_agg(col("row_id")).alias("rows"),
                    count(col("vix")).alias("vix_count"),
                    bit_or(col("vix")).alias("vix_set"),
                ],
            )?;
            df.clone()
                //.filter(col("vix_count").eq(lit(n as i32)))?
                .sort(vec![
                    col("chrom_id").sort(true, true),
                    col("start").sort(true, true),
                    col("end").sort(true, true),
                    col("kind").sort(true, true),
                    col("length").sort(true, true),
                    col("chrom2_id").sort(true, true),
                    col("end2").sort(true, true),
                    col("seq_hash").sort(true, true),
                ])?
                .show()
                .await?;
        }
    }
    Ok(())
}
