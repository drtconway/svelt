use std::{collections::HashMap, io::BufRead, rc::Rc};

use clap::{Parser, Subcommand};
use log::Log;
use noodles::vcf::{self, Header};
use svelt::{chroms::ChromSet, keys::VariantKey, vcf_reader::VcfReader};

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

            let mut sigs: HashMap<VariantKey, Vec<(usize,usize)>> = HashMap::new();
            for vix in 0..readers.len() {
                log::info!("reading {}", readers[vix].path);
                let reader: &mut VcfReader = &mut readers[vix];
                for (rn, rec) in reader.reader.records().enumerate() {
                    let rec = rec?;
                    let sig = VariantKey::new(&rec, &reader.header, reader.chroms.as_ref())?;
                    sigs.entry(sig).or_default().push((vix, rn));
                }
            }
            log::info!("scanned {} signatures", sigs.len());
        }
    }
    Ok(())
}
