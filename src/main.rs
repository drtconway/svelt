use clap::{Parser, Subcommand};
use svelt::{homology::index_features, merge::merge_vcfs, options::{CommonOptions, IndexingOptions, MergeOptions}};

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
        /// Write out the final merge table
        #[arg(long)]
        write_merge_table: Option<String>,

        /// INFO fields to drop (if they exist)
        #[arg(short, long, value_delimiter = ',')]
        unwanted_info: Vec<String>,

        /// Reference sequence. Required for some extended type of merging.
        #[arg(short, long)]
        reference: Option<String>,

        /// The output filename
        #[arg(short, long)]
        out: String,

        /// SV VCF files to merge
        #[arg(num_args(1..))]
        vcf: Vec<String>,

        #[command(flatten)]
        options: MergeOptions,

        #[command(flatten)]
        common: CommonOptions
    },

    /// Index a set of features for annotating homology
    #[command(arg_required_else_help = true)]
    IndexFeatures {
        /// FASTA file with feature sequences
        #[arg(short, long)]
        features: String,

        /// The output filename
        #[arg(short, long)]
        out: String,

        #[command(flatten)]
        options: IndexingOptions,

        #[command(flatten)]
        common: CommonOptions
    },
}

#[tokio::main]
async fn main() -> std::io::Result<()> {
    env_logger::init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Merge {
            out,
            vcf,
            unwanted_info,
            reference,
            write_merge_table,
            options,
            common
        } => {
            merge_vcfs(
                &out,
                &vcf,
                &unwanted_info,
                &reference,
                &write_merge_table,
                &options,
                &common
            )
            .await?;
        }
        Commands::IndexFeatures { out, features, options, common } => {
            index_features(&features, &out, &options, &common).await?;
        },
    }
    Ok(())
}
