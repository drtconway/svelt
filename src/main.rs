use std::rc::Rc;

use clap::{Parser, Subcommand};
use svelt::{
    features::FeatureIndex,
    homology::{cluster_sequences, find_similar},
    merge::merge_vcfs,
    options::{CommonOptions, IndexingOptions, MergeOptions, QueryOptions, make_session_context},
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
        /// The output filename
        #[arg(short, long)]
        out: String,

        /// SV VCF files to merge
        #[arg(num_args(1..))]
        vcf: Vec<String>,

        #[command(flatten)]
        options: MergeOptions,

        #[command(flatten)]
        common: CommonOptions,
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
        common: CommonOptions,
    },

    /// Identify whether a given sequence is similar to a previously indexed one.
    #[command(arg_required_else_help = true)]
    FindSimilar {
        /// Base name for previously sequenced features
        #[arg(short, long)]
        features: String,

        /// k-mer length
        #[arg(short, long, required = false, default_value = "11")]
        k: usize,

        #[command(flatten)]
        query: QueryOptions,

        #[command(flatten)]
        common: CommonOptions,
    },

    /// Group similar sequences, and select a representative from each group
    #[command(arg_required_else_help = true)]
    ClusterSequences {
        /// Name of FASTA file with sequences to cluster
        #[arg(short, long)]
        sequences: String,

        /// Name of FASTA file to write the output to
        #[arg(short, long)]
        out: String,

        /// k-mer length
        #[arg(short, long, required = false, default_value = "11")]
        k: usize,

        /// dot-product cutoff for clustering
        #[arg(short, long, required = false, default_value = "0.01")]
        cutoff: f64,

        #[command(flatten)]
        common: CommonOptions,
    },
}

#[tokio::main]
async fn main() -> std::io::Result<()> {
    env_logger::builder()
        .filter_level(log::LevelFilter::Info)
        .init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Merge {
            out,
            vcf,
            options,
            common,
        } => {
            let options = Rc::new(options);
            merge_vcfs(&out, &vcf, options, &common).await?;
        }
        Commands::IndexFeatures {
            out,
            features,
            options,
            common,
        } => {
            let ctx = make_session_context(&common);

            let idx = FeatureIndex::build(&features, &options, &ctx).await?;
            idx.save(&out).await?;
        }
        Commands::FindSimilar {
            features,
            query,
            k,
            common,
        } => {
            find_similar(&features, &query, k, &common).await?;
        }
        Commands::ClusterSequences {
            sequences,
            out,
            k,
            cutoff,
            common,
        } => {
            cluster_sequences(&sequences, &out, k, cutoff, &common).await?;
        }
    }
    Ok(())
}
