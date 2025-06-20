use clap::{ArgAction, Args};
use datafusion::prelude::{SessionConfig, SessionContext};

use crate::errors::SveltError;

/// Options controlling the merge process
#[derive(Debug, Args)]
pub struct MergeOptions {
    /// Allowed distance for merging of events
    #[arg(long, required = false, default_value = "25")]
    pub position_window: u32,

    /// Allowed end2 distance for merging BND events
    #[arg(long, required = false, default_value = "150")]
    pub end2_window: u32,

    /// Minimum length ratio (shorter/longer) for merging two events
    #[arg(long, required = false, default_value = "0.9")]
    pub length_ratio: f64,

    /// Write out the final merge table
    #[arg(long)]
    pub write_merge_table: Option<String>,

    /// INFO fields to drop (if they exist)
    #[arg(short, long, value_delimiter = ',')]
    pub unwanted_info: Vec<String>,

    /// Reference sequence. Required for some extended type of merging.
    #[arg(short, long)]
    pub reference: Option<String>,

    /// Force ALTs to be symbolic
    #[arg(long,
        action = ArgAction::Set,
        default_value_t = true,
        default_missing_value = "true",
        num_args = 0..=1,)]
    pub force_alt_tags: bool,

    /// Fill in correct reference bases
    #[arg(long,
        action = ArgAction::Set,
        default_value_t = true,
        default_missing_value = "true",
        num_args = 0..=1,)]
    pub fill_in_refs: bool,

    /// Allow breakend variants to be flipped
    #[arg(long,
        action = ArgAction::Set,
        default_value_t = true,
        default_missing_value = "true",
        num_args = 0..=1,)]
    pub allow_breakend_flipping: bool,
}

impl MergeOptions {
    /// Check merge options for mutual consistency
    pub fn check(&self) -> std::result::Result<(), SveltError> {
        if self.fill_in_refs && self.reference.is_none() {
            return Err(SveltError::OptionReferenceRequired(String::from("--fill-in-refs")));
        }
        if self.allow_breakend_flipping && self.reference.is_none() {
            return Err(SveltError::OptionReferenceRequired(String::from("--allow-breakend-flipping")));
        }
        Ok(())
    }
}

/// Options controlling feature indexing
#[derive(Debug, Args)]
pub struct IndexingOptions {
    /// k-mer length
     #[arg(short, long, required = false, default_value = "11")]
    pub k: usize,

    /// Regular expression for parsing names
     #[arg(long, required = false, default_value = "^(?<name>[[:word:]]+)[[:blank:]]+(?<class>[[:word:]]+)")]
    pub pattern: String,

    /// Reference to get the name (number for group number, or name for named group)
    #[arg(long, required = false, default_value = "name")]
    pub name: String,

    /// Reference to get the class (number for group number, or name for named group)
    #[arg(long, required = false, default_value = "class")]
    pub class: String,
}

/// Options for different kinds of query
#[derive(Debug, Args)]
#[group(required = true, multiple = false)]
pub struct QueryOptions {
    /// A literal query sequence given on the command line
    #[arg(short, long)]
    pub query: Option<String>,

    /// Read query sequences from a FASTA file
    #[arg(short, long)]
    pub query_file: Option<String>,

    /// Read query sequences from a VCF (either insersion sequences,
    /// or from SVELT_ALT_SEQ).
    #[arg(short, long)]
    pub vcf: Option<String>,

}

/// Options common to all commands
#[derive(Debug, Args)]
pub struct CommonOptions {
    /// Number of threads to use (0 to use all available)
    #[arg(long, required = false, default_value = "4")]
    pub threads: usize,
}

pub fn make_session_context(options: &CommonOptions) -> SessionContext {
    let mut cfg = SessionConfig::new();
    cfg.options_mut().execution.target_partitions = options.threads;

    SessionContext::new_with_config(cfg)
}