use clap::{ArgAction, Args};

/// Options controlling the merge process
#[derive(Debug, Args)]
pub struct Options {
    /// Allowed distance for merging of events
    #[arg(long, required = false, default_value = "25")]
    pub position_window: u32,

    /// Allowed end2 distance for merging BND events
    #[arg(long, required = false, default_value = "150")]
    pub end2_window: u32,

    /// Minimum length ratio (shorter/longer) for merging two events
    #[arg(long, required = false, default_value = "0.9")]
    pub length_ratio: f64,

    /// Force ALTs to be symbolic
    #[arg(long,
        action = ArgAction::Set,
        default_value_t = true,
        default_missing_value = "true",
        num_args = 0..=1,)]
    pub force_alt_tags: bool,

    /// Allow breakend variants to be flipped
    #[arg(long,
        action = ArgAction::Set,
        default_value_t = true,
        default_missing_value = "true",
        num_args = 0..=1,)]
    pub allow_breakend_flipping: bool,
}
