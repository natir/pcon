//! Command Line Interface declaration of project pcon

/* std use */

/* crate use */

/* project use */

/// Prompt COuNter, a short kmer counter.
#[derive(clap::Parser, std::fmt::Debug)]
#[clap(
    name = "pcon",
    version = "0.1",
    author = "Pierre Marijon <pierre@marijon.fr>"
)]
pub struct Command {
    /// Silence all output
    #[clap(short = 'q', long = "quiet")]
    quiet: bool,

    /// Verbose mode (-v, -vv, -vvv, etc)
    #[clap(short = 'v', long = "verbosity", parse(from_occurrences))]
    verbosity: usize,

    /// Timestamp (sec, ms, ns, none)
    #[clap(short = 'T', long = "timestamp")]
    ts: Option<stderrlog::Timestamp>,
}

impl Command {
    /// Get verbosity level
    pub fn verbosity(&self) -> usize {
        self.verbosity
    }

    /// Get quiet
    pub fn quiet(&self) -> bool {
        self.quiet
    }

    /// Get timestamp granularity
    pub fn timestamp(&self) -> stderrlog::Timestamp {
        self.ts.unwrap_or(stderrlog::Timestamp::Off)
    }
}
