//! Prompt COuNter, a short kmer counter.

#![warn(missing_docs)]

/* std use */

/* crate use */

use anyhow::Context as _;
use clap::Parser as _;

/* project use */
use pcon::cli;
use pcon::count;
use pcon::dump;
use pcon::error;
#[cfg(not(feature = "parallel"))]
use pcon::minicount;

fn main() -> error::Result<()> {
    // parse cli
    let params = cli::Command::parse();

    // Setup logger
    stderrlog::new()
        .module(module_path!())
        .quiet(params.quiet())
        .verbosity(params.verbosity())
        .timestamp(params.timestamp())
        .init()
        .context("stderrlog already create a logger")?;

    #[cfg(feature = "parallel")]
    rayon::ThreadPoolBuilder::new()
        .num_threads(params.threads())
        .build_global()?;

    match params.subcommand {
        cli::SubCommand::Count(params) => count::count(params),
        #[cfg(not(feature = "parallel"))]
        cli::SubCommand::MiniCount(params) => minicount::minicount(params),
        cli::SubCommand::Dump(params) => dump::dump(params),
    }
}
