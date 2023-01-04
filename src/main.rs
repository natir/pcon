//! Prompt COuNter, a short kmer counter.

#![warn(missing_docs)]

/* std use */

/* crate use */

use anyhow::Context as _;
use clap::Parser as _;

/* project use */
use pcon::cli;
use pcon::error;

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

    log::error!("Hello, word!");

    Ok(())
}
