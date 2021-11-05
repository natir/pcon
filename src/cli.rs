/*
Copyright (c) 2020 Pierre Marijon <pmarijon@mpi-inf.mpg.de>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

#[derive(clap::Parser, Debug)]
#[clap(
    version = "0.1",
    author = "Pierre Marijon <pmarijon@mpi-inf.mpg.de>",
    about = "Prompt COuNter is short kmer counter"
)]
pub struct Command {
    #[clap(subcommand)]
    pub subcmd: SubCommand,

    #[clap(
        short = 't',
        long = "threads",
        about = "Number of thread use by pcon to count, 0 use all avaible core, default value 0"
    )]
    pub threads: Option<usize>,

    #[clap(
        short = 'v',
        long = "verbosity",
        parse(from_occurrences),
        about = "verbosity level also control by environment variable PCON_LOG if flag is set PCON_LOG value is ignored"
    )]
    pub verbosity: i8,
}

#[derive(clap::Parser, Debug)]
pub enum SubCommand {
    Count(SubCommandCount),
    Dump(SubCommandDump),
}

#[derive(clap::Parser, Debug)]
#[clap(about = "Perform kmer count")]
pub struct SubCommandCount {
    #[clap(short = 'k', long = "kmer-size", about = "Size of kmer size")]
    pub kmer: u8,

    #[clap(short = 'i', long = "inputs", about = "Path to inputs")]
    pub inputs: Vec<String>,

    #[clap(
        short = 'o',
        long = "output",
        about = "Path where count are store in binary format"
    )]
    pub output: Option<String>,

    #[clap(
        short = 'b',
        long = "record_buffer",
        about = "Number of sequence record load in buffer, default 8192"
    )]
    pub record_buffer: Option<usize>,

    #[clap(
        short = 'a',
        long = "abundance",
        default_value = "0",
        about = "Minimal abundance"
    )]
    pub abundance: crate::counter::Count,

    #[clap(short = 'c', long = "csv", about = "Path where count is write in csv")]
    pub csv: Option<String>,

    #[clap(
        short = 's',
        long = "solid",
        about = "Path where count is write in solid format"
    )]
    pub solid: Option<String>,

    #[clap(
        short = 'S',
        long = "spectrum",
        about = "Path where kmer spectrum is write"
    )]
    pub spectrum: Option<String>,
}

#[derive(clap::Parser, Debug)]
#[clap(about = "Convert count in usable format")]
pub struct SubCommandDump {
    #[clap(short = 'i', long = "input", about = "Path to count file")]
    pub input: String,

    #[clap(
        short = 'a',
        long = "abundance",
        default_value = "0",
        about = "Minimal abundance"
    )]
    pub abundance: crate::counter::Count,

    #[clap(short = 'c', long = "csv", about = "Path where count is write in csv")]
    pub csv: Option<String>,

    #[clap(
        short = 's',
        long = "solid",
        about = "Path where count is write in solid format"
    )]
    pub solid: Option<String>,

    #[clap(
        short = 'S',
        long = "spectrum",
        about = "Path where kmer spectrum is write"
    )]
    pub spectrum: Option<String>,

    #[clap(short = 'b', long = "bin", about = "Path where count is write in bin")]
    pub bin: Option<String>,
}

use crate::error::{Cli, Error};
use Cli::*;

pub fn check_count_param(params: SubCommandCount) -> Result<SubCommandCount, Error> {
    if (params.kmer & 1) == 0 {
        return Err(Error::Cli(KMustBeOdd));
    }

    if params.kmer > 32 {
        return Err(Error::Cli(KMustBeLower32));
    }

    Ok(params)
}

pub fn check_dump_param(params: SubCommandDump) -> Result<SubCommandDump, Error> {
    if ![&params.csv, &params.solid, &params.spectrum, &params.bin]
        .iter()
        .any(|x| x.is_some())
    {
        return Err(Error::Cli(ADumpOptionMustBeSet));
    }

    Ok(params)
}

pub fn i82level(level: i8) -> Option<log::Level> {
    match level {
        std::i8::MIN..=0 => None,
        1 => Some(log::Level::Error),
        2 => Some(log::Level::Warn),
        3 => Some(log::Level::Info),
        4 => Some(log::Level::Debug),
        5..=std::i8::MAX => Some(log::Level::Trace),
    }
}
