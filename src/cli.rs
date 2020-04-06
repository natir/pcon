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

/* crate use */
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
#[structopt(
    version = "0.1",
    author = "Pierre Marijon <pmarijon@mpi-inf.mpg.de>",
    about = "Prompt COuNter is short kmer counter"
)]
pub struct Command {
    #[structopt(subcommand)]
    pub subcmd: SubCommand,
}

#[derive(StructOpt, Debug)]
pub enum SubCommand {
    Count(SubCommandCount),
    Dump(SubCommandDump),
}

#[derive(StructOpt, Debug)]
#[structopt(about = "Perform kmer count")]
pub struct SubCommandCount {
    #[structopt(short = "k", long = "kmer-size", help = "Size of kmer size")]
    pub kmer: u8,

    #[structopt(short = "i", long = "inputs", help = "Path to inputs")]
    pub inputs: Vec<String>,

    #[structopt(short = "o", long = "output", help = "Path where count are store")]
    pub output: String,
}

#[derive(StructOpt, Debug)]
#[structopt(about = "Convert count in usable format")]
pub struct SubCommandDump {
    #[structopt(short = "i", long = "input", help = "Path to count file")]
    pub input: String,

    #[structopt(short = "o", long = "output", help = "Path to output")]
    pub output: String,

    #[structopt(
        short = "a",
        long = "abundance",
        default_value = "0",
        help = "Minimal abundance"
    )]
    pub abundance: crate::counter::Count,

    #[structopt(
	short = "m",
	long = "mode",
	possible_values = &["csv", "solid", "spectrum"],
	help = "Write mode:
- csv: a two column csv, first for kmer sequence, second for occurence of this kmer
- solid: an binary format to store only if kmer is solid or not
- spectrum: a two column csv, first for count, second for occurence of this count
",
	parse(from_str),
    )]
    pub mode: DumpMode,
}

#[derive(Debug)]
pub enum DumpMode {
    Csv,
    Solid,
    Spectrum,
}

impl From<&str> for DumpMode {
    fn from(s: &str) -> Self {
        match s {
            "csv" => DumpMode::Csv,
            "solid" => DumpMode::Solid,
            "spectrum" => DumpMode::Spectrum,
            _ => DumpMode::Csv,
        }
    }
}

use crate::error::{Cli, Error};
use Cli::*;

pub fn check_count_param(params: SubCommandCount) -> Result<SubCommandCount, Error> {
    if (params.kmer & 1) == 0 {
        return Err(Error::Cli(KMustBeOdd));
    }

    if params.kmer < 32 {
        return Err(Error::Cli(KMustBeLower32));
    }

    Ok(params)
}

pub fn check_dump_param(params: SubCommandDump) -> Result<SubCommandDump, Error> {
    Ok(params)
}
