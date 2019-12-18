/*
Copyright (c) 2019 Pierre Marijon <pierre.marijon@inria.fr>

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
#![feature(stdsimd)]

/* crate declaration */
extern crate bio;
extern crate clap;
extern crate csv;

/* std use */

/* crate use */
use clap::{App, Arg, SubCommand};

/* project mod declaration */
pub mod bucketizer;
pub mod count;
pub mod counter;
pub mod dump;

fn main() {
    let matches = App::new("pcon")
        .version("0.1")
        .author("Pierre Marijon <pierre.marijon@inria.fr>")
        .about("Prompt COuNter is a short kmer counter")
        .subcommand(SubCommand::with_name("count")
                    .about("count kmer in fasta file")
                    .arg(
                        Arg::with_name("input")
                            .short("i")
                            .long("input")
                            .multiple(true)
                            .required(true)
                            .takes_value(true)
                            .help("sequence input in fasta format")
                    )
                    .arg(
                        Arg::with_name("output")
                            .short("o")
                            .long("output")
                            .required(true)
                            .takes_value(true)
                            .help("path where kmer count was write")
                    )
                    .arg(
                        Arg::with_name("kmer-size")
                            .short("k")
                            .long("kmer-size")
                            .required(true)
                            .takes_value(true)
                            .help("kmer size, if kmer size is even real value is equal to k-1, max value 31")
                    )
                    .arg(
                        Arg::with_name("nb-bit-count")
                            .short("n")
                            .long("nb-bit-count")
                            .default_value("4")
                            .takes_value(true)
                            .possible_values(&["4", "8"])
                            .help("number of bit use to count kmer")
                    )
        )
        .subcommand(SubCommand::with_name("dump")
                    .about("take binary file produce by count step and generate a csv with kmer")
                    .arg(
                        Arg::with_name("input")
                            .short("i")
                            .long("input")
                            .required(true)
                            .takes_value(true)
                            .help("binary file generate by count step"),
                    )
                    .arg(
                        Arg::with_name("output")
                            .short("o")
                            .long("output")
                            .required(true)
                            .takes_value(true)
                            .help("path where kmer count was write"),
                    )
                    .arg(
                        Arg::with_name("abundance-min")
                            .short("a")
                            .long("abudance-min")
                            .takes_value(true)
                            .default_value("1")
                            .help("write only kmer with abudance is higher than this parametre")
                    )
                    .arg(
                        Arg::with_name("mode")
                            .short("m")
                            .long("mode")
                            .takes_value(true)
                            .possible_values(&["csv", "solidity", "spectrum"])
                            .help("write only kmer with abudance is higher than this parametre")
                    )
        )
        .get_matches();

    if let Some(count_matches) = matches.subcommand_matches("count") {
        let k = normalize_size_of_k(
            count_matches
                .value_of("kmer-size")
                .unwrap()
                .parse::<u8>()
                .unwrap(),
        );

        let n = count_matches
            .value_of("nb-bit-count")
            .unwrap()
            .parse::<u8>()
            .unwrap();

        count::count(
            count_matches.values_of("input").unwrap().collect(),
            count_matches.value_of("output").unwrap(),
            k,
            n,
        );
    } else if let Some(dump_matches) = matches.subcommand_matches("dump") {
        let abundance = dump_matches
            .value_of("abundance-min")
            .unwrap()
            .parse::<u8>()
            .unwrap();

        let mode = dump::Mode::from(dump_matches.value_of("mode").unwrap());

        dump::dump(
            dump_matches.value_of("input").unwrap(),
            dump_matches.value_of("output").unwrap(),
            abundance,
            mode,
        );
    }
}

fn normalize_size_of_k(k: u8) -> u8 {
    if k > 31 {
        return 31;
    }

    k - (!k & 1)
}
