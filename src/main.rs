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

/* crate declaration */
extern crate bio;
extern crate clap;
extern crate csv;
#[macro_use]
extern crate itertools;

/* std use */

/* crate use */
use clap::{App, Arg, SubCommand};

/* project mod declaration */
mod convert;
mod count;
mod dump;
mod prepare;

fn main() {
    let matches = App::new("ssik")
        .version("0.1")
        .author("Pierre Marijon <pierre.marijon@inria.fr>")
        .about("Scorer for Stupidly Insufficiently long Kmer")
        .subcommand(SubCommand::with_name("count")
                    .about("count kmer in fasta file and generate a binary file with size 2^k-1 bytes")
                    .arg(
                        Arg::with_name("input")
                            .short("i")
                            .long("input")
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
                            .help("kmer size, if kmer size is even real value is equal to k-1, max value 63")
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
        )
        .subcommand(SubCommand::with_name("prepare")
                    .about("indicate for a k the memory usage of count")
                    .arg(
                        Arg::with_name("kmer-size")
                            .short("k")
                            .long("kmer-size")
                            .takes_value(true)
                            .help(
                                "kmer size you want use",
                            )
                    )
                    .arg(
                        Arg::with_name("minimizer-size")
                            .short("m")
                            .long("minimizer-size")
                            .takes_value(true)
                            .help(
                                "minimizer size you want use",
                            )
                    )
        )
        .get_matches();

    if let Some(count_matches) = matches.subcommand_matches("count") {
        let mut k = count_matches
            .value_of("kmer-size")
            .unwrap()
            .parse::<u8>()
            .unwrap();
        k -= !k & 1;
        if k > 63 {
            k = 63;
        }

        count::count(
            count_matches.value_of("input").unwrap(),
            count_matches.value_of("output").unwrap(),
            k,
            1,
        );
    } else if let Some(dump_matches) = matches.subcommand_matches("dump") {
        let abudance = dump_matches
            .value_of("abundance-min")
            .unwrap()
            .parse::<u8>()
            .unwrap();

        dump::dump(
            dump_matches.value_of("input").unwrap(),
            dump_matches.value_of("output").unwrap(),
            abudance,
        );
        
    } else if let Some(prepare_matches) = matches.subcommand_matches("prepare") {
        if prepare_matches.is_present("kmer-size") && prepare_matches.is_present("minimizer-size") {
            eprintln!("You can't use k and m in same time");
            return ();
        }

        if !prepare_matches.is_present("kmer-size") && !prepare_matches.is_present("minimizer-size")
        {
            eprintln!("You need use k or m");
            return ();
        }

        if prepare_matches.is_present("kmer-size") {
            let mut k = prepare_matches
                .value_of("kmer-size")
                .unwrap()
                .parse::<u8>()
                .unwrap();
            k -= !k & 1;
            if k > 63 {
                k = 63;
            }

            prepare::prepare(k, "k");

            return ();
        }

        if prepare_matches.is_present("minimizer-size") {
            let mut m = prepare_matches
                .value_of("minimizer-size")
                .unwrap()
                .parse::<u8>()
                .unwrap();
            m -= !m & 1;
            if m > 63 {
                m = 63;
            }

            prepare::prepare(m, "m");

            return ();
        }
    }
}
