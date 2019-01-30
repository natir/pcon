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

/* std use */
use std::io::Write;

/* crate use */
use clap::{App, Arg};

/* project mod declaration */
mod convert;

fn main() {
    let matches = App::new("ssik")
        .version("0.1")
        .author("Pierre Marijon <pierre.marijon@inria.fr>")
        .about("Scorer for Stupidly Insufficiently long Kmer")
        .arg(
            Arg::with_name("input")
                .short("i")
                .long("input")
                .takes_value(true)
                .help("sequence input in fasta format"),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .takes_value(true)
                .help("path where kmer count was write"),
        )
        .arg(
            Arg::with_name("kmer-size")
                .short("k")
                .long("kmer-size")
                .takes_value(true)
                .default_value("13")
                .help(
                    "kmer size, if kmer size is even real value is equal to k -= 1, max value 63",
                ),
        )
        .get_matches();

    // parse argument
    let mut k = matches
        .value_of("kmer-size")
        .unwrap()
        .parse::<u8>()
        .unwrap();
    k -= !k & 1;
    if k > 63 {
        k = 63;
    }

    let reader = bio::io::fasta::Reader::new(std::io::BufReader::new(
        std::fs::File::open(matches.value_of("input").unwrap()).unwrap(),
    ));

    // init counter
    let mut kmer2count: Vec<u8> = vec![0; 1 << (k * 2 - 1)];

    // count
    for result in reader.records() {
        let record = result.unwrap();

        for subseq in record.seq().windows(k as usize) {
            let hash = hash(subseq, k) as usize;
            if kmer2count[hash] != 255 {
                kmer2count[hash] += 1;
            }
        }
    }

    // write result
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(matches.value_of("output").unwrap()).unwrap(),
    );
    let mut writer = csv::WriterBuilder::new().from_writer(out);

    for i in 0..(1 << (k * 2 - 1)) {
        if kmer2count[i as usize] == 0 {
            continue;
        }
        writer
            .write_record(&[reverse_hash(i, k), kmer2count[i as usize].to_string()])
            .unwrap();
    }
}

fn reverse_hash(mut kmer: u128, k: u8) -> String {
    kmer <<= 1;

    if !convert::parity_even(kmer) {
        kmer = kmer + 1;
    }

    return convert::bit2seq(kmer, k);
}

fn hash(kmer: &[u8], k: u8) -> u128 {
    return convert::cannonical(convert::seq2bit(kmer), k) >> 1;
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn reverse_hash_() {
        // 100011110 -> TAGGC
        assert_eq!(reverse_hash(0b100011110, 5), "TAGGC");
    }

    #[test]
    fn hash_() {
        // TAGGC -> 100011110
        assert_eq!(hash(b"TAGGC", 5), 0b100011110);

        // GCCTA -> 110101100
        assert_eq!(hash(b"GCCTA", 5), 0b100011110);
    }
}
