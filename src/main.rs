/*
Copyright (c) 2018 Pierre Marijon <pierre.marijon@inria.fr>

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

extern crate bio;
extern crate clap;

use clap::{App, Arg};

fn main() {
    let matches = App::new("ssik")
        .version("0.1")
        .author("Pierre Marijon <pierre.marijon@inria.fr>")
        .about("Scorer for Stupidly Insufficiently long Kmer")
        .arg(Arg::with_name("input")
             .short("i")
             .long("input")
             .takes_value(true)
             .help("sequence input in fasta format")
        )
        .arg(Arg::with_name("output")
             .short("o")
             .long("output")
             .takes_value(true)
             .help("path where kmer count was write")
        )
        .arg(Arg::with_name("kmer-size")
             .short("k")
             .long("kmer-size")
             .takes_value(true)
             .default_value("13")
             .help("kmer size, if kmer size is even, k -= 1")
        )
        .get_matches();

    let reader = bio::io::fasta::Reader::new(
        std::io::BufReader::new(
            std::fs::File::open(matches.value_of("input").unwrap()
            ).unwrap()
        )
    );

    let mut k = matches.value_of("kmer-size").unwrap().parse::<u8>().unwrap();
    k -= !k & 1;

    let mut kmer2count: Vec<u8> = vec![0; 1 << (k * 2 - 1)];
    
    for result in reader.records() {
        let record = result.unwrap();
        
        for subseq in record.seq().windows(k as usize) {
            let hash = hash(subseq, k) as usize; 
            if kmer2count[hash] != 255 {
                kmer2count[hash] += 1;
            }
        }
    }

    let out = std::io::BufWriter::new(
        std::fs::File::create(
            matches.value_of("output").unwrap()
        ).unwrap()
    );
    let mut writer = csv::WriterBuilder::new().from_writer(out);
    
    for i in 0..(1 << (k * 2 - 1)) {
        if kmer2count[i as usize] == 0 {
            continue;
        }
        writer.write_record(&[hash2str(i, k), kmer2count[i as usize].to_string()]).unwrap();
    }
}

fn hash2str(mut kmer: u32, k: u8) -> String {
    kmer <<= 1;
    if !parity_even(kmer) {
        kmer = kmer + 1; 
    }

    let mut result = vec![0; k as usize];
    for i in 1..=k {
        let val = kmer & 0b11;

        if val == 0 {
            result[(k - i) as usize] = b'A';
        } else if val == 1 {
            result[(k - i) as usize] = b'C';
        } else if val == 2 {
            result[(k - i) as usize] = b'G';
        } else {
            result[(k - i) as usize] = b'T';
        }

        kmer >>= 2;
    }

    return String::from_utf8(result).unwrap();
}

fn hash(kmer: &[u8], k: u8) -> u32 {
    return cannonical(seq2bit(kmer), k) >> 1;
}

fn seq2bit(subseq: &[u8]) -> u32 {
  let mut kmer: u32 = 0;
  
  for n in subseq {
    kmer <<= 2;
    kmer |= (*n as u32 >> 1) & 0b11;
  }
  
  return kmer;
}

fn cannonical(kmer: u32, k: u8) -> u32 {
    if parity_even(kmer) {
        return kmer;
    } else {
        return revcomp(kmer, k);
    }
}

fn parity_even(kmer: u32) -> bool {
    return kmer.count_ones() % 2 == 0;
}

fn revcomp(kmer: u32, k: u8) -> u32 {
    return rev(comp(kmer), k);
}

fn comp(kmer: u32) -> u32 {
    return kmer ^ 0b10101010101010101010101010101010;
}

fn rev(kmer: u32, k: u8) -> u32 {
    let clean_move = 32 - k * 2;

    let mut reverse = reverse_2(kmer, k);
    reverse <<= clean_move;
    reverse >>= clean_move;
    
    return reverse;
}

fn reverse_2(mut kmer: u32, k: u8) -> u32 {
    let mut reversed: u32 = 0;
    
    for _ in 0..(k-1) {
        reversed = (reversed ^ (kmer & 0b11)) << 2;
        kmer >>= 2;
    }
  
    return reversed ^ (kmer & 0b11);
}
