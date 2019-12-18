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

/* std use */
use std::io::Write;

/* crate use */
use csv;

/* project use */
use crate::count;
use cocktail;

#[repr(C)]
#[derive(Clone, PartialEq)]
pub enum Mode {
    Csv,
    Solidity,
    Spectrum,
}

impl From<&str> for Mode {
    fn from(mode: &str) -> Self {
        match mode {
            "csv" => Mode::Csv,
            "solidity" => Mode::Solidity,
            "spectrum" => Mode::Spectrum,
            _ => Mode::Csv,
        }
    }
}

impl std::fmt::Display for Mode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match *self {
            Mode::Csv => write!(f, "csv"),
            Mode::Solidity => write!(f, "solidity"),
            Mode::Spectrum => write!(f, "spectrum"),
        }
    }
}

pub fn dump(input_path: &str, output_path: &str, abundance: u8, mode: Mode) {
    let mut reader = std::io::BufReader::new(std::fs::File::open(input_path).unwrap());

    let count = count::Count::deserialize(&mut reader);
    match mode {
        Mode::Csv => dump_csv(count, output_path, abundance),
        Mode::Solidity => dump_solidity(count, output_path, abundance),
        Mode::Spectrum => dump_spectrum(count, output_path),
    }
}

fn dump_csv(mut count: count::Count, output_path: &str, abundance: u8) {
    let out = std::io::BufWriter::new(std::fs::File::create(output_path).unwrap());
    let mut writer = csv::WriterBuilder::new().from_writer(out);

    for hash in 0..cocktail::kmer::get_kmer_space_size(count.get_k()) {
        let value = count.get_count(hash);

        if value >= abundance {
            let kmer = cocktail::kmer::kmer2seq(reverse_hash(hash), count.get_k());
            writer
                .write_record(&[kmer, value.to_string()])
                .expect("Error durring csv writting")
        }
    }
}

fn dump_solidity(count: count::Count, output_path: &str, abundance: u8) {
    let mut out = std::io::BufWriter::new(std::fs::File::create(output_path).unwrap());

    let bit_vec = count.generate_bitfield(abundance);
    out.write_all(&bit_vec.into_boxed_slice())
        .expect("Error durring solidity writting");
}

fn dump_spectrum(mut count: count::Count, output_path: &str) {
    let mut spectrum = vec![0; 256];

    for hash in 0..cocktail::kmer::get_kmer_space_size(count.get_k()) {
        let value = count.get_count(hash);

        spectrum[value as usize] += 1;
    }

    let out = std::io::BufWriter::new(std::fs::File::create(output_path).unwrap());
    let mut writer = csv::WriterBuilder::new().from_writer(out);
    writer
        .write_record(&["count", "nb_kmer"])
        .expect("Error durring spectrum writting");

    for (count, value) in spectrum.into_iter().enumerate() {
        writer
            .write_record(&[count.to_string(), value.to_string()])
            .expect("Error durring spectrum writting");
    }
}

pub fn reverse_hash(mut kmer: u64) -> u64 {
    kmer <<= 1;

    if !cocktail::kmer::parity_even(kmer) {
        return kmer + 1;
    }

    kmer
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn reverse_hash_() {
        assert_eq!(reverse_hash(0b100011110), 0b1000111101);
    }
}
