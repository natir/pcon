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

/* project use */
use crate::convert;
use crate::write;

pub fn minimizer(input_path: &str, output_path: &str, k: u8, m: u8, abundance_min: u8) -> () {
    let reader = bio::io::fasta::Reader::new(std::io::BufReader::new(
        std::fs::File::open(input_path).unwrap(),
    ));

    // init counter
    let mut minimizer2count: Vec<u8> = vec![0; 1 << (m * 2 - 1)];

    // count
    for result in reader.records() {
        let record = result.unwrap();

        for subseq in record.seq().windows(k as usize) {
            found_minimizer(subseq, &mut minimizer2count, m);
        }
    }

    write::write(minimizer2count, output_path, m, abundance_min);
}

fn found_minimizer(subseq: &[u8], mut minimizer2count: &mut Vec<u8>, m: u8) -> () {
    let mut mini: u64;
    let mut mini_hash = u64::max_value();

    for subk in subseq.windows(m as usize) {
        let kmer = hash(subk, m);
        let hash = revhash(kmer);
        
        if hash < mini_hash {
            mini = kmer;
            mini_hash = hash;
        }
    }

    add_in_counter(&mut minimizer2count, mini);
}

fn add_in_counter(minimizer2count: &mut Vec<u8>, mini: u64) -> () {
    minimizer2count[mini as usize] = minimizer2count[mini as usize].saturating_add(1);
}

fn hash(kmer: &[u8], k: u8) -> u64 {
    return convert::cannonical(convert::seq2bit(kmer), k) >> 1;
}

fn revhash(mut x: u64) -> u64 {
    x = ((x >> 32) ^ x).wrapping_mul(0xD6E8FEB86659FD93);
    x = ((x >> 32) ^ x).wrapping_mul(0xD6E8FEB86659FD93);
    x = (x >> 32) ^ x;

    return x;
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn hash_() {
        // TAGGC -> 100011110
        assert_eq!(unrevhash(hash(b"TAGGC", 5)), 0b100011110);

        // GCCTA -> 110101100
        assert_eq!(unrevhash(hash(b"GCCTA", 5)), 0b100011110);
    }
}