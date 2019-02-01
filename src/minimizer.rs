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

/* std use*/
use std::io::Write;

/* project use */
use crate::convert;

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

    write(minimizer2count, output_path, k, abundance_min);
}

fn found_minimizer(subseq: &[u8], mut minimizer2count: &mut Vec<u8>, m: u8) -> () {
    let mut mini = u64::max_value();
    
    for subk in subseq.windows(m as usize) {
        let hash = hash(subk, m);
        
        if hash < mini {
            mini = hash;
        }
    }
    
    add_in_counter(&mut minimizer2count, unrevhash(mini));
}

fn add_in_counter(kmer2count: &mut Vec<u8>, mini: u64) -> () {
    if kmer2count[mini as usize] != 255 {
        kmer2count[mini as usize] += 1;
    }
}

fn write(kmer2count: Vec<u8>, output_path: &str, k: u8, abundance_min: u8) -> () {
    let mut out = std::io::BufWriter::new(std::fs::File::create(output_path).unwrap());

    // write k in first bytes
    out.write(&[k])
        .expect("Error during write of count on disk");

    let mut last_write = 0;
    for (i, val) in kmer2count.iter().enumerate() {
        if val < &abundance_min {
            continue;
        }
        let dist = i - last_write;
        last_write = i;

        // write dist to last value and count of k
        out.write(&[dist as u8, *val])
            .expect("Error durring write count on disk");
    }
}

fn hash(kmer: &[u8], k: u8) -> u64 {
    return revhash(convert::cannonical(convert::seq2bit(kmer), k) >> 1);
}

fn revhash(mut x: u64) -> u64 {
    x = ((x >> 32) ^ x).wrapping_mul(0xD6E8FEB86659FD93);
    x = ((x >> 32) ^ x).wrapping_mul(0xD6E8FEB86659FD93);
    x = (x >> 32) ^ x;

    return x;
}

fn unrevhash(mut x: u64) -> u64 {
    x = ((x >> 32) ^ x).wrapping_mul(0xCFEE444D8B59A89B);
    x = ((x >> 32) ^ x).wrapping_mul(0xCFEE444D8B59A89B);
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
