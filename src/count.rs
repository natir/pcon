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
use crate::counter;
use crate::write;
use crate::bucketizer;

use crate::bucketizer::Bucket;

pub fn count(input_path: &str, output_path: &str, k: u8) -> () {
    let reader = bio::io::fasta::Reader::new(std::io::BufReader::new(std::fs::File::open(input_path).unwrap()));

    let mut counter: counter::ShortCounter = counter::ShortCounter::new(k);
    let mut bucketizer: bucketizer::Prefix<u8> = bucketizer::Prefix::new(&mut counter, k);
    
    for record in reader.records() {
        let line = record.unwrap();

        for subseq in line.seq().windows(k as usize) {
            bucketizer.add_kmer(convert::cannonical(convert::seq2bit(subseq), k));
        }
    }

    bucketizer.clean_all_buckets();

    write::write(&counter, output_path, k);
}

#[cfg(test)]
mod test {
    use super::*;
}
