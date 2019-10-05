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

pub fn minimizer(input_path: &str, output_path: &str, k: u8, m: u8) -> () {
    let reader = bio::io::fasta::Reader::new(std::io::BufReader::new(
        std::fs::File::open(input_path).unwrap(),
    ));

    // init counter
    let mut counter: counter::BasicCounter<u16> = counter::BasicCounter::<u16>::new(m);
    let bucketizer: bucketizer::Prefix<u16> = bucketizer::Prefix::new(&mut counter, m);

    minimizer_work::<u16, counter::BasicCounter<u16>, bucketizer::Prefix<u16> ,std::fs::File>(reader, bucketizer, k, m);

    write::write(&counter, output_path, k);
}

fn minimizer_work<'a, T, C: counter::Counter<T,  u64>, B: bucketizer::Bucket<'a, T>, R: std::io::Read>(
    reader: bio::io::fasta::Reader<std::io::BufReader<R>>,
    mut bucketizer: B,
    k: u8,
    m: u8,
) -> () {
    let pos_begin_last_minimizer: usize = ((k as i16) - (m as i16)) as usize;

    // count
    for result in reader.records() {
        let record = result.unwrap();
        let mut last_minimizer: u64 = 0;
        let mut last_mini_pos: i64 = 0;
        let mut last_mini_hash: i64 = 0;

        for (i, subseq) in record.seq().windows(k as usize).enumerate() {
            last_mini_pos -= 1;
            if last_mini_pos < 0 {
                // minimizer isn't in kmer
                let tmp = found_minimizer(subseq, m);
                last_minimizer = tmp.0;
                last_mini_pos = tmp.1;
                last_mini_hash = tmp.2;
            } else {
                // test if new subkmer is minimizer
                let subk = &subseq[pos_begin_last_minimizer..];
                let kmer = hash(subk, m);
                let kash = revhash(kmer);

                if kash < last_mini_hash {
                    last_minimizer = kmer;
                    last_mini_pos = i as i64;
                    last_mini_hash = kash;
                }
            }

            bucketizer.add_bit(last_minimizer);
        }
    }

    bucketizer.clean_all_buckets();
}

fn found_minimizer(subseq: &[u8], m: u8) -> (u64, i64, i64) {
    let mut mini: u64 = 0;
    let mut mash: i64 = std::i64::MAX;
    let mut pos: i64 = 0;

    for (i, subk) in subseq.windows(m as usize).enumerate() {
        let kmer = hash(subk, m);
        let kash = revhash(kmer);

        if kash < mash {
            mini = kmer;
            mash = kash;
            pos = i as i64;
        }
    }

    return (mini, pos, mash);
}

fn hash(kmer: &[u8], k: u8) -> u64 {
    return convert::cannonical(convert::seq2bit(kmer), k) >> 1;
}

fn revhash(mut x: u64) -> i64 {
    x = ((x >> 32) ^ x).wrapping_mul(0xD6E8FEB86659FD93);
    x = ((x >> 32) ^ x).wrapping_mul(0xD6E8FEB86659FD93);
    x = (x >> 32) ^ x;

    return x as i64;
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn hash_() {
        // TAGGC -> 100011110
        assert_eq!(hash(b"TAGGC", 5), 0b100011110);

        // GCCTA -> 110101100
        assert_eq!(hash(b"GCCTA", 5), 0b100011110);
    }

    #[test]
    fn minimizer_() {
        let file: &[u8] = b">1\nAAACCCTTTGGG";

        let reader = bio::io::fasta::Reader::new(std::io::BufReader::new(file));

        let mut counter: counter::BasicCounter<u8> = counter::BasicCounter::<u8>::new(3);
        let bucketizer: counter::Bucketizer<u8> = counter::Bucketizer::new(&mut counter, 3);

        minimizer_work::<u8, counter::BasicCounter<u8>, &[u8]>(reader, bucketizer, 5, 3);

        assert_eq!(
            counter.data,
            [
                2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 4
            ]
        );
    }
}
