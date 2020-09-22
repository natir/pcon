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

/* std use */
use std::io::Write;
use std::sync::atomic;

/* crate use */
use anyhow::Result;
//use anyhow::{anyhow, Context, Result};
use rayon::iter::ParallelBridge;
use rayon::prelude::*;

/* local use */
// use crate::error::IO::*;
// use crate::error::*;

pub type AtoCount = atomic::AtomicU8;
pub type Count = u8;

/// A counter of kmer based on cocktail crate 2bit conversion, canonicalisation and hashing.
/// If kmer occure more than 256 other occurence are ignored
#[derive(serde::Serialize, serde::Deserialize)]
pub struct Counter {
    pub k: u8,
    record_buffer_len: usize,
    count: Box<[AtoCount]>,
}

impl Counter {
    /// Create a new Counter for kmer size equal to k, record_buffer_len control the number of record read in same time
    pub fn new(k: u8, record_buffer_len: usize) -> Self {
        let tmp = vec![0u8; cocktail::kmer::get_hash_space_size(k) as usize];

        Self {
            k,
            record_buffer_len,
            count: unsafe {
                std::mem::transmute::<Box<[u8]>, Box<[AtoCount]>>(tmp.into_boxed_slice())
            },
        }
    }

    /// Read the given an instance of io::Read as a fasta format and count kmer init
    pub fn count_fasta<R>(&mut self, fasta: R)
    where
        R: std::io::Read,
    {
        let reader = bio::io::fasta::Reader::new(fasta);

        let mut iter = reader.records();
        let mut records = Vec::with_capacity(self.record_buffer_len);

        let mut end = false;
        loop {
            for _ in 0..self.record_buffer_len {
                if let Some(Ok(record)) = iter.next() {
                    records.push(record);
                } else {
                    end = true;
                    break;
                }
            }

            log::info!("Buffer len: {}", records.len());

            records.drain(..).par_bridge().for_each(|record| {
                if record.seq().len() > self.k as usize {
                    let tokenizer = cocktail::tokenizer::Canonical::new(record.seq(), self.k);

                    for canonical in tokenizer {
                        Counter::inc_canonic_ato(&self.count, canonical);
                    }
                }
            });

            records.clear();

            if end {
                break;
            }
        }
    }

    /// Increase the counter of a kmer
    pub fn inc(&mut self, kmer: u64) {
        self.inc_canonic(cocktail::kmer::canonical(kmer, self.k));
    }

    /// Increase the counter of a canonical kmer
    pub fn inc_canonic(&mut self, canonical: u64) {
        Counter::inc_canonic_ato(&self.count, canonical);
    }

    fn inc_canonic_ato(count: &[atomic::AtomicU8], canonical: u64) {
        let hash = (canonical >> 1) as usize;

        if count[hash].load(std::sync::atomic::Ordering::SeqCst) != std::u8::MAX {
            count[hash].fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        }
    }

    /// Get the counter of a kmer
    pub fn get(&self, kmer: u64) -> Count {
        self.get_canonic(cocktail::kmer::canonical(kmer, self.k))
    }

    /// Get the counter of a canonical kmer
    pub fn get_canonic(&self, canonical: u64) -> Count {
        let hash = (canonical >> 1) as usize;

        self.count[hash].load(atomic::Ordering::SeqCst)
    }

    pub(crate) fn get_raw_count(&self) -> &[AtoCount] {
        &self.count
    }

    /// Serialize counter in given [std::io::Write]
    pub fn serialize<W>(&self, w: W, min_abundance: u8) -> Result<()>
    where
        W: std::io::Write,
    {
        let mut writer = niffler::get_writer(
            Box::new(w),
            niffler::compression::Format::Gzip,
            niffler::compression::Level::Two,
        )?;

        writer.write_all(&[self.k])?;

        for (hash, count) in self.count.iter().enumerate() {
            let val = count.load(atomic::Ordering::SeqCst);
            if val > min_abundance {
                let bytes = unsafe { std::mem::transmute::<u64, [u8; 8]>(hash as u64) };

                writer.write_all(&bytes)?;
                writer.write_all(&[val])?;
            }
        }

        writer.flush()?;

        Ok(())
    }

    /// Deserialize counter for given [std::io::Read]
    pub fn deserialize<R>(r: R) -> Result<Self>
    where
        R: std::io::Read,
    {
        let mut reader = niffler::get_reader(Box::new(r))?.0;

        let mut buf_u8: [u8; 1] = [0; 1];
        let mut buf_u64: [u8; 8] = [0; 8];

        reader.read_exact(&mut buf_u8)?;

        let k = buf_u8[0];

        let mut tmp = vec![0u8; cocktail::kmer::get_hash_space_size(k) as usize];

        loop {
            if reader.read_exact(&mut buf_u64).is_err() {
                break;
            }
            if reader.read_exact(&mut buf_u8).is_err() {
                break;
            };
            unsafe { tmp[std::mem::transmute::<[u8; 8], u64>(buf_u64) as usize] = buf_u8[0] };
        }

        Ok(Self {
            k,
            record_buffer_len: 8192,
            count: unsafe {
                std::mem::transmute::<Box<[u8]>, Box<[AtoCount]>>(tmp.into_boxed_slice())
            },
        })
    }
}

#[cfg(test)]
mod tests {

    const FASTA_FILE: &[u8] = b">random_seq 0
GTTCTGCAAATTAGAACAGACAATACACTGGCAGGCGTTGCGTTGGGGGAGATCTTCCGTAACGAGCCGGCATTTGTAAGAAAGAGATTTCGAGTAAATG
>random_seq 1
AGGATAGAAGCTTAAGTACAAGATAATTCCCATAGAGGAAGGGTGGTATTACAGTGCCGCCTGTTGAAAGCCCCAATCCCGCTTCAATTGTTGAGCTCAG
";

    const FASTA_COUNT: &[u8] = &[
        31, 139, 8, 0, 0, 0, 0, 0, 0, 255, 52, 138, 59, 72, 16, 80, 24, 133, 59, 127, 106, 5, 189,
        232, 253, 24, 162, 41, 200, 33, 130, 202, 10, 26, 172, 104, 72, 133, 178, 208, 130, 162,
        130, 178, 162, 69, 40, 42, 162, 197, 158, 208, 214, 32, 136, 162, 136, 155, 14, 186, 248,
        0, 5, 81, 113, 20, 20, 5, 69, 81, 196, 193, 23, 34, 234, 230, 3, 212, 115, 240, 158, 111,
        248, 238, 199, 249, 111, 102, 198, 174, 29, 144, 41, 19, 100, 201, 4, 7, 101, 18, 135, 101,
        18, 199, 100, 18, 199, 101, 18, 39, 100, 130, 83, 50, 193, 105, 153, 224, 140, 76, 112, 86,
        38, 113, 78, 38, 56, 47, 147, 184, 32, 19, 100, 203, 4, 57, 50, 193, 77, 153, 32, 87, 38,
        184, 45, 19, 220, 151, 9, 242, 100, 130, 124, 153, 68, 129, 76, 226, 161, 76, 80, 40, 19,
        60, 150, 9, 138, 100, 130, 98, 153, 224, 37, 37, 240, 74, 38, 120, 45, 147, 120, 35, 19,
        148, 200, 4, 111, 101, 130, 15, 50, 65, 169, 76, 240, 89, 38, 248, 34, 19, 124, 149, 9,
        202, 100, 18, 63, 100, 18, 63, 101, 18, 191, 100, 130, 127, 50, 193, 127, 153, 160, 92, 38,
        81, 33, 19, 84, 203, 36, 106, 100, 18, 117, 50, 65, 189, 76, 208, 32, 147, 104, 148, 9, 90,
        101, 130, 54, 153, 160, 93, 38, 232, 144, 9, 122, 101, 130, 62, 153, 160, 95, 38, 24, 144,
        9, 6, 101, 130, 17, 153, 96, 84, 38, 24, 147, 9, 38, 100, 18, 51, 50, 193, 172, 76, 48, 39,
        19, 44, 200, 36, 22, 101, 130, 101, 153, 96, 93, 38, 1, 232, 225, 178, 219, 145, 225, 200,
        114, 236, 113, 236, 77, 17, 251, 188, 236, 119, 28, 112, 28, 114, 28, 113, 28, 117, 156,
        116, 92, 116, 100, 59, 46, 57, 46, 167, 136, 43, 94, 174, 58, 174, 57, 114, 28, 215, 29,
        55, 28, 183, 28, 185, 41, 226, 142, 151, 187, 142, 123, 142, 60, 71, 126, 138, 40, 240,
        242, 192, 81, 232, 120, 228, 40, 114, 20, 59, 158, 164, 136, 167, 94, 158, 57, 158, 59, 94,
        248, 79, 137, 227, 157, 79, 239, 29, 165, 62, 125, 244, 242, 201, 241, 205, 241, 221, 81,
        150, 34, 126, 123, 249, 227, 248, 235, 40, 119, 84, 58, 170, 82, 68, 141, 151, 90, 71, 189,
        163, 193, 209, 148, 34, 154, 29, 45, 62, 117, 58, 186, 28, 221, 142, 30, 199, 80, 138, 24,
        118, 140, 248, 52, 238, 152, 116, 76, 249, 207, 180, 151, 25, 199, 156, 99, 222, 177, 228,
        88, 113, 172, 58, 214, 28, 27, 142, 77, 199, 150, 99, 27, 0, 0, 255, 255, 3, 0, 101, 172,
        200, 158, 143, 5, 0, 0,
    ];

    #[test]
    fn count_fasta() {
        let mut counter = crate::counter::Counter::new(5, 8192);

        counter.count_fasta(FASTA_FILE);

        let mut outfile = Vec::new();
        counter.serialize(&mut outfile, 0).unwrap();

        assert_eq!(&outfile[..], &FASTA_COUNT[..]);
    }

    const FASTQ_FILE: &[u8] = b"@random_seq 0
CCAGTAGCTTGGTGTACCGACGCTGTAGAGTTACAGTCTCGCGTGGATATAAGCTACTATCGACAGCAGGGTACGTTGTGAGTAATCTAACGTCATCTCT
+
X-Ee--b`x6h6Yy,c7S`p_C*K~SAgs;LFManA`-oLL6!Pgi>`X'P~6np^M1jQ+xQc9.ZCTEn+Yy?5r+b|ta=EyHil%Z}>>(%y\\=IC
@random_seq 1
TCAAATTGGCCGCCGCACAGTGAACCCGGAACTAAACAAGCACCGCACCGTTTGGTACACTTGAACACCGTATAAATTCATGGTGTTTATAAGCCAATGG
+
<nS{ADz8'0V$`mP@o]1n7f6(!nq%CpN^Vq)EV,{gz:aQ`jSc/l&(ZYi3\\OFBM<Ee?%#:jdKF]bR#{5Yj\"[}B!AG2.QZU9xyU3;\"y
";

    const FASTQ_COUNT: &[u8] = &[
        5, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 2, 2, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0,
        1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 2, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 2, 0, 0,
        1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 2, 1, 1, 1, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 2, 1, 0, 0, 0, 1,
        0, 0, 0, 1, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 3, 0, 1, 2, 1, 0, 1, 0,
        0, 2, 0, 0, 2, 1, 0, 1, 0, 0, 1, 4, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0,
        1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
        0, 2, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 2, 1, 1, 0, 0, 1, 2, 1, 0, 0, 1, 0, 0, 1,
        2, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 2, 1, 1, 1, 2, 1, 1, 0, 2, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 3, 0, 0, 0, 0,
        0, 2, 3, 0, 0, 0, 0, 0, 1, 0, 0,
    ];

    #[test]
    #[ignore]
    fn count_fastq() {
        let mut counter = crate::counter::Counter::new(5, 8192);

        counter.count_fasta(FASTQ_FILE);

        let mut outfile = Vec::new();
        counter.serialize(&mut outfile, 0).unwrap();

        assert_eq!(&outfile[..], &FASTQ_COUNT[..]);
    }

    lazy_static::lazy_static! {
    static ref COUNTER: crate::counter::Counter = {
            let mut counter = crate::counter::Counter::new(5, 8192);

            for i in 0..cocktail::kmer::get_kmer_space_size(5) {
        counter.inc(i);
            }

            counter
    };
    }

    const ALLKMERSEEONE: &[u8] = &[
        31, 139, 8, 0, 0, 0, 0, 0, 0, 255, 52, 150, 211, 18, 32, 6, 0, 3, 123, 219, 94, 109, 219,
        182, 109, 219, 182, 109, 91, 103, 219, 182, 109, 219, 182, 109, 219, 106, 50, 211, 236,
        195, 102, 63, 33, 5, 247, 251, 31, 10, 216, 2, 108, 193, 254, 182, 224, 0, 91, 80, 208, 22,
        28, 104, 11, 14, 178, 5, 7, 219, 130, 67, 108, 193, 161, 182, 224, 48, 91, 112, 184, 45,
        56, 194, 22, 28, 105, 11, 142, 178, 5, 71, 219, 130, 99, 108, 193, 177, 182, 224, 56, 91,
        112, 188, 45, 56, 193, 22, 156, 104, 11, 78, 178, 5, 39, 219, 130, 83, 108, 193, 169, 182,
        224, 52, 91, 112, 186, 45, 56, 195, 22, 156, 105, 11, 206, 178, 5, 103, 75, 134, 115, 108,
        193, 185, 182, 224, 60, 91, 112, 190, 45, 184, 192, 22, 92, 104, 11, 46, 178, 5, 23, 219,
        130, 75, 108, 193, 165, 182, 224, 50, 91, 112, 185, 45, 184, 194, 22, 92, 105, 11, 174,
        178, 5, 87, 219, 130, 107, 108, 193, 181, 182, 224, 58, 91, 112, 189, 45, 184, 193, 22,
        220, 104, 11, 110, 178, 5, 55, 219, 130, 91, 108, 193, 173, 182, 224, 54, 91, 112, 187, 45,
        184, 195, 22, 220, 105, 11, 238, 178, 5, 119, 75, 134, 123, 108, 193, 189, 182, 224, 62,
        91, 112, 191, 45, 120, 192, 22, 60, 104, 11, 30, 178, 5, 15, 219, 130, 71, 108, 193, 163,
        182, 224, 49, 91, 240, 184, 45, 120, 194, 22, 60, 105, 11, 158, 178, 5, 79, 219, 130, 103,
        108, 193, 179, 182, 224, 57, 91, 240, 188, 45, 120, 193, 22, 188, 104, 11, 94, 178, 5, 47,
        219, 130, 87, 108, 193, 171, 182, 224, 53, 91, 240, 186, 45, 120, 195, 22, 188, 105, 11,
        222, 178, 5, 111, 75, 134, 119, 108, 193, 187, 182, 224, 61, 91, 240, 190, 45, 248, 192,
        22, 124, 104, 11, 62, 178, 5, 31, 219, 130, 79, 108, 193, 167, 182, 224, 51, 91, 240, 185,
        45, 248, 194, 22, 124, 105, 11, 190, 178, 5, 95, 219, 130, 111, 108, 193, 183, 182, 224,
        59, 91, 240, 189, 45, 248, 193, 22, 252, 104, 11, 126, 178, 5, 63, 219, 130, 95, 108, 193,
        175, 182, 224, 55, 91, 240, 187, 45, 248, 195, 22, 252, 105, 11, 254, 178, 5, 127, 75, 134,
        127, 108, 193, 191, 182, 160, 144, 45, 40, 108, 11, 138, 216, 130, 162, 182, 160, 152, 45,
        40, 110, 11, 74, 216, 130, 146, 182, 160, 148, 45, 40, 109, 11, 202, 216, 130, 178, 182,
        160, 156, 45, 40, 111, 11, 42, 216, 130, 138, 182, 160, 146, 45, 168, 108, 11, 170, 216,
        130, 170, 182, 160, 154, 45, 168, 110, 11, 106, 216, 130, 154, 182, 160, 150, 45, 168, 109,
        11, 234, 216, 130, 186, 182, 160, 158, 45, 168, 47, 25, 26, 216, 130, 134, 182, 160, 145,
        45, 104, 108, 11, 154, 216, 130, 166, 182, 160, 153, 45, 104, 110, 11, 90, 216, 130, 150,
        182, 160, 149, 45, 104, 109, 11, 218, 216, 130, 182, 182, 160, 157, 45, 104, 111, 11, 58,
        216, 130, 142, 182, 160, 147, 45, 232, 108, 11, 186, 216, 130, 174, 182, 160, 155, 45, 232,
        110, 11, 122, 216, 130, 158, 182, 160, 151, 45, 232, 109, 11, 250, 216, 130, 190, 182, 160,
        159, 45, 232, 47, 25, 6, 216, 130, 129, 182, 96, 144, 45, 24, 108, 11, 134, 216, 130, 161,
        182, 96, 152, 45, 24, 110, 11, 70, 216, 130, 145, 182, 96, 148, 45, 24, 109, 11, 198, 216,
        130, 177, 182, 96, 156, 45, 24, 111, 11, 38, 216, 130, 137, 182, 96, 146, 45, 152, 108, 11,
        166, 216, 130, 169, 182, 96, 154, 45, 152, 110, 11, 102, 216, 130, 153, 182, 96, 150, 45,
        152, 109, 11, 230, 216, 130, 185, 182, 96, 158, 45, 152, 47, 25, 22, 216, 130, 133, 182,
        96, 145, 45, 88, 108, 11, 150, 216, 130, 165, 182, 96, 153, 45, 88, 110, 11, 86, 216, 130,
        149, 182, 96, 149, 45, 88, 109, 11, 214, 216, 130, 181, 182, 96, 157, 45, 88, 111, 11, 54,
        216, 130, 141, 182, 96, 147, 45, 216, 108, 11, 182, 216, 130, 173, 182, 96, 155, 45, 216,
        110, 11, 118, 216, 130, 157, 182, 96, 151, 45, 216, 109, 11, 246, 216, 130, 189, 182, 96,
        159, 45, 216, 175, 128, 71, 81, 32, 145, 39, 196, 254, 137, 3, 18, 5, 19, 7, 38, 14, 74,
        28, 156, 56, 36, 113, 104, 226, 176, 196, 225, 137, 35, 18, 71, 38, 142, 74, 28, 157, 56,
        38, 113, 108, 226, 184, 196, 241, 137, 19, 18, 39, 38, 78, 74, 156, 156, 56, 37, 113, 106,
        226, 180, 196, 233, 137, 51, 18, 103, 38, 206, 74, 156, 157, 56, 39, 113, 110, 226, 188,
        196, 249, 137, 11, 18, 23, 38, 46, 74, 92, 156, 184, 36, 113, 105, 226, 178, 196, 229, 137,
        43, 18, 87, 38, 174, 74, 92, 157, 184, 38, 113, 109, 226, 186, 196, 245, 137, 27, 18, 55,
        38, 110, 74, 220, 156, 184, 37, 113, 107, 226, 182, 196, 237, 137, 59, 18, 119, 38, 238,
        74, 220, 157, 184, 39, 113, 111, 226, 190, 196, 253, 137, 7, 18, 15, 38, 30, 74, 60, 156,
        120, 36, 241, 104, 226, 177, 196, 227, 137, 39, 18, 79, 38, 158, 74, 60, 157, 120, 38, 241,
        108, 226, 185, 196, 243, 137, 23, 18, 47, 38, 94, 74, 188, 156, 120, 37, 241, 106, 226,
        181, 196, 235, 137, 55, 18, 111, 38, 222, 74, 188, 157, 120, 39, 241, 110, 226, 189, 196,
        251, 137, 15, 18, 31, 38, 62, 74, 124, 156, 248, 36, 241, 105, 226, 179, 196, 231, 137, 47,
        18, 95, 38, 190, 74, 124, 157, 248, 38, 241, 109, 226, 187, 196, 247, 137, 31, 18, 63, 38,
        126, 74, 252, 156, 248, 37, 241, 107, 226, 183, 196, 239, 137, 63, 18, 127, 38, 254, 74,
        252, 157, 248, 39, 241, 111, 162, 80, 162, 112, 162, 72, 162, 104, 162, 88, 162, 120, 162,
        68, 162, 100, 162, 84, 162, 116, 162, 76, 162, 108, 162, 92, 162, 124, 162, 66, 162, 98,
        162, 82, 162, 114, 162, 74, 162, 106, 162, 90, 162, 122, 162, 70, 162, 102, 162, 86, 162,
        118, 162, 78, 162, 110, 162, 94, 162, 126, 162, 65, 162, 97, 162, 81, 162, 113, 162, 73,
        162, 105, 162, 89, 162, 121, 162, 69, 162, 101, 162, 85, 162, 117, 162, 77, 162, 109, 162,
        93, 162, 125, 162, 67, 162, 99, 162, 83, 162, 115, 162, 75, 162, 107, 162, 91, 162, 123,
        162, 71, 162, 103, 162, 87, 162, 119, 162, 79, 162, 111, 162, 95, 162, 127, 98, 64, 98, 96,
        98, 80, 98, 112, 98, 72, 98, 104, 98, 88, 98, 120, 98, 68, 98, 100, 98, 84, 98, 116, 98,
        76, 98, 108, 98, 92, 98, 124, 98, 66, 98, 98, 98, 82, 98, 114, 98, 74, 98, 106, 98, 90, 98,
        122, 98, 70, 98, 102, 98, 86, 98, 118, 98, 78, 98, 110, 98, 94, 98, 126, 98, 65, 98, 97,
        98, 81, 98, 113, 98, 73, 98, 105, 98, 89, 98, 121, 98, 69, 98, 101, 98, 85, 98, 117, 98,
        77, 98, 109, 98, 93, 98, 125, 98, 67, 98, 99, 98, 83, 98, 115, 98, 75, 98, 107, 98, 91, 98,
        123, 98, 71, 98, 103, 98, 87, 98, 119, 98, 79, 98, 111, 98, 95, 226, 63, 0, 0, 0, 255, 255,
        3, 0, 196, 128, 33, 253, 1, 18, 0, 0,
    ];

    #[test]
    fn serialize() {
        let mut outfile = Vec::new();

        let counter = &COUNTER;
        counter.serialize(&mut outfile, 0).unwrap();

        assert_eq!(&outfile[..], &ALLKMERSEEONE[..]);
    }

    #[test]
    fn deserialize() {
        let counter = crate::counter::Counter::deserialize(&ALLKMERSEEONE[..]).unwrap();

        assert_eq!(counter.k, COUNTER.k);

        for (a, b) in counter
            .get_raw_count()
            .iter()
            .zip(COUNTER.get_raw_count().iter())
        {
            assert_eq!(
                a.load(std::sync::atomic::Ordering::SeqCst),
                b.load(std::sync::atomic::Ordering::SeqCst)
            );
        }
    }
}
