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
    pub fn serialize<W>(&self, w: W) -> Result<()>
    where
        W: std::io::Write,
    {
        let mut writer = niffler::get_writer(
            Box::new(w),
            niffler::compression::Format::Gzip,
            niffler::compression::Level::One,
        )?;

        writer.write(&[self.k])?;

        for (hash, count) in self.count.iter().enumerate() {
            let val = count.load(atomic::Ordering::SeqCst);
            if val > 0 {
                let bytes = unsafe { std::mem::transmute::<u64, [u8; 8]>(hash as u64) };

                writer.write(&bytes)?;
                writer.write(&[val])?;
            }
        }

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
        31, 139, 8, 0, 0, 0, 0, 0, 4, 255, 53, 208, 73, 40, 39, 0, 28, 197, 113, 191, 223, 96, 102,
        106, 182, 102, 95, 14, 211, 156, 166, 198, 97, 82, 246, 114, 176, 228, 96, 41, 91, 150, 34,
        148, 53, 23, 69, 72, 46, 246, 114, 115, 80, 34, 146, 27, 7, 46, 150, 162, 132, 28, 21, 81,
        68, 36, 7, 91, 18, 110, 150, 194, 225, 255, 190, 239, 242, 62, 215, 247, 66, 130, 131, 2,
        177, 144, 64, 7, 89, 40, 120, 39, 248, 7, 240, 25, 124, 1, 95, 5, 251, 14, 126, 128, 159,
        224, 151, 224, 191, 5, 251, 35, 248, 95, 193, 194, 64, 52, 136, 3, 9, 32, 17, 164, 130, 52,
        144, 46, 120, 6, 200, 18, 44, 27, 228, 130, 60, 144, 15, 74, 65, 25, 40, 23, 188, 66, 176,
        74, 80, 5, 106, 65, 29, 104, 4, 77, 160, 25, 180, 9, 222, 14, 58, 64, 167, 96, 189, 160,
        15, 244, 11, 62, 32, 216, 176, 224, 35, 96, 76, 176, 113, 48, 33, 248, 164, 96, 179, 96,
        14, 204, 131, 5, 176, 10, 214, 192, 58, 216, 0, 155, 96, 7, 236, 130, 61, 112, 32, 248,
        137, 96, 167, 224, 12, 92, 8, 126, 41, 216, 53, 184, 23, 220, 44, 32, 123, 1, 130, 65, 40,
        120, 9, 94, 9, 254, 90, 176, 55, 224, 45, 120, 15, 62, 130, 79, 224, 27, 248, 7, 194, 192,
        127, 16, 46, 120, 132, 96, 145, 32, 10, 68, 131, 24, 16, 11, 226, 65, 130, 224, 73, 130,
        37, 131, 20, 144, 6, 210, 5, 207, 16, 44, 211, 116, 75, 54, 200, 1, 121, 32, 31, 20, 8, 94,
        40, 88, 17, 40, 6, 37, 130, 87, 130, 106, 193, 106, 64, 157, 224, 245, 130, 53, 128, 22,
        208, 10, 218, 4, 239, 18, 172, 27, 244, 128, 126, 48, 8, 134, 4, 31, 17, 108, 212, 180,
        116, 28, 76, 128, 41, 193, 167, 193, 140, 96, 139, 96, 9, 44, 131, 21, 176, 37, 248, 54,
        216, 17, 108, 31, 28, 130, 35, 193, 143, 5, 59, 1, 103, 224, 28, 92, 129, 27, 112, 11, 238,
        192, 3, 120, 4, 79, 224, 25, 101, 172, 200, 158, 143, 5, 0, 0,
    ];

    #[test]
    fn count_fasta() {
        let mut counter = crate::counter::Counter::new(5, 8192);

        counter.count_fasta(FASTA_FILE);

        let mut outfile = Vec::new();
        counter.serialize(&mut outfile).unwrap();

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
    #[should_panic("fastq isn't support now")]
    fn count_fastq() {
        let mut counter = crate::counter::Counter::new(5, 8192);

        counter.count_fasta(FASTQ_FILE);

        let mut outfile = Vec::new();
        counter.serialize(&mut outfile).unwrap();

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
        31, 139, 8, 0, 0, 0, 0, 0, 4, 255, 53, 214, 195, 178, 32, 6, 0, 4, 192, 228, 37, 89, 219,
        182, 109, 219, 182, 109, 219, 182, 109, 219, 182, 109, 219, 182, 109, 35, 151, 237, 57,
        117, 77, 205, 7, 204, 127, 127, 253, 73, 192, 223, 16, 0, 255, 192, 191, 96, 29, 16, 72,
        19, 24, 130, 64, 80, 8, 6, 193, 33, 4, 132, 132, 80, 16, 26, 194, 64, 88, 8, 7, 225, 33, 2,
        68, 132, 72, 16, 25, 162, 64, 84, 136, 6, 209, 33, 6, 196, 132, 88, 16, 27, 226, 64, 92,
        136, 7, 241, 33, 1, 36, 132, 68, 144, 24, 146, 64, 82, 72, 6, 201, 33, 5, 164, 132, 84,
        144, 26, 210, 64, 90, 72, 7, 233, 33, 3, 100, 132, 76, 144, 25, 178, 64, 86, 200, 6, 217,
        33, 7, 228, 132, 92, 144, 27, 242, 64, 94, 200, 7, 249, 161, 0, 20, 132, 66, 80, 24, 138,
        64, 81, 40, 6, 197, 161, 4, 148, 132, 82, 80, 26, 202, 64, 89, 40, 7, 229, 161, 2, 84, 132,
        74, 80, 25, 170, 64, 85, 168, 6, 213, 161, 6, 212, 132, 90, 80, 27, 234, 64, 93, 168, 7,
        245, 161, 1, 52, 132, 70, 208, 24, 154, 64, 83, 104, 6, 205, 161, 5, 180, 132, 86, 208, 26,
        218, 64, 91, 104, 7, 237, 161, 3, 116, 132, 78, 208, 25, 186, 64, 87, 232, 6, 221, 161, 7,
        244, 132, 94, 208, 27, 250, 64, 95, 232, 7, 253, 97, 0, 12, 132, 65, 48, 24, 134, 192, 80,
        24, 6, 195, 97, 4, 140, 132, 81, 48, 26, 198, 192, 88, 24, 7, 227, 97, 2, 76, 132, 73, 48,
        25, 166, 192, 84, 152, 6, 211, 97, 6, 204, 132, 89, 48, 27, 230, 192, 92, 152, 7, 243, 97,
        1, 44, 132, 69, 176, 24, 150, 192, 82, 88, 6, 203, 97, 5, 172, 132, 85, 176, 26, 214, 192,
        90, 88, 7, 235, 97, 3, 108, 132, 77, 176, 25, 182, 192, 86, 216, 6, 219, 97, 7, 236, 132,
        93, 176, 27, 246, 192, 94, 216, 7, 251, 225, 0, 28, 132, 67, 112, 24, 142, 192, 81, 56, 6,
        199, 225, 4, 156, 132, 83, 112, 26, 206, 192, 89, 56, 7, 231, 225, 2, 92, 132, 75, 112, 25,
        174, 192, 85, 184, 6, 215, 225, 6, 220, 132, 91, 112, 27, 238, 192, 93, 184, 7, 247, 225,
        1, 60, 132, 71, 240, 24, 158, 192, 83, 120, 6, 207, 225, 5, 188, 132, 87, 240, 26, 222,
        192, 91, 120, 7, 239, 225, 3, 124, 132, 79, 240, 25, 190, 192, 87, 248, 6, 223, 225, 7,
        252, 132, 95, 240, 27, 254, 250, 243, 123, 2, 254, 134, 0, 248, 7, 254, 133, 255, 32, 16,
        4, 134, 32, 16, 20, 130, 65, 112, 8, 1, 33, 33, 20, 132, 134, 48, 16, 22, 194, 65, 120,
        136, 0, 17, 33, 18, 68, 134, 40, 16, 21, 162, 65, 116, 136, 1, 49, 33, 22, 196, 134, 56,
        16, 23, 226, 65, 124, 72, 0, 9, 33, 17, 36, 134, 36, 144, 20, 146, 65, 114, 72, 1, 41, 33,
        21, 164, 134, 52, 144, 22, 210, 65, 122, 200, 0, 25, 33, 19, 100, 134, 44, 144, 21, 178,
        65, 118, 200, 1, 57, 33, 23, 228, 134, 60, 144, 23, 242, 65, 126, 40, 0, 5, 161, 16, 20,
        134, 34, 80, 20, 138, 65, 113, 40, 1, 37, 161, 20, 148, 134, 50, 80, 22, 202, 65, 121, 168,
        0, 21, 161, 18, 84, 134, 42, 80, 21, 170, 65, 117, 168, 1, 53, 161, 22, 212, 134, 58, 80,
        23, 234, 65, 125, 104, 0, 13, 161, 17, 52, 134, 38, 208, 20, 154, 65, 115, 104, 1, 45, 161,
        21, 180, 134, 54, 208, 22, 218, 65, 123, 232, 0, 29, 161, 19, 116, 134, 46, 208, 21, 186,
        65, 119, 232, 1, 61, 161, 23, 244, 134, 62, 208, 23, 250, 65, 127, 24, 0, 3, 97, 16, 12,
        134, 33, 48, 20, 134, 193, 112, 24, 1, 35, 97, 20, 140, 134, 49, 48, 22, 198, 193, 120,
        152, 0, 19, 97, 18, 76, 134, 41, 48, 21, 166, 193, 116, 152, 1, 51, 97, 22, 204, 134, 57,
        48, 23, 230, 193, 124, 88, 0, 11, 97, 17, 44, 134, 37, 176, 20, 150, 193, 114, 88, 1, 43,
        97, 21, 172, 134, 53, 176, 22, 214, 193, 122, 216, 0, 27, 97, 19, 108, 134, 45, 176, 21,
        182, 193, 118, 216, 1, 59, 97, 23, 236, 134, 61, 176, 23, 246, 193, 126, 56, 0, 7, 225, 16,
        28, 134, 35, 112, 20, 142, 193, 113, 56, 1, 39, 225, 20, 156, 134, 51, 112, 22, 206, 193,
        121, 184, 0, 23, 225, 18, 92, 134, 43, 112, 21, 174, 193, 117, 184, 1, 55, 225, 22, 220,
        134, 59, 112, 23, 238, 193, 125, 120, 0, 15, 225, 17, 60, 134, 39, 240, 20, 158, 193, 115,
        120, 1, 47, 225, 21, 188, 134, 55, 240, 22, 222, 193, 123, 248, 0, 31, 225, 19, 124, 134,
        47, 240, 21, 190, 193, 119, 248, 1, 63, 225, 23, 252, 134, 255, 1, 196, 128, 33, 253, 1,
        18, 0, 0,
    ];

    #[test]
    fn serialize() {
        let mut outfile = Vec::new();

        let counter = &COUNTER;
        counter.serialize(&mut outfile).unwrap();

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
