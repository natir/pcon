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
use anyhow::{anyhow, Context, Result};
use rayon::iter::ParallelBridge;
use rayon::prelude::*;

/* local use */
use crate::error::IO::*;
use crate::error::*;

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
        let tmp = vec![0u64; cocktail::kmer::get_hash_space_size(k) as usize];

        Self {
            k,
            record_buffer_len,
            count: unsafe {
                std::mem::transmute::<Box<[u64]>, Box<[AtoCount]>>(tmp.into_boxed_slice())
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
    pub fn serialize<W>(&self, writer: W) -> Result<()>
    where
        W: std::io::Write,
    {
        bincode::serialize_into(
            niffler::get_writer(
                Box::new(writer),
                niffler::compression::Format::Gzip,
                niffler::compression::Level::One,
            )?,
            self,
        )
        .with_context(|| Error::IO(ErrorDurringWrite))
        .with_context(|| anyhow!("Error durring serialize counter"))
    }

    /// Deserialize counter for given [std::io::Read]
    pub fn deserialize<R>(reader: R) -> Result<Self>
    where
        R: std::io::Read,
    {
        bincode::deserialize_from(niffler::get_reader(Box::new(reader))?.0)
            .with_context(|| Error::IO(ErrorDurringRead))
            .with_context(|| anyhow!("Error durring deserialize counter"))
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
        31, 139, 8, 0, 0, 0, 0, 0, 4, 255, 93, 142, 9, 110, 4, 49, 8, 4, 187, 144, 242, 143, 252,
        255, 149, 169, 6, 175, 34, 45, 3, 166, 15, 240, 248, 39, 191, 185, 152, 107, 61, 129, 54,
        107, 82, 121, 134, 0, 133, 172, 144, 144, 108, 229, 133, 220, 44, 225, 245, 192, 236, 6,
        241, 163, 86, 171, 246, 64, 121, 75, 239, 212, 37, 133, 51, 115, 248, 206, 29, 240, 158,
        101, 35, 88, 193, 193, 247, 28, 145, 217, 251, 90, 66, 115, 103, 175, 43, 46, 35, 162, 116,
        93, 61, 18, 34, 67, 148, 79, 60, 19, 77, 179, 143, 112, 133, 136, 147, 216, 243, 21, 186,
        81, 119, 141, 141, 68, 150, 222, 202, 130, 176, 119, 4, 52, 14, 163, 155, 198, 172, 122,
        232, 3, 9, 79, 8, 44, 44, 221, 90, 198, 36, 248, 229, 226, 70, 212, 146, 254, 231, 137, 1,
        14, 254, 159, 206, 220, 132, 14, 145, 57, 212, 140, 65, 240, 52, 93, 19, 193, 31, 33, 156,
        6, 64, 17, 2, 0, 0,
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
        31, 139, 8, 0, 0, 0, 0, 0, 4, 255, 237, 208, 49, 13, 0, 0, 12, 2, 65, 66, 82, 31, 245, 175,
        146, 1, 92, 192, 45, 191, 255, 225, 97, 76, 56, 205, 7, 4, 100, 224, 155, 94, 17, 2, 0, 0,
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
