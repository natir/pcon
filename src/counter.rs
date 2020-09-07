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
        let mut tmp = Vec::with_capacity(cocktail::kmer::get_hash_space_size(k) as usize);
        for _ in 0..cocktail::kmer::get_hash_space_size(k) {
            tmp.push(AtoCount::new(0));
        }

        Self {
            k,
            record_buffer_len,
            count: tmp.into_boxed_slice(),
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

            log::debug!("Buffer len: {}", records.len());

            records.par_iter().for_each(|record| {
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
        let mut old = count[hash].load(atomic::Ordering::SeqCst);

        if old == 255 {
            return;
        }

        while count[hash]
            .compare_exchange(
                old,
                old + 1,
                atomic::Ordering::SeqCst,
                atomic::Ordering::Acquire,
            )
            .is_err()
        {
            old = count[hash].load(atomic::Ordering::SeqCst);

            if old == 255 {
                return;
            }
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

    pub(crate) fn get_raw_count(&self) -> &Box<[AtoCount]> {
        &self.count
    }

    /// Serialize counter in given [std::io::Write]
    pub fn serialize<W>(&self, writer: W) -> Result<()>
    where
        W: std::io::Write,
    {
        bincode::serialize_into(writer, self)
            .with_context(|| Error::IO(ErrorDurringWrite))
            .with_context(|| anyhow!("Error durring serialize counter"))
    }

    /// Deserialize counter for given [std::io::Read]
    pub fn deserialize<R>(reader: R) -> Result<Self>
    where
        R: std::io::Read,
    {
        bincode::deserialize_from(reader)
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
        5, 0, 32, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 2, 0, 2, 0, 0, 0, 2, 2, 1, 0, 1, 1, 1, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 2, 2, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 2, 1, 1,
        1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 1,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 2, 2, 0, 0, 0, 1,
        0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 2, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 2, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0,
        0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 2, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 2, 0, 0, 1, 1, 1, 1, 1,
        1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 2, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1,
        1, 2, 1, 0, 0, 1, 1, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 1, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 2, 0,
        0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 2, 0, 0, 2, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 2, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0,
        0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1,
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
        5, 0, 32, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
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
