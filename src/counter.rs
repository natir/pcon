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

use anyhow::{anyhow, Context, Result};

use crate::error::IO::*;
use crate::error::*;

pub type Count = u8;

/// A counter of kmer based on cocktail crate 2bit conversion, canonicalisation and hashing.
/// If kmer occure more than 256 other occurence are ignored
#[derive(serde::Serialize, serde::Deserialize)]
pub struct Counter {
    pub k: u8,
    count: Box<[Count]>,
}

impl Counter {
    /// Create a new Counter for kmer size equal to k
    pub fn new(k: u8) -> Self {
        Self {
            k,
            count: vec![0; cocktail::kmer::get_hash_space_size(k) as usize].into_boxed_slice(),
        }
    }

    /// Read the given an instance of io::Read as a fasta format and count kmer init
    pub fn count_fasta<R>(&mut self, fasta: R)
    where
        R: std::io::Read,
    {
        let reader = bio::io::fasta::Reader::new(fasta);

        let mut iter = reader.records();
        while let Some(Ok(record)) = iter.next() {
            let tokenizer = cocktail::tokenizer::Canonical::new(record.seq(), self.k);

            for canonical in tokenizer {
                self.canonical_inc(canonical);
            }
        }
    }

    /// Read the given an instance of io::Read as a fastq format and count kmer init
    pub fn count_fastq<R>(&mut self, fasta: R)
    where
        R: std::io::Read,
    {
        let reader = bio::io::fastq::Reader::new(fasta);

        let mut iter = reader.records();
        while let Some(Ok(record)) = iter.next() {
            let tokenizer = cocktail::tokenizer::Tokenizer::new(record.seq(), self.k);

            for kmer in tokenizer {
                self.inc(kmer);
            }
        }
    }

    /// Increase the counter of a kmer
    pub fn inc(&mut self, kmer: u64) {
        self.canonical_inc(cocktail::kmer::canonical(kmer, self.k));
    }

    fn canonical_inc(&mut self, canonical: u64) {
        let hash = (canonical >> 1) as usize;

        self.count[hash] = self.count[hash].saturating_add(1);
    }

    /// Get the counter of a kmer
    pub fn get(&self, kmer: u64) -> Count {
        let cano = cocktail::kmer::canonical(kmer, self.k);
        let hash = (cano >> 1) as usize;

        self.count[hash]
    }

    pub(crate) fn get_raw_count(&self) -> Box<[Count]> {
        self.count.clone()
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
        5, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 2,
        2, 1, 0, 1, 1, 1, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 0, 0, 0,
        1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 2, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 1, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 2, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 2,
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        1, 1, 0, 0, 2, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1,
        0, 1, 1, 2, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 2, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0,
        2, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 2, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 2, 1, 0, 0, 1, 1, 0,
        2, 0, 0, 0, 0, 0, 0, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1,
        0, 0, 0, 0, 2, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
        2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2,
        1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 2, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0,
        0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1,
    ];

    #[test]
    fn count_fasta() {
        let mut counter = crate::counter::Counter::new(5);

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
    fn count_fastq() {
        let mut counter = crate::counter::Counter::new(5);

        counter.count_fastq(FASTQ_FILE);

        let mut outfile = Vec::new();
        counter.serialize(&mut outfile).unwrap();

        assert_eq!(&outfile[..], &FASTQ_COUNT[..]);
    }

    lazy_static::lazy_static! {
    static ref COUNTER: crate::counter::Counter = {
            let mut counter = crate::counter::Counter::new(5);

            for i in 0..cocktail::kmer::get_kmer_space_size(5) {
        counter.inc(i);
            }

            counter
    };
    }

    const ALLKMERSEEONE: &[u8] = &[
        5, 0, 2, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
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
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
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
        assert_eq!(counter.get_raw_count(), COUNTER.get_raw_count());
    }
}
