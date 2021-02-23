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

/* crate use */
use anyhow::{anyhow, Context, Result};
use bitvec::prelude::*;
use byteorder::{ReadBytesExt, WriteBytesExt};

/* local use */
use crate::error::IO::*;
use crate::error::*;
use crate::*;

/// A struct to store if a kmer is Solid or not. Only kmer with abundance upper than a threshold is solid
pub struct Solid {
    pub k: u8,
    solid: BitBox<Lsb0, u8>,
}

impl Solid {
    /// Create a new Solid for kmer size equal to `k`
    pub fn new(k: u8) -> Self {
        Self {
            k,
            solid: bitbox![Lsb0, u8; 0; cocktail::kmer::get_hash_space_size(k) as usize],
        }
    }

    /// Create a new Solid with count in `counter` only kmer upper than `abundance` are solid
    pub fn from_counter(counter: &counter::Counter, abundance: counter::Count) -> Self {
        let counts = unsafe {
            &(*(counter.get_raw_count() as *const [counter::AtoCount] as *const [counter::Count]))
        };

        let mut solid = bitbox![Lsb0, u8; 0; counts.len()];

        unsafe {
            for (index, count) in (*(counter.get_raw_count() as *const [counter::AtoCount]
                as *const [counter::Count]))
                .iter()
                .enumerate()
            {
                if *count > abundance {
                    solid.set(index, true);
                }
            }
        }

        Self {
            k: counter.k,
            solid,
        }
    }

    /// Solidity status of `kmer` is set to `value`
    pub fn set(&mut self, kmer: u64, value: bool) {
        self.set_canonic(cocktail::kmer::canonical(kmer, self.k), value);
    }

    /// Solidity status of a canonical`kmer` is set to `value`
    pub fn set_canonic(&mut self, canonical: u64, value: bool) {
        let hash = (canonical >> 1) as usize;

        if let Some(mut v) = self.solid.get_mut(hash) {
            *v = value;
        }
    }

    /// Get the solidity status of `kmer`
    pub fn get(&self, kmer: u64) -> bool {
        self.get_canonic(cocktail::kmer::canonical(kmer, self.k))
    }

    /// Get the solidity status of a canonical `kmer`
    pub fn get_canonic(&self, canonical: u64) -> bool {
        let hash = (canonical >> 1) as usize;

        self.solid[hash]
    }

    #[allow(dead_code)]
    pub(crate) fn get_raw_solid(&self) -> &BitBox<Lsb0, u8> {
        &self.solid
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

        writer
            .write_u8(self.k)
            .with_context(|| Error::IO(ErrorDurringWrite))
            .with_context(|| anyhow!("Error durring serialize solid"))?;

        writer
            .write_all(self.solid.as_slice())
            .with_context(|| Error::IO(ErrorDurringWrite))
            .with_context(|| anyhow!("Error durring serialize solid"))?;

        Ok(())
    }

    /// Deserialize counter for given [std::io::Read]
    pub fn deserialize<R>(r: R) -> Result<Self>
    where
        R: std::io::Read,
    {
        let mut reader = niffler::get_reader(Box::new(r))?.0;

        let k = reader
            .read_u8()
            .with_context(|| Error::IO(ErrorDurringRead))
            .with_context(|| anyhow!("Error durring deserialize solid"))?;

        // >> 3 <-> divide by 8
        let mut tmp =
            vec![0u8; (cocktail::kmer::get_hash_space_size(k) >> 3) as usize].into_boxed_slice();

        reader
            .read_exact(&mut tmp)
            .with_context(|| Error::IO(ErrorDurringRead))
            .with_context(|| anyhow!("Error durring deserialize solid"))?;

        Ok(Self {
            k,
            solid: BitBox::from_boxed_slice(tmp),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const FASTA_FILE: &[u8] = b">random_seq 0
GTTCTGCAAATTAGAACAGACAATACACTGGCAGGCGTTGCGTTGGGGGAGATCTTCCGTAACGAGCCGGCATTTGTAAGAAAGAGATTTCGAGTAAATG
>random_seq 1
AGGATAGAAGCTTAAGTACAAGATAATTCCCATAGAGGAAGGGTGGTATTACAGTGCCGCCTGTTGAAAGCCCCAATCCCGCTTCAATTGTTGAGCTCAG
";

    fn get_solid() -> solid::Solid {
        let mut counter = crate::counter::Counter::new(5);

        counter.count_fasta(FASTA_FILE, 1);

        solid::Solid::from_counter(&counter, 0)
    }

    const SOLID: &[u8] = &[
        112, 64, 113, 143, 130, 8, 128, 4, 6, 60, 214, 0, 243, 8, 193, 1, 30, 4, 34, 97, 4, 70,
        192, 12, 16, 144, 133, 38, 192, 41, 1, 4, 218, 179, 140, 0, 0, 140, 242, 35, 90, 56, 205,
        179, 64, 3, 25, 20, 226, 0, 32, 76, 1, 134, 48, 64, 7, 0, 200, 144, 98, 131, 2, 203,
    ];

    #[test]
    fn presence() {
        let solid = get_solid();

        assert_eq!(solid.get_raw_solid().as_slice(), SOLID);
    }

    const SOLID_SET: &[u8] = &[
        112, 64, 113, 143, 130, 8, 128, 4, 6, 52, 214, 0, 243, 8, 193, 1, 30, 4, 2, 97, 4, 70, 192,
        12, 16, 144, 133, 36, 192, 41, 1, 4, 218, 179, 140, 0, 0, 140, 242, 35, 90, 56, 205, 179,
        64, 3, 25, 20, 226, 0, 32, 76, 1, 134, 48, 64, 7, 0, 192, 144, 98, 131, 2, 203,
    ];

    #[test]
    fn set_value() {
        let mut solid = get_solid();

        solid.set(cocktail::kmer::seq2bit(b"GTTCT"), false);
        solid.set(cocktail::kmer::seq2bit(b"AAATG"), false);
        solid.set(cocktail::kmer::seq2bit(b"AGGAT"), false);
        solid.set(cocktail::kmer::seq2bit(b"CTCAG"), false);

        assert_eq!(solid.get_raw_solid().as_slice(), SOLID_SET);
    }

    const FASTA_SOLID: &[u8] = &[
        31, 139, 8, 0, 0, 0, 0, 0, 4, 255, 1, 65, 0, 190, 255, 5, 112, 64, 113, 143, 130, 8, 128,
        4, 6, 60, 214, 0, 243, 8, 193, 1, 30, 4, 34, 97, 4, 70, 192, 12, 16, 144, 133, 38, 192, 41,
        1, 4, 218, 179, 140, 0, 0, 140, 242, 35, 90, 56, 205, 179, 64, 3, 25, 20, 226, 0, 32, 76,
        1, 134, 48, 64, 7, 0, 200, 144, 98, 131, 2, 203, 186, 210, 139, 120, 65, 0, 0, 0,
    ];

    #[test]
    fn serialize() {
        let mut outfile = Vec::new();

        let solid = get_solid();
        solid.serialize(&mut outfile).unwrap();

        assert_eq!(&outfile[..], &FASTA_SOLID[..]);
    }

    #[test]
    fn deserialize() {
        let solid = crate::solid::Solid::deserialize(&FASTA_SOLID[..]).unwrap();

        assert_eq!(solid.k, 5);
        assert_eq!(solid.get_raw_solid().as_slice(), SOLID);
    }
}
