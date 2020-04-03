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

use bitvec::prelude::*;

use crate::*;

use crate::error::IO::*;
use crate::error::*;

//#[repr(C)]
#[derive(serde::Serialize, serde::Deserialize)]
pub struct Solid {
    pub k: u8,
    solid: BitVec,
}

impl Solid {
    pub fn new(k: u8) -> Self {
	Self {
	    k,
	    solid: bitvec![Lsb0; 0; cocktail::kmer::get_hash_space_size(k) as usize]
	}
    }
    
    pub fn from_counter(counter: &counter::Counter, abundance: u8) -> Self {
	let count = counter.get_raw_count();

        let mut solid = bitvec![Lsb0; 0; count.len()];

        for (index, val) in count.iter().enumerate() {
            if val > &abundance {
                if let Some(mut v) = solid.get_mut(index) {
                    *v = true;
                }
            }
        }

	Self {
	    k: counter.k,
	    solid,
	}
    }
    
    pub fn set(&mut self, kmer: u64, value: bool) {
	let cano = cocktail::kmer::cannonical(kmer, self.k);
        let hash = (cano >> 1) as usize;

        if let Some(mut v) = self.solid.get_mut(hash) {
            *v = value;
        }
    }

    pub fn get(&self, kmer: u64) -> bool {
        let cano = cocktail::kmer::cannonical(kmer, self.k);
        let hash = (cano >> 1) as usize;
	
        self.solid[hash]
    }

    pub fn serialize<W>(&self, writer: W) -> Result<()>
    where
        W: std::io::Write,
    {
        bincode::serialize_into(writer, &self)
	    .with_context(|| {
		Error::IO(ErrorDurringWrite)
            })
	    .with_context(|| {
		anyhow!("Error durring serialize counter")
	    })	    
    }

    pub fn deserialize<R>(reader: R) -> Result<Self>
    where
        R: std::io::Read,
    {
        bincode::deserialize_from(reader)
	    .with_context(|| {
		Error::IO(ErrorDurringRead)
	    })
	    .with_context(|| {
		anyhow!("Error durring deserialize counter")
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

    lazy_static::lazy_static! {
	static ref SOLID: std::sync::Mutex<crate::solid::Solid> = {
	    let mut counter = crate::counter::Counter::new(5);

	    counter.count_fasta(FASTA_FILE);

            std::sync::Mutex::new(crate::solid::Solid::from_counter(&counter, 0))
	};
    }
    
    #[test]
    fn presence() {
	let solid = SOLID.lock().unwrap();

	assert_eq!(solid.get(cocktail::kmer::seq2bit(b"GTTCT")), true);
	assert_eq!(solid.get(cocktail::kmer::seq2bit(b"AAATG")), true);
	assert_eq!(solid.get(cocktail::kmer::seq2bit(b"AGGAT")), true);
	assert_eq!(solid.get(cocktail::kmer::seq2bit(b"CTCAG")), true);
    }

    #[test]
    fn set_value() {
	let mut solid = SOLID.lock().unwrap();

	solid.set(cocktail::kmer::seq2bit(b"GTTCT"), false);
	solid.set(cocktail::kmer::seq2bit(b"AAATG"), false);
	solid.set(cocktail::kmer::seq2bit(b"AGGAT"), false);
	solid.set(cocktail::kmer::seq2bit(b"CTCAG"), false);	
	
	assert_eq!(solid.get(cocktail::kmer::seq2bit(b"GTTCT")), false);
	assert_eq!(solid.get(cocktail::kmer::seq2bit(b"AAATG")), false);
	assert_eq!(solid.get(cocktail::kmer::seq2bit(b"AGGAT")), false);
	assert_eq!(solid.get(cocktail::kmer::seq2bit(b"CTCAG")), false);
    }

    const FASTA_SOLID: &[u8] = &[5, 0, 0, 2, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 112, 64, 113, 143, 130, 8, 128, 4, 6, 60, 214, 0, 243, 8, 193, 1, 30, 4, 34, 97, 4, 70, 192, 12, 16, 144, 133, 38, 192, 41, 1, 4, 218, 179, 140, 0, 0, 140, 242, 35, 90, 56, 205, 179, 64, 3, 25, 20, 226, 0, 32, 76, 1, 134, 48, 64, 7, 0, 200, 144, 98, 131, 2, 203];
    
    #[test]
    fn serialize() {
	let mut outfile = Vec::new();

	let solid = SOLID.lock().unwrap();
	solid.serialize(&mut outfile).unwrap();

	assert_eq!(&outfile[..], &FASTA_SOLID[..]);
    }

    #[test]
    fn deserialize() {
	let solid = crate::solid::Solid::deserialize(&FASTA_SOLID[..]).unwrap();
	
	assert_eq!(solid.k, 5);
	assert_eq!(solid.get(cocktail::kmer::seq2bit(b"GTTCT")), true);
	assert_eq!(solid.get(cocktail::kmer::seq2bit(b"AAATG")), true);
	assert_eq!(solid.get(cocktail::kmer::seq2bit(b"AGGAT")), true);
	assert_eq!(solid.get(cocktail::kmer::seq2bit(b"CTCAG")), true);
    }
}
