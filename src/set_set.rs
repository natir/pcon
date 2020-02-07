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

/* standard use */

/* crates use */
use bv::BitVec;

/* project use */
use crate::counter;

pub struct SetOfSet {
    data: std::collections::HashMap<u64, std::collections::HashSet<u64>>,
    k: u8,
    m: u8,
    mask: u64,
}

impl SetOfSet {
    pub fn new(counter: &Box<dyn counter::Counter<u8, u64>>, k: u8, m: u8, abundance_min: u8) -> SetOfSet {
	let mut set_set = std::collections::HashMap::new();
	let mask = (1 << (m * 2 + 1)) - 1;
	
        for kmer in 0..cocktail::kmer::get_kmer_space_size(k) {
	    let hash = cocktail::kmer::remove_first_bit(kmer);
            if counter.get(hash) >= abundance_min {
		let mini = SetOfSet::get_minimizer(kmer, mask, k, m);
		set_set.entry(mini).or_insert(std::collections::HashSet::new()).insert(kmer);
            }
        }

	SetOfSet {
	    data: set_set,
	    k: k,
	    m: m,
	    mask: mask,
	}
    }

    pub fn from_bitfield(exist: &bv::BitVec<u8>, k: u8, m: u8) -> SetOfSet {
	let mask = (1 << (m * 2 + 1)) - 1;
	let mut set_set = std::collections::HashMap::new();
	
	for val in 0..cocktail::kmer::get_kmer_space_size(k) {
	    let kmer = cocktail::kmer::revcomp(val, k);
	    let hash = cocktail::kmer::remove_first_bit(kmer);
	    
	    if exist.get(hash) {
		let mini = SetOfSet::get_minimizer(kmer, mask, k, m);
		set_set.entry(mini).or_insert(std::collections::HashSet::new()).insert(kmer);
	    }
	}

	SetOfSet {
	    data: set_set,
	    k: k,
	    m: m,
	    mask: mask,
	}
    }

    fn get_minimizer(kmer: u64, mask: u64, k: u8, m: u8) -> u64 {
	let mut mini = xorshift64(kmer & mask);
	for shift in 0..(k - m + 1) {
	    let sub = xorshift64((kmer >> (shift * 2)) & mask);
	    if sub < mini {
		mini = sub;
	    }
	}
	mini
    }

    pub fn contains(&self, kmer: u64) -> bool {
	let mini = SetOfSet::get_minimizer(kmer, self.mask, self.k, self.m);
	if let Some(set) = self.data.get(&mini) {
	    set.contains(&kmer)
	} else {
	    false
	}
    }
}


fn xorshift64(mut x: u64) -> u64
{
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
 
    x
}
