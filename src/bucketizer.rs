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

use crate::counter;
use cocktail;

const BUCKET_SIZE: usize = 1 << 10;

pub struct NoTemporalArray {
    pos: usize,
    data: Box<[u64; BUCKET_SIZE]>,
}

impl NoTemporalArray {
    pub fn new() -> Self {
        NoTemporalArray {
            data: Box::new(unsafe { std::mem::MaybeUninit::uninit().assume_init() }),
            pos: 0,
        }
    }

    pub fn push(&mut self, val: u64) {
        self.data[self.pos] = val;
        self.pos += 1
    }

    fn len(&self) -> usize {
        self.pos
    }
}

impl<'a> std::iter::IntoIterator for &'a NoTemporalArray {
    type Item = &'a u64;
    type IntoIter = std::slice::Iter<'a, u64>;

    fn into_iter(self) -> std::slice::Iter<'a, u64> {
        (&(self.data)[..self.pos]).iter()
    }
}

pub trait Bucket<T> {
    fn add_kmer(&mut self, kmer: u64) -> ();
    fn clean_all_buckets(&mut self) -> ();
    fn clean_bucket(&mut self, prefix: usize) -> ();
    fn counter(&self) -> &Box<dyn counter::Counter<T, u64>>;
}

pub struct Prefix<T> {
    counter: Box<dyn counter::Counter<T, u64>>,
    buckets: Vec<NoTemporalArray>,
    k: u8,
    bucket_size: usize,
    prefix_mask: u64,
}

impl<T> Prefix<T> {
    pub fn new(counter: Box<dyn counter::Counter<T, u64>>, k: u8) -> Self {
        Prefix {
            counter,
            buckets: (0..nb_bucket(k)).map(|_| NoTemporalArray::new()).collect(),
            k,
            bucket_size: BUCKET_SIZE,
            prefix_mask: mask_prefix(k) as u64,
        }
    }

    fn get_prefix(&self, hash: u64) -> usize {
        ((self.prefix_mask & hash) >> nb_bit(self.k)) as usize
    }
}

impl<T> Bucket<T> for Prefix<T> {
    #[inline(always)]
    fn add_kmer(&mut self, kmer: u64) {
        let hash = cocktail::kmer::remove_first_bit(kmer);
        let prefix: usize = self.get_prefix(hash);

        self.buckets[prefix].push(hash);

        if self.buckets[prefix].len() == self.bucket_size {
            self.clean_bucket(prefix);
        }
    }

    fn clean_all_buckets(&mut self) {
        for prefix in 0..nb_bucket(self.k) {
            self.clean_bucket(prefix);
        }
    }

    fn clean_bucket(&mut self, prefix: usize) {
        self.counter.incs(&self.buckets[prefix]);
        self.buckets[prefix].pos = 0;
    }

    fn counter(&self) -> &Box<dyn counter::Counter<T, u64>> {
        &self.counter
    }
}

fn nb_bucket(k: u8) -> usize {
    1 << prefix_size(k)
}

pub fn nb_bit(k: u8) -> usize {
    (k * 2 - 1) as usize
}

fn prefix_size(k: u8) -> usize {
    nb_bit(k) / 4
}

fn suffix_size(k: u8) -> usize {
    nb_bit(k) - prefix_size(k)
}

fn mask_prefix(k: u8) -> usize {
    ((1 << prefix_size(k)) - 1) << suffix_size(k)
}

#[cfg(test)]
mod test {
    use super::*;

    mod bucket {
        use super::*;

        #[test]
        fn no_temporal_array() {
            let mut array = NoTemporalArray::new();

            assert_eq!(array.len(), 0);

            for i in 0..BUCKET_SIZE {
                array.push(i as u64);
            }

            assert_eq!(array.len(), BUCKET_SIZE);
        }

        #[test]
        #[should_panic]
        fn no_temporal_array_overflow() {
            let mut array = NoTemporalArray::new();

            for i in 0..(BUCKET_SIZE + 1) {
                array.push(i as u64);
            }
        }

        #[test]
        fn add_kmer() {
            let mut buckets = Prefix::new(Box::new(counter::BasicCounter::<u16>::new(5)), 5);

            buckets.add_kmer(0);

            assert_eq!(buckets.counter().get(0), 0); // with buckets value isn't directly add in counter

            for _ in 0..BUCKET_SIZE {
                buckets.add_kmer(0);
            }

            assert_eq!(buckets.counter().get(0), 1024); // we reach max buckets size the buckket is clean

            buckets.clean_all_buckets();

            assert_eq!(buckets.counter().get(0), 1025); // we clean all buckets and last kmer is add in counter
        }
    }

    #[test]
    fn nb_bit_() -> () {
        assert_eq!(nb_bit(5), 9);
    }

    #[test]
    fn prefix_size_() -> () {
        assert_eq!(prefix_size(5), 2);
    }

    #[test]
    fn suffix_size_() -> () {
        assert_eq!(suffix_size(5), 7);
    }

    #[test]
    fn nb_bucket_() -> () {
        assert_eq!(nb_bucket(5), 4);
        assert_eq!(nb_bucket(7), 8);
    }

    #[test]
    fn mask_prefix_() -> () {
        assert_eq!(mask_prefix(5), 0b110000000);
    }
}
