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

trait Inc<T> {
    fn inc(&mut self, val: &mut T);
}

struct IncUnsigned;

macro_rules! impl_incunsigned {
    ($type:ty) => {
        impl Inc<$type> for IncUnsigned {
            fn inc(&mut self, val: &mut $type) {
                *val = val.saturating_add(1);
            }
        }
    };
}

impl_incunsigned!(u8);
impl_incunsigned!(u16);

pub trait Counter<CounterType, KmerType> {
    fn incs(&mut self, bucket: &NoTemporalArray);
    fn get(&self, kmer: KmerType) -> CounterType;
}

pub struct BasicCounter<T> {
    incrementor: IncUnsigned,
    pub data: Vec<T>,
}

macro_rules! impl_basiccounter {
    ($type:ty) => {
        impl BasicCounter<$type>
        where
            $type: std::clone::Clone,
        {
            pub fn new(k: u8) -> Self {
                BasicCounter {
                    incrementor: IncUnsigned {},
                    data: vec![0; 1 << nb_bit(k)],
                }
            }
        }

        impl Counter<$type, u64> for BasicCounter<$type> {
            fn incs(&mut self, bucket: &NoTemporalArray) {
                for i in bucket {
                    self.incrementor.inc(&mut self.data[*i as usize]);
                }
            }

            fn get(&self, kmer: u64) -> $type {
                self.data[kmer as usize]
            }
        }
    };
}

impl_basiccounter!(u8);
impl_basiccounter!(u16);

pub struct ShortCounter {
    pub data: Vec<u8>,
}

impl ShortCounter {
    pub fn new(k: u8) -> Self {
        ShortCounter {
            data: vec![0; (1 << nb_bit(k)) / 2],
        }
    }
}

impl Counter<u8, u64> for ShortCounter {
    fn incs(&mut self, bucket: &NoTemporalArray) {
        for i in bucket {
            let weight = i & 1 != 0;
            let key: usize = *i as usize >> 1;
            
            match weight {
                true => match self.data[key] & 0b11110000 == 240 {
                    true => (),
                    false => self.data[key] += 16,
                },
                false => match self.data[key] & 0b00001111 == 15 {
                    true => (),
                    false => self.data[key] += 1,
                }
            }
        }
    }

    fn get(&self, kmer: u64) -> u8 {
        let weight = kmer & 1 != 0;
        let key: usize = kmer as usize >> 1;
            
        return match weight {
            true => self.data[key] >> 4,
            false => self.data[key] & 0b00001111,
        };
    }
}

const BUCKET_SIZE: usize = 1 << 12;

pub struct NoTemporalArray {
    pos: usize,
    data: Box<[u64; BUCKET_SIZE]>,
}

impl NoTemporalArray {
    fn new() -> Self {
        NoTemporalArray {
            data: Box::new(unsafe { std::mem::MaybeUninit::uninit().assume_init() }),
            pos: 0,
        }
    }

    fn push(&mut self, val: u64) -> () {
        unsafe {
            core::arch::x86_64::_mm_stream_pi(
                self.data.as_mut_ptr().add(self.pos) as *mut std::arch::x86_64::__m64,
                std::mem::transmute(val),
            )
        }
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

pub struct Bucketizer<'a, T> {
    pub counter: &'a mut Counter<T, u64>,
    buckets: Vec<NoTemporalArray>,
    k: u8,
    bucket_size: usize,
}

impl<'a, T> Bucketizer<'a, T> {
    pub fn new(counter: &'a mut Counter<T, u64>, k: u8) -> Self {
        Bucketizer {
            counter: counter,
            buckets: (0..nb_bucket(k)).map(|_| NoTemporalArray::new()).collect(),
            k: k,
            bucket_size: BUCKET_SIZE,
        }
    }

    #[inline(always)]
    pub fn add_bit(&mut self, hash: u64) -> () {
        let prefix: usize = get_prefix(self.k, hash);

        self.buckets[prefix].push(hash);

        if self.buckets[prefix].len() == self.bucket_size {
            self.clean_bucket(prefix);
        }
    }

    pub fn clean_all_buckets(&mut self) -> () {
        for prefix in 0..nb_bucket(self.k) {
            self.clean_bucket(prefix);
        }
    }

    fn clean_bucket(&mut self, prefix: usize) -> () {
        self.counter.incs(&self.buckets[prefix]);
        self.buckets[prefix].pos = 0;
    }
}

fn get_prefix(k: u8, hash: u64) -> usize {
    (mask_prefix(k) & hash as usize) >> suffix_size(k)
}

fn get_suffix(k: u8, hash: u64) -> usize {
    let mov = 64 - suffix_size(k);
    ((hash as usize) << mov) >> mov
}

fn nb_bit(k: u8) -> usize {
    (k * 2 - 1) as usize
}

fn prefix_size(k: u8) -> usize {
    nb_bit(k) / 2
}

fn suffix_size(k: u8) -> usize {
    nb_bit(k) - prefix_size(k)
}

fn nb_bucket(k: u8) -> usize {
    1 << prefix_size(k)
}

fn mask_prefix(k: u8) -> usize {
    ((1 << prefix_size(k)) - 1) << suffix_size(k)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn nb_bit_() -> () {
        assert_eq!(nb_bit(5), 9);
    }

    #[test]
    fn prefix_size_() -> () {
        assert_eq!(prefix_size(5), 4);
    }

    #[test]
    fn suffix_size_() -> () {
        assert_eq!(suffix_size(5), 5);
    }

    #[test]
    fn nb_bucket_() -> () {
        assert_eq!(nb_bucket(5), 16);
        assert_eq!(nb_bucket(7), 64);
    }

    #[test]
    fn mask_prefix_() -> () {
        assert_eq!(mask_prefix(5), 0b111100000);
    }

    #[test]
    fn get_prefix_() -> () {
        let kmer = 0b111111111;
        assert_eq!(get_prefix(5, kmer), 0b1111);
    }

    #[test]
    fn get_suffix_() -> () {
        let kmer: u64 = 0b111111111;

        assert_eq!(get_suffix(5, kmer), 0b11111);
    }

    use crate::convert;
    #[test]
    fn short_counter() -> () {
        let kmer1 = convert::cannonical(convert::seq2bit("ACTGG".as_bytes()), 5) >> 1;
        let kmer2 = convert::cannonical(convert::seq2bit("ACTGA".as_bytes()), 5) >> 1;

        let mut count = ShortCounter::new(5);

        let mut bucket = NoTemporalArray::new();
        for _ in 0..2 {
            bucket.push(kmer2);
        }
        count.incs(&bucket);
        assert_eq!(count.data[27], 2);

        bucket = NoTemporalArray::new();
        for _ in 0..2 {
            bucket.push(kmer1);
        }
        count.incs(&bucket);
        assert_eq!(count.data[27], 34);

        assert_eq!(count.get(kmer1), 2);
        assert_eq!(count.get(kmer2), 2);

        bucket = NoTemporalArray::new();
        for _ in 0..15 {
            bucket.push(kmer1);
        }
        count.incs(&bucket);
        assert_eq!(count.data[27], 242);

        bucket = NoTemporalArray::new();
        for _ in 0..15 {
            bucket.push(kmer2);
        }
        count.incs(&bucket);
        assert_eq!(count.data[27], 255);

        assert_eq!(count.get(kmer1), 15);
        assert_eq!(count.get(kmer2), 15);
    }
}
