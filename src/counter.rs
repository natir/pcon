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

use crate::bucketizer;
use cocktail;

pub trait Counter<CounterType, KmerType> {
    fn set_data(&mut self, data: Box<[CounterType]>);

    fn inc(&mut self, hash: KmerType);

    fn incs(&mut self, bucket: &bucketizer::NoTemporalArray);

    fn get(&self, hash: KmerType) -> CounterType;

    fn data(&self) -> &Box<[CounterType]>;

    fn nb_bit(&self) -> u8;

    fn clean(&mut self) -> ();
}

pub struct BasicCounter<T> {
    pub data: Box<[T]>,
}

macro_rules! impl_basiccounter {
    ($type:ty) => {
        impl BasicCounter<$type>
        where
            $type: std::clone::Clone,
        {
            pub fn new(k: u8) -> Self {
                BasicCounter {
                    data: vec![0; 1 << bucketizer::nb_bit(k)].into_boxed_slice(),
                }
            }
        }

        impl Counter<$type, u64> for BasicCounter<$type> {
            fn set_data(&mut self, data: Box<[$type]>) {
                self.data = data;
            }

            fn inc(&mut self, kmer: u64) {
                self.data[kmer as usize] = self.data[kmer as usize].saturating_add(1);
            }

            fn incs(&mut self, bucket: &bucketizer::NoTemporalArray) {
                for i in bucket {
                    self.inc(*i);
                }
            }

            fn get(&self, hash: u64) -> $type {
                self.data[hash as usize]
            }

            fn data(&self) -> &Box<[$type]> {
                &self.data
            }

            fn nb_bit(&self) -> u8 {
                std::mem::size_of::<$type>() as u8 * 8
            }

            fn clean(&mut self) -> () {
                for elt in self.data.iter_mut() {
                    *elt = 0;
                }
            }
        }
    };
}

impl_basiccounter!(u8);
impl_basiccounter!(u16);

pub struct ShortCounter {
    pub data: Box<[u8]>,
}

impl ShortCounter {
    pub fn new(k: u8) -> Self {
        ShortCounter {
            data: vec![0; 1usize << (bucketizer::nb_bit(k) - 1)].into_boxed_slice(),
        }
    }
}

impl Counter<u8, u64> for ShortCounter {
    fn set_data(&mut self, data: Box<[u8]>) {
        self.data = data;
    }

    fn inc(&mut self, hash: u64) {
        let key: usize = cocktail::kmer::remove_first_bit(hash) as usize;

        match cocktail::kmer::get_first_bit(hash) {
            true => match self.data[key] & 0b1111_0000 == 240 {
                true => (),
                false => self.data[key] += 16,
            },
            false => match self.data[key] & 0b0000_1111 == 15 {
                true => (),
                false => self.data[key] += 1,
            },
        }
    }

    fn incs(&mut self, bucket: &bucketizer::NoTemporalArray) {
        for i in bucket {
            self.inc(*i);
        }
    }

    fn get(&self, hash: u64) -> u8 {
        let key: usize = cocktail::kmer::remove_first_bit(hash) as usize;

        match cocktail::kmer::get_first_bit(hash) {
            true => self.data[key] >> 4,
            false => self.data[key] & 0b0000_1111,
        }
    }

    fn data(&self) -> &Box<[u8]> {
        &self.data
    }

    fn nb_bit(&self) -> u8 {
        4
    }

    fn clean(&mut self) {
        for elt in self.data.iter_mut() {
            *elt = 0;
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use cocktail;
    #[test]
    fn short_counter() -> () {
        let kmer1 = cocktail::kmer::cannonical(cocktail::kmer::seq2bit("ACTGG".as_bytes()), 5) >> 1;
        let kmer2 = cocktail::kmer::cannonical(cocktail::kmer::seq2bit("ACTGA".as_bytes()), 5) >> 1;

        let mut count = ShortCounter::new(5);

        let mut bucket = bucketizer::NoTemporalArray::new();
        for _ in 0..2 {
            bucket.push(kmer2);
        }
        count.incs(&bucket);
        assert_eq!(count.data[27], 2);

        bucket = bucketizer::NoTemporalArray::new();
        for _ in 0..2 {
            bucket.push(kmer1);
        }
        count.incs(&bucket);
        assert_eq!(count.data[27], 34);

        assert_eq!(count.get(kmer1), 2);
        assert_eq!(count.get(kmer2), 2);

        bucket = bucketizer::NoTemporalArray::new();
        for _ in 0..15 {
            bucket.push(kmer1);
        }
        count.incs(&bucket);
        assert_eq!(count.data[27], 242);

        bucket = bucketizer::NoTemporalArray::new();
        for _ in 0..15 {
            bucket.push(kmer2);
        }
        count.incs(&bucket);
        assert_eq!(count.data[27], 255);

        assert_eq!(count.get(kmer1), 15);
        assert_eq!(count.get(kmer2), 15);
    }

    #[test]
    fn basic_counter() -> () {
        let kmer1 = cocktail::kmer::cannonical(cocktail::kmer::seq2bit("ACTGG".as_bytes()), 5) >> 1;
        let kmer2 = cocktail::kmer::cannonical(cocktail::kmer::seq2bit("ACTGA".as_bytes()), 5) >> 1;

        let mut count = BasicCounter::<u8>::new(5);

        let mut bucket = bucketizer::NoTemporalArray::new();
        for _ in 0..2 {
            bucket.push(kmer2);
        }
        count.incs(&bucket);
        assert_eq!(count.data[54], 2);

        bucket = bucketizer::NoTemporalArray::new();
        for _ in 0..2 {
            bucket.push(kmer1);
        }
        count.incs(&bucket);
        assert_eq!(count.data[55], 2);

        assert_eq!(count.get(kmer1), 2);
        assert_eq!(count.get(kmer2), 2);

        bucket = bucketizer::NoTemporalArray::new();
        for _ in 0..255 {
            bucket.push(kmer1);
        }
        count.incs(&bucket);
        assert_eq!(count.data[54], 2);
        assert_eq!(count.data[55], 255);

        bucket = bucketizer::NoTemporalArray::new();
        for _ in 0..255 {
            bucket.push(kmer2);
        }
        count.incs(&bucket);
        assert_eq!(count.data[55], 255);

        assert_eq!(count.get(kmer1), 255);
        assert_eq!(count.get(kmer2), 255);
    }
}
