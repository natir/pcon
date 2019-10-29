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

use crate::convert;
use crate::bucketizer;

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
    fn inc(&mut self, kmer: KmerType);

    fn incs(&mut self, bucket: &bucketizer::NoTemporalArray);
    
    fn get(&self, kmer: KmerType) -> CounterType;

    fn data(&self) -> &Box<[CounterType]>;

    fn nb_bit(&self) -> u8;
}

pub struct BasicCounter<T> {
    incrementor: IncUnsigned,
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
                    incrementor: IncUnsigned {},
                    data: vec![0; 1 << bucketizer::nb_bit(k)].into_boxed_slice(),
                }
            }
        }

        impl Counter<$type, u64> for BasicCounter<$type> {
            fn inc(&mut self, kmer: u64) {
                self.incrementor.inc(&mut self.data[kmer as usize]);
            }        

            fn incs(&mut self, bucket: &bucketizer::NoTemporalArray){
                for i in bucket {
                    self.inc(*i);
                }
            }

            fn get(&self, kmer: u64) -> $type {
                self.data[kmer as usize]
            }

            fn data(&self) -> &Box<[$type]> {
                &self.data
            }

            fn nb_bit(&self) -> u8 {
                std::mem::size_of::<$type>() as u8 * 8
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
    fn inc (&mut self, kmer: u64) {
        let key: usize = convert::remove_first_bit(kmer) as usize;
        
        match convert::get_first_bit(kmer) {
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

    fn incs(&mut self, bucket: &bucketizer::NoTemporalArray){
        for i in bucket {
            self.inc(*i);
        }
    }
    
    fn get(&self, kmer: u64) -> u8 {
        let key: usize = kmer as usize >> 1;
            
        return match convert::get_first_bit(kmer) {
            true => self.data[key] >> 4,
            false => self.data[key] & 0b00001111,
        };
    }

    fn data(&self) -> &Box<[u8]> {
        &self.data
    }

    fn nb_bit(&self) -> u8 {
        4
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::convert;
    #[test]
    fn short_counter() -> () {
        let kmer1 = convert::cannonical(convert::seq2bit("ACTGG".as_bytes()), 5) >> 1;
        let kmer2 = convert::cannonical(convert::seq2bit("ACTGA".as_bytes()), 5) >> 1;

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
}
