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
    fn inc(&mut self, val :&mut T);
}

struct IncUnsigned;

impl Inc<u8> for IncUnsigned {
    fn inc(&mut self, val: &mut u8) {
        *val = val.saturating_add(1);
    }
}

impl Inc<u16> for IncUnsigned {
    fn inc(&mut self, val: &mut u16) {
        *val = val.saturating_add(1);
    }
}

pub trait Counter<CounterType, BucketId, KmerType> {
    fn incs(&mut self, bucket_id: BucketId, bucket: Vec<KmerType>);
    fn get(&self, kmer: KmerType) -> CounterType;
}

pub struct MultiThreadCounter<T> {
    incrementor: IncUnsigned,
    data: Vec<std::sync::Arc<std::sync::Mutex<Vec<T>>>>,
    k: u8,
}

macro_rules! impl_multithreadcounter {
    ($type:ty) => (        
        impl MultiThreadCounter<$type> where $type: std::clone::Clone {
            pub fn new(k: u8) -> Self {
                let mut data = Vec::with_capacity(nb_bucket(k));

                for i in 0..nb_bucket(k) {
                    data.push(std::sync::Arc::new(std::sync::Mutex::new(vec![0; (1 << suffix_size(k))])));
                }
                
                MultiThreadCounter  {
                    incrementor: IncUnsigned {},
                    data: data,
                    k: k,
                }
            }
        }
    )
}

macro_rules! impl_counter_for_multithreadcounter {
    ($type:ty) => (
        impl Counter<$type, u64, u64> for MultiThreadCounter<$type> {
            fn incs(&mut self, bucket_id: u64, bucket: Vec<u64>) {
                let ptr = self.data[bucket_id as usize].clone();
                let mut data = ptr.lock().unwrap();

                for hash in bucket {
                    let suffix = get_suffix(self.k, hash);
                    self.incrementor.inc(&mut data[suffix as usize]);
                }
            }
            
            fn get(&self, kmer: u64) -> $type {
                let prefix = get_prefix(self.k, kmer);
                let suffix = get_suffix(self.k, kmer);
                self.data[prefix as usize].lock().unwrap()[suffix as usize]
            }
        }
    )
}

impl_multithreadcounter!(u8);
impl_multithreadcounter!(u16);

impl_counter_for_multithreadcounter!(u8);
impl_counter_for_multithreadcounter!(u16);

pub struct BasicCounter<T> {
    incrementor: IncUnsigned,
    pub data: Vec<T>,
}

macro_rules! impl_basiccounter {
    ($type:ty) => (        
        impl BasicCounter<$type> where $type: std::clone::Clone {
            pub fn new(k: u8) -> Self {
                BasicCounter  {
                    incrementor: IncUnsigned {},
                    data: vec![0; 1 << nb_bit(k)],
                }
            }
        }
    )
}

macro_rules! impl_counter_for_basiccounter {
    ($type:ty) => (
        impl Counter<$type, u64, u64> for BasicCounter<$type> {
            fn incs(&mut self, _bucket_id: u64, bucket: Vec<u64>) {
                for i in bucket {
                    self.incrementor.inc(&mut self.data[i as usize]);
                }
            }
            
            fn get(&self, kmer: u64) -> $type {
                self.data[kmer as usize]
            }
        }
    )
}

impl_basiccounter!(u8);
impl_basiccounter!(u16);
impl_counter_for_basiccounter!(u8);
impl_counter_for_basiccounter!(u16);

pub struct Bucketizer<'a, T> {
    pub counter: &'a mut Counter<T, u64, u64>,
    buckets: Vec<Vec<u64>>,
    k: u8,
    bucket_size: usize,
}

impl<'a, T> Bucketizer<'a, T> {
    pub fn new(counter: &'a mut Counter<T, u64, u64>, k: u8) -> Self {
        let bucket_size: usize = 1024;

        Bucketizer {
            counter: counter,
            buckets: vec![Vec::with_capacity(bucket_size); nb_bucket(k)],
            k: k,
            bucket_size: bucket_size,
        }
    }

    pub fn add_bit(&mut self, hash: u64) -> () {
        let prefix: usize = get_prefix(self.k, hash);
        
        self.buckets[prefix].push(hash);

        if self.buckets[prefix].len() == self.bucket_size {
            self.clean_bucket(prefix, self.bucket_size);
        }
    }
    
    pub fn clean_all_buckets(&mut self) -> () {
        for prefix in 0..nb_bucket(self.k) {
            self.clean_bucket(prefix, 0);
        }
    }

    fn clean_bucket(&mut self, prefix: usize, new_size: usize) -> () {
        self.counter.incs(prefix as u64,
            std::mem::replace(&mut self.buckets[prefix], Vec::with_capacity(new_size))
        );
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
}
