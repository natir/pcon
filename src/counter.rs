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
    fn incs(&mut self, bucket_id: BucketId, bucket: &NoTemporalArray);
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
    };
}

macro_rules! impl_counter_for_basiccounter {
    ($type:ty) => {
        impl Counter<$type, u64, u64> for BasicCounter<$type> {
            fn incs(&mut self, _bucket_id: u64, bucket: &NoTemporalArray) {
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
impl_counter_for_basiccounter!(u8);
impl_counter_for_basiccounter!(u16);

const BUCKET_SIZE: usize = (2 << 20) / 8;

//#[derive(Copy, Clone)]
pub struct NoTemporalArray {
    pos: usize,
    data: Box<[u64; BUCKET_SIZE]>,
}

impl NoTemporalArray {
    fn new() -> Self {
        NoTemporalArray {
            data: Box::new(unsafe { std::mem::uninitialized() }),
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
    pub counter: &'a mut Counter<T, u64, u64>,
    buckets: Vec<NoTemporalArray>,
    k: u8,
    bucket_size: usize,
}

impl<'a, T> Bucketizer<'a, T> {
    pub fn new(counter: &'a mut Counter<T, u64, u64>, k: u8) -> Self {
        Bucketizer {
            counter: counter,
            buckets: (0..nb_bucket(k)).map(|x| NoTemporalArray::new()).collect(),
            k: k,
            bucket_size: BUCKET_SIZE,
        }
    }

    #[inline(always)]
    pub fn add_bit(&mut self, hash: u64) -> () {
        let prefix: usize =
            (self.mask_prefix() & hash as usize) >> (self.nb_bit() - self.mask_size());

        self.buckets[prefix].push(hash);

        if self.buckets[prefix].len() == self.bucket_size {
            self.clean_bucket(prefix, self.bucket_size);
        }
    }

    pub fn clean_all_buckets(&mut self) -> () {
        //         println!("cleanall");
        for prefix in 0..self.nb_bucket() {
            self.clean_bucket(prefix, 0);
        }
    }

    //#[inline(never)]
    fn clean_bucket(&mut self, prefix: usize, new_size: usize) -> () {
        //         println!("clean: {} {}", prefix, self.buckets[prefix].pos);
        //         self.counter.incs(prefix as u64,
        //                           std::mem::replace(&mut self.buckets[prefix], NoTemporalArray::new())
        //         );
        self.counter.incs(prefix as u64, &self.buckets[prefix]);
        self.buckets[prefix].pos = 0;
    }

    fn nb_bit(&self) -> usize {
        nb_bit(self.k)
    }

    fn mask_size(&self) -> usize {
        match mask_size(self.k) > 12 {
            true => 12,
            false => mask_size(self.k),
        }
    }

    fn nb_bucket(&self) -> usize {
        nb_bucket(self.k)
    }

    fn mask_prefix(&self) -> usize {
        mask_prefix(self.k)
    }
}

fn nb_bit(k: u8) -> usize {
    (k * 2 - 1) as usize
}

fn mask_size(k: u8) -> usize {
    //nb_bit(k) / 2
    nb_bit(k) - 16
}

fn nb_bucket(k: u8) -> usize {
    let n = 1 << mask_size(k);
    println!("nbucket: {} {}", n, 2 * (1 << nb_bit(k)) / n);
    n
}

fn mask_prefix(k: u8) -> usize {
    ((1 << mask_size(k)) - 1) << (nb_bit(k) - mask_size(k))
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn nb_bit_() -> () {
        assert_eq!(nb_bit(5), 9);
    }

    #[test]
    fn mask_size_() -> () {
        assert_eq!(mask_size(5), 4);
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
}
