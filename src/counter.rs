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


pub trait Counter {
    fn inc(&mut self, pos: usize) -> ();
    fn get(&self, pos: usize) -> u8;
}

pub struct VecCounter {
    inner: Vec<u8>,
}

impl VecCounter {
    pub fn new(k: u8) -> Self {
        VecCounter {
            inner: vec![0; 1 << (k * 2 - 1)]
        }
    }
}

impl Counter for VecCounter {
    fn inc(&mut self, pos: usize) -> () {
        self.inner[pos] = self.inner[pos].saturating_add(1);
    }

    fn get(&self, pos: usize) -> u8 {
        self.inner[pos]
    }
}

pub struct Bucketizer<'a, C: Counter> {
    pub counter: &'a mut C,
    buckets: Vec<Vec<u64>>,
    k: u8,
    bucket_size: usize,
}

impl<'a, C: Counter> Bucketizer<'a, C> {
    pub fn new(counter: &'a mut C, k: u8) -> Self {
        let bucket_size: usize = 1024;

        Bucketizer {
            counter: counter,
            buckets: vec![Vec::with_capacity(bucket_size); nb_bucket(k)],
            k: k,
            bucket_size: bucket_size,
        }
    }

    pub fn add_bit(&mut self, hash: u64) -> () {
        let prefix: usize = (self.mask_prefix() & hash as usize) >> (self.nb_bit() - self.mask_size());
        
        self.buckets[prefix].push(hash);

        if self.buckets[prefix].len() == self.bucket_size {
            self.clean_bucket(prefix, self.bucket_size);
        }
    }

    pub fn clean_all_buckets(&mut self) -> () {
        for prefix in 0..self.nb_bucket() {
            self.clean_bucket(prefix, 0);
        }
    }

    fn clean_bucket(&mut self, prefix: usize, new_size: usize) -> () {
        for hash in std::mem::replace(&mut self.buckets[prefix], Vec::with_capacity(new_size)) {
            self.counter.inc(hash as usize);
        }
    }

    fn nb_bit(&self) -> usize {
        nb_bit(self.k)
    }

    fn mask_size(&self) -> usize {
        match mask_size(self.k) > 12 { true => 12, false => mask_size(self.k)}
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
    nb_bit(k) / 2
}

fn nb_bucket(k: u8) -> usize {
    1 << mask_size(k)
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
    }

    #[test]
    fn mask_prefix_() -> () {
        assert_eq!(mask_prefix(5), 0b111100000);
    }

}
