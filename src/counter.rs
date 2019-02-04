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

// use std::thread;
// use std::sync::{Arc, Mutex};

pub struct Counter {
    pub counts: Vec<u8>,
    buckets: Vec<Vec<u64>>,
    mask_size: u8,
    bucket_size: usize,
    nb_bit: u8,
}

impl Counter {
    pub fn new(k: u8) -> Self {
        let nb_bit = k * 2 - 1;
        let mask_size = nb_bit / 2;
        let nb_bucket = 1 << mask_size;
        let bucket_size: usize = 4096;

        Counter {
            counts: vec![0; 1 << nb_bit],
            buckets: vec![Vec::with_capacity(bucket_size); nb_bucket],
            mask_size: mask_size,
            bucket_size: bucket_size,
            nb_bit: nb_bit,
        }
    }

    pub fn add_bit(&mut self, hash: u64) -> () {
        let mask_prefix = ((1 << (self.mask_size + 1)) - 1) << (self.nb_bit - self.mask_size);
        let prefix: usize = (mask_prefix & hash as usize) >> (self.nb_bit - self.mask_size);

        let bucket: &mut Vec<u64> = &mut self.buckets[prefix];
        bucket.push(hash);

        if bucket.len() == self.bucket_size {
            for hash in bucket.drain(..).map(|x| x as usize) {
                self.counts[hash] = self.counts[hash].saturating_add(1);
            }
        }
    }

    pub fn clean_all_buckets(&mut self) -> () {
        for mut bucket in self.buckets.drain(..) {
            for hash in bucket.drain(..).map(|x| x as usize) {
                self.counts[hash] = self.counts[hash].saturating_add(1);
            }
        }
    }
}
