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

/* project use */
use crate::counter;
use crate::io::Mode;

/* std use*/
use std::io::Write;

fn write_header(out: &mut std::io::BufWriter<std::fs::File>, k: u8, mode: &Mode, count_nb_bit: u8) {
    if mode != &Mode::Numpy {
        out.write(&[k, mode.into(), count_nb_bit]).expect("Error when try to write header.");
    }
}


pub trait AbstractWriter<T> {
    fn write(&self, count: &counter::Counter<T, u64, u64>, output_path: &str, k: u8, mode: Mode) -> () {

        let mut out = std::io::BufWriter::new(std::fs::File::create(output_path).unwrap());
        
        write_header(&mut out, k, &mode, std::mem::size_of::<T>() as u8);
        
        match mode {
            Mode::AllCounts => (self.write_all_counts(count, &mut out, k)),
            Mode::Counts => (self.write_counts(count, &mut out, k)),
            Mode::KmerCounts => (self.write_kmer_counts(count, &mut out, k)),
            Mode::Numpy => (self.write_numpy(count, &mut out, k)),
        }
    }

    fn write_values<W: std::io::Write>(&self, out: &mut W, val: &[u8]) -> () {
        out.write(val).expect("Error durring write count on disk");
    }
    
    fn write_value<W: std::io::Write>(&self, out: &mut W, val: T) -> ();
    fn write_all_counts(&self, count: &counter::Counter<T, u64, u64>, out: &mut std::io::BufWriter<std::fs::File>, k: u8) -> ();
    fn write_counts(&self, count: &counter::Counter<T, u64, u64>, out: &mut std::io::BufWriter<std::fs::File>, k: u8) -> ();
    fn write_kmer_counts(&self, count: &counter::Counter<T, u64, u64>, out: &mut std::io::BufWriter<std::fs::File>, k: u8) -> ();
    fn write_numpy(&self, count: &counter::Counter<T, u64, u64>, out: &mut std::io::BufWriter<std::fs::File>, k: u8) -> ();
}

pub struct Writer<T> {
    phantom: std::marker::PhantomData<T>,
}

impl<T> Writer<T> {
    pub fn new() -> Self {
        Writer {
            phantom: std::marker::PhantomData,
        }
    }
}

impl AbstractWriter<u8> for Writer<u8> {
    fn write_value<W: std::io::Write>(&self, out: &mut W, val: u8) -> () {
        self.write_values(out, &[val]);
    }

    fn write_all_counts(&self, count: &counter::Counter<u8, u64, u64>, out: &mut std::io::BufWriter<std::fs::File>, k: u8) -> () {
        for i in 0..(1 << (k * 2 - 1)) {
            self.write_value(out, count.get(i));
        }
    }
    
    fn write_counts(&self, count: &counter::Counter<u8, u64, u64>, out: &mut std::io::BufWriter<std::fs::File>, k: u8) -> () {
        let mut last_write = 0;
        for i in 0..(1 << (k * 2 - 1)) {
            let val = count.get(i);
            
            if val == 0 {
                continue;
            }

            let mut dist = i - last_write;
            last_write = i;

            if dist > u8::max_value() as u64 {
                // dist overflow u8 we need write a some kmer with 0 count
                let n = dist / u8::max_value() as u64;
                for _ in 0..n {
                    self.write_values(out, &[u8::max_value(), count.get(i)]);
                }

                dist %= 255;
            }

            // write dist to last value and count of k
            self.write_values(out, &[dist as u8, val]);
        }
    }
        
    fn write_kmer_counts(&self, count: &counter::Counter<u8, u64, u64>, out: &mut std::io::BufWriter<std::fs::File>, k: u8) -> () {
        for i in 0..(1 << (k * 2 - 1)) {
            let val = count.get(i);
            
            if val == 0 {
                continue;
            }

            let k: [u8; 8] = unsafe{ std::mem::transmute(i.to_be()) };

            self.write_values(out, &k);
            self.write_value(out, val);
        }
    }

    fn write_numpy(&self, count: &counter::Counter<u8, u64, u64>, out: &mut std::io::BufWriter<std::fs::File>, k: u8) -> () {
        out.write(&[147, b'N', b'U', b'M', b'P', b'Y', 1, 0]).expect("Error durring numpy magic number write");

        let nb_value = 1 << (k * 2 - 1);
        let header = format!("{{'descr': '|u1', 'fortran_order': False, 'shape': (1, {}), }}", nb_value);

        out.write(&[b'v', 0]).expect("Error durring numpy magic number write");
        out.write(header.as_bytes()).expect("Error durring numpy magic number write");
        out.write(&vec![b' '; 128 - 9 - header.len() - 1]).expect("Error durring numpy magic number write");

        for i in 0..(1 << (k * 2 - 1)) {
            out.write(&[count.get(i)]).expect("Error durring numpy value write");
        }
    }
}


impl AbstractWriter<u16> for Writer<u16> {
    fn write_value<W: std::io::Write>(&self, out: &mut W, v: u16) -> () {
        let val: [u8; 2] = unsafe{ std::mem::transmute(v.to_be()) };
        self.write_values(out, &val);
    }

    fn write_all_counts(&self, count: &counter::Counter<u16, u64, u64>, out: &mut std::io::BufWriter<std::fs::File>, k: u8) -> () {
        for i in 0..(1 << (k * 2 - 1)) {
            self.write_value(out, count.get(i));
        }
    }
    
    fn write_counts(&self, count: &counter::Counter<u16, u64, u64>, out: &mut std::io::BufWriter<std::fs::File>, k: u8) -> () {
        let mut last_write = 0;
        for i in 0..(1 << (k * 2 - 1)) {
            let val = count.get(i);
            
            if val == 0 {
                continue;
            }

            let mut dist = i - last_write;
            last_write = i;

            if dist > u16::max_value() as u64 {
                // dist overflow u8 we need write a some kmer with 0 count
                let n = dist / u16::max_value() as u64;
                for _ in 0..n {
                    self.write_value(out, u16::max_value());
                    self.write_value(out, 0u16);
                }

                dist %= u16::max_value() as u64;
            }
            
            // write dist to last value and count of k
            self.write_value(out, dist as u16);
            self.write_value(out, val);
        }
    }
        
    fn write_kmer_counts(&self, count: &counter::Counter<u16, u64, u64>, out: &mut std::io::BufWriter<std::fs::File>, k: u8) -> () {
        for i in 0..((1 << (k * 2 - 1)) as u64 ) {
            let v = count.get(i);
            
            if v == 0 {
                continue;
            }

            let kmer: [u8; 8] = unsafe{ std::mem::transmute(i.to_be()) };

            self.write_values(out, &kmer);
            self.write_value(out, v);
        }
    }    

    fn write_numpy(&self, count: &counter::Counter<u16, u64, u64>, out: &mut std::io::BufWriter<std::fs::File>, k: u8) -> () {
        out.write(&[147, b'N', b'U', b'M', b'P', b'Y', 1, 0]).expect("Error durring numpy magic number write");

        let nb_value = 1 << (k * 2 - 1);
        let header = format!("{{'descr': '>u2', 'fortran_order': False, 'shape': (1, {}), }}", nb_value);

        out.write(&[b'v', 0]).expect("Error durring numpy magic number write");
        out.write(header.as_bytes()).expect("Error durring numpy magic number write");
        out.write(&vec![b' '; 128 - 9 - header.len() - 1]).expect("Error durring numpy magic number write");

        for i in 0..(1 << (k * 2 - 1)) {
            self.write_value(out, count.get(i));
        }
    }
}
