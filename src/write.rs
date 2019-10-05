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

/* std use*/
use std::io::Write;

pub fn write<T>(count: &dyn counter::Counter<T, u64>, output_path: &str, k: u8)
where
    T: std::marker::Copy,
    u64: std::convert::From<T>,
    Writer: AbstractWriter<T>,
{
    AbstractWriter::<T>::run(&Writer {}, count, output_path, k)
}

fn write_header(out: &mut std::io::BufWriter<std::fs::File>, k: u8, count_nb_bit: u8) {
    out.write(&[k, count_nb_bit])
            .expect("Error when try to write header.");
}

pub trait AbstractWriter<T>
where
    T: Into<u64> + Copy,
{
    fn run(
        &self,
        count: &dyn counter::Counter<T, u64>,
        output_path: &str,
        k: u8,
    ) -> () {
        let mut out = std::io::BufWriter::new(std::fs::File::create(output_path).unwrap());

        write_header(&mut out, k, std::mem::size_of::<T>() as u8);

        self.write_all_counts(count, &mut out, k); 
    }

    fn write_values<W: std::io::Write>(&self, out: &mut W, val: &[u8]) -> () {
        out.write(val).expect("Error durring write count on disk");
    }

    fn write_all_counts(
        &self,
        count: &dyn counter::Counter<T, u64>,
        out: &mut std::io::BufWriter<std::fs::File>,
        k: u8,
    ) -> () {
        for i in 0..(1 << (k * 2 - 1)) {
            self.write_value(out, count.get(i));
        }
    }
    
    fn write_value<W: std::io::Write>(&self, out: &mut W, val: T) -> ();
    fn max_value_count(&self) -> T;
    fn u64_to_type(&self, val: u64) -> T;
    fn zero(&self) -> T;
}

pub struct Writer {}

impl AbstractWriter<u8> for Writer {
    fn write_value<W: std::io::Write>(&self, out: &mut W, val: u8) -> () {
        AbstractWriter::<u8>::write_values(self, out, &[val]);
    }

    fn max_value_count(&self) -> u8 {
        u8::max_value()
    }

    fn u64_to_type(&self, val: u64) -> u8 {
        val as u8
    }

    fn zero(&self) -> u8 {
        0
    }
}

impl AbstractWriter<u16> for Writer {
    fn write_value<W: std::io::Write>(&self, out: &mut W, v: u16) -> () {
        let val: [u8; 2] = unsafe { std::mem::transmute(v.to_be()) };
        AbstractWriter::<u16>::write_values(self, out, &val);
    }

    fn max_value_count(&self) -> u16 {
        u16::max_value()
    }

    fn u64_to_type(&self, val: u64) -> u16 {
        val as u16
    }

    fn zero(&self) -> u16 {
        0
    }
}
