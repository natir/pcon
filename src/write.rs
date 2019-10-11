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
use std::io::Write;

/* std use*/

pub trait AbstractWrite<W: std::io::Write, C: counter::Counter<u8, u64>>  {
    fn do_it(out: &mut W, count: &C, k: u8) -> ();
}

pub struct Ssik;

impl AbstractWrite<std::io::BufWriter<std::fs::File>, counter::BasicCounter<u8>> for Ssik {
    fn do_it(out: &mut std::io::BufWriter<std::fs::File>, count: &counter::BasicCounter<u8>, k: u8) -> () {

        write_header(out, k, 8);

        out.write(&count.data).expect("Error durring data write");
    }
}

impl AbstractWrite<std::io::BufWriter<std::fs::File>, counter::ShortCounter> for Ssik {
    fn do_it(out: &mut std::io::BufWriter<std::fs::File>, count: &counter::ShortCounter, k: u8) -> () {

        write_header(out, k, 4);

        out.write(&count.data).expect("Error durring data write");
    }
}

fn write_header<W: std::io::Write>(out: &mut W, k: u8, count_nb_bit: u8) {
    out.write(&[k, count_nb_bit])
            .expect("Error durring write header");
}
