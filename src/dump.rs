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

/* std use */
use std::io::Read;

/* crate use */
use csv;

/* project use */
use crate::convert;

pub fn dump<T>(input_path: &str, output_path: &str, abundance: T) -> ()
where
    T: Into<u64> + Copy + std::string::ToString,
    Dumper: AbstractDump<T>,
{
    AbstractDump::<T>::run(&Dumper {}, input_path, output_path, abundance);
}

pub trait AbstractDump<T>
where
    T: Into<u64> + Copy + std::string::ToString,
{
    fn run(&self, input_path: &str, output_path: &str, abundance: T) -> () {
        let mut reader = std::io::BufReader::new(std::fs::File::open(input_path).unwrap());
        let out = std::io::BufWriter::new(std::fs::File::create(output_path).unwrap());
        let writer = csv::WriterBuilder::new().from_writer(out);

        let (k, nb_bytes) = read_header(&mut reader);

        self.dump_all_counts(writer, reader, k, abundance);
    }

    fn dump_all_counts(
        &self,
        mut writer: csv::Writer<std::io::BufWriter<std::fs::File>>,
        mut reader: std::io::BufReader<std::fs::File>,
        k: u8,
        abundance: T,
    ) -> () {
        let mut read_buf = vec![0; self.all_counts_bufsize()];

        let mut hash = 0;
        while reader.read_exact(&mut read_buf).is_ok() {
            let count = self.all_counts_read_buffer(&read_buf);

            if count.into() >= abundance.into() {
                let kmer = reverse_hash(hash as u64, k);

                writer.write_record(&[kmer, count.to_string()]).unwrap();
            }

            hash += 1;
        }
    }

    fn all_counts_bufsize(&self) -> usize;
    fn all_counts_read_buffer(&self, buffer: &[u8]) -> T;
}

pub struct Dumper {}

impl AbstractDump<u8> for Dumper {
    fn all_counts_bufsize(&self) -> usize {
        1
    }

    fn all_counts_read_buffer(&self, buffer: &[u8]) -> u8 {
        buffer[0]
    }
}

impl AbstractDump<u16> for Dumper {
    fn all_counts_bufsize(&self) -> usize {
        2
    }

    fn all_counts_read_buffer(&self, buffer: &[u8]) -> u16 {
        u16::from_be_bytes([buffer[0], buffer[1]])
    }
}

fn read_header(reader: &mut std::io::BufReader<std::fs::File>) -> (u8, u8) {
    let header_buff: &mut [u8] = &mut [0; 2];
    reader
        .read_exact(header_buff)
        .expect("Error when try to read the header.");

    return (header_buff[0], header_buff[2]);
}

fn reverse_hash(mut kmer: u64, k: u8) -> String {
    kmer <<= 1;

    if !convert::parity_even(kmer) {
        kmer = kmer + 1;
    }

    return convert::bit2seq(kmer, k);
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn reverse_hash_() {
        // 100011110 -> TAGGC
        assert_eq!(reverse_hash(0b100011110, 5), "TAGGC");
    }
}
