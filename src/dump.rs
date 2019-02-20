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
use std::io::Seek;

/* crate use */
use csv;

/* project use */
use crate::convert;
use crate::io::Mode;

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

        let (k, mode, nb_bytes) = read_header(&mut reader);

        match nb_bytes {
            2 => match mode {
                Mode::AllCounts => self.dump_all_counts(writer, reader, k, abundance),
                Mode::Counts => self.dump_counts(writer, reader, k, abundance),
                Mode::KmerCounts => self.dump_kmer_counts(writer, reader, k, abundance),
                Mode::Numpy => self.dump_numpy(writer, reader, k, abundance),
            },
            _ => match mode {
                Mode::AllCounts => self.dump_all_counts(writer, reader, k, abundance),
                Mode::Counts => self.dump_counts(writer, reader, k, abundance),
                Mode::KmerCounts => self.dump_kmer_counts(writer, reader, k, abundance),
                Mode::Numpy => self.dump_numpy(writer, reader, k, abundance),
            },
        }
    }

    fn dump_counts(
        &self,
        mut writer: csv::Writer<std::io::BufWriter<std::fs::File>>,
        mut reader: std::io::BufReader<std::fs::File>,
        k: u8,
        abundance: T,
    ) -> () {
        let mut i: u64 = 0;

        let mut read_buf = vec![0; self.counts_bufsize() as usize];
        while reader.read_exact(&mut read_buf).is_ok() {
            let (dist, val) = self.counts_read_buffer(&read_buf);

            i += dist as u64;
            if val.into() >= abundance.into() {
                writer
                    .write_record(&[reverse_hash(i, k), val.to_string()])
                    .unwrap();
            }
        }
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

    fn dump_kmer_counts(
        &self,
        mut writer: csv::Writer<std::io::BufWriter<std::fs::File>>,
        mut reader: std::io::BufReader<std::fs::File>,
        k: u8,
        abundance: T,
    ) -> () {
        let mut read_buf = vec![0; self.kmer_counts_bufsize()];

        while reader.read_exact(&mut read_buf).is_ok() {
            let (kmer, count): (u64, T) = self.kmer_counts_read_buffer(&read_buf);

            if count.into() >= abundance.into() {
                writer
                    .write_record(&[reverse_hash(kmer as u64, k), count.to_string()])
                    .unwrap();
            }
        }
    }

    fn dump_numpy(
        &self,
        writer: csv::Writer<std::io::BufWriter<std::fs::File>>,
        mut reader: std::io::BufReader<std::fs::File>,
        k: u8,
        abundance: T,
    ) -> () {
        reader
            .seek(std::io::SeekFrom::Start(128))
            .expect("Error durring numpy reading");
        self.dump_all_counts(writer, reader, k, abundance);
    }

    fn counts_bufsize(&self) -> usize;
    fn counts_read_buffer(&self, buffer: &[u8]) -> (u64, T);
    fn all_counts_bufsize(&self) -> usize;
    fn all_counts_read_buffer(&self, buffer: &[u8]) -> T;
    fn kmer_counts_bufsize(&self) -> usize;
    fn kmer_counts_read_buffer(&self, buffer: &[u8]) -> (u64, T);
}

pub struct Dumper {}

impl AbstractDump<u8> for Dumper {
    fn counts_bufsize(&self) -> usize {
        2
    }

    fn counts_read_buffer(&self, buffer: &[u8]) -> (u64, u8) {
        (buffer[0] as u64, buffer[1])
    }

    fn all_counts_bufsize(&self) -> usize {
        1
    }

    fn all_counts_read_buffer(&self, buffer: &[u8]) -> u8 {
        buffer[0]
    }

    fn kmer_counts_bufsize(&self) -> usize {
        9
    }

    fn kmer_counts_read_buffer(&self, buffer: &[u8]) -> (u64, u8) {
        (
            u64::from_be_bytes([
                buffer[0], buffer[1], buffer[2], buffer[3], buffer[4], buffer[5], buffer[6],
                buffer[7],
            ]),
            buffer[8],
        )
    }
}

impl AbstractDump<u16> for Dumper {
    fn counts_bufsize(&self) -> usize {
        4
    }

    fn counts_read_buffer(&self, buffer: &[u8]) -> (u64, u16) {
        (
            u64::from_be_bytes([0, 0, 0, 0, 0, 0, buffer[0], buffer[1]]),
            u16::from_be_bytes([buffer[2], buffer[3]]),
        )
    }

    fn all_counts_bufsize(&self) -> usize {
        2
    }

    fn all_counts_read_buffer(&self, buffer: &[u8]) -> u16 {
        u16::from_be_bytes([buffer[0], buffer[1]])
    }

    fn kmer_counts_bufsize(&self) -> usize {
        10
    }

    fn kmer_counts_read_buffer(&self, buffer: &[u8]) -> (u64, u16) {
        (
            u64::from_be_bytes([
                buffer[0], buffer[1], buffer[2], buffer[3], buffer[4], buffer[5], buffer[6],
                buffer[7],
            ]),
            u16::from_be_bytes([buffer[8], buffer[9]]),
        )
    }
}

fn read_header(reader: &mut std::io::BufReader<std::fs::File>) -> (u8, Mode, u8) {
    let header_buff: &mut [u8] = &mut [0; 3];
    reader
        .read_exact(header_buff)
        .expect("Error when try to read the header.");

    if header_buff == &[147, b'N', b'U'] {
        reader
            .seek(std::io::SeekFrom::Start(22))
            .expect("Error when try move in count file");
        let type_buff: &mut [u8] = &mut [0; 2];
        reader
            .read_exact(type_buff)
            .expect("Error when try to read nb_bytes");

        let nb_bytes = if type_buff == &[b'u', b'1'] {
            1
        } else if type_buff == &[b'u', b'2'] {
            2
        } else {
            1
        };

        let k =
            (f64::log2((reader.get_ref().metadata().unwrap().len() - 128) as f64) as u8 + 1) / 2;

        reader
            .seek(std::io::SeekFrom::Start(0))
            .expect("Error when try move in count file");
        return (k, Mode::Numpy, nb_bytes);
    }

    return (header_buff[0], Mode::from(header_buff[1]), header_buff[2]);
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
