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
use std::io::Write;

/* crate use */
use csv;

/* project use */
use crate::convert;
use crate::io::read;

#[repr(C)]
#[derive(Clone, PartialEq)]
pub enum Mode {
    Csv,
    Exist,
}

impl From<&str> for Mode {
    fn from(mode: &str) -> Self {
        return match mode {
            "csv" => Mode::Csv,
            "exist" => Mode::Exist,
            _ => Mode::Csv,
        };
    }
}

impl std::fmt::Display for Mode{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match *self {
            Mode::Csv => write!(f, "csv"),
            Mode::Exist => write!(f, "exist"),
        }
    }
}

pub fn dump(input_path: &str, output_path: &str, abundance: u8, mode: Mode) -> () {
    let mut reader = std::io::BufReader::new(std::fs::File::open(input_path).unwrap());

    let (k, nb_bit) = read::read_header(&mut reader);

    
    match (nb_bit, mode) {
        (8, Mode::Csv) => dump_count_u8_csv(reader, output_path, k, abundance),
        (8, Mode::Exist) => dump_count_u8_exist(reader, output_path, abundance),
        (4, Mode::Csv) => dump_count_u4_csv(reader, output_path, k, abundance),
        (4, Mode::Exist) => dump_count_u4_exist(reader, output_path, abundance),
        _ => panic!("This file wasn't a ssik count output"),
    };
}


fn dump_count_u8_csv(
    mut reader: std::io::BufReader<std::fs::File>,
    output_path: &str,
    k: u8,
    abundance: u8
) -> () {
    let out = std::io::BufWriter::new(std::fs::File::create(output_path).unwrap());
    let mut writer = csv::WriterBuilder::new().from_writer(out);
    let mut read_buf = vec![0; 1];

    let mut hash = 0;
    while reader.read_exact(&mut read_buf).is_ok() {
        let count = read_buf[0];

        if count >= abundance {
            let kmer = reverse_hash(hash as u64, k);

            writer.write_record(&[kmer, count.to_string()]).unwrap();
        }

        hash += 1;
    }
}

fn dump_count_u4_csv(
    mut reader: std::io::BufReader<std::fs::File>,
    output_path: &str,
    k: u8,
    abundance: u8
) -> () {
    let out = std::io::BufWriter::new(std::fs::File::create(output_path).unwrap());
    let mut writer = csv::WriterBuilder::new().from_writer(out);
    let mut read_buf = vec![0; 1];

    let mut hash = 0;
    while reader.read_exact(&mut read_buf).is_ok() {
        let data = read_buf[0];
        let mut count = data & 0b1111;
        
        if count >= abundance {
            let kmer = reverse_hash(hash as u64, k);

            writer.write_record(&[kmer, count.to_string()]).unwrap();
        }

        hash += 1;

        count = (data & 0b11110000) >> 4;

        if count >= abundance {
            let kmer = reverse_hash(hash as u64, k);

            writer.write_record(&[kmer, count.to_string()]).unwrap();
        }

        hash += 1;
    }
}


fn dump_count_u8_exist(
    mut reader: std::io::BufReader<std::fs::File>,
    output_path: &str,
    abundance: u8
) -> () {
    let mut out = std::io::BufWriter::new(std::fs::File::create(output_path).unwrap());
    let mut read_buf = vec![0; 1];
    let mut write_buf: u8 = 0;
    let mut write_buf_size: u8 = 0;
    
    while reader.read_exact(&mut read_buf).is_ok() {
        write_buf = populate_buf(write_buf, read_buf[0], abundance);        
        write_buf_size += 1;
        
        if write_buf_size == 8 {
            out.write(&[write_buf]).expect("Error durring write bitfield");
            write_buf_size = 0;
        }
    }
}

fn dump_count_u4_exist(
    mut reader: std::io::BufReader<std::fs::File>,
    output_path: &str,
    abundance: u8
) -> () {
    let mut out = std::io::BufWriter::new(std::fs::File::create(output_path).unwrap());
    let mut read_buf = vec![0; 1];
    let mut write_buf: u8 = 0;
    let mut write_buf_size: u8 = 0;

    while reader.read_exact(&mut read_buf).is_ok() {
        let data = read_buf[0];

        write_buf = populate_buf(write_buf, data & 0b1111, abundance);
        write_buf = populate_buf(write_buf, data & 0b11110000, abundance);
        
        write_buf_size += 2;
        
        if write_buf_size == 8 {
            out.write(&[write_buf]).expect("Error durring write bitfield");
            write_buf_size = 0;
        }
    }
}

fn populate_buf(buf: u8, count: u8, abundance: u8) -> u8 {
    if count >= abundance {
        return (buf ^ 1) << 1;
    } else {
        return buf << 1;
    }
}

pub fn reverse_hash(mut kmer: u64, k: u8) -> String {
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
