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
use itertools::Itertools;

/* project use */
use crate::convert;
use crate::io::Mode;

fn read_header(reader: &mut std::io::BufReader<std::fs::File>) -> (u8, Mode, u8) {
    let header_buff: &mut [u8] = &mut [0; 3];
    reader.read_exact(header_buff).expect("Error when try to read the header.");
    
    if header_buff == &[147, b'N', b'U'] {
        reader.seek(std::io::SeekFrom::Start(22)).expect("Error when try move in count file");
        let type_buff: &mut [u8] = &mut[0; 2];
        reader.read_exact(type_buff).expect("Error when try to read nb_bytes");
        
        let nb_bytes = if type_buff == &[b'u', b'1'] {
            1
        } else if type_buff == &[b'u', b'2'] {
            2
        } else {
            1
        };
        
        let k = ((f64::log2((reader.get_ref().metadata().unwrap().len() - 128) as f64) as u8 + 1) / 2);

        reader.seek(std::io::SeekFrom::Start(0)).expect("Error when try move in count file");
        return (k, Mode::Numpy, nb_bytes)
    }

    return (header_buff[0], Mode::from(header_buff[1]), header_buff[2])
}

pub fn dump(input_path: &str, output_path: &str, abundance: u16) -> () {
    let mut reader = std::io::BufReader::new(std::fs::File::open(input_path).unwrap());
    let out = std::io::BufWriter::new(std::fs::File::create(output_path).unwrap());
    let writer = csv::WriterBuilder::new().from_writer(out);

    let (k, mode, nb_bytes) = read_header(&mut reader);

    match nb_bytes {
        2 => match mode {
            Mode::AllCounts => dump_all_counts_u16(writer, reader, k, abundance),
            Mode::Counts => dump_counts_u16(writer, reader, k, abundance),
            Mode::KmerCounts => dump_kmer_counts_u16(writer, reader, k, abundance),
            Mode::Numpy => dump_numpy_u16(writer, reader, k, abundance),
        },
        _ => match mode {
            Mode::AllCounts => dump_all_counts_u8(writer, reader, k, abundance),
            Mode::Counts => dump_counts_u8(writer, reader, k, abundance),
            Mode::KmerCounts => dump_kmer_counts_u8(writer, reader, k, abundance),
            Mode::Numpy => dump_numpy_u8(writer, reader, k, abundance),
        },
    }
}

pub fn dump_all_counts_u8(mut writer: csv::Writer<std::io::BufWriter<std::fs::File>>, reader: std::io::BufReader<std::fs::File>, k: u8, abundance: u16) -> () {
    for (hash, c) in reader.bytes().enumerate() {
        let count = c.unwrap();
        
        if count as u16 >= abundance {
            let kmer = reverse_hash(hash as u64, k);

            writer.write_record(&[kmer, count.to_string()]).unwrap();
        }
    }
}

pub fn dump_all_counts_u16(mut writer: csv::Writer<std::io::BufWriter<std::fs::File>>, mut reader: std::io::BufReader<std::fs::File>, k: u8, abundance: u16) -> () {
    let mut read_buf = [0; 2];

    let mut hash = 0;
    while reader.read_exact(&mut read_buf).is_ok() {
        let count = u16::from_be_bytes(read_buf);
        
        if count >= abundance {    
            let kmer = reverse_hash(hash as u64, k);

            writer.write_record(&[kmer, count.to_string()]).unwrap();
        }

        hash += 1;
    }
}

pub fn dump_counts_u8(mut writer: csv::Writer<std::io::BufWriter<std::fs::File>>, mut reader: std::io::BufReader<std::fs::File>, k: u8, abundance: u16) -> () {
    let mut i: u64 = 0;

    let mut read_buf = [0; 2];
    while reader.read_exact(&mut read_buf).is_ok() {
        let dist = read_buf[0];
        let val = read_buf[1];

        i += dist as u64;
        if val as u16 >= abundance {
            writer
                .write_record(&[reverse_hash(i as u64, k), val.to_string()])
                .unwrap();
        }
    }
}

pub fn dump_counts_u16(mut writer: csv::Writer<std::io::BufWriter<std::fs::File>>, mut reader: std::io::BufReader<std::fs::File>, k: u8, abundance: u16) -> () {
    let mut i: u64 = 0;

    let mut read_buf = [0; 4];
    while reader.read_exact(&mut read_buf).is_ok() {
        let dist = u16::from_be_bytes([read_buf[0], read_buf[1]]);
        let val = u16::from_be_bytes([read_buf[2], read_buf[3]]);

        i += dist as u64;
        if val >= abundance {
            writer
                .write_record(&[reverse_hash(i as u64, k), val.to_string()])
                .unwrap();
        }
    }
}

pub fn dump_kmer_counts_u8(mut writer: csv::Writer<std::io::BufWriter<std::fs::File>>, mut reader: std::io::BufReader<std::fs::File>, k: u8, abundance: u16) -> () {
    let mut read_buf = [0; 9];

    while reader.read_exact(&mut read_buf).is_ok() {
        let kmer: u64 = u64::from_be_bytes([read_buf[0], read_buf[1], read_buf[2], read_buf[3], read_buf[4], read_buf[5], read_buf[6], read_buf[7]]);
        let count = read_buf[8];

        if count as u16 >= abundance {
            writer
                .write_record(&[reverse_hash(kmer as u64, k), count.to_string()])
                .unwrap();
        }
    }
}

pub fn dump_kmer_counts_u16(mut writer: csv::Writer<std::io::BufWriter<std::fs::File>>, mut reader: std::io::BufReader<std::fs::File>, k: u8, abundance: u16) -> () {
    let mut read_buf = [0; 10];
    
    while reader.read_exact(&mut read_buf).is_ok() {
        let kmer: u64 = u64::from_be_bytes([read_buf[0], read_buf[1], read_buf[2], read_buf[3], read_buf[4], read_buf[5], read_buf[6], read_buf[7]]);
        let count = u16::from_be_bytes([read_buf[8], read_buf[9]]);

        if count >= abundance {
            writer
                .write_record(&[reverse_hash(kmer as u64, k), count.to_string()])
                .unwrap();
        }
    }
}

pub fn dump_numpy_u8(mut writer: csv::Writer<std::io::BufWriter<std::fs::File>>, mut reader: std::io::BufReader<std::fs::File>, k: u8, abundance: u16) -> () {
    reader.seek(std::io::SeekFrom::Start(128)).expect("Error durring numpy reading");

    let mut read_buf = [0; 1];
    let mut hash = 0;

    while reader.read_exact(&mut read_buf).is_ok() {
        let count = read_buf[0];
        
        if count as u16 >= abundance {    
            let kmer = reverse_hash(hash as u64, k);

            writer.write_record(&[kmer, count.to_string()]).unwrap();
        }

        hash += 1;
    }
}

pub fn dump_numpy_u16(mut writer: csv::Writer<std::io::BufWriter<std::fs::File>>, mut reader: std::io::BufReader<std::fs::File>, k: u8, abundance: u16) -> () {
    reader.seek(std::io::SeekFrom::Start(128)).expect("Error durring numpy reading");

    let mut read_buf = [0; 2];
    let mut hash = 0;

    while reader.read_exact(&mut read_buf).is_ok() {
        let count = u16::from_be_bytes(read_buf);
        
        if count >= abundance {    
            let kmer = reverse_hash(hash as u64, k);

            writer.write_record(&[kmer, count.to_string()]).unwrap();
        }

        hash += 1;
    }
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
