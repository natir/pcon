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

pub fn dump(input_path: &str, output_path: &str, abundance: u8) -> () {
    let mut reader = std::io::BufReader::new(std::fs::File::open(input_path).unwrap());
    let out = std::io::BufWriter::new(std::fs::File::create(output_path).unwrap());
    let writer = csv::WriterBuilder::new().from_writer(out);

    let numpy_buff: &mut [u8] = &mut [0; 6];
    reader.read_exact(numpy_buff).unwrap();

    if &numpy_buff == &[147, b'N', b'U', b'M', b'P', b'Y'] {
        let k = (f64::log2((std::fs::metadata(input_path).unwrap().len() - 128) as f64) as u8 + 1) / 2;

        return dump_numpy(writer, reader, abundance, k);
    }
    
    reader.seek(std::io::SeekFrom::Start(0)).unwrap();
    
    let header_buff: &mut [u8] = &mut [0; 2]; 
    reader
        .read_exact(header_buff)
        .expect("Error durring read count on disk");
    let k = header_buff[0];
    let mode = header_buff[1];
    
    match mode {
        0 => dump_all_counts(writer, reader, k, abundance),
        1 => dump_counts(writer, reader, k, abundance),
        2 => dump_kmer_counts(writer, reader, k, abundance),
        _ => dump_counts(writer, reader, k, abundance),
    }


}

pub fn dump_all_counts(mut writer: csv::Writer<std::io::BufWriter<std::fs::File>>, reader: std::io::BufReader<std::fs::File>, k: u8, abundance: u8) -> () {
    for (hash, c) in reader.bytes().enumerate() {
        let count = c.unwrap();
        
        if count >= abundance {    
            let kmer = reverse_hash(hash as u64, k);

            writer.write_record(&[kmer, count.to_string()]).unwrap();
        }
    }
}

pub fn dump_counts(mut writer: csv::Writer<std::io::BufWriter<std::fs::File>>, reader: std::io::BufReader<std::fs::File>, k: u8, abundance: u8) -> () {
    let mut i: u64 = 0;
    for mut chunks in reader.bytes().chunks(2).into_iter() {
        let dist = chunks.next().unwrap().unwrap();
        let val = chunks.next().unwrap().unwrap();

        i += dist as u64;
        if val >= abundance {
            writer
                .write_record(&[reverse_hash(i as u64, k), val.to_string()])
                .unwrap();
        }
    }
}


pub fn dump_kmer_counts(mut writer: csv::Writer<std::io::BufWriter<std::fs::File>>, mut reader: std::io::BufReader<std::fs::File>, k: u8, abundance: u8) -> () {
    let mut read_buf = [0; 9];
    
    while reader.read_exact(&mut read_buf).is_ok() {
        let kmer: u64 = u64::from_be_bytes([read_buf[0], read_buf[1], read_buf[2], read_buf[3], read_buf[4], read_buf[5], read_buf[6], read_buf[7]]);
        let count = read_buf[8];

        if count >= abundance {
            writer
                .write_record(&[reverse_hash(kmer as u64, k), count.to_string()])
                .unwrap();
        }
    }
    
    
}

pub fn dump_numpy(mut writer: csv::Writer<std::io::BufWriter<std::fs::File>>, mut reader: std::io::BufReader<std::fs::File>, abundance: u8, k: u8) -> () {
    reader.seek(std::io::SeekFrom::Start(128)).expect("Error durring numpy reading");

    for (hash, c) in reader.bytes().enumerate() {
        let count = c.unwrap();
        
        if count >= abundance {    
            let kmer = reverse_hash(hash as u64, k);

            writer.write_record(&[kmer, count.to_string()]).unwrap();
        }
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
