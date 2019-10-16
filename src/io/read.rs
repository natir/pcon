/*
Copyright (c) 2019 Pierre Marijon <pmarijon@mmci.uni-saarland.de>

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

use std::io::Read;

use crate::convert;

pub fn read_header(reader: &mut std::io::BufReader<std::fs::File>) -> (u8, u8) {
    let header_buff: &mut [u8] = &mut [0; 2];

    reader
        .read_exact(header_buff)
        .expect("Error when try to read the header.");

    return (header_buff[0], header_buff[1]);
}

pub fn read(path: &str) -> (u8, u8, Vec<u8>) {
    let mut reader = std::io::BufReader::new(std::fs::File::open(path).unwrap());

    let (k, nb_bit) = read_header(&mut reader);

    let mut data: Vec<u8> = vec![0u8; get_data_size(k, nb_bit) as usize];

    reader.read_exact(&mut data).expect("Error durring read data");

    return (k, nb_bit, data);
}

pub fn get_count(count: &[u8], hash: u64, nb_bit: u8) -> u8 {
    return match nb_bit {
        8 => get_count_8bit(count, hash),
        4 => get_count_4bit(count, hash),
        _ => panic!("Number bit isn't valid check your ssik count file."),
    };
}

fn get_count_8bit(count: &[u8], hash: u64) -> u8 {
    return count[hash as usize];
}

fn get_count_4bit(count: &[u8], hash: u64) -> u8 {
    return match convert::get_first_bit(hash) {
        true => get_count_8bit(count, convert::remove_first_bit(hash)) & 0b11110000,
        false => get_count_8bit(count, convert::remove_first_bit(hash)) & 0b1111,
    };
}

pub fn get_kmer_space_size(k: u8) -> u64 {
    return 1 << (k * 2 - 1);
}

pub fn get_data_size(k: u8, nb_bit: u8) -> u64 {
    return match nb_bit {
        8 => get_kmer_space_size(k),
        4 => get_kmer_space_size(k) / 2,
        _ => panic!("Number bit isn't valid check your ssik count file."),
    };
}
