/*
Copyright (c) 2019 Pierre Marijon <pmarijon@mmci.univ-saarland.de>

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

/* project mod declaration */
pub mod convert;
pub mod count;
pub mod counter;
pub mod dump;
pub mod io;
pub mod bucketizer;
mod lookup_table;

use std::io::Read;

#[no_mangle]
pub extern "C" fn pcon_hash(subseq: *const std::os::raw::c_uchar, k: u8) -> u64 {
    return convert::hash(unsafe{std::slice::from_raw_parts(subseq, k as usize)}, k);
}

#[no_mangle]
pub extern "C" fn pcon_revhash(kmer: u64, k: u8) -> *const std::os::raw::c_char {
    return std::ffi::CString::new(dump::reverse_hash(kmer, k)).unwrap().into_raw();
}

#[no_mangle]
pub extern "C" fn pcon_get_first_bit(kmer: u64) -> bool {
    return convert::get_first_bit(kmer);
}

#[no_mangle]
pub extern "C" fn pcon_cannonical(kmer: u64, k: u8) -> u64 {
    return convert::cannonical(kmer, k);
}

#[no_mangle]
pub extern "C" fn pcon_nuc2bit(nuc: u8) -> u64 {
    return convert::nuc2bit(nuc);
}

#[no_mangle]
pub extern "C" fn pcon_seq2bit(subseq: *const std::os::raw::c_uchar, len: usize) -> u64 {
    return convert::seq2bit(unsafe{std::slice::from_raw_parts(subseq, len)});
}

#[no_mangle]
pub extern "C" fn pcon_read_count(path: *const std::os::raw::c_char, data: *mut std::os::raw::c_uchar) -> *const std::os::raw::c_uchar {
    let input_path = unsafe{std::ffi::CStr::from_ptr(path).to_str().unwrap()};

    let mut reader = std::io::BufReader::new(std::fs::File::open(input_path).unwrap());

    let (k, nb_bit) = io::read::read_header(&mut reader);

    let data: &mut [u8] = unsafe{std::slice::from_raw_parts_mut(data, pcon_get_data_size(k, nb_bit) as usize)};

    reader.read_exact(data).expect("Error durring read data");

    return data.as_ptr();
}

#[no_mangle]
pub extern "C" fn pcon_get_header(subseq: *const std::os::raw::c_char, k: &mut u8, nb_bit: &mut u8) -> () {
    let input_path = unsafe{std::ffi::CStr::from_ptr(subseq).to_str().unwrap()};

    let mut reader = std::io::BufReader::new(std::fs::File::open(input_path).unwrap());

    let (tmp_k, tmp_nb_bit) = io::read::read_header(&mut reader);

    *k = tmp_k;
    *nb_bit = tmp_nb_bit;
} 

#[no_mangle]
pub extern "C" fn pcon_get_count(count: *const std::os::raw::c_uchar, hash: u64, k: u8, nb_bit: u8) -> u8 {
    return io::read::get_count(unsafe{std::slice::from_raw_parts(count, io::read::get_data_size(k, nb_bit) as usize)},
                               hash,
                               nb_bit);
}

#[no_mangle]
pub extern "C" fn pcon_get_kmer_space_size(k: u8) -> u64 {
    return io::read::get_kmer_space_size(k);
}

#[no_mangle]
pub extern "C" fn pcon_get_data_size(k: u8, nb_bit: u8) -> u64 {
    return io::read::get_data_size(k, nb_bit);
}
