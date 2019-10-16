/* project mod declaration */
mod convert;
mod count;
mod counter;
mod dump;
mod io;
//mod minimizer;
mod write;
mod bucketizer;
mod lookup_table;

use std::io::Read;

#[no_mangle]
pub extern "C" fn ssik_hash(subseq: *const std::os::raw::c_uchar, k: u8) -> u64 {
    return convert::hash(unsafe{std::slice::from_raw_parts(subseq, k as usize)}, k);
}

#[no_mangle]
pub extern "C" fn ssik_revhash(mut kmer: u64, k: u8) -> *const std::os::raw::c_char {
    return std::ffi::CString::new(dump::reverse_hash(kmer, k)).unwrap().into_raw();
}

#[no_mangle]
pub extern "C" fn ssik_get_first_bit(kmer: u64) -> bool {
    return convert::get_first_bit(kmer);
}

#[no_mangle]
pub extern "C" fn ssik_cannonical(kmer: u64, k: u8) -> u64 {
    return convert::cannonical(kmer, k);
}

#[no_mangle]
pub extern "C" fn ssik_nuc2bit(nuc: u8) -> u64 {
    return convert::nuc2bit(nuc);
}

#[no_mangle]
pub extern "C" fn ssik_seq2bit(subseq: *const std::os::raw::c_uchar, len: usize) -> u64 {
    return convert::seq2bit(unsafe{std::slice::from_raw_parts(subseq, len)});
}

#[no_mangle]
pub extern "C" fn ssik_read_count(path: *const std::os::raw::c_char, data: *mut std::os::raw::c_uchar) -> *const std::os::raw::c_uchar {
    let input_path = unsafe{std::ffi::CStr::from_ptr(path).to_str().unwrap()};

    let mut reader = std::io::BufReader::new(std::fs::File::open(input_path).unwrap());

    let (k, nb_bit) = dump::read_header(&mut reader);

    let mut data: &mut [u8] = unsafe{std::slice::from_raw_parts_mut(data, ssik_get_data_size(k, nb_bit) as usize)};

    reader.read_exact(data).expect("Error durring read data");

    return data.as_ptr();
}

#[no_mangle]
pub extern "C" fn ssik_get_header(subseq: *const std::os::raw::c_char, mut k: &mut u8, nb_bit: &mut u8) -> () {
    let input_path = unsafe{std::ffi::CStr::from_ptr(subseq).to_str().unwrap()};

    let mut reader = std::io::BufReader::new(std::fs::File::open(input_path).unwrap());

    let (tmp_k, tmp_nb_bit) = dump::read_header(&mut reader);

    *k = tmp_k;
    *nb_bit = tmp_nb_bit;
} 

#[no_mangle]
pub extern "C" fn ssik_get_count(count: *const std::os::raw::c_uchar, hash: u64, nb_bit: u8) -> u8 {
    return match nb_bit {
        8 => get_count_8bit(count, hash),
        4 => get_count_4bit(count, hash),
        _ => panic!("Number bit isn't valid check your ssik count file."),
    };
}

fn get_count_8bit(count: *const std::os::raw::c_uchar, hash: u64) -> u8 {
    return unsafe{*count.offset(hash as isize)};
}

fn get_count_4bit(count: *const std::os::raw::c_uchar, hash: u64) -> u8 {
    return match convert::get_first_bit(hash) {
        true => get_count_8bit(count, convert::remove_first_bit(hash)) & 0b11110000,
        false => get_count_8bit(count, convert::remove_first_bit(hash)) & 0b1111,
    };
}

#[no_mangle]
pub extern "C" fn ssik_get_kmer_space_size(k: u8) -> u64 {
    return 1 << (k * 2 - 1);
}

#[no_mangle]
pub extern "C" fn ssik_get_data_size(k: u8, nb_bit: u8) -> u64 {
    return match nb_bit {
        8 => ssik_get_kmer_space_size(k),
        4 => ssik_get_kmer_space_size(k) / 2,
        _ => panic!("Number bit isn't valid check your ssik count file."),
    };
}
