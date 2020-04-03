/*
Copyright (c) 2020 Pierre Marijon <pmarijon@mpi-inf.mpg.de>

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

use crate::error;
use error::*;

use crate::counter;
use crate::solid;
use crate::dump;

/* Helper */
fn cstr2string(c_path: *const std::os::raw::c_char) -> String {
    unsafe { String::from_utf8_unchecked(std::ffi::CStr::from_ptr(c_path).to_bytes().to_owned()) }
}

fn reader_from_c_path(c_path: *const std::os::raw::c_char) -> Result<std::io::BufReader<Box<dyn std::io::Read>>, error::IO> {
    let path = cstr2string(c_path);

    let file = match std::fs::File::open(&path) {
	Ok(f) => f,
	Err(_) => return Err(IO::CantOpenFile)
    };

    let niffler = match niffler::get_reader(Box::new(file)) {
	Ok(f) => f,
	Err(_) => return Err(IO::ErrorDurringRead)
    };

    Ok(std::io::BufReader::new(niffler.0))
}

fn writer_from_c_path(c_path: *const std::os::raw::c_char) -> Result<std::io::BufWriter<Box<dyn std::io::Write>>, error::IO> {
    let path = cstr2string(c_path);

    let file = match std::fs::File::create(&path) {
	Ok(f) => f,
	Err(_) => return Err(IO::CantOpenFile)
    };

    Ok(std::io::BufWriter::new(Box::new(file)))
}

/* Error section */
#[no_mangle]
pub extern fn pcon_error_new() -> *mut error::IO {
    Box::into_raw(Box::new(error::IO::None))
}

#[no_mangle]
pub extern fn pcon_error_free(error: *mut error::IO) {
    if error.is_null() {
        return;
    }
    
    unsafe { Box::from_raw(error) };
}

/* Counter section */
#[no_mangle]
pub extern fn pcon_counter_new(k: u8) -> *mut counter::Counter {
    Box::into_raw(Box::new(counter::Counter::new(k)))
}

#[no_mangle]
pub extern fn pcon_counter_free(counter: *mut counter::Counter) {
    if counter.is_null() {
        return;
    }
    
    unsafe { Box::from_raw(counter) };
}

#[no_mangle]
pub extern fn pcon_counter_count_fasta(counter: &mut counter::Counter, c_path: *const std::os::raw::c_char, io_error: &mut error::IO) {
    let reader = reader_from_c_path(c_path);

    match reader {
	Ok(r) => counter.count_fasta(r),
	Err(e) => *io_error = e,
    }
}

#[no_mangle]
pub extern fn pcon_counter_count_fastq(counter: &mut counter::Counter, c_path: *const std::os::raw::c_char, io_error: &mut error::IO) {
    let reader = reader_from_c_path(c_path);

    match reader {
	Ok(r) => counter.count_fastq(r),
	Err(e) => *io_error = e,
    }
}

#[no_mangle]
pub extern fn pcon_counter_inc(counter: &mut counter::Counter, kmer: u64) {
    counter.inc(kmer);
}

#[no_mangle]
pub extern fn pcon_counter_get(counter: &counter::Counter, kmer: u64) -> counter::Count {
    counter.get(kmer)
}

#[no_mangle]
pub extern fn pcon_serialize_counter(counter: &counter::Counter, c_path: *const std::os::raw::c_char, io_error: &mut error::IO) {
    let writer = writer_from_c_path(c_path);

    match writer {
	Ok(w) => match counter.serialize(w) {
	    Ok(_) => (),
	    Err(_) => *io_error = IO::ErrorDurringWrite,
	},
	Err(e) => *io_error = e,
    }
}

#[no_mangle]
pub extern fn pcon_deserialize_counter(counter: &mut counter::Counter, c_path: *const std::os::raw::c_char, io_error: &mut error::IO) {
    let reader = reader_from_c_path(c_path);

    match reader {
	Ok(r) => match counter::Counter::deserialize(r) {
	    Ok(c) => *counter = c,
	    Err(_) => *io_error = IO::ErrorDurringRead,
	},
	Err(e) => *io_error = e,
    }
} 

/* solid section */
#[no_mangle]
pub extern fn pcon_solid_new(k: u8) -> *mut solid::Solid {
    Box::into_raw(Box::new(solid::Solid::new(k)))
}

#[no_mangle]
pub extern fn pcon_solid_from_counter(counter: &counter::Counter, abundance: u8) -> *mut solid::Solid {
    Box::into_raw(Box::new(solid::Solid::from_counter(counter, abundance)))
}

#[no_mangle]
pub extern fn pcon_solid_free(solid: *mut solid::Solid) {
    if solid.is_null() {
        return;
    }
    
    unsafe { Box::from_raw(solid) };
}

#[no_mangle]
pub extern fn pcon_solid_set(solid: &mut solid::Solid, kmer: u64, value: bool) {
    solid.set(kmer, value);
}

#[no_mangle]
pub extern fn pcon_solid_get(solid: &mut solid::Solid, kmer: u64) -> bool {
    solid.get(kmer)
}

#[no_mangle]
pub extern fn pcon_serialize_solid(solid: &solid::Solid, c_path: *const std::os::raw::c_char, io_error: &mut error::IO) {
    let writer = writer_from_c_path(c_path);

    match writer {
	Ok(w) => match solid.serialize(w) {
	    Ok(_) => (),
	    Err(_) => *io_error = IO::ErrorDurringWrite,
	},
	Err(e) => *io_error = e,
    }
}

#[no_mangle]
pub extern fn pcon_deserialize_solid(solid: &mut solid::Solid, c_path: *const std::os::raw::c_char, io_error: &mut error::IO) {
    let reader = reader_from_c_path(c_path);

    match reader {
	Ok(r) => match solid::Solid::deserialize(r) {
	    Ok(s) => *solid = s,
	    Err(_) => *io_error = IO::ErrorDurringRead,
	},
	Err(e) => *io_error = e,
    }
} 

/* Dump section */
#[no_mangle]
pub extern fn pcon_dump_csv(counter: &counter::Counter, abundance: u8, c_path: *const std::os::raw::c_char, io_error: &mut error::IO) {
    let writer = writer_from_c_path(c_path);

    match writer {
	Ok(w) => match dump::csv(w, counter, abundance) {
	    Ok(_) => (),
	    Err(_) => *io_error = IO::ErrorDurringWrite,
	},
	Err(e) => *io_error = e,
    }
}

#[no_mangle]
pub extern fn pcon_dump_solid(counter: &counter::Counter, abundance: u8, c_path: *const std::os::raw::c_char, io_error: &mut error::IO) {
    let writer = writer_from_c_path(c_path);

    match writer {
	Ok(w) => match dump::solid(w, counter, abundance) {
	    Ok(_) => (),
	    Err(_) => *io_error = IO::ErrorDurringWrite,
	},
	Err(e) => *io_error = e,
    }
}

#[no_mangle]
pub extern fn pcon_dump_spectrum(counter: &counter::Counter, c_path: *const std::os::raw::c_char, io_error: &mut error::IO) {
    let writer = writer_from_c_path(c_path);

    match writer {
	Ok(w) => match dump::spectrum(w, counter) {
	    Ok(_) => (),
	    Err(_) => *io_error = IO::ErrorDurringWrite,
	},
	Err(e) => *io_error = e,
    }
}

