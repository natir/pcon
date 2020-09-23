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
use crate::dump;
use crate::solid;

/* Helper */
fn cstr2string(c_path: *const std::os::raw::c_char) -> String {
    unsafe { String::from_utf8_unchecked(std::ffi::CStr::from_ptr(c_path).to_bytes().to_owned()) }
}

fn reader_from_c_path(
    c_path: *const std::os::raw::c_char,
) -> Result<std::io::BufReader<Box<dyn std::io::Read>>, error::IO> {
    let path = cstr2string(c_path);

    let file = match std::fs::File::open(&path) {
        Ok(f) => f,
        Err(_) => return Err(IO::CantOpenFile),
    };

    let niffler = match niffler::get_reader(Box::new(file)) {
        Ok(f) => f,
        Err(_) => return Err(IO::ErrorDurringRead),
    };

    Ok(std::io::BufReader::new(niffler.0))
}

fn writer_from_c_path(
    c_path: *const std::os::raw::c_char,
) -> Result<std::io::BufWriter<Box<dyn std::io::Write>>, error::IO> {
    let path = cstr2string(c_path);

    let file = match std::fs::File::create(&path) {
        Ok(f) => f,
        Err(_) => return Err(IO::CantOpenFile),
    };

    Ok(std::io::BufWriter::new(Box::new(file)))
}

/* Error section */
/// Create a new pcon io error it's init to no error, see [error::IO]. In python corresponding string error is emit.
#[no_mangle]
pub extern "C" fn pcon_error_new() -> *mut error::IO {
    Box::into_raw(Box::new(error::IO::NoError))
}

/// Free a pcon io error
///
/// # Safety
/// It's safe
#[no_mangle]
pub unsafe extern "C" fn pcon_error_free(error: *mut error::IO) {
    if error.is_null() {
        return;
    }

    let boxed = Box::from_raw(error);

    drop(boxed);
}

/* Counter section */
/// Create a new Counter. In python binding Counter is an object, new is the default constructor.
/// See [counter::Counter::new].
#[no_mangle]
pub extern "C" fn pcon_counter_new(k: u8) -> *mut counter::Counter {
    Box::into_raw(Box::new(counter::Counter::new(k)))
}

/// Free a Counter. In Python use del on Counter object.
///
/// # Safety
/// It's safe
#[no_mangle]
pub unsafe extern "C" fn pcon_counter_free(counter: *mut counter::Counter) {
    if counter.is_null() {
        return;
    }

    let boxed = Box::from_raw(counter);

    drop(boxed);
}

/// Perform count of kmer in fasta file in path, this file can be compress in gzip, bzip2, xz.
/// You must check value of `io_error` is equal to NoError before use `counter`.
///
/// In Python it's count_fasta method of Counter object.
/// See [counter::Counter::count_fasta].
#[no_mangle]
pub extern "C" fn pcon_counter_count_fasta(
    counter: &mut counter::Counter,
    c_path: *const std::os::raw::c_char,
    read_buffer_len: usize,
    io_error: &mut error::IO,
) {
    let reader = reader_from_c_path(c_path);

    match reader {
        Ok(r) => counter.count_fasta(r, read_buffer_len),
        Err(e) => *io_error = e,
    }
}

/// Increase the count of `kmer`
///
/// In Python it's inc method of Counter object.
/// See [counter::Counter::inc].
#[no_mangle]
pub extern "C" fn pcon_counter_inc(counter: &mut counter::Counter, kmer: u64) {
    counter.inc(kmer);
}

/// Increase the count of a canonical `kmer`
///
/// In Python it's inc_canonic method of Counter object.
/// See [counter::Counter::inc_canonic].
#[no_mangle]
pub extern "C" fn pcon_counter_inc_canonic(counter: &mut counter::Counter, kmer: u64) {
    counter.inc_canonic(kmer);
}

/// Get the count of value `kmer`
///
/// In Python it's get method of Counter object.
/// See [counter::Counter::get].
#[no_mangle]
pub extern "C" fn pcon_counter_get(counter: &counter::Counter, kmer: u64) -> counter::Count {
    counter.get(kmer)
}

/// Get the count of value a canonical `kmer`
///
/// In Python it's get_canonic method of Counter object.
/// See [counter::Counter::get_canonic].
#[no_mangle]
pub extern "C" fn pcon_counter_get_canonic(
    counter: &counter::Counter,
    kmer: u64,
) -> counter::Count {
    counter.get_canonic(kmer)
}

/// Serialize Counter in path of file
/// You must check value of `io_error` is equal to NoError before use `counter`
///
/// In Python it's serialize method of Counter object.
/// See [counter::Counter::serialize].
#[no_mangle]
pub extern "C" fn pcon_serialize_counter(
    counter: &counter::Counter,
    c_path: *const std::os::raw::c_char,
    min_abundance: u8,
    io_error: &mut error::IO,
) {
    let writer = writer_from_c_path(c_path);

    match writer {
        Ok(w) => match counter.serialize(w, min_abundance) {
            Ok(_) => (),
            Err(_) => *io_error = IO::ErrorDurringWrite,
        },
        Err(e) => *io_error = e,
    }
}

/// Deserialize Counter from `c_path` in `counter`
/// You must check value of `io_error` is equal to NoError before use `counter`
///
/// In Python it's deserialize class method of Counter.
/// See [counter::Counter::deserialize].
#[no_mangle]
pub extern "C" fn pcon_deserialize_counter(
    counter: &mut counter::Counter,
    c_path: *const std::os::raw::c_char,
    io_error: &mut error::IO,
) {
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
/// Create a new Solid. In python binding Solid is an object, new is the default constructor.
/// See [solid::Solid::new]
#[no_mangle]
pub extern "C" fn pcon_solid_new(k: u8) -> *mut solid::Solid {
    Box::into_raw(Box::new(solid::Solid::new(k)))
}

/// Create a new Solid from value in Counter
/// In python binding, this is a Solid class method from_counter.
/// See [solid::Solid::from_counter].
#[no_mangle]
pub extern "C" fn pcon_solid_from_counter(
    counter: &counter::Counter,
    abundance: counter::Count,
) -> *mut solid::Solid {
    Box::into_raw(Box::new(solid::Solid::from_counter(counter, abundance)))
}

/// Free a Solid. In Python use del on Solid object.
///
/// # Safety
/// It's safe
#[no_mangle]
pub unsafe extern "C" fn pcon_solid_free(solid: *mut solid::Solid) {
    if solid.is_null() {
        return;
    }

    let boxed = Box::from_raw(solid);

    drop(boxed);
}

/// Set the solidity status of `kmer` to `value`
///
/// In Python it's set method of Solid object.
/// See [solid::Solid::set].
#[no_mangle]
pub extern "C" fn pcon_solid_set(solid: &mut solid::Solid, kmer: u64, value: bool) {
    solid.set(kmer, value);
}

/// Set the solidity status of a canonical `kmer` to `value`
///
/// In Python it's set_canonic method of Solid object.
/// See [solid::Solid::set_canonic].
#[no_mangle]
pub extern "C" fn pcon_solid_set_canonic(solid: &mut solid::Solid, kmer: u64, value: bool) {
    solid.set_canonic(kmer, value);
}

/// Get the solidity status of `kmer`
///
/// In Python it's get method of Solid object.
/// See [solid::Solid::get].
#[no_mangle]
pub extern "C" fn pcon_solid_get(solid: &mut solid::Solid, kmer: u64) -> bool {
    solid.get(kmer)
}

/// Get the solidity status of a canonical `kmer`
///
/// In Python it's get_canonic method of Solid object.
/// See [solid::Solid::get_canonic].
#[no_mangle]
pub extern "C" fn pcon_solid_get_canonic(solid: &mut solid::Solid, kmer: u64) -> bool {
    solid.get_canonic(kmer)
}

/// Serialize Solid in path of file
/// You must check value of `io_error` is equal to NoError before use `solid`
///
/// In Python it's serialize method of Solid object.
/// See [solid::Solid::serialize].
#[no_mangle]
pub extern "C" fn pcon_serialize_solid(
    solid: &solid::Solid,
    c_path: *const std::os::raw::c_char,
    io_error: &mut error::IO,
) {
    let writer = writer_from_c_path(c_path);

    match writer {
        Ok(w) => match solid.serialize(w) {
            Ok(_) => (),
            Err(_) => *io_error = IO::ErrorDurringWrite,
        },
        Err(e) => *io_error = e,
    }
}

/// Deserialize Solid from `c_path` in `counter`
/// You must check value of `io_error` is equal to NoError before use `solid`
///
/// In Python it's deserialize class method of solid.
/// See [solid::Solid::deserialize].
#[no_mangle]
pub extern "C" fn pcon_deserialize_solid(
    solid: &mut solid::Solid,
    c_path: *const std::os::raw::c_char,
    io_error: &mut error::IO,
) {
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
/// See [dump::csv].
/// You must check value of `io_error` is equal to NoError to be sure no problem occure durring write
///
/// In Python it's csv function of dump module.
#[no_mangle]
pub extern "C" fn pcon_dump_csv(
    counter: &counter::Counter,
    abundance: counter::Count,
    c_path: *const std::os::raw::c_char,
    io_error: &mut error::IO,
) {
    let writer = writer_from_c_path(c_path);

    match writer {
        Ok(w) => match dump::csv(w, counter, abundance) {
            Ok(_) => (),
            Err(_) => *io_error = IO::ErrorDurringWrite,
        },
        Err(e) => *io_error = e,
    }
}

/// See [dump::solid()].
/// You must check value of `io_error` is equal to NoError to be sure no problem occure durring write
///
/// In Python it's solid function of dump module.
#[no_mangle]
pub extern "C" fn pcon_dump_solid(
    counter: &counter::Counter,
    abundance: counter::Count,
    c_path: *const std::os::raw::c_char,
    io_error: &mut error::IO,
) {
    let writer = writer_from_c_path(c_path);

    match writer {
        Ok(w) => match dump::solid(w, counter, abundance) {
            Ok(_) => (),
            Err(_) => *io_error = IO::ErrorDurringWrite,
        },
        Err(e) => *io_error = e,
    }
}

/// See [dump::spectrum].
/// You must check value of `io_error` is equal to NoError to be sure no problem occure durring write
///
/// In Python it's spectrum function of dump module.
#[no_mangle]
pub extern "C" fn pcon_dump_spectrum(
    counter: &counter::Counter,
    c_path: *const std::os::raw::c_char,
    io_error: &mut error::IO,
) {
    let writer = writer_from_c_path(c_path);

    match writer {
        Ok(w) => match dump::spectrum(w, counter) {
            Ok(_) => (),
            Err(_) => *io_error = IO::ErrorDurringWrite,
        },
        Err(e) => *io_error = e,
    }
}

/// See [set_count_nb_threads]
#[no_mangle]
pub extern "C" fn pcon_set_nb_threads(nb_threads: usize) {
    crate::set_nb_threads(nb_threads);
}
