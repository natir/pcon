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

pub fn dump(input_path: &str, output_path: &str, abundance: u8) -> () {
    let reader = std::io::BufReader::new(
        std::fs::File::create(input_path).unwrap()
    );

    let k = (std::fs::metadata(input_path).unwrap().len() as f64).log2() as u8;
    
    let out = std::io::BufWriter::new(
        std::fs::File::create(output_path).unwrap(),
    );
    
    let mut writer = csv::WriterBuilder::new().from_writer(out);
    
    for (i, v) in reader.bytes().enumerate() {
        let val = v.unwrap();
        
        if val < abundance {
            writer
                .write_record(&[reverse_hash(i as u128, k), val.to_string()])
                .unwrap();
        }
    }
}

fn reverse_hash(mut kmer: u128, k: u8) -> String {
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
