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

/* standard use */

/* project use */
use crate::convert;
use crate::counter;
use crate::bucketizer;
use crate::io::write;

use crate::bucketizer::Bucket;

pub fn count(multi_input_path: Vec<&str>, output_path: &str, k: u8, m: u8, n: u8) -> () {

    let mut counter: Box<dyn counter::Counter<u8, u64>> = match n {
        4 => Box::new(counter::ShortCounter::new(k)),
        8 => Box::new(counter::BasicCounter::<u8>::new(k)),
        _ => Box::new(counter::ShortCounter::new(k)),
    };
    
    for input_path in multi_input_path {    
        let reader = bio::io::fasta::Reader::new(std::io::BufReader::new(std::fs::File::open(input_path).unwrap()));
        perform_count(reader, &mut counter, k, m);
    }
        
    let mut out = std::io::BufWriter::new(std::fs::File::create(output_path).unwrap());

    write::do_it(&mut out, counter, k);
}

fn perform_count<R: std::io::Read>(reader: bio::io::fasta::Reader<R>, counter: &mut Box<dyn counter::Counter<u8, u64>>, k: u8, m: u8)  {

    let mut bucketizer: bucketizer::Prefix<u8> = bucketizer::Prefix::new(counter, k);

    // let minimizer_mask = generate_mask(m * 2); required for minimizer search
    let kmer_mask = generate_mask(k * 2);

    for record in reader.records() {
        let result = record.unwrap();
        let line = result.seq();

        if line.len() < k as usize { continue; } // if line is lower than k we have nothing to count
        
        let mut kmer = convert::seq2bit(&line[0..k as usize]);
        bucketizer.add_kmer(convert::cannonical(kmer, k));
        
        if line.len() == k as usize { continue; } // if line is equal to k we count just the first kmer
        
        // found minimizer of first kmer
        //let (mut minimizer, mut mini_hash, mut mini_pos) = found_minimizer(kmer, k, m, minimizer_mask);

        for nuc in &line[1..(line.len() - (k as usize))] {

            // create space for new nuc; clean old nuc; add new nuc
            kmer = ((kmer << 2) & kmer_mask) | convert::nuc2bit(*nuc);

            /* Minimizer search
            mini_pos += 2; 
            
            if mini_pos > (k * 2 - m * 2) as i8 { // last minimizer go away we need to found the new minimizer
                let tmp = found_minimizer(kmer, k, m, minimizer_mask);
                minimizer = tmp.0; mini_hash = tmp.1; mini_pos = tmp.2;
            } else {
                if let Some(tmp) = update_minimizer(kmer, minimizer_mask, mini_hash) { // check if new nuc create a lowest minimizer
                    minimizer = tmp.0; mini_hash = tmp.1; mini_pos = tmp.2;
                }
            }
            */

            // add kmer in count structure
            bucketizer.add_kmer(convert::cannonical(kmer, k));
        }
    }

    bucketizer.clean_all_buckets();
}

fn generate_mask(size: u8) -> u64 {
    return (1 << size) - 1;
}

fn minimizer_mask(mask: u64, offset: usize) -> u64 {
    return mask << offset;
}

fn update_minimizer(kmer: u64, mask: u64, hash: i64) -> Option<(u64, i64, i8)> {
    let new_mini = kmer & mask;
    let new_hash = revhash(new_mini);
    
    if new_hash < hash {
        return Some((new_mini, new_hash, 0));
    }

    return None
}

fn found_minimizer(kmer: u64, k: u8, m: u8, mask: u64) -> (u64, i64, i8) {
    let mut minimizer = kmer & minimizer_mask(mask, 0);
    let mut minimizer_hash = revhash(minimizer);
    let mut index = 0;

    for i in (2..=(k as usize * 2 - m as usize * 2)).step_by(2) {
        let subk = (kmer & minimizer_mask(mask, i)) >> i;
        let hash = revhash(subk);
        if minimizer_hash > hash {
            minimizer = subk;
            minimizer_hash = hash;
            index = i;
        }
    }

    return (minimizer, minimizer_hash, index as i8);
}

fn revhash(mut x: u64) -> i64 {
    x = ((x >> 32) ^ x).wrapping_mul(0xD6E8FEB86659FD93);
    x = ((x >> 32) ^ x).wrapping_mul(0xD6E8FEB86659FD93);
    x = (x >> 32) ^ x;

    return x as i64;
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::counter::Counter;
    
    #[test]
    fn global() -> () {
        let mut counter: Box<dyn counter::Counter<u8, u64>> = Box::new(counter::ShortCounter::new(5));
        
        let mut canonical = bio::io::fasta::Reader::new(std::io::Cursor::new(">1
ACTAG
>2
CTAGT
"));

        perform_count(canonical, &mut counter, 5, 1);
        
        assert_eq!(counter.get(convert::hash(b"ACTAG", 5)), 2);
    }
}
