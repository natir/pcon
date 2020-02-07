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

/* crates use */
use bv;
use cocktail;

/* project use */
use crate::bucketizer;
use crate::counter;
use crate::set_set;

use crate::bucketizer::Bucket;

pub struct Count {
    k: u8,
    n: u8,
    buckets: bucketizer::Prefix<u8>,
}

impl Count {
    pub fn new(k: u8, n: u8) -> Self {
        let counter: Box<dyn counter::Counter<u8, u64>> = match n {
            4 => Box::new(counter::ShortCounter::new(k)),
            8 => Box::new(counter::BasicCounter::<u8>::new(k)),
            _ => Box::new(counter::ShortCounter::new(k)),
        };

        Count {
            k,
            n,
            buckets: bucketizer::Prefix::new(counter, k),
        }
    }

    pub fn deserialize<R>(reader: &mut R) -> Self
    where
        R: std::io::BufRead,
    {
        let (k, n) = Count::read_header(reader);

        let mut counter: Box<dyn counter::Counter<u8, u64>> = match n {
            4 => Box::new(counter::ShortCounter::new(k)),
            8 => Box::new(counter::BasicCounter::<u8>::new(k)),
            _ => Box::new(counter::ShortCounter::new(k)),
        };

        let mut data: Vec<u8> = vec![0u8; get_data_size(k, n) as usize];
        reader
            .read_exact(&mut data)
            .expect("Error durring read data");

        counter.set_data(data.into_boxed_slice());

        Count {
            k,
            n,
            buckets: bucketizer::Prefix::new(counter, k),
        }
    }

    pub fn serialize<W>(&self, out: &mut W)
    where
        W: std::io::Write,
    {
        self.write_header(out);

        let mut cumulative_bytes_write = 0;

        while cumulative_bytes_write < self.buckets.counter().data().len() {
            cumulative_bytes_write += out
                .write(&self.buckets.counter().data()[cumulative_bytes_write..])
                .expect("Error durring data write");
        }
    }

    pub fn add_sequence(&mut self, seq: &[u8]) {
        if seq.len() < self.k as usize {
            return;
        } // if line is lower than k we have nothing to count

        let k = self.k;
        for kmer in
            cocktail::tokenizer::Tokenizer::new(seq, k).map(|x| cocktail::kmer::cannonical(x, k))
        {
            self.buckets.add_kmer(kmer);
        }
    }

    pub fn clean_buckets(&mut self) {
        self.buckets.clean_all_buckets();
    }

    pub fn get_count(&mut self, kmer: u64) -> u8 {
        self.buckets.counter().get(kmer)
    }

    pub fn generate_bitfield(&self, abundance_min: u8) -> bv::BitVec<u8> {
        let mut bit_vec = bv::BitVec::new_fill(false, cocktail::kmer::get_hash_space_size(self.k));

        for kmer in 0..cocktail::kmer::get_hash_space_size(self.k) {
            if self.buckets.counter().get(kmer) >= abundance_min {
                bit_vec.set(kmer, true);
            }
        }

        bit_vec
    }

    pub fn generate_set_of_set(&self, abundance_min: u8, m: u8) -> set_set::SetOfSet {
	set_set::SetOfSet::new(self.buckets.counter(), self.k, m, abundance_min)
    }
    
    pub fn get_k(&self) -> u8 {
        self.k
    }

    pub fn get_n(&self) -> u8 {
        self.n
    }

    fn write_header<W>(&self, out: &mut W)
    where
        W: std::io::Write,
    {
        out.write(&[self.k, self.n])
            .expect("Error durring write header");
    }

    fn read_header<R>(reader: &mut R) -> (u8, u8)
    where
        R: std::io::BufRead,
    {
        let header_buff: &mut [u8] = &mut [0; 2];

        reader
            .read_exact(header_buff)
            .expect("Error when try to read the header.");

        (header_buff[0], header_buff[1])
    }
}

pub(crate) fn count(multi_input_path: Vec<&str>, output_path: &str, k: u8, n: u8) {
    let mut count = Count::new(k, n);

    for input_path in multi_input_path {
        let reader = bio::io::fasta::Reader::new(std::io::BufReader::new(
            std::fs::File::open(input_path).unwrap(),
        ));

        for record in reader.records() {
            let result = record.expect("Error durring read parsing");

            count.add_sequence(result.seq());
        }
    }

    count.clean_buckets();

    let mut out = std::io::BufWriter::new(std::fs::File::create(output_path).unwrap());

    count.serialize(&mut out);
}

pub(crate) fn get_data_size(k: u8, nb_bit: u8) -> u64 {
    match nb_bit {
        8 => cocktail::kmer::get_hash_space_size(k),
        4 => cocktail::kmer::get_hash_space_size(k) / 2,
        _ => panic!("Number bit isn't valid check your pcon count file."),
    }
}


#[cfg(test)]
mod test {
    use super::*;

    mod object {
        use super::*;

        #[test]
        fn check_size_of_counter() {
            let k5_n4 = Count::new(5, 4);
            assert_eq!(k5_n4.buckets.counter().data().len(), 256);

            let k5_n8 = Count::new(5, 8);
            assert_eq!(k5_n8.buckets.counter().data().len(), 512);

            let k5_n6 = Count::new(5, 6);
            assert_eq!(k5_n6.buckets.counter().data().len(), 256);
        }

        #[test]
        fn deserlialize() {
            let mut input = std::io::Cursor::new([0x01, 0x08, 0xff, 0x01]);

            let mut k1_n8 = Count::deserialize(&mut input);

            assert_eq!(k1_n8.get_k(), 1);
            assert_eq!(k1_n8.get_n(), 8);

            assert_eq!(k1_n8.get_count(0), 255);
            assert_eq!(k1_n8.get_count(1), 1);

            let mut input2 = std::io::Cursor::new([0x01, 0x04, 0x2f]);

            let mut k1_n4 = Count::deserialize(&mut input2);

            assert_eq!(k1_n4.get_k(), 1);
            assert_eq!(k1_n4.get_n(), 4);

            assert_eq!(k1_n4.get_count(0), 15);
            assert_eq!(k1_n4.get_count(1), 2);
        }

        #[test]
        fn serlialize() {
            let mut output = std::io::Cursor::new(Vec::new());

            let k1_n4 = Count::new(1, 4);
            k1_n4.serialize(&mut output);
            assert_eq!(output.into_inner(), vec![0x01, 0x04, 0x00]);

            output = std::io::Cursor::new(Vec::new());

            let k1_n8 = Count::new(1, 8);
            k1_n8.serialize(&mut output);
            assert_eq!(output.into_inner(), vec![0x01, 0x08, 0x00, 0x00]);
        }

        #[test]
        fn add_seq() {
            let mut k1_n8 = Count::new(1, 8);
            let seq = b"AAAAA";

            k1_n8.add_sequence(seq);
            assert_eq!(k1_n8.get_count(0), 0);

            k1_n8.clean_buckets();
            assert_eq!(k1_n8.get_count(0), 5);
        }

        #[test]
        fn bitfield() {
            let mut k1_n8 = Count::new(1, 8);
            let seq = b"AAAAA";

            k1_n8.add_sequence(seq);
            k1_n8.clean_buckets();

            assert_eq!(k1_n8.generate_bitfield(0).as_slice(), [true, true]);
            assert_eq!(k1_n8.generate_bitfield(1).as_slice(), [true, false]);
            assert_eq!(k1_n8.generate_bitfield(6).as_slice(), [false, false]);
        }

	#[test]
        fn set_set() {
            let mut k5_n8 = Count::new(5, 8);
            let seq = b"AGGAAGCTAC";

            k5_n8.add_sequence(seq);
            k5_n8.clean_buckets();

	    let set_of_set = k5_n8.generate_set_of_set(1, 3);

	    for seq in &[b"AGGAA", b"GGAAG", b"GAAGC", b"AAGCT", b"GCTAC"] {
		let kmer = cocktail::kmer::cannonical(cocktail::kmer::seq2bit(*seq), 5);
		assert_eq!(set_of_set.contains(kmer), true);
	    }

	    assert_eq!(set_of_set.contains(0), false); // 0 is kmer AAAAA
	}
    }

    #[test]
    fn good_data_size() {
        assert_eq!(get_data_size(5, 8), 512);
        assert_eq!(get_data_size(5, 4), 256);
    }

    #[test]
    #[should_panic]
    fn data_size_panic() {
        get_data_size(1, 2);
    }
}
