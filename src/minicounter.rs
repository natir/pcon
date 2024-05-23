//! Generic struct of minicounter and implementation for many type

/* std use */

/* crate use */

/* project use */
use crate::counter;
use crate::error;

/// A counter of kmer, count only if minimizer is present more than a threshold.
#[derive(Clone, Eq, PartialEq, Debug, Default)]
pub struct MiniCounter<T> {
    k: u64,
    threshold: T,
    mini_count: counter::Counter<T>,
    kmer_count: rustc_hash::FxHashMap<Vec<u8>, T>,
}

/**************************/
/* generic implementation */
/**************************/
impl<T> MiniCounter<T> {
    /// Get value of k
    pub fn k(&self) -> u64 {
        self.k
    }

    /// Get value of m
    pub fn m(&self) -> u8 {
        self.mini_count.k()
    }

    /// Get mini_count at one index
    pub fn get_mini(&self, index: usize) -> &T {
        self.mini_count.get_raw(index)
    }

    /// Get mini count data
    pub fn mini_raw(&self) -> &[T] {
        self.mini_count.raw()
    }

    /// Get kmer count data
    pub fn kmer_raw(&self) -> &rustc_hash::FxHashMap<Vec<u8>, T> {
        &self.kmer_count
    }
}

/*****************************/
/* sequential implementation */
/*****************************/
macro_rules! impl_sequential (
    ($type:ty) => {
	impl MiniCounter<$type> {
	    /// Create a new MiniCounter with kmer size equal to k and minimizer size equal to m
	    pub fn new(k: u64, m: u8, threshold: $type) -> Self {
		Self {
		    k,
		    threshold,
		    mini_count: counter::Counter::<$type>::new(m),
		    kmer_count: rustc_hash::FxHashMap::default(),
		}
	    }

	    /// Perform count on fasta input
	    pub fn count_fasta(&mut self, fasta: Box<dyn std::io::BufRead>, _record_buffer: u64) {
		let mut reader = noodles::fasta::Reader::new(fasta);
		let mut records = reader.records();

		while let Some(Ok(record)) = records.next() {
		    if record.sequence().len() > self.k() as usize {
			let minimizer = cocktail::tokenizer::MiniBstr::new(
			    record.sequence().as_ref(),
			    self.k(),
			    self.m(),
			);

			let mut prev_mini = None;
			for (kmer, minimizer) in minimizer {
			    if prev_mini != Some(minimizer) {
				counter::Counter::<$type>::inc(self.mini_count.raw_mut(), (minimizer >> 1) as usize);
			    }
			    Self::inc(
				&mut self.mini_count,
				&mut self.kmer_count,
				minimizer as usize,
				kmer.to_ascii_uppercase(),
				self.threshold,
			    );

			    prev_mini = Some(minimizer);
			}
		    }
		}
	    }

	    #[cfg(feature = "fastq")]
	    /// Perform count on fastq input
	    pub fn count_fastq(&mut self, fastq: Box<dyn std::io::BufRead>, _record_buffer: u64) {
		let mut reader = noodles::fastq::Reader::new(fastq);
		let mut records = reader.records();

		while let Some(Ok(record)) = records.next() {
		    if record.sequence().len() >= self.k() as usize {
			let minimizer = cocktail::tokenizer::MiniBstr::new(
			    record.sequence().as_ref(),
			    self.k(),
			    self.m(),
			);

			let mut prev_mini = None;
			for (kmer, minimizer) in minimizer {
			    if prev_mini != Some(minimizer) {
				counter::Counter::<$type>::inc(self.mini_count.raw_mut(), (minimizer >> 1) as usize);
			    }

			    Self::inc(
				&mut self.mini_count,
				&mut self.kmer_count,
				minimizer as usize,
				kmer.to_ascii_uppercase(),
				self.threshold,
			    );

			    prev_mini = Some(minimizer);
			}
		    }
		}
	    }


	    /// Increment value at index
	    pub(crate) fn inc(
		mini_count: &mut counter::Counter<$type>,
		kmer_count: &mut rustc_hash::FxHashMap<Vec<u8>, $type>,
		minimizer: usize,
		kmer: Vec<u8>,
		threshold: $type,
	    ) {
		if mini_count.get(minimizer as u64) > threshold {
		    kmer_count
			.entry(kmer)
			.and_modify(|c| *c += 1)
			.or_insert(1);
		}
	    }

	    /// Get count of a kmer
	    pub fn get(&self, kmer: &[u8]) -> &$type {
		self.kmer_count.get(kmer).unwrap_or(&0)
	    }

	    /// Write minicounter result in csv
	    pub fn serialize<W>(&self, abundance: $type, mut output: W) -> error::Result<()>
	    where
		W: std::io::Write,
	    {
		for (kmer, count) in self.kmer_count.iter() {
		    if count > &abundance {
			writeln!(
			    output,
			    "{},{}",
			    unsafe { String::from_utf8_unchecked(kmer.to_vec()) },
			count
			)?;
		    }
		}
		Ok(())
	    }
	}
    }
);

impl_sequential!(u8);
impl_sequential!(u16);
impl_sequential!(u32);
impl_sequential!(u64);
impl_sequential!(u128);

#[cfg(test)]
mod tests {
    use super::*;

    const FASTA_FILE: &[u8] = b">random_seq 0
GTTCTGCAAATTAGAACAGACAATACACTGGCAGGCGTTGCGTTGGGGGAGATCttCCGTAACGAGCCGGCATTTGTAAGAAAGAGATTTCGAGTAAATG
>random_seq 1
AGGATAGAAGCTTAAGTACAAGATAATTCCCATAGAGGAAGGGTGGTAttACAGTGCCGCCTGTTGAAAGCCCCAATCCCGCTTCAATTGTTGAGCTCAG
";

    macro_rules! fasta_count {
        ($type:ty, $name:ident) => {
            const $name: &[$type] = &[
                0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 1,
                0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0,
                0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                1, 0, 0, 1, 0, 0, 0, 0,
            ];
        };
    }

    fasta_count!(u8, FASTA_COUNT_U8);
    fasta_count!(u16, FASTA_COUNT_U16);
    fasta_count!(u32, FASTA_COUNT_U32);
    fasta_count!(u64, FASTA_COUNT_U64);
    fasta_count!(u128, FASTA_COUNT_U128);

    macro_rules! sequential {
        ($type:ty, $name:ident, $truth:ident) => {
            #[test]
            fn $name() {
                let mut mini_count = MiniCounter::<$type>::new(10, 5, 1);

                mini_count.count_fasta(Box::new(FASTA_FILE), 1);

                assert_eq!(mini_count.mini_raw(), $truth);

                let mut result = vec![];
                for (kmer, count) in mini_count.kmer_raw().iter() {
                    result.push((kmer.to_vec(), *count));
                }
                result.sort();

                assert_eq!(
                    result,
                    vec![
                        (b"AAGATAATTC".to_vec(), 1),
                        (b"AGATAATTCC".to_vec(), 1),
                        (b"ATAATTCCCA".to_vec(), 1),
                        (b"ATTACAGTGC".to_vec(), 1),
                        (b"GATAATTCCC".to_vec(), 1),
                        (b"GGTATTACAG".to_vec(), 1),
                        (b"GTATTACAGT".to_vec(), 1),
                        (b"TAATTCCCAT".to_vec(), 1),
                        (b"TATTACAGTG".to_vec(), 1),
                        (b"TGGTATTACA".to_vec(), 1),
                        (b"TTACAGTGCC".to_vec(), 1),
                    ]
                )
            }
        };
    }

    sequential!(u8, sequential_u8, FASTA_COUNT_U8);
    sequential!(u16, sequential_u16, FASTA_COUNT_U16);
    sequential!(u32, sequential_u32, FASTA_COUNT_U32);
    sequential!(u64, sequential_u64, FASTA_COUNT_U64);
    sequential!(u128, sequential_u128, FASTA_COUNT_U128);
}
