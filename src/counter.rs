//! Run count command

/* std use */

/* crate use */
#[cfg(feature = "parallel")]
use rayon::prelude::*;

/* project use */

/// A counter of kmer based on cocktail crate 2bit conversion, canonicalisation and hashing.
/// Implement only for u8, std::sync::atomic::AtomicU8
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default)]
pub struct Counter<T> {
    k: u8,
    count: Box<[T]>,
}

/**************************/
/* generic implementation */
/**************************/
impl<T> Counter<T> {
    /// Get value of k
    pub fn k(&self) -> u8 {
        self.k
    }

    /// Get count at on index
    pub fn get_raw(&self, index: usize) -> &T {
        &self.count[index]
    }

    /// Get raw data
    pub fn raw(&self) -> &[T] {
        &self.count
    }
}

macro_rules! impl_sequential (
    ($type:ty, $fill:expr) => {
	impl Counter<$type> {
	    /// Create a new kmer Counter with kmer size equal to k
	    pub fn new(k: u8) -> Self {
		let tmp = vec![$fill; cocktail::kmer::get_hash_space_size(k) as usize];

		Self {
		    k,
		    count: tmp.into_boxed_slice(),
		}
	    }

	    /// Perform count on fasta input
	    pub fn count_fasta(&mut self, fasta: Box<dyn std::io::BufRead>, _record_buffer: u64) {
		let mut reader = noodles::fasta::Reader::new(fasta);
		let mut records = reader.records();

		while let Some(Ok(record)) = records.next() {
		    if record.sequence().len() >= self.k() as usize {
			let kmerizer =
			    cocktail::tokenizer::Canonical::new(record.sequence().as_ref(), self.k());

			for canonical in kmerizer {
			    Self::inc(&mut self.count, (canonical >> 1) as usize);
			}
		    }
		}
	    }

	    /// Increment value at index
	    fn inc(count: &mut [$type], index: usize) {
		count[index] = count[index].saturating_add(1);
	    }

	    /// Get count of a kmer
	    pub fn get(&self, kmer: u64) -> $type {
		self.get_canonic(cocktail::kmer::canonical(kmer, self.k))
	    }

	    /// Get the counter of a canonical kmer
	    fn get_canonic(&self, canonical: u64) -> $type {
		self.count[(canonical >> 1) as usize]
	    }
	}
    }
);

impl_sequential!(u8, 0u8);
impl_sequential!(u16, 0u16);
impl_sequential!(u32, 0u32);
impl_sequential!(u64, 0u64);
impl_sequential!(u128, 0u128);

macro_rules! impl_rayon (
    ($type:ty, $fill:expr, $out_type:ty, $max:expr) => {
	impl Counter<$type> {
	    /// Create a new kmer Counter with kmer size equal to k
	    pub fn new(k: u8) -> Self {
		let tmp = vec![$fill; cocktail::kmer::get_hash_space_size(k) as usize];

		Self {
		    k,
		    count: unsafe {
			std::mem::transmute::<Box<[$out_type]>, Box<[$type]>>(
			    tmp.into_boxed_slice(),
			)
		    },
		}
	    }

	    /// Perform count on fasta input
	    pub fn count_fasta(&mut self, fasta: Box<dyn std::io::BufRead>, record_buffer: u64) {
		let mut reader = noodles::fasta::Reader::new(fasta);
		let mut iter = reader.records();
		let mut records = Vec::with_capacity(record_buffer as usize);

		let mut end = true;
		while end {
		    log::info!("Start populate buffer");
		    end = populate_buffer(&mut iter, &mut records, record_buffer);
		    log::info!("End populate buffer {}", records.len());

		    records.par_iter().for_each(|record| {
			if record.sequence().len() >= self.k as usize {
			    let tokenizer =
				cocktail::tokenizer::Canonical::new(record.sequence().as_ref(), self.k);

			    for canonical in tokenizer {
				Self::inc(&self.count, (canonical >> 1) as usize);
			    }
			}
		    });
		}
	    }

	    /// Increment value at index
	    fn inc(count: &[$type], index: usize) {
		if count[index].load(std::sync::atomic::Ordering::SeqCst) != $max {
		    count[index].fetch_add(1, std::sync::atomic::Ordering::SeqCst);
		}
	    }

	    /// Get count of a kmer
	    pub fn get(&self, kmer: u64) -> $out_type {
		self.get_canonic(cocktail::kmer::canonical(kmer, self.k))
	    }

	    /// Get the counter of a canonical kmer
	    pub fn get_canonic(&self, canonical: u64) -> $out_type {
		self.count[(canonical >> 1) as usize].load(std::sync::atomic::Ordering::SeqCst)
	    }
	}

    }
);

#[cfg(feature = "parallel")]
impl_rayon!(std::sync::atomic::AtomicU8, 0u8, u8, u8::MAX);
#[cfg(feature = "parallel")]
impl_rayon!(std::sync::atomic::AtomicU16, 0u16, u16, u16::MAX);
#[cfg(feature = "parallel")]
impl_rayon!(std::sync::atomic::AtomicU32, 0u32, u32, u32::MAX);
#[cfg(feature = "parallel")]
impl_rayon!(std::sync::atomic::AtomicU64, 0u64, u64, u64::MAX);

/// Populate record buffer with content of iterator
fn populate_buffer(
    iter: &mut noodles::fasta::reader::Records<'_, Box<dyn std::io::BufRead>>,
    records: &mut Vec<noodles::fasta::Record>,
    record_buffer: u64,
) -> bool {
    records.clear();

    for i in 0..record_buffer {
        if let Some(Ok(record)) = iter.next() {
            records.push(record);
        } else {
            records.truncate(i as usize);
            return false;
        }
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;

    const FASTA_FILE: &[u8] = b">random_seq 0
GTTCTGCAAATTAGAACAGACAATACACTGGCAGGCGTTGCGTTGGGGGAGATCTTCCGTAACGAGCCGGCATTTGTAAGAAAGAGATTTCGAGTAAATG
>random_seq 1
AGGATAGAAGCTTAAGTACAAGATAATTCCCATAGAGGAAGGGTGGTATTACAGTGCCGCCTGTTGAAAGCCCCAATCCCGCTTCAATTGTTGAGCTCAG
";

    const FASTA_COUNT: &[u8] = &[
        0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 2, 2, 1, 0, 1, 1, 1, 2, 0, 0,
        0, 1, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 2, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2,
        0, 0, 1, 0, 0, 0, 0, 2, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0,
        1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 2, 0, 1, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 2, 1, 0, 0, 1, 1,
        0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 1, 0, 2, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 1, 1, 0, 1, 0, 0, 0,
        0, 1, 2, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 2, 1, 0, 0, 1, 1, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0,
        1, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1,
        1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 2, 0, 0, 1, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 1, 0, 0, 0, 0, 1, 0, 0, 1,
        0, 2, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0,
        1, 1,
    ];

    #[test]
    fn sequential() {
        let mut counter = Counter::<u8>::new(5);

        counter.count_fasta(Box::new(FASTA_FILE), 1);

        assert_eq!(*counter.get_raw(14), 2);

        assert_eq!(counter.get(cocktail::kmer::seq2bit(b"GTTCT")), 2);

        assert_eq!(
            counter.get(cocktail::kmer::canonical(
                cocktail::kmer::seq2bit(b"GTTCT"),
                5
            )),
            2
        );

        assert_eq!(&counter.raw()[..], &FASTA_COUNT[..]);
    }

    #[cfg(feature = "parallel")]
    #[test]
    fn parallel() {
        let mut counter = Counter::<std::sync::atomic::AtomicU8>::new(5);

        counter.count_fasta(Box::new(FASTA_FILE), 1);

        assert_eq!(
            counter
                .get_raw(14)
                .load(std::sync::atomic::Ordering::Relaxed),
            2
        );

        assert_eq!(counter.get(cocktail::kmer::seq2bit(b"GTTCT")), 2);

        assert_eq!(
            counter.get(cocktail::kmer::canonical(
                cocktail::kmer::seq2bit(b"GTTCT"),
                5
            )),
            2
        );

        assert_eq!(
            &unsafe { std::mem::transmute::<&[std::sync::atomic::AtomicU8], &[u8]>(counter.raw()) }
                [..],
            &FASTA_COUNT[..]
        );
    }
}
