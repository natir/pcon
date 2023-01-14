//! Run count command

/* std use */

use std::io::Read as _;

/* crate use */

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/* project use */
use crate::error;
use crate::serialize;

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

    #[cfg(test)]
    /// Get raw data
    pub(crate) fn raw_mut(&mut self) -> &mut [T] {
        &mut self.count
    }

    /// Convert counter in serializer
    pub fn serialize(self) -> serialize::Serialize<T> {
        serialize::Serialize::new(self)
    }
}

/*****************************/
/* sequential implementation */
/*****************************/
macro_rules! impl_sequential (
    ($type:ty, $convert:expr) => {
	impl Counter<$type> {
	    /// Create a new kmer Counter with kmer size equal to k
	    pub fn new(k: u8) -> Self {
		let tmp = vec![0u8; cocktail::kmer::get_hash_space_size(k) as usize * std::mem::size_of::<$type>()];

		Self {
		    k,
		    count: $convert(tmp.into_boxed_slice()),
		}
	    }

	    /// Create a new kmer by read a file
	    pub fn from_stream<R>(input: R) -> error::Result<Self>
		where R: std::io::Read
	    {
		let (k, data) = read_pcon_format::<R, $type>(input)?;

		#[allow(clippy::useless_transmute)]
		Ok(Self {
		    k,
		    count: $convert(data.into_boxed_slice()),
		})
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
	    pub(crate) fn inc(count: &mut [$type], index: usize) {
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

impl_sequential!(u8, |x| x);
impl_sequential!(u16, transmute::<u16>);
impl_sequential!(u32, transmute::<u32>);
impl_sequential!(u64, transmute::<u64>);
impl_sequential!(u128, transmute::<u128>);

/***************************/
/* parallel implementation */
/***************************/
#[cfg(feature = "parallel")]
macro_rules! impl_atomic (
    ($type:ty, $out_type:ty, $max:expr, $convert:expr) => {
	impl Counter<$type> {
	    /// Create a new kmer Counter with kmer size equal to k
	    pub fn new(k: u8) -> Self {
		let data = vec![0u8; cocktail::kmer::get_hash_space_size(k) as usize * std::mem::size_of::<$type>()];

		Self {
		    k,
		    count: $convert(data.into_boxed_slice()),
		}
	    }

	    /// Create a new kmer by read a file
	    pub fn from_stream<R>(input: R) -> error::Result<Self>
		where R: std::io::Read
	    {
		let (k, data) = read_pcon_format::<R, $type>(input)?;

		Ok(Self {
		    k,
		    count: $convert(data.into_boxed_slice()),
		})
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
	    pub(crate) fn inc(count: &[$type], index: usize) {
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

	    /// Get raw data in not atomic type
	    pub fn raw_noatomic(&self) -> &[$out_type] {
		unsafe {std::mem::transmute::<&[$type], &[$out_type]>(&self.count)}
	    }
	}

    }
);

#[cfg(feature = "parallel")]
impl_atomic!(
    std::sync::atomic::AtomicU8,
    u8,
    u8::MAX,
    transmute::<std::sync::atomic::AtomicU8>
);
#[cfg(feature = "parallel")]
impl_atomic!(
    std::sync::atomic::AtomicU16,
    u16,
    u16::MAX,
    transmute::<std::sync::atomic::AtomicU16>
);
#[cfg(feature = "parallel")]
impl_atomic!(
    std::sync::atomic::AtomicU32,
    u32,
    u32::MAX,
    transmute::<std::sync::atomic::AtomicU32>
);
#[cfg(feature = "parallel")]
impl_atomic!(
    std::sync::atomic::AtomicU64,
    u64,
    u64::MAX,
    transmute::<std::sync::atomic::AtomicU64>
);

/*******************/
/* utils function */
/******************/
/// Read data from pcon file
fn read_pcon_format<R, T>(mut reader: R) -> error::Result<(u8, Vec<u8>)>
where
    R: std::io::Read,
{
    let mut read_buffer = [0u8; 2];
    reader.read_exact(&mut read_buffer)?;
    let k = read_buffer[0];

    if std::mem::size_of::<T>() != read_buffer[1] as usize {
        return Err(error::Error::TypeNotMatch.into());
    }

    let mut deflate = flate2::read::MultiGzDecoder::new(reader);
    let mut tmp =
        vec![0u8; cocktail::kmer::get_hash_space_size(k) as usize * std::mem::size_of::<T>()];

    deflate.read_exact(&mut tmp)?;

    Ok((k, tmp))
}

/// Transmute from u8
fn transmute<T>(data: Box<[u8]>) -> Box<[T]>
where
    T: std::marker::Sized,
{
    unsafe { std::mem::transmute::<Box<[u8]>, Box<[T]>>(data) }
}

#[cfg(feature = "parallel")]
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

    #[test]
    fn sequential_from_stream() -> error::Result<()> {
        let mut file = vec![];

        let mut counter = Counter::<u8>::new(5);
        counter.count_fasta(Box::new(FASTA_FILE), 1);

        let serialize = counter.clone().serialize();
        serialize.pcon(std::io::Cursor::new(&mut file))?;

        let second_counter = Counter::<u8>::from_stream(&file[..])?;

        assert_eq!(counter.raw(), second_counter.raw());

        Ok(())
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

    #[cfg(feature = "parallel")]
    #[test]
    fn parallel_from_stream() -> error::Result<()> {
        let mut file = vec![];

        let mut counter = Counter::<std::sync::atomic::AtomicU8>::new(5);
        counter.count_fasta(Box::new(FASTA_FILE), 1);

        let mut truth = counter.raw_noatomic().to_vec();

        let serialize = counter.serialize();
        serialize.pcon(std::io::Cursor::new(&mut file))?;

        let second_counter = Counter::<std::sync::atomic::AtomicU8>::from_stream(&file[..])?;

        assert_eq!(truth, second_counter.raw_noatomic());

        Ok(())
    }
}
