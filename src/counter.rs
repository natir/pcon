//! Generic struct of counter and implementation for many type

/* std use */

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
    ($type:ty, $init:expr, $read:expr) => {
	impl Counter<$type> {
	    /// Create a new kmer Counter with kmer size equal to k
	    pub fn new(k: u8) -> Self {
		let data: Box<[$type]> = $init(k, 0 as $type);
		Self {
		    k,
		    count: data,
		}
	    }

	    /// Create a new kmer by read a file
	    pub fn from_stream<R>(mut input: R) -> error::Result<Self>
		where R: std::io::Read
	    {
		let mut read_buffer = [0u8; 2];
		input.read_exact(&mut read_buffer)?;
		let k = read_buffer[0];

		if std::mem::size_of::<$type>() != read_buffer[1] as usize {
		    return Err(error::Error::TypeNotMatch.into());
		}

		let mut deflate = flate2::read::MultiGzDecoder::new(input);
		let mut data = $init(k, 0 as $type);

		$read(&mut deflate, &mut data)?;

		Ok(Self {
		    k,
		    count: data,
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

	    /// Perform count on fastq input
	    pub fn count_fastq(&mut self, fastq: Box<dyn std::io::BufRead>, _record_buffer: u64) {
		let mut reader = noodles::fastq::Reader::new(fastq);
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

impl_sequential!(u8, init_data, std::io::Read::read_exact);
impl_sequential!(
    u16,
    init_data,
    byteorder::ReadBytesExt::read_u16_into::<crate::ByteOrder>
);
impl_sequential!(
    u32,
    init_data,
    byteorder::ReadBytesExt::read_u32_into::<crate::ByteOrder>
);
impl_sequential!(
    u64,
    init_data,
    byteorder::ReadBytesExt::read_u64_into::<crate::ByteOrder>
);
impl_sequential!(
    u128,
    init_data,
    byteorder::ReadBytesExt::read_u128_into::<crate::ByteOrder>
);

/***************************/
/* parallel implementation */
/***************************/
#[cfg(feature = "parallel")]
macro_rules! impl_atomic (
    ($type:ty, $out_type:ty, $max:expr, $init:expr, $read:expr) => {
	impl Counter<$type> {
	    /// Create a new kmer Counter with kmer size equal to k
	    pub fn new(k: u8) -> Self {
		Self {
		    k,
		    count: transmute_box($init(k, 0 as $out_type)),
		}
	    }

	    /// Create a new kmer by read a file
	    pub fn from_stream<R>(mut input: R) -> error::Result<Self>
		where R: std::io::Read
	    {
		let mut read_buffer = [0u8; 2];
		input.read_exact(&mut read_buffer)?;
		let k = read_buffer[0];

		if std::mem::size_of::<$type>() != read_buffer[1] as usize {
		    return Err(error::Error::TypeNotMatch.into());
		}

		let mut deflate = flate2::read::MultiGzDecoder::new(input);
		let mut data = $init(k, 0 as $out_type);

		$read(&mut deflate, &mut data)?;

		Ok(Self {
		    k,
		    count: transmute_box(data),
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


	    /// Perform count on fastq input
	    pub fn count_fastq(&mut self, fastq: Box<dyn std::io::BufRead>, record_buffer: u64) {
		let mut reader = noodles::fastq::Reader::new(fastq);
		let mut iter = reader.records();
		let mut records = Vec::with_capacity(record_buffer as usize);

		let mut end = true;
		while end {
		    log::info!("Start populate buffer");
		    end = populate_bufferq(&mut iter, &mut records, record_buffer);
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

	    /// Get raw data in no atomic type
	    pub fn raw_noatomic(&self) -> &[$out_type] {
		transmute(&self.count)
	    }
	}

    }
);

#[cfg(feature = "parallel")]
impl_atomic!(
    std::sync::atomic::AtomicU8,
    u8,
    u8::MAX,
    init_data,
    std::io::Read::read_exact
);
#[cfg(feature = "parallel")]
impl_atomic!(
    std::sync::atomic::AtomicU16,
    u16,
    u16::MAX,
    init_data,
    byteorder::ReadBytesExt::read_u16_into::<crate::ByteOrder>
);
#[cfg(feature = "parallel")]
impl_atomic!(
    std::sync::atomic::AtomicU32,
    u32,
    u32::MAX,
    init_data,
    byteorder::ReadBytesExt::read_u32_into::<crate::ByteOrder>
);
#[cfg(feature = "parallel")]
impl_atomic!(
    std::sync::atomic::AtomicU64,
    u64,
    u64::MAX,
    init_data,
    byteorder::ReadBytesExt::read_u64_into::<crate::ByteOrder>
);

/*******************/
/* utils function */
/******************/
/// Initialize counter
fn init_data<T>(k: u8, value: T) -> Box<[T]>
where
    T: std::marker::Sized + std::clone::Clone,
{
    vec![value; cocktail::kmer::get_hash_space_size(k) as usize].into_boxed_slice()
}

#[cfg(feature = "parallel")]
/// Perform transmutation on box
fn transmute<I, O>(data: &[I]) -> &[O]
where
    I: std::marker::Sized,
    O: std::marker::Sized,
{
    unsafe { std::mem::transmute::<&[I], &[O]>(data) }
}

#[cfg(feature = "parallel")]
/// Perform transmutation on box
fn transmute_box<I, O>(data: Box<[I]>) -> Box<[O]>
where
    I: std::marker::Sized,
    O: std::marker::Sized,
{
    unsafe { std::mem::transmute::<Box<[I]>, Box<[O]>>(data) }
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

#[cfg(feature = "parallel")]
/// Populate record buffer with content of iterator
fn populate_bufferq(
    iter: &mut noodles::fastq::reader::Records<'_, Box<dyn std::io::BufRead>>,
    records: &mut Vec<noodles::fastq::Record>,
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

    #[cfg(feature = "parallel")]
    /// Perform transmutation on box
    fn transmute<I, O>(data: &[I]) -> &[O]
    where
        I: std::marker::Sized,
        O: std::marker::Sized,
    {
        unsafe { std::mem::transmute::<&[I], &[O]>(data) }
    }

    const FASTA_FILE: &[u8] = b">random_seq 0
GTTCTGCAAATTAGAACAGACAATACACTGGCAGGCGTTGCGTTGGGGGAGATCTTCCGTAACGAGCCGGCATTTGTAAGAAAGAGATTTCGAGTAAATG
>random_seq 1
AGGATAGAAGCTTAAGTACAAGATAATTCCCATAGAGGAAGGGTGGTATTACAGTGCCGCCTGTTGAAAGCCCCAATCCCGCTTCAATTGTTGAGCTCAG
";

    macro_rules! fasta_count {
        ($type:ty, $name:ident) => {
            const $name: &[$type] = &[
                0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 2, 2, 1, 0, 1, 1, 1, 2,
                0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 0, 0, 0, 1, 1, 0,
                1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 2, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0,
                1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 1, 0, 0, 0, 0, 0, 1, 0,
                0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 2, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                0, 1, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 2, 0, 0,
                0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 2, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,
                0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 2, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 2, 0, 0,
                1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 2, 1, 0, 0,
                1, 0, 1, 1, 0, 0, 1, 1, 2, 1, 0, 0, 1, 1, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 1, 1, 0, 0,
                0, 0, 0, 0, 2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 1, 1,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 2, 0, 0, 1, 0, 1, 0, 0, 0,
                0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 1, 0, 0, 0, 0,
                1, 0, 0, 1, 0, 2, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
                1, 1, 0, 1, 0, 0, 1, 1,
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
                let mut counter = Counter::<$type>::new(5);

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

                assert_eq!(&counter.raw()[..], &$truth[..]);
            }
        };
    }

    sequential!(u8, sequential_u8, FASTA_COUNT_U8);
    sequential!(u16, sequential_u16, FASTA_COUNT_U16);
    sequential!(u32, sequential_u32, FASTA_COUNT_U32);
    sequential!(u64, sequential_u64, FASTA_COUNT_U64);
    sequential!(u128, sequential_u128, FASTA_COUNT_U128);

    macro_rules! sequential_serialize {
        ($type:ty, $failled_type:ty, $name:ident, $failled_name:ident) => {
            #[test]
            fn $name() -> error::Result<()> {
                let mut file = vec![];

                let mut counter = Counter::<$type>::new(5);
                counter.count_fasta(Box::new(FASTA_FILE), 1);

                let serialize = counter.clone().serialize();
                serialize.pcon(std::io::Cursor::new(&mut file))?;

                let second_counter = Counter::<$type>::from_stream(&file[..])?;

                assert_eq!(counter.raw(), second_counter.raw());

                Ok(())
            }

            #[test]
            fn $failled_name() -> error::Result<()> {
                let mut file = vec![];

                let mut counter = Counter::<$type>::new(5);
                counter.count_fasta(Box::new(FASTA_FILE), 1);

                let serialize = counter.serialize();
                serialize.pcon(std::io::Cursor::new(&mut file))?;

                let second_counter = Counter::<$failled_type>::from_stream(&file[..]);

                assert!(second_counter.is_err());

                Ok(())
            }
        };
    }

    sequential_serialize!(
        u8,
        u16,
        sequential_serialize_u8,
        failled_sequential_serialize_u8
    );
    sequential_serialize!(
        u16,
        u8,
        sequential_serialize_u16,
        failled_sequential_serialize_u16
    );
    sequential_serialize!(
        u32,
        u8,
        sequential_serialize_u32,
        failled_sequential_serialize_u32
    );
    sequential_serialize!(
        u64,
        u8,
        sequential_serialize_u64,
        failled_sequential_serialize_u64
    );
    sequential_serialize!(
        u128,
        u8,
        sequential_serialize_u128,
        failled_sequential_serialize_u128
    );

    #[cfg(feature = "parallel")]
    macro_rules! parallel {
        ($type:ty, $out_type:ty, $name:ident, $truth:ident) => {
            #[test]
            fn $name() {
                let mut counter = Counter::<$type>::new(5);

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

                assert_eq!(counter.raw_noatomic(), &$truth[..]);
            }
        };
    }

    #[cfg(feature = "parallel")]
    parallel!(std::sync::atomic::AtomicU8, u8, parallel_u8, FASTA_COUNT_U8);
    #[cfg(feature = "parallel")]
    parallel!(
        std::sync::atomic::AtomicU16,
        u16,
        parallel_u16,
        FASTA_COUNT_U16
    );
    #[cfg(feature = "parallel")]
    parallel!(
        std::sync::atomic::AtomicU32,
        u32,
        parallel_u32,
        FASTA_COUNT_U32
    );
    #[cfg(feature = "parallel")]
    parallel!(
        std::sync::atomic::AtomicU64,
        u64,
        parallel_u64,
        FASTA_COUNT_U64
    );

    #[cfg(feature = "parallel")]
    macro_rules! parallel_serialize {
        ($type:ty, $out_type:ty, $failled_type:ty, $name:ident, $failled_name:ident) => {
            #[test]
            fn $name() -> error::Result<()> {
                let mut file = vec![];

                {
                    let mut counter = Counter::<$type>::new(5);
                    counter.count_fasta(Box::new(FASTA_FILE), 1);

                    let serialize = counter.serialize();
                    serialize.pcon(std::io::Cursor::new(&mut file))?;
                }

                let mut counter = Counter::<$type>::new(5);
                counter.count_fasta(Box::new(FASTA_FILE), 1);

                let second_counter = Counter::<$type>::from_stream(&file[..])?;

                assert_eq!(
                    &transmute::<$type, $out_type>(counter.raw())[..],
                    &transmute::<$type, $out_type>(second_counter.raw())[..],
                );

                Ok(())
            }

            #[test]
            fn $failled_name() -> error::Result<()> {
                let mut file = vec![];

                let mut counter = Counter::<$type>::new(5);
                counter.count_fasta(Box::new(FASTA_FILE), 1);

                let serialize = counter.serialize();
                serialize.pcon(std::io::Cursor::new(&mut file))?;

                let second_counter = Counter::<$failled_type>::from_stream(&file[..]);

                assert!(second_counter.is_err());

                Ok(())
            }
        };
    }

    #[cfg(feature = "parallel")]
    parallel_serialize!(
        std::sync::atomic::AtomicU8,
        u8,
        std::sync::atomic::AtomicU16,
        parallel_serialize_u8,
        failled_parallel_serialize_u8
    );
    #[cfg(feature = "parallel")]
    parallel_serialize!(
        std::sync::atomic::AtomicU16,
        u16,
        std::sync::atomic::AtomicU8,
        parallel_serialize_u16,
        failled_parallel_serialize_u16
    );
    #[cfg(feature = "parallel")]
    parallel_serialize!(
        std::sync::atomic::AtomicU32,
        u32,
        std::sync::atomic::AtomicU8,
        parallel_serialize_u32,
        failled_parallel_serialize_u32
    );
    #[cfg(feature = "parallel")]
    parallel_serialize!(
        std::sync::atomic::AtomicU64,
        u64,
        std::sync::atomic::AtomicU8,
        parallel_serialize_u64,
        failled_parallel_serialize_u64
    );
}
