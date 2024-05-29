//! Generic struct of counter and implementation for many type

/* std use */

/* crate use */

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/* project use */
use crate::error;
use crate::serialize;
use crate::utils;

/// A counter of kmer based on cocktail crate 2bit conversion, canonicalisation and hashing.
/// Implement only for u8, std::sync::atomic::AtomicU8
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default)]
pub struct Counter<T> {
    k: u8,
    pub(crate) count: Box<[T]>,
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

    #[allow(dead_code)]
    /// Get raw data mut
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

	    #[cfg(feature = "fastq")]
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

impl_sequential!(u8, utils::init_data, std::io::Read::read_exact);
impl_sequential!(
    u16,
    utils::init_data,
    byteorder::ReadBytesExt::read_u16_into::<crate::ByteOrder>
);
impl_sequential!(
    u32,
    utils::init_data,
    byteorder::ReadBytesExt::read_u32_into::<crate::ByteOrder>
);
impl_sequential!(
    u64,
    utils::init_data,
    byteorder::ReadBytesExt::read_u64_into::<crate::ByteOrder>
);
impl_sequential!(
    u128,
    utils::init_data,
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
		    count: utils::transmute_box($init(k, 0 as $out_type)),
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
		    count: utils::transmute_box(data),
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
		    end = utils::populate_buffer(&mut iter, &mut records, record_buffer);
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

	    #[cfg(feature = "fastq")]
	    /// Perform count on fastq input
	    pub fn count_fastq(&mut self, fastq: Box<dyn std::io::BufRead>, record_buffer: u64) {
		let mut reader = noodles::fastq::Reader::new(fastq);
		let mut iter = reader.records();
		let mut records = Vec::with_capacity(record_buffer as usize);

		let mut end = true;
		while end {
		    log::info!("Start populate buffer");
		    end = utils::populate_bufferq(&mut iter, &mut records, record_buffer);
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
		utils::transmute(&self.count)
	    }
	}

    }
);

#[cfg(feature = "parallel")]
impl_atomic!(
    std::sync::atomic::AtomicU8,
    u8,
    u8::MAX,
    utils::init_data,
    std::io::Read::read_exact
);
#[cfg(feature = "parallel")]
impl_atomic!(
    std::sync::atomic::AtomicU16,
    u16,
    u16::MAX,
    utils::init_data,
    byteorder::ReadBytesExt::read_u16_into::<crate::ByteOrder>
);
#[cfg(feature = "parallel")]
impl_atomic!(
    std::sync::atomic::AtomicU32,
    u32,
    u32::MAX,
    utils::init_data,
    byteorder::ReadBytesExt::read_u32_into::<crate::ByteOrder>
);
#[cfg(feature = "parallel")]
impl_atomic!(
    std::sync::atomic::AtomicU64,
    u64,
    u64::MAX,
    utils::init_data,
    byteorder::ReadBytesExt::read_u64_into::<crate::ByteOrder>
);

#[cfg(test)]
mod tests {
    use super::*;

    const FASTA_FILE: &[u8] = b">random_seq 0
GTTCTGCAAATTAGAACAGACAATACACTGGCAGGCGTTGCGTTGGGGGAGATCTTCCGTAACGAGCCGGCATTTGTAAGAAAGAGATTTCGAGTAAATG
>random_seq 1
AGGATAGAAGCTTAAGTACAAGATAATTCCCATAGAGGAAGGGTGGTATTACAGTGCCGCCTGTTGAAAGCCCCAATCCCGCTTCAATTGTTGAGCTCAG
";

    #[cfg(feature = "fastq")]
    const FASTQ_FILE: &[u8] = b"@random_seq 0
GTTCTGCAAATTAGAACAGACAATACACTGGCAGGCGTTGCGTTGGGGGAGATCttCCGTAACGAGCCGGCATTTGTAAGAAAGAGATTTCGAGTAAATG
+
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
@random_seq 1
AGGATAGAAGCTTAAGTACAAGATAATTCCCATAGAGGAAGGGTGGTAttACAGTGCCGCCTGTTGAAAGCCCCAATCCCGCTTCAATTGTTGAGCTCAG
+
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
";

    macro_rules! truth_count {
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

    truth_count!(u8, TRUTH_COUNT_U8);
    truth_count!(u16, TRUTH_COUNT_U16);
    truth_count!(u32, TRUTH_COUNT_U32);
    truth_count!(u64, TRUTH_COUNT_U64);
    truth_count!(u128, TRUTH_COUNT_U128);

    macro_rules! sequential_fasta {
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

    sequential_fasta!(u8, sequential_fasta_u8, TRUTH_COUNT_U8);
    sequential_fasta!(u16, sequential_fasta_u16, TRUTH_COUNT_U16);
    sequential_fasta!(u32, sequential_fasta_u32, TRUTH_COUNT_U32);
    sequential_fasta!(u64, sequential_fasta_u64, TRUTH_COUNT_U64);
    sequential_fasta!(u128, sequential_fasta_u128, TRUTH_COUNT_U128);

    #[cfg(feature = "fastq")]
    macro_rules! sequential_fastq {
        ($type:ty, $name:ident, $truth:ident) => {
            #[test]
            fn $name() {
                let mut counter = Counter::<$type>::new(5);

                counter.count_fastq(Box::new(FASTQ_FILE), 1);

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

    #[cfg(feature = "fastq")]
    sequential_fastq!(u8, sequential_fastq_u8, TRUTH_COUNT_U8);
    #[cfg(feature = "fastq")]
    sequential_fastq!(u16, sequential_fastq_u16, TRUTH_COUNT_U16);
    #[cfg(feature = "fastq")]
    sequential_fastq!(u32, sequential_fastq_u32, TRUTH_COUNT_U32);
    #[cfg(feature = "fastq")]
    sequential_fastq!(u64, sequential_fastq_u64, TRUTH_COUNT_U64);
    #[cfg(feature = "fastq")]
    sequential_fastq!(u128, sequential_fastq_u128, TRUTH_COUNT_U128);

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
    macro_rules! parallel_fasta {
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
    parallel_fasta!(
        std::sync::atomic::AtomicU8,
        u8,
        parallel_fasta_u8,
        TRUTH_COUNT_U8
    );
    #[cfg(feature = "parallel")]
    parallel_fasta!(
        std::sync::atomic::AtomicU16,
        u16,
        parallel_fasta_u16,
        TRUTH_COUNT_U16
    );
    #[cfg(feature = "parallel")]
    parallel_fasta!(
        std::sync::atomic::AtomicU32,
        u32,
        parallel_fasta_u32,
        TRUTH_COUNT_U32
    );
    #[cfg(feature = "parallel")]
    parallel_fasta!(
        std::sync::atomic::AtomicU64,
        u64,
        parallel_fasta_u64,
        TRUTH_COUNT_U64
    );

    #[cfg(feature = "parallel")]
    macro_rules! parallel_fastq {
        ($type:ty, $out_type:ty, $name:ident, $truth:ident) => {
            #[test]
            fn $name() {
                let mut counter = Counter::<$type>::new(5);

                counter.count_fastq(Box::new(FASTQ_FILE), 1);

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
    parallel_fastq!(
        std::sync::atomic::AtomicU8,
        u8,
        parallel_fastq_u8,
        TRUTH_COUNT_U8
    );
    #[cfg(feature = "parallel")]
    parallel_fastq!(
        std::sync::atomic::AtomicU16,
        u16,
        parallel_fastq_u16,
        TRUTH_COUNT_U16
    );
    #[cfg(feature = "parallel")]
    parallel_fastq!(
        std::sync::atomic::AtomicU32,
        u32,
        parallel_fastq_u32,
        TRUTH_COUNT_U32
    );
    #[cfg(feature = "parallel")]
    parallel_fastq!(
        std::sync::atomic::AtomicU64,
        u64,
        parallel_fastq_u64,
        TRUTH_COUNT_U64
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
                    &utils::transmute::<$type, $out_type>(counter.raw())[..],
                    &utils::transmute::<$type, $out_type>(second_counter.raw())[..],
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
