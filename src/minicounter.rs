//! Generic struct of minicounter and implementation for many type

/* std use */

/* crate use */

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/* project use */
use crate::counter;
use crate::error;
#[cfg(feature = "parallel")]
use crate::utils;

/// A counter of kmer, count only if minimizer is present more than a threshold.
#[derive(Clone, Eq, PartialEq, Debug, Default)]
pub struct MiniCounter<T, U> {
    k: u64,
    threshold: U,
    mini_count: counter::Counter<T>,
    kmer_count: rustc_hash::FxHashMap<Vec<u8>, U>,
}

/**************************/
/* generic implementation */
/**************************/
impl<T, U> MiniCounter<T, U>
where
    U: std::fmt::Display + std::cmp::PartialOrd,
{
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
    pub fn kmer_raw(&self) -> &rustc_hash::FxHashMap<Vec<u8>, U> {
        &self.kmer_count
    }

    /// Write minicounter result in csv
    pub fn serialize<W>(&self, abundance: U, mut output: W) -> error::Result<()>
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

/*****************************/
/* sequential implementation */
/*****************************/
macro_rules! impl_sequential (
    ($type:ty) => {
	impl MiniCounter<$type, $type> {
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
			.and_modify(|c| *c = c.saturating_add(1))
			.or_insert(1);
		}
	    }

	    /// Get count of a kmer
	    pub fn get(&self, kmer: &[u8]) -> &$type {
		self.kmer_count.get(kmer).unwrap_or(&0)
	    }
	}
    }
);

impl_sequential!(u8);
impl_sequential!(u16);
impl_sequential!(u32);
impl_sequential!(u64);
impl_sequential!(u128);

#[cfg(feature = "parallel")]
macro_rules! impl_atomic(
    ($type:ty, $out_type:ty, $max:expr) => {
	impl MiniCounter<$type, $out_type> {
	    /// Create a new kmer MiniCounter with kmer size equal to k and minimizer equal to m
	    pub fn new(k: u64, m: u8, threshold: $out_type) -> Self {
		Self {
		    k,
		    threshold,
		    mini_count: counter::Counter::<$type>::new(m),
		    kmer_count: rustc_hash::FxHashMap::default(),
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
		    end = utils::populate_buffer(&mut iter, &mut records, record_buffer);
		    log::info!("End populate buffer {}", records.len());

		    let local = records.par_iter().map(|record| {
			let mut values = std::collections::HashMap::new();

			if record.sequence().len() >= self.k as usize {
			    let minimizer = cocktail::tokenizer::MiniBstr::new(
				record.sequence().as_ref(),
				self.k(),
				self.m(),
			    );

			    let mut prev_mini = None;
			    for (kmer, mini) in minimizer {
				if prev_mini != Some(mini) {
				    Self::mini_inc(&self.mini_count.count, (mini >> 1) as usize);
				}

				let normalize = kmer.to_ascii_uppercase();
				if self.mini_count.get(mini as u64) > self.threshold {

				    values.entry(normalize).and_modify(|c: &mut $out_type| *c = c.saturating_add(1)).or_insert(1);
				}
				prev_mini = Some(mini);
			    }
			}

			values
		    }).reduce(|| std::collections::HashMap::new(), |mut a, b| {
			for (k, v) in b.into_iter() {
			    a.entry(k.to_vec()).and_modify(|c: &mut $out_type| *c = c.saturating_add(v)).or_insert(v);
			}

			a
		    });

		    for (k, v) in local.into_iter() {
			self.kmer_count.entry(k.to_vec()).and_modify(|c: &mut $out_type| *c = c.saturating_add(v)).or_insert(v);
		    }
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

		    let local = records.par_iter().map(|record| {
			let mut values = std::collections::HashMap::new();

			if record.sequence().len() >= self.k as usize {
			    let minimizer = cocktail::tokenizer::MiniBstr::new(
				record.sequence().as_ref(),
				self.k(),
				self.m(),
			    );

			    let mut prev_mini = None;
			    for (kmer, mini) in minimizer {
				if prev_mini != Some(mini) {
				    Self::mini_inc(&self.mini_count.count, (mini >> 1) as usize);
				}

				let normalize = kmer.to_ascii_uppercase();
				if self.mini_count.get(mini as u64) > self.threshold {
				    values.entry(normalize).and_modify(|c: &mut $out_type| *c = c.saturating_add(1)).or_insert(1);
				}
				prev_mini = Some(mini);
			    }
			}

			values
		    }).reduce(|| std::collections::HashMap::new(), |mut a, b| {
			for (k, v) in b.into_iter() {
			    a.entry(k).and_modify(|c: &mut $out_type| *c = c.saturating_add(v)).or_insert(v);
			}

			a
		    });

		    for (k, v) in local.into_iter() {
			self.kmer_count.entry(k).and_modify(|c: &mut $out_type| *c = c.saturating_add(v)).or_insert(v);
		    }
		}
	    }

	    /// Increment value at index
	    pub(crate) fn mini_inc(count: &[$type], index: usize) {
		if count[index].load(std::sync::atomic::Ordering::SeqCst) != $max {
		    count[index].fetch_add(1, std::sync::atomic::Ordering::SeqCst);
		}
	    }

	    /// Get count of a kmer
	    pub fn get(&self, kmer: &[u8]) -> $out_type {
		*self.kmer_count.get(kmer).unwrap_or(&0)
	    }
	}
    }
);

#[cfg(feature = "parallel")]
impl_atomic!(std::sync::atomic::AtomicU8, u8, u8::MAX);
#[cfg(feature = "parallel")]
impl_atomic!(std::sync::atomic::AtomicU16, u16, u16::MAX);
#[cfg(feature = "parallel")]
impl_atomic!(std::sync::atomic::AtomicU32, u32, u32::MAX);
#[cfg(feature = "parallel")]
impl_atomic!(std::sync::atomic::AtomicU64, u64, u64::MAX);

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(feature = "parallel")]
    use crate::utils;

    const FASTA_FILE: &[u8] = b">random_seq 0
GTTCTGCAAATTAGAACAGACAATACACTGGCAGGCGTTGCGTTGGGGGAGATCttCCGTAACGAGCCGGCATTTGTAAGAAAGAGATTTCGAGTAAATG
>random_seq 1
AGGATAGAAGCTTAAGTACAAGATAATTCCCATAGAGGAAGGGTGGTAttACAGTGCCGCCTGTTGAAAGCCCCAATCCCGCTTCAATTGTTGAGCTCAG
>random_seq 2
AAGATAATTC";

    #[cfg(feature = "fastq")]
    const FASTQ_FILE: &[u8] = b"@random_seq 0
GTTCTGCAAATTAGAACAGACAATACACTGGCAGGCGTTGCGTTGGGGGAGATCttCCGTAACGAGCCGGCATTTGTAAGAAAGAGATTTCGAGTAAATG
+
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
@random_seq 1
AGGATAGAAGCTTAAGTACAAGATAATTCCCATAGAGGAAGGGTGGTAttACAGTGCCGCCTGTTGAAAGCCCCAATCCCGCTTCAATTGTTGAGCTCAG
+
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
@random_seq 2
AAGATAATTC
+
+
!!!!!!!!!!
";

    macro_rules! truth_count {
        ($type:ty, $name:ident) => {
            const $name: &[$type] = &[
                0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 1, 1,
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

    truth_count!(u8, TRUTH_COUNT_U8);
    truth_count!(u16, TRUTH_COUNT_U16);
    truth_count!(u32, TRUTH_COUNT_U32);
    truth_count!(u64, TRUTH_COUNT_U64);
    truth_count!(u128, TRUTH_COUNT_U128);

    macro_rules! sequential_fasta {
        ($type:ty, $name:ident, $truth:ident) => {
            #[test]
            fn $name() {
                let mut mini_count = MiniCounter::<$type, $type>::new(10, 5, 1);

                mini_count.count_fasta(Box::new(FASTA_FILE), 1);

                assert_eq!(mini_count.mini_raw(), $truth);

                let mut result = vec![];
                for (kmer, count) in mini_count.kmer_raw().iter() {
                    println!(
                        "b\"{}\".to_vec(), {}",
                        String::from_utf8(kmer.to_vec()).unwrap(),
                        count
                    );
                    result.push((kmer.to_vec(), *count));
                }
                result.sort();

                assert_eq!(
                    result,
                    vec![
                        (b"AAGATAATTC".to_vec(), 2),
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
                let mut mini_count = MiniCounter::<$type, $type>::new(10, 5, 1);

                mini_count.count_fastq(Box::new(FASTQ_FILE), 1);

                assert_eq!(mini_count.mini_raw(), $truth);

                let mut result = vec![];
                for (kmer, count) in mini_count.kmer_raw().iter() {
                    println!(
                        "b\"{}\".to_vec(), {}",
                        String::from_utf8(kmer.to_vec()).unwrap(),
                        count
                    );
                    result.push((kmer.to_vec(), *count));
                }
                result.sort();

                assert_eq!(
                    result,
                    vec![
                        (b"AAGATAATTC".to_vec(), 2),
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

    #[cfg(feature = "parallel")]
    macro_rules! parallel_fasta {
        ($type:ty, $out_type:ty, $name:ident, $truth:ident) => {
            #[test]
            fn $name() {
                let mut mini_count = MiniCounter::<$type, $out_type>::new(10, 5, 1);

                mini_count.count_fasta(Box::new(FASTA_FILE), 1);

                assert_eq!(
                    utils::transmute::<$type, $out_type>(mini_count.mini_raw()),
                    $truth
                );

                let mut result: Vec<(Vec<u8>, $out_type)> = vec![];
                for (kmer, count) in mini_count.kmer_raw().iter() {
                    result.push((kmer.to_vec(), *count));
                }
                result.sort();

                assert_eq!(
                    result,
                    vec![
                        (b"AAGATAATTC".to_vec(), 2),
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

    #[cfg(all(feature = "parallel", feature = "fastq"))]
    macro_rules! parallel_fastq {
        ($type:ty, $out_type:ty, $name:ident, $truth:ident) => {
            #[test]
            fn $name() {
                let mut mini_count = MiniCounter::<$type, $out_type>::new(10, 5, 1);

                mini_count.count_fastq(Box::new(FASTQ_FILE), 1);

                assert_eq!(
                    utils::transmute::<$type, $out_type>(mini_count.mini_raw()),
                    $truth
                );

                let mut result: Vec<(Vec<u8>, $out_type)> = vec![];
                for (kmer, count) in mini_count.kmer_raw().iter() {
                    result.push((kmer.to_vec(), *count));
                }
                result.sort();

                assert_eq!(
                    result,
                    vec![
                        (b"AAGATAATTC".to_vec(), 2),
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

    #[cfg(all(feature = "parallel", feature = "fastq"))]
    parallel_fastq!(
        std::sync::atomic::AtomicU8,
        u8,
        parallel_fastq_u8,
        TRUTH_COUNT_U8
    );
    #[cfg(all(feature = "parallel", feature = "fastq"))]
    parallel_fastq!(
        std::sync::atomic::AtomicU16,
        u16,
        parallel_fastq_u16,
        TRUTH_COUNT_U16
    );
    #[cfg(all(feature = "parallel", feature = "fastq"))]
    parallel_fastq!(
        std::sync::atomic::AtomicU32,
        u32,
        parallel_fastq_u32,
        TRUTH_COUNT_U32
    );
    #[cfg(all(feature = "parallel", feature = "fastq"))]
    parallel_fastq!(
        std::sync::atomic::AtomicU64,
        u64,
        parallel_fastq_u64,
        TRUTH_COUNT_U64
    );
}
