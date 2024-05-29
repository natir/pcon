//! Tools to serialize a Counter

/* std use */

use std::io::Write as _;

/* crate use */
use byteorder::WriteBytesExt as _;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[cfg(feature = "kff")]
use kff;

/* project use */
use crate::counter;
use crate::error;
use crate::solid;
#[cfg(feature = "parallel")]
use crate::utils;

/// Struct to serialize counter
pub struct Serialize<T> {
    counter: counter::Counter<T>,
}

impl<T> Serialize<T> {
    /// Create a new Serialize from a Counter
    pub fn new(counter: counter::Counter<T>) -> Self {
        Self { counter }
    }
}

macro_rules! impl_sequential {
    ($type:ty) => {
        impl Serialize<$type> {
            /// Write counter in pcon format
            pub fn pcon<W>(&self, mut output: W) -> error::Result<()>
            where
                W: std::io::Write,
            {
                output.write_all(&[self.counter.k(), std::mem::size_of::<$type>() as u8])?;

                // Magic number choose empirically
                let chunk_size = (1 << 21) / std::mem::size_of::<$type>();

                let compress_block: Vec<error::Result<Vec<u8>>> = self
                    .counter
                    .raw()
                    .chunks(chunk_size)
                    .map(|input_buffer| {
                        input_buffer
                            .iter()
                            .map(|x| x.to_le_bytes())
                            .flatten()
                            .collect::<Vec<u8>>()
                    })
                    .map(|input_buffer| {
                        let mut output_buffer = Vec::with_capacity(1 << 25);

                        {
                            let mut encoder = flate2::write::GzEncoder::new(
                                &mut output_buffer,
                                flate2::Compression::fast(),
                            );
                            encoder.write_all(&input_buffer)?;
                        }

                        Ok(output_buffer)
                    })
                    .collect();

                for result in compress_block {
                    output.write_all(&result?)?;
                }

                Ok(())
            }

            /// Write kmer count in csv format
            pub fn csv<W>(&self, abundance: $type, mut output: W) -> error::Result<()>
            where
                W: std::io::Write,
            {
                let counts = self.counter.raw();

                for (hash, value) in counts.iter().enumerate() {
                    let kmer = if cocktail::kmer::parity_even(hash as u64) {
                        cocktail::kmer::kmer2seq((hash as u64) << 1, self.counter.k())
                    } else {
                        cocktail::kmer::kmer2seq(((hash as u64) << 1) ^ 0b1, self.counter.k())
                    };

                    if value > &abundance {
                        writeln!(output, "{},{}", kmer, value)?;
                    }
                }

                Ok(())
            }

            /// Convert counter in solid and write it
            ///
            /// The first bytes contains the size of k the rest of the file are a
            /// bitfield of absence for each kmer
            pub fn solid<W>(&self, abundance: $type, output: W) -> error::Result<()>
            where
                W: std::io::Write,
            {
                let solid =
                    solid::Solid::from_count(self.counter.k(), self.counter.raw(), abundance);

                let mut writer = niffler::get_writer(
                    Box::new(output),
                    niffler::compression::Format::Gzip,
                    niffler::compression::Level::One,
                )?;

                writer.write_u8(solid.k())?;

                writer.write_all(solid.get_raw_solid().as_raw_slice())?;

                Ok(())
            }

            #[cfg(feature = "kff")]
            /// Write kmer count in kff format
            pub fn kff<W>(&self, abundance: $type, output: W) -> error::Result<()>
            where
                W: std::io::Write,
            {
                let header = kff::section::Header::new(
                    1,
                    0,
                    0b00011110,
                    true,
                    true,
                    b"producer: pcon".to_vec(),
                )?;
                let mut writer = kff::Kff::write(output, header)?;
                let mut values = kff::section::Values::default();
                values.insert("k".to_string(), self.counter.k() as u64);
                values.insert("ordered".to_string(), true as u64);
                values.insert("max".to_string(), <$type>::MAX as u64);
                values.insert("data_size".to_string(), std::mem::size_of::<$type>() as u64);

                writer.write_values(values.clone())?;

                let mut kmers = vec![];

                for (hash, value) in self.counter.raw().iter().enumerate() {
                    if value > &abundance {
                        let kmer = if cocktail::kmer::parity_even(hash as u64) {
                            ((hash << 1) | 0b0) as u64
                        } else {
                            ((hash << 1) | 0b1) as u64
                        };

                        kmers.push(kff::section::Block::new(
                            self.counter.k() as u64,
                            std::mem::size_of::<crate::CountType>(),
                            kff::Kmer::new(
                                bitvec::boxed::BitBox::<u8, bitvec::order::Msb0>::from_boxed_slice(
                                    Box::new(kmer.to_be_bytes()),
                                ),
                                value.to_be_bytes().to_vec(),
                            ),
                            0,
                        ));
                    }
                }

                writer.write_raw(kff::section::Raw::new(&values)?, kmers)?;

                writer.finalize()?;

                Ok(())
            }
        }
    };
}

impl_sequential!(u8);
impl_sequential!(u16);
impl_sequential!(u32);
impl_sequential!(u64);
impl_sequential!(u128);

#[cfg(feature = "parallel")]
macro_rules! impl_atomic {
    ($type:ty, $out_type:ty) => {
        impl Serialize<$type> {
            /// Write counter in pcon format
            pub fn pcon<W>(&self, mut output: W) -> error::Result<()>
            where
                W: std::io::Write,
            {
                output.write_all(&[self.counter.k(), std::mem::size_of::<$type>() as u8])?;

                // Magic number choose empirically
                let chunk_size = (1 << 21) / std::mem::size_of::<$type>();

                let count = utils::transmute::<$type, $out_type>(self.counter.raw());

                let compress_block: Vec<error::Result<Vec<u8>>> = count
                    .par_chunks(chunk_size)
                    .map(|input_buffer| {
                        input_buffer
                            .iter()
                            .map(|x| x.to_le_bytes())
                            .flatten()
                            .collect::<Vec<u8>>()
                    })
                    .map(|input_buffer| {
                        let mut output_buffer = Vec::with_capacity(1 << 25);

                        {
                            let mut encoder = flate2::write::GzEncoder::new(
                                &mut output_buffer,
                                flate2::Compression::fast(),
                            );
                            encoder.write_all(&input_buffer)?;
                        }

                        Ok(output_buffer)
                    })
                    .collect();

                for result in compress_block {
                    output.write_all(&result?)?;
                }

                Ok(())
            }

            /// Write kmer in csv format
            pub fn csv<W>(&self, abundance: $out_type, mut output: W) -> error::Result<()>
            where
                W: std::io::Write,
            {
                let counts = utils::transmute::<$type, $out_type>(self.counter.raw());

                for (hash, value) in counts.iter().enumerate() {
                    let kmer = if cocktail::kmer::parity_even(hash as u64) {
                        cocktail::kmer::kmer2seq((hash as u64) << 1, self.counter.k())
                    } else {
                        cocktail::kmer::kmer2seq(((hash as u64) << 1) ^ 0b1, self.counter.k())
                    };

                    if value > &abundance {
                        writeln!(output, "{},{}", kmer, value)?;
                    }
                }

                Ok(())
            }

            /// Convert counter in solid and write it
            ///
            /// The first bytes contains the size of k the rest of the file and a
            /// bitfield of absence for each kmer
            pub fn solid<W>(&self, abundance: $out_type, output: W) -> error::Result<()>
            where
                W: std::io::Write,
            {
                let solid = solid::Solid::from_count(
                    self.counter.k(),
                    utils::transmute::<$type, $out_type>(self.counter.raw()),
                    abundance,
                );

                let mut writer = niffler::get_writer(
                    Box::new(output),
                    niffler::compression::Format::Gzip,
                    niffler::compression::Level::One,
                )?;

                writer.write_all(&[solid.k()])?;

                writer.write_all(solid.get_raw_solid().as_raw_slice())?;

                Ok(())
            }

            #[cfg(feature = "kff")]
            /// Write kmer count in kff format
            pub fn kff<W>(&self, abundance: $out_type, output: W) -> error::Result<()>
            where
                W: std::io::Write,
            {
                let header = kff::section::Header::new(
                    1,
                    0,
                    0b00011110,
                    true,
                    true,
                    b"producer: pcon".to_vec(),
                )?;
                let mut writer = kff::Kff::write(output, header)?;
                let mut values = kff::section::Values::default();
                values.insert("k".to_string(), self.counter.k() as u64);
                values.insert("ordered".to_string(), true as u64);
                values.insert("max".to_string(), <$out_type>::MAX as u64);
                values.insert(
                    "data_size".to_string(),
                    std::mem::size_of::<$out_type>() as u64,
                );

                writer.write_values(values.clone())?;

                let mut kmers = vec![];

                let counts = utils::transmute::<$type, $out_type>(self.counter.raw());
                for (hash, value) in counts.iter().enumerate() {
                    if value > &abundance {
                        let kmer = if cocktail::kmer::parity_even(hash as u64) {
                            ((hash << 1) | 0b0) as u64
                        } else {
                            ((hash << 1) | 0b1) as u64
                        };

                        kmers.push(kff::section::Block::new(
                            self.counter.k() as u64,
                            std::mem::size_of::<crate::CountType>(),
                            kff::Kmer::new(
                                bitvec::boxed::BitBox::<u8, bitvec::order::Msb0>::from_boxed_slice(
                                    Box::new(kmer.to_be_bytes()),
                                ),
                                value.to_be_bytes().to_vec(),
                            ),
                            0,
                        ));
                    }
                }

                writer.write_raw(kff::section::Raw::new(&values)?, kmers)?;

                writer.finalize()?;

                Ok(())
            }
        }
    };
}

#[cfg(feature = "parallel")]
impl_atomic!(std::sync::atomic::AtomicU8, u8);
#[cfg(feature = "parallel")]
impl_atomic!(std::sync::atomic::AtomicU16, u16);
#[cfg(feature = "parallel")]
impl_atomic!(std::sync::atomic::AtomicU32, u32);
#[cfg(feature = "parallel")]
impl_atomic!(std::sync::atomic::AtomicU64, u64);

#[cfg(test)]
mod tests {
    use super::*;

    fn generate_counter() -> counter::Counter<u8> {
        let mut counter: counter::Counter<u8> = counter::Counter::<u8>::new(5);

        for i in 0..cocktail::kmer::get_kmer_space_size(5) {
            counter::Counter::<u8>::inc(
                counter.raw_mut(),
                (cocktail::kmer::canonical(i, 5) >> 1) as usize,
            );
        }

        counter::Counter::<u8>::inc(counter.raw_mut(), 0);

        counter
    }

    #[cfg(feature = "parallel")]
    fn generate_atomic_counter() -> counter::Counter<std::sync::atomic::AtomicU8> {
        let mut counter: counter::Counter<std::sync::atomic::AtomicU8> =
            counter::Counter::<std::sync::atomic::AtomicU8>::new(5);

        for i in 0..cocktail::kmer::get_kmer_space_size(5) {
            counter::Counter::<std::sync::atomic::AtomicU8>::inc(
                counter.raw_mut(),
                (cocktail::kmer::canonical(i, 5) >> 1) as usize,
            );
        }

        counter::Counter::<std::sync::atomic::AtomicU8>::inc(counter.raw_mut(), 0);

        counter
    }

    const PCON_ABUNDANCE: &[u8] = &[
        5, 1, 31, 139, 8, 0, 0, 0, 0, 0, 4, 255, 237, 208, 1, 13, 0, 0, 0, 130, 176, 77, 251, 119,
        38, 8, 60, 194, 191, 152, 7, 0, 94, 201, 71, 192, 0, 2, 0, 0,
    ];

    #[test]
    fn pcon() -> error::Result<()> {
        let mut outfile = Vec::new();
        let counter = generate_counter();
        let serialize = counter.serialize();

        serialize.pcon(&mut outfile)?;
        assert_eq!(&outfile[..], &PCON_ABUNDANCE[..]);

        Ok(())
    }

    #[cfg(feature = "parallel")]
    #[test]
    fn atomic_pcon() -> error::Result<()> {
        let mut outfile = Vec::new();
        let counter = generate_atomic_counter();
        let serialize = counter.serialize();

        serialize.pcon(&mut outfile)?;
        assert_eq!(&outfile[..], &PCON_ABUNDANCE[..]);

        Ok(())
    }

    const CSV_ABUNDANCE_MIN_1: &[u8] = b"AAAAA,3\nAAAAG,2\nAAACC,2\nAAACT,2\nAAATC,2\nAAATT,2\nAAAGA,2\nAAAGG,2\nAACAC,2\nAACAT,2\nAACCA,2\nAACCG,2\nAACTA,2\nAACTG,2\nAACGC,2\nAACGT,2\nAATAC,2\nAATAT,2\nAATCA,2\nAATCG,2\nAATTA,2\nAATTG,2\nAATGC,2\nAATGT,2\nAAGAA,2\nAAGAG,2\nAAGCC,2\nAAGCT,2\nAAGTC,2\nAAGTT,2\nAAGGA,2\nAAGGG,2\nACAAC,2\nACAAT,2\nACACA,2\nACACG,2\nACATA,2\nACATG,2\nACAGC,2\nACAGT,2\nACCAA,2\nACCAG,2\nACCCC,2\nACCCT,2\nACCTC,2\nACCTT,2\nACCGA,2\nACCGG,2\nACTAA,2\nACTAG,2\nACTCC,2\nACTCT,2\nACTTC,2\nACTTT,2\nACTGA,2\nACTGG,2\nACGAC,2\nACGAT,2\nACGCA,2\nACGCG,2\nACGTA,2\nACGTG,2\nACGGC,2\nACGGT,2\nATAAC,2\nATAAT,2\nATACA,2\nATACG,2\nATATA,2\nATATG,2\nATAGC,2\nATAGT,2\nATCAA,2\nATCAG,2\nATCCC,2\nATCCT,2\nATCTC,2\nATCTT,2\nATCGA,2\nATCGG,2\nATTAA,2\nATTAG,2\nATTCC,2\nATTCT,2\nATTTC,2\nATTTT,2\nATTGA,2\nATTGG,2\nATGAC,2\nATGAT,2\nATGCA,2\nATGCG,2\nATGTA,2\nATGTG,2\nATGGC,2\nATGGT,2\nAGAAA,2\nAGAAG,2\nAGACC,2\nAGACT,2\nAGATC,2\nAGATT,2\nAGAGA,2\nAGAGG,2\nAGCAC,2\nAGCAT,2\nAGCCA,2\nAGCCG,2\nAGCTA,2\nAGCTG,2\nAGCGC,2\nAGCGT,2\nAGTAC,2\nAGTAT,2\nAGTCA,2\nAGTCG,2\nAGTTA,2\nAGTTG,2\nAGTGC,2\nAGTGT,2\nAGGAA,2\nAGGAG,2\nAGGCC,2\nAGGCT,2\nAGGTC,2\nAGGTT,2\nAGGGA,2\nAGGGG,2\nCAAAC,2\nCAAAT,2\nCAACA,2\nCAACG,2\nCAATA,2\nCAATG,2\nCAAGC,2\nCAAGT,2\nCACAA,2\nCACAG,2\nCACCC,2\nCACCT,2\nCACTC,2\nCACTT,2\nCACGA,2\nCACGG,2\nCATAA,2\nCATAG,2\nCATCC,2\nCATCT,2\nCATTC,2\nCATTT,2\nCATGA,2\nCATGG,2\nCAGAC,2\nCAGAT,2\nCAGCA,2\nCAGCG,2\nCAGTA,2\nCAGTG,2\nCAGGC,2\nCAGGT,2\nCCAAA,2\nCCAAG,2\nCCACC,2\nCCACT,2\nCCATC,2\nCCATT,2\nCCAGA,2\nCCAGG,2\nCCCAC,2\nCCCAT,2\nCCCCA,2\nCCCCG,2\nCCCTA,2\nCCCTG,2\nCCCGC,2\nCCCGT,2\nCCTAC,2\nCCTAT,2\nCCTCA,2\nCCTCG,2\nCCTTA,2\nCCTTG,2\nCCTGC,2\nCCTGT,2\nCCGAA,2\nCCGAG,2\nCCGCC,2\nCCGCT,2\nCCGTC,2\nCCGTT,2\nCCGGA,2\nCCGGG,2\nCTAAA,2\nCTAAG,2\nCTACC,2\nCTACT,2\nCTATC,2\nCTATT,2\nCTAGA,2\nCTAGG,2\nCTCAC,2\nCTCAT,2\nCTCCA,2\nCTCCG,2\nCTCTA,2\nCTCTG,2\nCTCGC,2\nCTCGT,2\nCTTAC,2\nCTTAT,2\nCTTCA,2\nCTTCG,2\nCTTTA,2\nCTTTG,2\nCTTGC,2\nCTTGT,2\nCTGAA,2\nCTGAG,2\nCTGCC,2\nCTGCT,2\nCTGTC,2\nCTGTT,2\nCTGGA,2\nCTGGG,2\nCGAAC,2\nCGAAT,2\nCGACA,2\nCGACG,2\nCGATA,2\nCGATG,2\nCGAGC,2\nCGAGT,2\nCGCAA,2\nCGCAG,2\nCGCCC,2\nCGCCT,2\nCGCTC,2\nCGCTT,2\nCGCGA,2\nCGCGG,2\nCGTAA,2\nCGTAG,2\nCGTCC,2\nCGTCT,2\nCGTTC,2\nCGTTT,2\nCGTGA,2\nCGTGG,2\nCGGAC,2\nCGGAT,2\nCGGCA,2\nCGGCG,2\nCGGTA,2\nCGGTG,2\nCGGGC,2\nCGGGT,2\nTAAAC,2\nTAAAT,2\nTAACA,2\nTAACG,2\nTAATA,2\nTAATG,2\nTAAGC,2\nTAAGT,2\nTACAA,2\nTACAG,2\nTACCC,2\nTACCT,2\nTACTC,2\nTACTT,2\nTACGA,2\nTACGG,2\nTATAA,2\nTATAG,2\nTATCC,2\nTATCT,2\nTATTC,2\nTATTT,2\nTATGA,2\nTATGG,2\nTAGAC,2\nTAGAT,2\nTAGCA,2\nTAGCG,2\nTAGTA,2\nTAGTG,2\nTAGGC,2\nTAGGT,2\nTCAAA,2\nTCAAG,2\nTCACC,2\nTCACT,2\nTCATC,2\nTCATT,2\nTCAGA,2\nTCAGG,2\nTCCAC,2\nTCCAT,2\nTCCCA,2\nTCCCG,2\nTCCTA,2\nTCCTG,2\nTCCGC,2\nTCCGT,2\nTCTAC,2\nTCTAT,2\nTCTCA,2\nTCTCG,2\nTCTTA,2\nTCTTG,2\nTCTGC,2\nTCTGT,2\nTCGAA,2\nTCGAG,2\nTCGCC,2\nTCGCT,2\nTCGTC,2\nTCGTT,2\nTCGGA,2\nTCGGG,2\nTTAAA,2\nTTAAG,2\nTTACC,2\nTTACT,2\nTTATC,2\nTTATT,2\nTTAGA,2\nTTAGG,2\nTTCAC,2\nTTCAT,2\nTTCCA,2\nTTCCG,2\nTTCTA,2\nTTCTG,2\nTTCGC,2\nTTCGT,2\nTTTAC,2\nTTTAT,2\nTTTCA,2\nTTTCG,2\nTTTTA,2\nTTTTG,2\nTTTGC,2\nTTTGT,2\nTTGAA,2\nTTGAG,2\nTTGCC,2\nTTGCT,2\nTTGTC,2\nTTGTT,2\nTTGGA,2\nTTGGG,2\nTGAAC,2\nTGAAT,2\nTGACA,2\nTGACG,2\nTGATA,2\nTGATG,2\nTGAGC,2\nTGAGT,2\nTGCAA,2\nTGCAG,2\nTGCCC,2\nTGCCT,2\nTGCTC,2\nTGCTT,2\nTGCGA,2\nTGCGG,2\nTGTAA,2\nTGTAG,2\nTGTCC,2\nTGTCT,2\nTGTTC,2\nTGTTT,2\nTGTGA,2\nTGTGG,2\nTGGAC,2\nTGGAT,2\nTGGCA,2\nTGGCG,2\nTGGTA,2\nTGGTG,2\nTGGGC,2\nTGGGT,2\nGAAAA,2\nGAAAG,2\nGAACC,2\nGAACT,2\nGAATC,2\nGAATT,2\nGAAGA,2\nGAAGG,2\nGACAC,2\nGACAT,2\nGACCA,2\nGACCG,2\nGACTA,2\nGACTG,2\nGACGC,2\nGACGT,2\nGATAC,2\nGATAT,2\nGATCA,2\nGATCG,2\nGATTA,2\nGATTG,2\nGATGC,2\nGATGT,2\nGAGAA,2\nGAGAG,2\nGAGCC,2\nGAGCT,2\nGAGTC,2\nGAGTT,2\nGAGGA,2\nGAGGG,2\nGCAAC,2\nGCAAT,2\nGCACA,2\nGCACG,2\nGCATA,2\nGCATG,2\nGCAGC,2\nGCAGT,2\nGCCAA,2\nGCCAG,2\nGCCCC,2\nGCCCT,2\nGCCTC,2\nGCCTT,2\nGCCGA,2\nGCCGG,2\nGCTAA,2\nGCTAG,2\nGCTCC,2\nGCTCT,2\nGCTTC,2\nGCTTT,2\nGCTGA,2\nGCTGG,2\nGCGAC,2\nGCGAT,2\nGCGCA,2\nGCGCG,2\nGCGTA,2\nGCGTG,2\nGCGGC,2\nGCGGT,2\nGTAAC,2\nGTAAT,2\nGTACA,2\nGTACG,2\nGTATA,2\nGTATG,2\nGTAGC,2\nGTAGT,2\nGTCAA,2\nGTCAG,2\nGTCCC,2\nGTCCT,2\nGTCTC,2\nGTCTT,2\nGTCGA,2\nGTCGG,2\nGTTAA,2\nGTTAG,2\nGTTCC,2\nGTTCT,2\nGTTTC,2\nGTTTT,2\nGTTGA,2\nGTTGG,2\nGTGAC,2\nGTGAT,2\nGTGCA,2\nGTGCG,2\nGTGTA,2\nGTGTG,2\nGTGGC,2\nGTGGT,2\nGGAAA,2\nGGAAG,2\nGGACC,2\nGGACT,2\nGGATC,2\nGGATT,2\nGGAGA,2\nGGAGG,2\nGGCAC,2\nGGCAT,2\nGGCCA,2\nGGCCG,2\nGGCTA,2\nGGCTG,2\nGGCGC,2\nGGCGT,2\nGGTAC,2\nGGTAT,2\nGGTCA,2\nGGTCG,2\nGGTTA,2\nGGTTG,2\nGGTGC,2\nGGTGT,2\nGGGAA,2\nGGGAG,2\nGGGCC,2\nGGGCT,2\nGGGTC,2\nGGGTT,2\nGGGGA,2\nGGGGG,2\n";

    const CSV_ABUNDANCE_MIN_2: &[u8] = b"AAAAA,3\n";

    #[test]
    fn csv() -> error::Result<()> {
        let mut outfile = Vec::new();
        let counter = generate_counter();
        let serialize = counter.serialize();

        serialize.csv(1, &mut outfile)?;
        assert_eq!(&outfile[..], &CSV_ABUNDANCE_MIN_1[..]);

        outfile.clear();

        serialize.csv(2, &mut outfile)?;
        assert_eq!(&outfile[..], &CSV_ABUNDANCE_MIN_2[..]);

        Ok(())
    }

    #[cfg(feature = "parallel")]
    #[test]
    fn atomic_csv() -> error::Result<()> {
        let mut outfile = Vec::new();
        let counter = generate_atomic_counter();
        let serialize = counter.serialize();

        serialize.csv(1, &mut outfile)?;
        assert_eq!(&outfile[..], &CSV_ABUNDANCE_MIN_1[..]);

        outfile.clear();

        serialize.csv(2, &mut outfile)?;
        assert_eq!(&outfile[..], &CSV_ABUNDANCE_MIN_2[..]);

        Ok(())
    }

    const SOLID_ABUNDANCE_MIN_1: &[u8] = &[
        31, 139, 8, 0, 0, 0, 0, 0, 4, 255, 165, 192, 49, 1, 0, 0, 0, 64, 176, 75, 255, 200, 132,
        48, 156, 2, 70, 0, 241, 137, 65, 0, 0, 0,
    ];

    const SOLID_ABUNDANCE_MIN_2: &[u8] = &[
        31, 139, 8, 0, 0, 0, 0, 0, 4, 255, 165, 192, 49, 1, 0, 0, 0, 130, 48, 61, 232, 95, 153, 16,
        140, 175, 17, 95, 201, 40, 124, 65, 0, 0, 0,
    ];

    #[test]
    fn solid() -> error::Result<()> {
        let mut outfile = Vec::new();
        let counter = generate_counter();
        let serialize = counter.serialize();

        serialize.solid(1, &mut outfile)?;
        assert_eq!(&outfile[..], &SOLID_ABUNDANCE_MIN_1[..]);

        outfile.clear();

        serialize.solid(2, &mut outfile)?;
        assert_eq!(&outfile[..], &SOLID_ABUNDANCE_MIN_2[..]);

        Ok(())
    }

    #[cfg(feature = "parallel")]
    #[test]
    fn atomic_solid() -> error::Result<()> {
        let mut outfile = Vec::new();
        let counter = generate_atomic_counter();
        let serialize = counter.serialize();

        serialize.solid(1, &mut outfile)?;
        assert_eq!(&outfile[..], &SOLID_ABUNDANCE_MIN_1[..]);

        outfile.clear();

        serialize.solid(2, &mut outfile)?;
        assert_eq!(&outfile[..], &SOLID_ABUNDANCE_MIN_2[..]);

        Ok(())
    }

    #[cfg(feature = "kff")]
    const KFF_ABUNDANCE_MIN_1: &[u8] = &[
        75, 70, 70, 1, 0, 30, 1, 1, 0, 0, 0, 14, 112, 114, 111, 100, 117, 99, 101, 114, 58, 32,
        112, 99, 111, 110, 0, 118, 0, 0, 0, 0, 0, 0, 0, 4, 111, 114, 100, 101, 114, 101, 100, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 100, 97, 116, 97, 95, 115, 105, 122, 101, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        107, 0, 0, 0, 0, 0, 0, 0, 0, 5, 109, 97, 120, 0, 0, 0, 0, 0, 0, 0, 0, 255, 114, 0, 0, 0, 0,
        0, 0, 2, 0, 28, 0, 0, 0, 0, 0, 0, 0, 0, 3, 28, 0, 0, 0, 0, 0, 0, 0, 3, 2, 28, 0, 0, 0, 0,
        0, 0, 0, 5, 2, 28, 0, 0, 0, 0, 0, 0, 0, 6, 2, 28, 0, 0, 0, 0, 0, 0, 0, 9, 2, 28, 0, 0, 0,
        0, 0, 0, 0, 10, 2, 28, 0, 0, 0, 0, 0, 0, 0, 12, 2, 28, 0, 0, 0, 0, 0, 0, 0, 15, 2, 28, 0,
        0, 0, 0, 0, 0, 0, 17, 2, 28, 0, 0, 0, 0, 0, 0, 0, 18, 2, 28, 0, 0, 0, 0, 0, 0, 0, 20, 2,
        28, 0, 0, 0, 0, 0, 0, 0, 23, 2, 28, 0, 0, 0, 0, 0, 0, 0, 24, 2, 28, 0, 0, 0, 0, 0, 0, 0,
        27, 2, 28, 0, 0, 0, 0, 0, 0, 0, 29, 2, 28, 0, 0, 0, 0, 0, 0, 0, 30, 2, 28, 0, 0, 0, 0, 0,
        0, 0, 33, 2, 28, 0, 0, 0, 0, 0, 0, 0, 34, 2, 28, 0, 0, 0, 0, 0, 0, 0, 36, 2, 28, 0, 0, 0,
        0, 0, 0, 0, 39, 2, 28, 0, 0, 0, 0, 0, 0, 0, 40, 2, 28, 0, 0, 0, 0, 0, 0, 0, 43, 2, 28, 0,
        0, 0, 0, 0, 0, 0, 45, 2, 28, 0, 0, 0, 0, 0, 0, 0, 46, 2, 28, 0, 0, 0, 0, 0, 0, 0, 48, 2,
        28, 0, 0, 0, 0, 0, 0, 0, 51, 2, 28, 0, 0, 0, 0, 0, 0, 0, 53, 2, 28, 0, 0, 0, 0, 0, 0, 0,
        54, 2, 28, 0, 0, 0, 0, 0, 0, 0, 57, 2, 28, 0, 0, 0, 0, 0, 0, 0, 58, 2, 28, 0, 0, 0, 0, 0,
        0, 0, 60, 2, 28, 0, 0, 0, 0, 0, 0, 0, 63, 2, 28, 0, 0, 0, 0, 0, 0, 0, 65, 2, 28, 0, 0, 0,
        0, 0, 0, 0, 66, 2, 28, 0, 0, 0, 0, 0, 0, 0, 68, 2, 28, 0, 0, 0, 0, 0, 0, 0, 71, 2, 28, 0,
        0, 0, 0, 0, 0, 0, 72, 2, 28, 0, 0, 0, 0, 0, 0, 0, 75, 2, 28, 0, 0, 0, 0, 0, 0, 0, 77, 2,
        28, 0, 0, 0, 0, 0, 0, 0, 78, 2, 28, 0, 0, 0, 0, 0, 0, 0, 80, 2, 28, 0, 0, 0, 0, 0, 0, 0,
        83, 2, 28, 0, 0, 0, 0, 0, 0, 0, 85, 2, 28, 0, 0, 0, 0, 0, 0, 0, 86, 2, 28, 0, 0, 0, 0, 0,
        0, 0, 89, 2, 28, 0, 0, 0, 0, 0, 0, 0, 90, 2, 28, 0, 0, 0, 0, 0, 0, 0, 92, 2, 28, 0, 0, 0,
        0, 0, 0, 0, 95, 2, 28, 0, 0, 0, 0, 0, 0, 0, 96, 2, 28, 0, 0, 0, 0, 0, 0, 0, 99, 2, 28, 0,
        0, 0, 0, 0, 0, 0, 101, 2, 28, 0, 0, 0, 0, 0, 0, 0, 102, 2, 28, 0, 0, 0, 0, 0, 0, 0, 105, 2,
        28, 0, 0, 0, 0, 0, 0, 0, 106, 2, 28, 0, 0, 0, 0, 0, 0, 0, 108, 2, 28, 0, 0, 0, 0, 0, 0, 0,
        111, 2, 28, 0, 0, 0, 0, 0, 0, 0, 113, 2, 28, 0, 0, 0, 0, 0, 0, 0, 114, 2, 28, 0, 0, 0, 0,
        0, 0, 0, 116, 2, 28, 0, 0, 0, 0, 0, 0, 0, 119, 2, 28, 0, 0, 0, 0, 0, 0, 0, 120, 2, 28, 0,
        0, 0, 0, 0, 0, 0, 123, 2, 28, 0, 0, 0, 0, 0, 0, 0, 125, 2, 28, 0, 0, 0, 0, 0, 0, 0, 126, 2,
        28, 0, 0, 0, 0, 0, 0, 0, 129, 2, 28, 0, 0, 0, 0, 0, 0, 0, 130, 2, 28, 0, 0, 0, 0, 0, 0, 0,
        132, 2, 28, 0, 0, 0, 0, 0, 0, 0, 135, 2, 28, 0, 0, 0, 0, 0, 0, 0, 136, 2, 28, 0, 0, 0, 0,
        0, 0, 0, 139, 2, 28, 0, 0, 0, 0, 0, 0, 0, 141, 2, 28, 0, 0, 0, 0, 0, 0, 0, 142, 2, 28, 0,
        0, 0, 0, 0, 0, 0, 144, 2, 28, 0, 0, 0, 0, 0, 0, 0, 147, 2, 28, 0, 0, 0, 0, 0, 0, 0, 149, 2,
        28, 0, 0, 0, 0, 0, 0, 0, 150, 2, 28, 0, 0, 0, 0, 0, 0, 0, 153, 2, 28, 0, 0, 0, 0, 0, 0, 0,
        154, 2, 28, 0, 0, 0, 0, 0, 0, 0, 156, 2, 28, 0, 0, 0, 0, 0, 0, 0, 159, 2, 28, 0, 0, 0, 0,
        0, 0, 0, 160, 2, 28, 0, 0, 0, 0, 0, 0, 0, 163, 2, 28, 0, 0, 0, 0, 0, 0, 0, 165, 2, 28, 0,
        0, 0, 0, 0, 0, 0, 166, 2, 28, 0, 0, 0, 0, 0, 0, 0, 169, 2, 28, 0, 0, 0, 0, 0, 0, 0, 170, 2,
        28, 0, 0, 0, 0, 0, 0, 0, 172, 2, 28, 0, 0, 0, 0, 0, 0, 0, 175, 2, 28, 0, 0, 0, 0, 0, 0, 0,
        177, 2, 28, 0, 0, 0, 0, 0, 0, 0, 178, 2, 28, 0, 0, 0, 0, 0, 0, 0, 180, 2, 28, 0, 0, 0, 0,
        0, 0, 0, 183, 2, 28, 0, 0, 0, 0, 0, 0, 0, 184, 2, 28, 0, 0, 0, 0, 0, 0, 0, 187, 2, 28, 0,
        0, 0, 0, 0, 0, 0, 189, 2, 28, 0, 0, 0, 0, 0, 0, 0, 190, 2, 28, 0, 0, 0, 0, 0, 0, 0, 192, 2,
        28, 0, 0, 0, 0, 0, 0, 0, 195, 2, 28, 0, 0, 0, 0, 0, 0, 0, 197, 2, 28, 0, 0, 0, 0, 0, 0, 0,
        198, 2, 28, 0, 0, 0, 0, 0, 0, 0, 201, 2, 28, 0, 0, 0, 0, 0, 0, 0, 202, 2, 28, 0, 0, 0, 0,
        0, 0, 0, 204, 2, 28, 0, 0, 0, 0, 0, 0, 0, 207, 2, 28, 0, 0, 0, 0, 0, 0, 0, 209, 2, 28, 0,
        0, 0, 0, 0, 0, 0, 210, 2, 28, 0, 0, 0, 0, 0, 0, 0, 212, 2, 28, 0, 0, 0, 0, 0, 0, 0, 215, 2,
        28, 0, 0, 0, 0, 0, 0, 0, 216, 2, 28, 0, 0, 0, 0, 0, 0, 0, 219, 2, 28, 0, 0, 0, 0, 0, 0, 0,
        221, 2, 28, 0, 0, 0, 0, 0, 0, 0, 222, 2, 28, 0, 0, 0, 0, 0, 0, 0, 225, 2, 28, 0, 0, 0, 0,
        0, 0, 0, 226, 2, 28, 0, 0, 0, 0, 0, 0, 0, 228, 2, 28, 0, 0, 0, 0, 0, 0, 0, 231, 2, 28, 0,
        0, 0, 0, 0, 0, 0, 232, 2, 28, 0, 0, 0, 0, 0, 0, 0, 235, 2, 28, 0, 0, 0, 0, 0, 0, 0, 237, 2,
        28, 0, 0, 0, 0, 0, 0, 0, 238, 2, 28, 0, 0, 0, 0, 0, 0, 0, 240, 2, 28, 0, 0, 0, 0, 0, 0, 0,
        243, 2, 28, 0, 0, 0, 0, 0, 0, 0, 245, 2, 28, 0, 0, 0, 0, 0, 0, 0, 246, 2, 28, 0, 0, 0, 0,
        0, 0, 0, 249, 2, 28, 0, 0, 0, 0, 0, 0, 0, 250, 2, 28, 0, 0, 0, 0, 0, 0, 0, 252, 2, 28, 0,
        0, 0, 0, 0, 0, 0, 255, 2, 28, 0, 0, 0, 0, 0, 0, 1, 1, 2, 28, 0, 0, 0, 0, 0, 0, 1, 2, 2, 28,
        0, 0, 0, 0, 0, 0, 1, 4, 2, 28, 0, 0, 0, 0, 0, 0, 1, 7, 2, 28, 0, 0, 0, 0, 0, 0, 1, 8, 2,
        28, 0, 0, 0, 0, 0, 0, 1, 11, 2, 28, 0, 0, 0, 0, 0, 0, 1, 13, 2, 28, 0, 0, 0, 0, 0, 0, 1,
        14, 2, 28, 0, 0, 0, 0, 0, 0, 1, 16, 2, 28, 0, 0, 0, 0, 0, 0, 1, 19, 2, 28, 0, 0, 0, 0, 0,
        0, 1, 21, 2, 28, 0, 0, 0, 0, 0, 0, 1, 22, 2, 28, 0, 0, 0, 0, 0, 0, 1, 25, 2, 28, 0, 0, 0,
        0, 0, 0, 1, 26, 2, 28, 0, 0, 0, 0, 0, 0, 1, 28, 2, 28, 0, 0, 0, 0, 0, 0, 1, 31, 2, 28, 0,
        0, 0, 0, 0, 0, 1, 32, 2, 28, 0, 0, 0, 0, 0, 0, 1, 35, 2, 28, 0, 0, 0, 0, 0, 0, 1, 37, 2,
        28, 0, 0, 0, 0, 0, 0, 1, 38, 2, 28, 0, 0, 0, 0, 0, 0, 1, 41, 2, 28, 0, 0, 0, 0, 0, 0, 1,
        42, 2, 28, 0, 0, 0, 0, 0, 0, 1, 44, 2, 28, 0, 0, 0, 0, 0, 0, 1, 47, 2, 28, 0, 0, 0, 0, 0,
        0, 1, 49, 2, 28, 0, 0, 0, 0, 0, 0, 1, 50, 2, 28, 0, 0, 0, 0, 0, 0, 1, 52, 2, 28, 0, 0, 0,
        0, 0, 0, 1, 55, 2, 28, 0, 0, 0, 0, 0, 0, 1, 56, 2, 28, 0, 0, 0, 0, 0, 0, 1, 59, 2, 28, 0,
        0, 0, 0, 0, 0, 1, 61, 2, 28, 0, 0, 0, 0, 0, 0, 1, 62, 2, 28, 0, 0, 0, 0, 0, 0, 1, 64, 2,
        28, 0, 0, 0, 0, 0, 0, 1, 67, 2, 28, 0, 0, 0, 0, 0, 0, 1, 69, 2, 28, 0, 0, 0, 0, 0, 0, 1,
        70, 2, 28, 0, 0, 0, 0, 0, 0, 1, 73, 2, 28, 0, 0, 0, 0, 0, 0, 1, 74, 2, 28, 0, 0, 0, 0, 0,
        0, 1, 76, 2, 28, 0, 0, 0, 0, 0, 0, 1, 79, 2, 28, 0, 0, 0, 0, 0, 0, 1, 81, 2, 28, 0, 0, 0,
        0, 0, 0, 1, 82, 2, 28, 0, 0, 0, 0, 0, 0, 1, 84, 2, 28, 0, 0, 0, 0, 0, 0, 1, 87, 2, 28, 0,
        0, 0, 0, 0, 0, 1, 88, 2, 28, 0, 0, 0, 0, 0, 0, 1, 91, 2, 28, 0, 0, 0, 0, 0, 0, 1, 93, 2,
        28, 0, 0, 0, 0, 0, 0, 1, 94, 2, 28, 0, 0, 0, 0, 0, 0, 1, 97, 2, 28, 0, 0, 0, 0, 0, 0, 1,
        98, 2, 28, 0, 0, 0, 0, 0, 0, 1, 100, 2, 28, 0, 0, 0, 0, 0, 0, 1, 103, 2, 28, 0, 0, 0, 0, 0,
        0, 1, 104, 2, 28, 0, 0, 0, 0, 0, 0, 1, 107, 2, 28, 0, 0, 0, 0, 0, 0, 1, 109, 2, 28, 0, 0,
        0, 0, 0, 0, 1, 110, 2, 28, 0, 0, 0, 0, 0, 0, 1, 112, 2, 28, 0, 0, 0, 0, 0, 0, 1, 115, 2,
        28, 0, 0, 0, 0, 0, 0, 1, 117, 2, 28, 0, 0, 0, 0, 0, 0, 1, 118, 2, 28, 0, 0, 0, 0, 0, 0, 1,
        121, 2, 28, 0, 0, 0, 0, 0, 0, 1, 122, 2, 28, 0, 0, 0, 0, 0, 0, 1, 124, 2, 28, 0, 0, 0, 0,
        0, 0, 1, 127, 2, 28, 0, 0, 0, 0, 0, 0, 1, 128, 2, 28, 0, 0, 0, 0, 0, 0, 1, 131, 2, 28, 0,
        0, 0, 0, 0, 0, 1, 133, 2, 28, 0, 0, 0, 0, 0, 0, 1, 134, 2, 28, 0, 0, 0, 0, 0, 0, 1, 137, 2,
        28, 0, 0, 0, 0, 0, 0, 1, 138, 2, 28, 0, 0, 0, 0, 0, 0, 1, 140, 2, 28, 0, 0, 0, 0, 0, 0, 1,
        143, 2, 28, 0, 0, 0, 0, 0, 0, 1, 145, 2, 28, 0, 0, 0, 0, 0, 0, 1, 146, 2, 28, 0, 0, 0, 0,
        0, 0, 1, 148, 2, 28, 0, 0, 0, 0, 0, 0, 1, 151, 2, 28, 0, 0, 0, 0, 0, 0, 1, 152, 2, 28, 0,
        0, 0, 0, 0, 0, 1, 155, 2, 28, 0, 0, 0, 0, 0, 0, 1, 157, 2, 28, 0, 0, 0, 0, 0, 0, 1, 158, 2,
        28, 0, 0, 0, 0, 0, 0, 1, 161, 2, 28, 0, 0, 0, 0, 0, 0, 1, 162, 2, 28, 0, 0, 0, 0, 0, 0, 1,
        164, 2, 28, 0, 0, 0, 0, 0, 0, 1, 167, 2, 28, 0, 0, 0, 0, 0, 0, 1, 168, 2, 28, 0, 0, 0, 0,
        0, 0, 1, 171, 2, 28, 0, 0, 0, 0, 0, 0, 1, 173, 2, 28, 0, 0, 0, 0, 0, 0, 1, 174, 2, 28, 0,
        0, 0, 0, 0, 0, 1, 176, 2, 28, 0, 0, 0, 0, 0, 0, 1, 179, 2, 28, 0, 0, 0, 0, 0, 0, 1, 181, 2,
        28, 0, 0, 0, 0, 0, 0, 1, 182, 2, 28, 0, 0, 0, 0, 0, 0, 1, 185, 2, 28, 0, 0, 0, 0, 0, 0, 1,
        186, 2, 28, 0, 0, 0, 0, 0, 0, 1, 188, 2, 28, 0, 0, 0, 0, 0, 0, 1, 191, 2, 28, 0, 0, 0, 0,
        0, 0, 1, 193, 2, 28, 0, 0, 0, 0, 0, 0, 1, 194, 2, 28, 0, 0, 0, 0, 0, 0, 1, 196, 2, 28, 0,
        0, 0, 0, 0, 0, 1, 199, 2, 28, 0, 0, 0, 0, 0, 0, 1, 200, 2, 28, 0, 0, 0, 0, 0, 0, 1, 203, 2,
        28, 0, 0, 0, 0, 0, 0, 1, 205, 2, 28, 0, 0, 0, 0, 0, 0, 1, 206, 2, 28, 0, 0, 0, 0, 0, 0, 1,
        208, 2, 28, 0, 0, 0, 0, 0, 0, 1, 211, 2, 28, 0, 0, 0, 0, 0, 0, 1, 213, 2, 28, 0, 0, 0, 0,
        0, 0, 1, 214, 2, 28, 0, 0, 0, 0, 0, 0, 1, 217, 2, 28, 0, 0, 0, 0, 0, 0, 1, 218, 2, 28, 0,
        0, 0, 0, 0, 0, 1, 220, 2, 28, 0, 0, 0, 0, 0, 0, 1, 223, 2, 28, 0, 0, 0, 0, 0, 0, 1, 224, 2,
        28, 0, 0, 0, 0, 0, 0, 1, 227, 2, 28, 0, 0, 0, 0, 0, 0, 1, 229, 2, 28, 0, 0, 0, 0, 0, 0, 1,
        230, 2, 28, 0, 0, 0, 0, 0, 0, 1, 233, 2, 28, 0, 0, 0, 0, 0, 0, 1, 234, 2, 28, 0, 0, 0, 0,
        0, 0, 1, 236, 2, 28, 0, 0, 0, 0, 0, 0, 1, 239, 2, 28, 0, 0, 0, 0, 0, 0, 1, 241, 2, 28, 0,
        0, 0, 0, 0, 0, 1, 242, 2, 28, 0, 0, 0, 0, 0, 0, 1, 244, 2, 28, 0, 0, 0, 0, 0, 0, 1, 247, 2,
        28, 0, 0, 0, 0, 0, 0, 1, 248, 2, 28, 0, 0, 0, 0, 0, 0, 1, 251, 2, 28, 0, 0, 0, 0, 0, 0, 1,
        253, 2, 28, 0, 0, 0, 0, 0, 0, 1, 254, 2, 28, 0, 0, 0, 0, 0, 0, 2, 1, 2, 28, 0, 0, 0, 0, 0,
        0, 2, 2, 2, 28, 0, 0, 0, 0, 0, 0, 2, 4, 2, 28, 0, 0, 0, 0, 0, 0, 2, 7, 2, 28, 0, 0, 0, 0,
        0, 0, 2, 8, 2, 28, 0, 0, 0, 0, 0, 0, 2, 11, 2, 28, 0, 0, 0, 0, 0, 0, 2, 13, 2, 28, 0, 0, 0,
        0, 0, 0, 2, 14, 2, 28, 0, 0, 0, 0, 0, 0, 2, 16, 2, 28, 0, 0, 0, 0, 0, 0, 2, 19, 2, 28, 0,
        0, 0, 0, 0, 0, 2, 21, 2, 28, 0, 0, 0, 0, 0, 0, 2, 22, 2, 28, 0, 0, 0, 0, 0, 0, 2, 25, 2,
        28, 0, 0, 0, 0, 0, 0, 2, 26, 2, 28, 0, 0, 0, 0, 0, 0, 2, 28, 2, 28, 0, 0, 0, 0, 0, 0, 2,
        31, 2, 28, 0, 0, 0, 0, 0, 0, 2, 32, 2, 28, 0, 0, 0, 0, 0, 0, 2, 35, 2, 28, 0, 0, 0, 0, 0,
        0, 2, 37, 2, 28, 0, 0, 0, 0, 0, 0, 2, 38, 2, 28, 0, 0, 0, 0, 0, 0, 2, 41, 2, 28, 0, 0, 0,
        0, 0, 0, 2, 42, 2, 28, 0, 0, 0, 0, 0, 0, 2, 44, 2, 28, 0, 0, 0, 0, 0, 0, 2, 47, 2, 28, 0,
        0, 0, 0, 0, 0, 2, 49, 2, 28, 0, 0, 0, 0, 0, 0, 2, 50, 2, 28, 0, 0, 0, 0, 0, 0, 2, 52, 2,
        28, 0, 0, 0, 0, 0, 0, 2, 55, 2, 28, 0, 0, 0, 0, 0, 0, 2, 56, 2, 28, 0, 0, 0, 0, 0, 0, 2,
        59, 2, 28, 0, 0, 0, 0, 0, 0, 2, 61, 2, 28, 0, 0, 0, 0, 0, 0, 2, 62, 2, 28, 0, 0, 0, 0, 0,
        0, 2, 64, 2, 28, 0, 0, 0, 0, 0, 0, 2, 67, 2, 28, 0, 0, 0, 0, 0, 0, 2, 69, 2, 28, 0, 0, 0,
        0, 0, 0, 2, 70, 2, 28, 0, 0, 0, 0, 0, 0, 2, 73, 2, 28, 0, 0, 0, 0, 0, 0, 2, 74, 2, 28, 0,
        0, 0, 0, 0, 0, 2, 76, 2, 28, 0, 0, 0, 0, 0, 0, 2, 79, 2, 28, 0, 0, 0, 0, 0, 0, 2, 81, 2,
        28, 0, 0, 0, 0, 0, 0, 2, 82, 2, 28, 0, 0, 0, 0, 0, 0, 2, 84, 2, 28, 0, 0, 0, 0, 0, 0, 2,
        87, 2, 28, 0, 0, 0, 0, 0, 0, 2, 88, 2, 28, 0, 0, 0, 0, 0, 0, 2, 91, 2, 28, 0, 0, 0, 0, 0,
        0, 2, 93, 2, 28, 0, 0, 0, 0, 0, 0, 2, 94, 2, 28, 0, 0, 0, 0, 0, 0, 2, 97, 2, 28, 0, 0, 0,
        0, 0, 0, 2, 98, 2, 28, 0, 0, 0, 0, 0, 0, 2, 100, 2, 28, 0, 0, 0, 0, 0, 0, 2, 103, 2, 28, 0,
        0, 0, 0, 0, 0, 2, 104, 2, 28, 0, 0, 0, 0, 0, 0, 2, 107, 2, 28, 0, 0, 0, 0, 0, 0, 2, 109, 2,
        28, 0, 0, 0, 0, 0, 0, 2, 110, 2, 28, 0, 0, 0, 0, 0, 0, 2, 112, 2, 28, 0, 0, 0, 0, 0, 0, 2,
        115, 2, 28, 0, 0, 0, 0, 0, 0, 2, 117, 2, 28, 0, 0, 0, 0, 0, 0, 2, 118, 2, 28, 0, 0, 0, 0,
        0, 0, 2, 121, 2, 28, 0, 0, 0, 0, 0, 0, 2, 122, 2, 28, 0, 0, 0, 0, 0, 0, 2, 124, 2, 28, 0,
        0, 0, 0, 0, 0, 2, 127, 2, 28, 0, 0, 0, 0, 0, 0, 2, 128, 2, 28, 0, 0, 0, 0, 0, 0, 2, 131, 2,
        28, 0, 0, 0, 0, 0, 0, 2, 133, 2, 28, 0, 0, 0, 0, 0, 0, 2, 134, 2, 28, 0, 0, 0, 0, 0, 0, 2,
        137, 2, 28, 0, 0, 0, 0, 0, 0, 2, 138, 2, 28, 0, 0, 0, 0, 0, 0, 2, 140, 2, 28, 0, 0, 0, 0,
        0, 0, 2, 143, 2, 28, 0, 0, 0, 0, 0, 0, 2, 145, 2, 28, 0, 0, 0, 0, 0, 0, 2, 146, 2, 28, 0,
        0, 0, 0, 0, 0, 2, 148, 2, 28, 0, 0, 0, 0, 0, 0, 2, 151, 2, 28, 0, 0, 0, 0, 0, 0, 2, 152, 2,
        28, 0, 0, 0, 0, 0, 0, 2, 155, 2, 28, 0, 0, 0, 0, 0, 0, 2, 157, 2, 28, 0, 0, 0, 0, 0, 0, 2,
        158, 2, 28, 0, 0, 0, 0, 0, 0, 2, 161, 2, 28, 0, 0, 0, 0, 0, 0, 2, 162, 2, 28, 0, 0, 0, 0,
        0, 0, 2, 164, 2, 28, 0, 0, 0, 0, 0, 0, 2, 167, 2, 28, 0, 0, 0, 0, 0, 0, 2, 168, 2, 28, 0,
        0, 0, 0, 0, 0, 2, 171, 2, 28, 0, 0, 0, 0, 0, 0, 2, 173, 2, 28, 0, 0, 0, 0, 0, 0, 2, 174, 2,
        28, 0, 0, 0, 0, 0, 0, 2, 176, 2, 28, 0, 0, 0, 0, 0, 0, 2, 179, 2, 28, 0, 0, 0, 0, 0, 0, 2,
        181, 2, 28, 0, 0, 0, 0, 0, 0, 2, 182, 2, 28, 0, 0, 0, 0, 0, 0, 2, 185, 2, 28, 0, 0, 0, 0,
        0, 0, 2, 186, 2, 28, 0, 0, 0, 0, 0, 0, 2, 188, 2, 28, 0, 0, 0, 0, 0, 0, 2, 191, 2, 28, 0,
        0, 0, 0, 0, 0, 2, 193, 2, 28, 0, 0, 0, 0, 0, 0, 2, 194, 2, 28, 0, 0, 0, 0, 0, 0, 2, 196, 2,
        28, 0, 0, 0, 0, 0, 0, 2, 199, 2, 28, 0, 0, 0, 0, 0, 0, 2, 200, 2, 28, 0, 0, 0, 0, 0, 0, 2,
        203, 2, 28, 0, 0, 0, 0, 0, 0, 2, 205, 2, 28, 0, 0, 0, 0, 0, 0, 2, 206, 2, 28, 0, 0, 0, 0,
        0, 0, 2, 208, 2, 28, 0, 0, 0, 0, 0, 0, 2, 211, 2, 28, 0, 0, 0, 0, 0, 0, 2, 213, 2, 28, 0,
        0, 0, 0, 0, 0, 2, 214, 2, 28, 0, 0, 0, 0, 0, 0, 2, 217, 2, 28, 0, 0, 0, 0, 0, 0, 2, 218, 2,
        28, 0, 0, 0, 0, 0, 0, 2, 220, 2, 28, 0, 0, 0, 0, 0, 0, 2, 223, 2, 28, 0, 0, 0, 0, 0, 0, 2,
        224, 2, 28, 0, 0, 0, 0, 0, 0, 2, 227, 2, 28, 0, 0, 0, 0, 0, 0, 2, 229, 2, 28, 0, 0, 0, 0,
        0, 0, 2, 230, 2, 28, 0, 0, 0, 0, 0, 0, 2, 233, 2, 28, 0, 0, 0, 0, 0, 0, 2, 234, 2, 28, 0,
        0, 0, 0, 0, 0, 2, 236, 2, 28, 0, 0, 0, 0, 0, 0, 2, 239, 2, 28, 0, 0, 0, 0, 0, 0, 2, 241, 2,
        28, 0, 0, 0, 0, 0, 0, 2, 242, 2, 28, 0, 0, 0, 0, 0, 0, 2, 244, 2, 28, 0, 0, 0, 0, 0, 0, 2,
        247, 2, 28, 0, 0, 0, 0, 0, 0, 2, 248, 2, 28, 0, 0, 0, 0, 0, 0, 2, 251, 2, 28, 0, 0, 0, 0,
        0, 0, 2, 253, 2, 28, 0, 0, 0, 0, 0, 0, 2, 254, 2, 28, 0, 0, 0, 0, 0, 0, 3, 0, 2, 28, 0, 0,
        0, 0, 0, 0, 3, 3, 2, 28, 0, 0, 0, 0, 0, 0, 3, 5, 2, 28, 0, 0, 0, 0, 0, 0, 3, 6, 2, 28, 0,
        0, 0, 0, 0, 0, 3, 9, 2, 28, 0, 0, 0, 0, 0, 0, 3, 10, 2, 28, 0, 0, 0, 0, 0, 0, 3, 12, 2, 28,
        0, 0, 0, 0, 0, 0, 3, 15, 2, 28, 0, 0, 0, 0, 0, 0, 3, 17, 2, 28, 0, 0, 0, 0, 0, 0, 3, 18, 2,
        28, 0, 0, 0, 0, 0, 0, 3, 20, 2, 28, 0, 0, 0, 0, 0, 0, 3, 23, 2, 28, 0, 0, 0, 0, 0, 0, 3,
        24, 2, 28, 0, 0, 0, 0, 0, 0, 3, 27, 2, 28, 0, 0, 0, 0, 0, 0, 3, 29, 2, 28, 0, 0, 0, 0, 0,
        0, 3, 30, 2, 28, 0, 0, 0, 0, 0, 0, 3, 33, 2, 28, 0, 0, 0, 0, 0, 0, 3, 34, 2, 28, 0, 0, 0,
        0, 0, 0, 3, 36, 2, 28, 0, 0, 0, 0, 0, 0, 3, 39, 2, 28, 0, 0, 0, 0, 0, 0, 3, 40, 2, 28, 0,
        0, 0, 0, 0, 0, 3, 43, 2, 28, 0, 0, 0, 0, 0, 0, 3, 45, 2, 28, 0, 0, 0, 0, 0, 0, 3, 46, 2,
        28, 0, 0, 0, 0, 0, 0, 3, 48, 2, 28, 0, 0, 0, 0, 0, 0, 3, 51, 2, 28, 0, 0, 0, 0, 0, 0, 3,
        53, 2, 28, 0, 0, 0, 0, 0, 0, 3, 54, 2, 28, 0, 0, 0, 0, 0, 0, 3, 57, 2, 28, 0, 0, 0, 0, 0,
        0, 3, 58, 2, 28, 0, 0, 0, 0, 0, 0, 3, 60, 2, 28, 0, 0, 0, 0, 0, 0, 3, 63, 2, 28, 0, 0, 0,
        0, 0, 0, 3, 65, 2, 28, 0, 0, 0, 0, 0, 0, 3, 66, 2, 28, 0, 0, 0, 0, 0, 0, 3, 68, 2, 28, 0,
        0, 0, 0, 0, 0, 3, 71, 2, 28, 0, 0, 0, 0, 0, 0, 3, 72, 2, 28, 0, 0, 0, 0, 0, 0, 3, 75, 2,
        28, 0, 0, 0, 0, 0, 0, 3, 77, 2, 28, 0, 0, 0, 0, 0, 0, 3, 78, 2, 28, 0, 0, 0, 0, 0, 0, 3,
        80, 2, 28, 0, 0, 0, 0, 0, 0, 3, 83, 2, 28, 0, 0, 0, 0, 0, 0, 3, 85, 2, 28, 0, 0, 0, 0, 0,
        0, 3, 86, 2, 28, 0, 0, 0, 0, 0, 0, 3, 89, 2, 28, 0, 0, 0, 0, 0, 0, 3, 90, 2, 28, 0, 0, 0,
        0, 0, 0, 3, 92, 2, 28, 0, 0, 0, 0, 0, 0, 3, 95, 2, 28, 0, 0, 0, 0, 0, 0, 3, 96, 2, 28, 0,
        0, 0, 0, 0, 0, 3, 99, 2, 28, 0, 0, 0, 0, 0, 0, 3, 101, 2, 28, 0, 0, 0, 0, 0, 0, 3, 102, 2,
        28, 0, 0, 0, 0, 0, 0, 3, 105, 2, 28, 0, 0, 0, 0, 0, 0, 3, 106, 2, 28, 0, 0, 0, 0, 0, 0, 3,
        108, 2, 28, 0, 0, 0, 0, 0, 0, 3, 111, 2, 28, 0, 0, 0, 0, 0, 0, 3, 113, 2, 28, 0, 0, 0, 0,
        0, 0, 3, 114, 2, 28, 0, 0, 0, 0, 0, 0, 3, 116, 2, 28, 0, 0, 0, 0, 0, 0, 3, 119, 2, 28, 0,
        0, 0, 0, 0, 0, 3, 120, 2, 28, 0, 0, 0, 0, 0, 0, 3, 123, 2, 28, 0, 0, 0, 0, 0, 0, 3, 125, 2,
        28, 0, 0, 0, 0, 0, 0, 3, 126, 2, 28, 0, 0, 0, 0, 0, 0, 3, 129, 2, 28, 0, 0, 0, 0, 0, 0, 3,
        130, 2, 28, 0, 0, 0, 0, 0, 0, 3, 132, 2, 28, 0, 0, 0, 0, 0, 0, 3, 135, 2, 28, 0, 0, 0, 0,
        0, 0, 3, 136, 2, 28, 0, 0, 0, 0, 0, 0, 3, 139, 2, 28, 0, 0, 0, 0, 0, 0, 3, 141, 2, 28, 0,
        0, 0, 0, 0, 0, 3, 142, 2, 28, 0, 0, 0, 0, 0, 0, 3, 144, 2, 28, 0, 0, 0, 0, 0, 0, 3, 147, 2,
        28, 0, 0, 0, 0, 0, 0, 3, 149, 2, 28, 0, 0, 0, 0, 0, 0, 3, 150, 2, 28, 0, 0, 0, 0, 0, 0, 3,
        153, 2, 28, 0, 0, 0, 0, 0, 0, 3, 154, 2, 28, 0, 0, 0, 0, 0, 0, 3, 156, 2, 28, 0, 0, 0, 0,
        0, 0, 3, 159, 2, 28, 0, 0, 0, 0, 0, 0, 3, 160, 2, 28, 0, 0, 0, 0, 0, 0, 3, 163, 2, 28, 0,
        0, 0, 0, 0, 0, 3, 165, 2, 28, 0, 0, 0, 0, 0, 0, 3, 166, 2, 28, 0, 0, 0, 0, 0, 0, 3, 169, 2,
        28, 0, 0, 0, 0, 0, 0, 3, 170, 2, 28, 0, 0, 0, 0, 0, 0, 3, 172, 2, 28, 0, 0, 0, 0, 0, 0, 3,
        175, 2, 28, 0, 0, 0, 0, 0, 0, 3, 177, 2, 28, 0, 0, 0, 0, 0, 0, 3, 178, 2, 28, 0, 0, 0, 0,
        0, 0, 3, 180, 2, 28, 0, 0, 0, 0, 0, 0, 3, 183, 2, 28, 0, 0, 0, 0, 0, 0, 3, 184, 2, 28, 0,
        0, 0, 0, 0, 0, 3, 187, 2, 28, 0, 0, 0, 0, 0, 0, 3, 189, 2, 28, 0, 0, 0, 0, 0, 0, 3, 190, 2,
        28, 0, 0, 0, 0, 0, 0, 3, 192, 2, 28, 0, 0, 0, 0, 0, 0, 3, 195, 2, 28, 0, 0, 0, 0, 0, 0, 3,
        197, 2, 28, 0, 0, 0, 0, 0, 0, 3, 198, 2, 28, 0, 0, 0, 0, 0, 0, 3, 201, 2, 28, 0, 0, 0, 0,
        0, 0, 3, 202, 2, 28, 0, 0, 0, 0, 0, 0, 3, 204, 2, 28, 0, 0, 0, 0, 0, 0, 3, 207, 2, 28, 0,
        0, 0, 0, 0, 0, 3, 209, 2, 28, 0, 0, 0, 0, 0, 0, 3, 210, 2, 28, 0, 0, 0, 0, 0, 0, 3, 212, 2,
        28, 0, 0, 0, 0, 0, 0, 3, 215, 2, 28, 0, 0, 0, 0, 0, 0, 3, 216, 2, 28, 0, 0, 0, 0, 0, 0, 3,
        219, 2, 28, 0, 0, 0, 0, 0, 0, 3, 221, 2, 28, 0, 0, 0, 0, 0, 0, 3, 222, 2, 28, 0, 0, 0, 0,
        0, 0, 3, 225, 2, 28, 0, 0, 0, 0, 0, 0, 3, 226, 2, 28, 0, 0, 0, 0, 0, 0, 3, 228, 2, 28, 0,
        0, 0, 0, 0, 0, 3, 231, 2, 28, 0, 0, 0, 0, 0, 0, 3, 232, 2, 28, 0, 0, 0, 0, 0, 0, 3, 235, 2,
        28, 0, 0, 0, 0, 0, 0, 3, 237, 2, 28, 0, 0, 0, 0, 0, 0, 3, 238, 2, 28, 0, 0, 0, 0, 0, 0, 3,
        240, 2, 28, 0, 0, 0, 0, 0, 0, 3, 243, 2, 28, 0, 0, 0, 0, 0, 0, 3, 245, 2, 28, 0, 0, 0, 0,
        0, 0, 3, 246, 2, 28, 0, 0, 0, 0, 0, 0, 3, 249, 2, 28, 0, 0, 0, 0, 0, 0, 3, 250, 2, 28, 0,
        0, 0, 0, 0, 0, 3, 252, 2, 28, 0, 0, 0, 0, 0, 0, 3, 255, 2, 75, 70, 70,
    ];

    #[cfg(feature = "kff")]
    const KFF_ABUNDANCE_MIN_2: &[u8] = &[
        75, 70, 70, 1, 0, 30, 1, 1, 0, 0, 0, 14, 112, 114, 111, 100, 117, 99, 101, 114, 58, 32,
        112, 99, 111, 110, 0, 118, 0, 0, 0, 0, 0, 0, 0, 4, 111, 114, 100, 101, 114, 101, 100, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 100, 97, 116, 97, 95, 115, 105, 122, 101, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        107, 0, 0, 0, 0, 0, 0, 0, 0, 5, 109, 97, 120, 0, 0, 0, 0, 0, 0, 0, 0, 255, 114, 0, 0, 0, 0,
        0, 0, 0, 1, 28, 0, 0, 0, 0, 0, 0, 0, 0, 3, 75, 70, 70,
    ];

    #[cfg(feature = "kff")]
    #[test]
    fn kff() -> error::Result<()> {
        let mut outfile = Vec::new();
        let counter = generate_counter();
        let serialize = counter.serialize();

        serialize.kff(1, &mut outfile)?;
        assert_eq!(&outfile, KFF_ABUNDANCE_MIN_1);

        outfile.clear();

        serialize.kff(2, &mut outfile)?;
        assert_eq!(&outfile[..], &KFF_ABUNDANCE_MIN_2[..]);

        Ok(())
    }

    #[cfg(all(feature = "kff", feature = "parallel"))]
    #[test]
    fn atomic_kff() -> error::Result<()> {
        let mut outfile = Vec::new();
        let counter = generate_atomic_counter();
        let serialize = counter.serialize();

        serialize.kff(1, &mut outfile)?;
        assert_eq!(&outfile, KFF_ABUNDANCE_MIN_1);

        outfile.clear();

        serialize.kff(2, &mut outfile)?;
        assert_eq!(&outfile[..], &KFF_ABUNDANCE_MIN_2[..]);

        Ok(())
    }
}
