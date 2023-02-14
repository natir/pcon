//! Tools to serialize a Counter

/* std use */

use std::io::Write as _;

/* crate use */
use byteorder::WriteBytesExt as _;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/* project use */
use crate::counter;
use crate::error;
use crate::solid;

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

            /// Write kmer in csv format
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
            /// The first bytes contains the size of k the rest of the file and a
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

                let count =
                    unsafe { std::mem::transmute::<&[$type], &[$out_type]>(self.counter.raw()) };

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
                let counts =
                    unsafe { std::mem::transmute::<&[$type], &[$out_type]>(self.counter.raw()) };

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
                    unsafe { std::mem::transmute::<&[$type], &[$out_type]>(self.counter.raw()) },
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

    const CSV_ABUNDANCE_MIN_2: &[u8] = &[65, 65, 65, 65, 65, 44, 51, 10];

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
}
