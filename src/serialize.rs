//! Serialize function

/* std use */

use std::io::Write as _;

/* crate use */

/* project use */
use crate::counter;
use crate::error;
use crate::solid;

/// Struct to serialize counter
pub struct Serialize<T> {
    counter: counter::Counter<T>,
}

impl<T> Serialize<T> {
    /// Create a new Serialize
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

                writer.write_all(&[solid.k()])?;

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
