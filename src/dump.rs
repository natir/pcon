/*
Copyright (c) 2020 Pierre Marijon <pmarijon@mpi-inf.mpg.de>

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

/* crate use */
use anyhow::{anyhow, Context, Result};
use rayon::iter::ParallelBridge;
use rayon::prelude::*;

/* local use */
use cocktail::*;

use crate::*;

use crate::error::IO::*;
use crate::error::*;

pub fn dump(params: cli::SubCommandDump) -> Result<()> {
    let params = cli::check_dump_param(params)?;

    log::info!("Start of read of count");
    let reader = std::fs::File::open(&params.input)
        .with_context(|| Error::IO(CantOpenFile))
        .with_context(|| anyhow!("File {}", params.input.clone()))?;

    let counter = counter::Counter::deserialize(reader)
        .with_context(|| Error::IO(ErrorDurringRead))
        .with_context(|| anyhow!("File {}", params.input.clone()))?;

    dump_worker(
        counter,
        None,
        params.csv,
        params.solid,
        params.spectrum,
        params.abundance,
    )?;

    log::info!("End of read of count");

    Ok(())
}

pub(crate) fn dump_worker(
    counter: counter::Counter,
    bin_path: Option<String>,
    csv_path: Option<String>,
    solid_path: Option<String>,
    spectrum_path: Option<String>,
    abundance: counter::Count,
) -> Result<()> {
    if let Some(output) = bin_path.clone() {
        log::info!("Start of dump count data in binary");
        let writer = std::fs::File::create(&output)
            .with_context(|| Error::IO(CantCreateFile))
            .with_context(|| anyhow!("File {}", output.clone()))?;

        counter
            .serialize(writer)
            .with_context(|| Error::IO(ErrorDurringWrite))
            .with_context(|| anyhow!("In file {}", output.clone()))?;

        log::info!("End of dump count data in binary");
    }

    if let Some(output) = csv_path.clone() {
        log::info!("Start of dump count data in csv");

        let writer = std::io::BufWriter::new(
            std::fs::File::create(&output)
                .with_context(|| Error::IO(CantOpenFile))
                .with_context(|| anyhow!("File {}", output.clone()))?,
        );
        csv(writer, &counter, abundance)
            .with_context(|| Error::IO(ErrorDurringWrite))
            .with_context(|| anyhow!("File {} in csv format", output))?;

        log::info!("End of dump count data in csv");
    }

    if let Some(output) = solid_path.clone() {
        log::info!("Start of dump count data in solid format");

        let writer = std::io::BufWriter::new(
            std::fs::File::create(&output)
                .with_context(|| Error::IO(CantOpenFile))
                .with_context(|| anyhow!("File {}", output.clone()))?,
        );
        solid(writer, &counter, abundance)
            .with_context(|| Error::IO(ErrorDurringWrite))
            .with_context(|| anyhow!("File {} in solid format", output))?;

        log::info!("End of dump count data in solid format");
    }

    if let Some(output) = spectrum_path.clone() {
        log::info!("Start of dump count data in spectrum");

        let writer = std::io::BufWriter::new(
            std::fs::File::create(&output)
                .with_context(|| Error::IO(CantOpenFile))
                .with_context(|| anyhow!("File {}", output.clone()))?,
        );
        spectrum(writer, &counter)
            .with_context(|| Error::IO(ErrorDurringWrite))
            .with_context(|| anyhow!("File {} in spectrum format", output))?;

        log::info!("End of dump count data in spectrum");
    }

    Ok(())
}

/// Write in the given instance of io::Write the count in `counter` in binary format.
pub fn binary<W>(writer: W, counter: &counter::Counter) -> Result<()>
where
    W: std::io::Write,
{
    counter
        .serialize(writer)
        .with_context(|| Error::IO(ErrorDurringWrite))?;

    Ok(())
}

/// Write in the given instance of io::Write the count in `counter` in csv format.
/// Only count upper than `abundance` is write.
pub fn csv<W>(mut writer: W, counter: &counter::Counter, abundance: counter::Count) -> Result<()>
where
    W: std::io::Write,
{
    for (hash, count) in counter.get_raw_count().iter().enumerate() {
        let kmer = if cocktail::kmer::parity_even(hash as u64) {
            kmer::kmer2seq((hash as u64) << 1, counter.k)
        } else {
            kmer::kmer2seq(((hash as u64) << 1) ^ 0b1, counter.k)
        };

        let c = count.load(std::sync::atomic::Ordering::SeqCst);
        if c > abundance {
            writeln!(writer, "{},{}", kmer, c).with_context(|| Error::IO(ErrorDurringWrite))?;
        }
    }

    Ok(())
}

/// Serialize in the given instance of io::Write an instance of [solid::Solid] build from counts in `counter` upper than `abundance`.
pub fn solid<W>(writer: W, counter: &counter::Counter, abundance: counter::Count) -> Result<()>
where
    W: std::io::Write,
{
    let solid = solid::Solid::from_counter(&counter, abundance);

    solid
        .serialize(writer)
        .with_context(|| Error::IO(ErrorDurringWrite))?;

    Ok(())
}

/// Write in the given instance of io::Write the kmer spectrum from counts in `counter`
pub fn spectrum<W>(mut writer: W, counter: &counter::Counter) -> Result<()>
where
    W: std::io::Write,
{
    for (i, nb) in compute_spectrum(counter).iter().enumerate() {
        writeln!(writer, "{},{}", i, nb).with_context(|| Error::IO(ErrorDurringWrite))?;
    }

    Ok(())
}

pub fn compute_spectrum(counter: &counter::Counter) -> Box<[u64]> {
    let spectrum: Box<[std::sync::atomic::AtomicU64]> = (0..256)
        .map(|_| std::sync::atomic::AtomicU64::new(0))
        .collect::<Box<[std::sync::atomic::AtomicU64]>>();

    counter.get_raw_count().iter().par_bridge().for_each(|x| {
        let index = x.load(std::sync::atomic::Ordering::SeqCst) as usize;

        if spectrum[index].load(std::sync::atomic::Ordering::SeqCst) != std::u64::MAX {
            spectrum[x.load(std::sync::atomic::Ordering::SeqCst) as usize]
                .fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        }
    });

    spectrum
        .iter()
        .map(|x| x.load(std::sync::atomic::Ordering::SeqCst) as u64)
        .collect()
}

#[cfg(test)]
mod tests {
    lazy_static::lazy_static! {
    static ref COUNTER: crate::counter::Counter = {
            let mut counter = crate::counter::Counter::new(5, 1);

            for i in 0..cocktail::kmer::get_kmer_space_size(5) {
        counter.inc(i);
            }

            counter.inc(0);

            counter
    };
    }

    const CSV_ABUNDANCE_MIN_1: &[u8] = &[
        65, 65, 65, 65, 65, 44, 51, 10, 65, 65, 65, 65, 71, 44, 50, 10, 65, 65, 65, 67, 67, 44, 50,
        10, 65, 65, 65, 67, 84, 44, 50, 10, 65, 65, 65, 84, 67, 44, 50, 10, 65, 65, 65, 84, 84, 44,
        50, 10, 65, 65, 65, 71, 65, 44, 50, 10, 65, 65, 65, 71, 71, 44, 50, 10, 65, 65, 67, 65, 67,
        44, 50, 10, 65, 65, 67, 65, 84, 44, 50, 10, 65, 65, 67, 67, 65, 44, 50, 10, 65, 65, 67, 67,
        71, 44, 50, 10, 65, 65, 67, 84, 65, 44, 50, 10, 65, 65, 67, 84, 71, 44, 50, 10, 65, 65, 67,
        71, 67, 44, 50, 10, 65, 65, 67, 71, 84, 44, 50, 10, 65, 65, 84, 65, 67, 44, 50, 10, 65, 65,
        84, 65, 84, 44, 50, 10, 65, 65, 84, 67, 65, 44, 50, 10, 65, 65, 84, 67, 71, 44, 50, 10, 65,
        65, 84, 84, 65, 44, 50, 10, 65, 65, 84, 84, 71, 44, 50, 10, 65, 65, 84, 71, 67, 44, 50, 10,
        65, 65, 84, 71, 84, 44, 50, 10, 65, 65, 71, 65, 65, 44, 50, 10, 65, 65, 71, 65, 71, 44, 50,
        10, 65, 65, 71, 67, 67, 44, 50, 10, 65, 65, 71, 67, 84, 44, 50, 10, 65, 65, 71, 84, 67, 44,
        50, 10, 65, 65, 71, 84, 84, 44, 50, 10, 65, 65, 71, 71, 65, 44, 50, 10, 65, 65, 71, 71, 71,
        44, 50, 10, 65, 67, 65, 65, 67, 44, 50, 10, 65, 67, 65, 65, 84, 44, 50, 10, 65, 67, 65, 67,
        65, 44, 50, 10, 65, 67, 65, 67, 71, 44, 50, 10, 65, 67, 65, 84, 65, 44, 50, 10, 65, 67, 65,
        84, 71, 44, 50, 10, 65, 67, 65, 71, 67, 44, 50, 10, 65, 67, 65, 71, 84, 44, 50, 10, 65, 67,
        67, 65, 65, 44, 50, 10, 65, 67, 67, 65, 71, 44, 50, 10, 65, 67, 67, 67, 67, 44, 50, 10, 65,
        67, 67, 67, 84, 44, 50, 10, 65, 67, 67, 84, 67, 44, 50, 10, 65, 67, 67, 84, 84, 44, 50, 10,
        65, 67, 67, 71, 65, 44, 50, 10, 65, 67, 67, 71, 71, 44, 50, 10, 65, 67, 84, 65, 65, 44, 50,
        10, 65, 67, 84, 65, 71, 44, 50, 10, 65, 67, 84, 67, 67, 44, 50, 10, 65, 67, 84, 67, 84, 44,
        50, 10, 65, 67, 84, 84, 67, 44, 50, 10, 65, 67, 84, 84, 84, 44, 50, 10, 65, 67, 84, 71, 65,
        44, 50, 10, 65, 67, 84, 71, 71, 44, 50, 10, 65, 67, 71, 65, 67, 44, 50, 10, 65, 67, 71, 65,
        84, 44, 50, 10, 65, 67, 71, 67, 65, 44, 50, 10, 65, 67, 71, 67, 71, 44, 50, 10, 65, 67, 71,
        84, 65, 44, 50, 10, 65, 67, 71, 84, 71, 44, 50, 10, 65, 67, 71, 71, 67, 44, 50, 10, 65, 67,
        71, 71, 84, 44, 50, 10, 65, 84, 65, 65, 67, 44, 50, 10, 65, 84, 65, 65, 84, 44, 50, 10, 65,
        84, 65, 67, 65, 44, 50, 10, 65, 84, 65, 67, 71, 44, 50, 10, 65, 84, 65, 84, 65, 44, 50, 10,
        65, 84, 65, 84, 71, 44, 50, 10, 65, 84, 65, 71, 67, 44, 50, 10, 65, 84, 65, 71, 84, 44, 50,
        10, 65, 84, 67, 65, 65, 44, 50, 10, 65, 84, 67, 65, 71, 44, 50, 10, 65, 84, 67, 67, 67, 44,
        50, 10, 65, 84, 67, 67, 84, 44, 50, 10, 65, 84, 67, 84, 67, 44, 50, 10, 65, 84, 67, 84, 84,
        44, 50, 10, 65, 84, 67, 71, 65, 44, 50, 10, 65, 84, 67, 71, 71, 44, 50, 10, 65, 84, 84, 65,
        65, 44, 50, 10, 65, 84, 84, 65, 71, 44, 50, 10, 65, 84, 84, 67, 67, 44, 50, 10, 65, 84, 84,
        67, 84, 44, 50, 10, 65, 84, 84, 84, 67, 44, 50, 10, 65, 84, 84, 84, 84, 44, 50, 10, 65, 84,
        84, 71, 65, 44, 50, 10, 65, 84, 84, 71, 71, 44, 50, 10, 65, 84, 71, 65, 67, 44, 50, 10, 65,
        84, 71, 65, 84, 44, 50, 10, 65, 84, 71, 67, 65, 44, 50, 10, 65, 84, 71, 67, 71, 44, 50, 10,
        65, 84, 71, 84, 65, 44, 50, 10, 65, 84, 71, 84, 71, 44, 50, 10, 65, 84, 71, 71, 67, 44, 50,
        10, 65, 84, 71, 71, 84, 44, 50, 10, 65, 71, 65, 65, 65, 44, 50, 10, 65, 71, 65, 65, 71, 44,
        50, 10, 65, 71, 65, 67, 67, 44, 50, 10, 65, 71, 65, 67, 84, 44, 50, 10, 65, 71, 65, 84, 67,
        44, 50, 10, 65, 71, 65, 84, 84, 44, 50, 10, 65, 71, 65, 71, 65, 44, 50, 10, 65, 71, 65, 71,
        71, 44, 50, 10, 65, 71, 67, 65, 67, 44, 50, 10, 65, 71, 67, 65, 84, 44, 50, 10, 65, 71, 67,
        67, 65, 44, 50, 10, 65, 71, 67, 67, 71, 44, 50, 10, 65, 71, 67, 84, 65, 44, 50, 10, 65, 71,
        67, 84, 71, 44, 50, 10, 65, 71, 67, 71, 67, 44, 50, 10, 65, 71, 67, 71, 84, 44, 50, 10, 65,
        71, 84, 65, 67, 44, 50, 10, 65, 71, 84, 65, 84, 44, 50, 10, 65, 71, 84, 67, 65, 44, 50, 10,
        65, 71, 84, 67, 71, 44, 50, 10, 65, 71, 84, 84, 65, 44, 50, 10, 65, 71, 84, 84, 71, 44, 50,
        10, 65, 71, 84, 71, 67, 44, 50, 10, 65, 71, 84, 71, 84, 44, 50, 10, 65, 71, 71, 65, 65, 44,
        50, 10, 65, 71, 71, 65, 71, 44, 50, 10, 65, 71, 71, 67, 67, 44, 50, 10, 65, 71, 71, 67, 84,
        44, 50, 10, 65, 71, 71, 84, 67, 44, 50, 10, 65, 71, 71, 84, 84, 44, 50, 10, 65, 71, 71, 71,
        65, 44, 50, 10, 65, 71, 71, 71, 71, 44, 50, 10, 67, 65, 65, 65, 67, 44, 50, 10, 67, 65, 65,
        65, 84, 44, 50, 10, 67, 65, 65, 67, 65, 44, 50, 10, 67, 65, 65, 67, 71, 44, 50, 10, 67, 65,
        65, 84, 65, 44, 50, 10, 67, 65, 65, 84, 71, 44, 50, 10, 67, 65, 65, 71, 67, 44, 50, 10, 67,
        65, 65, 71, 84, 44, 50, 10, 67, 65, 67, 65, 65, 44, 50, 10, 67, 65, 67, 65, 71, 44, 50, 10,
        67, 65, 67, 67, 67, 44, 50, 10, 67, 65, 67, 67, 84, 44, 50, 10, 67, 65, 67, 84, 67, 44, 50,
        10, 67, 65, 67, 84, 84, 44, 50, 10, 67, 65, 67, 71, 65, 44, 50, 10, 67, 65, 67, 71, 71, 44,
        50, 10, 67, 65, 84, 65, 65, 44, 50, 10, 67, 65, 84, 65, 71, 44, 50, 10, 67, 65, 84, 67, 67,
        44, 50, 10, 67, 65, 84, 67, 84, 44, 50, 10, 67, 65, 84, 84, 67, 44, 50, 10, 67, 65, 84, 84,
        84, 44, 50, 10, 67, 65, 84, 71, 65, 44, 50, 10, 67, 65, 84, 71, 71, 44, 50, 10, 67, 65, 71,
        65, 67, 44, 50, 10, 67, 65, 71, 65, 84, 44, 50, 10, 67, 65, 71, 67, 65, 44, 50, 10, 67, 65,
        71, 67, 71, 44, 50, 10, 67, 65, 71, 84, 65, 44, 50, 10, 67, 65, 71, 84, 71, 44, 50, 10, 67,
        65, 71, 71, 67, 44, 50, 10, 67, 65, 71, 71, 84, 44, 50, 10, 67, 67, 65, 65, 65, 44, 50, 10,
        67, 67, 65, 65, 71, 44, 50, 10, 67, 67, 65, 67, 67, 44, 50, 10, 67, 67, 65, 67, 84, 44, 50,
        10, 67, 67, 65, 84, 67, 44, 50, 10, 67, 67, 65, 84, 84, 44, 50, 10, 67, 67, 65, 71, 65, 44,
        50, 10, 67, 67, 65, 71, 71, 44, 50, 10, 67, 67, 67, 65, 67, 44, 50, 10, 67, 67, 67, 65, 84,
        44, 50, 10, 67, 67, 67, 67, 65, 44, 50, 10, 67, 67, 67, 67, 71, 44, 50, 10, 67, 67, 67, 84,
        65, 44, 50, 10, 67, 67, 67, 84, 71, 44, 50, 10, 67, 67, 67, 71, 67, 44, 50, 10, 67, 67, 67,
        71, 84, 44, 50, 10, 67, 67, 84, 65, 67, 44, 50, 10, 67, 67, 84, 65, 84, 44, 50, 10, 67, 67,
        84, 67, 65, 44, 50, 10, 67, 67, 84, 67, 71, 44, 50, 10, 67, 67, 84, 84, 65, 44, 50, 10, 67,
        67, 84, 84, 71, 44, 50, 10, 67, 67, 84, 71, 67, 44, 50, 10, 67, 67, 84, 71, 84, 44, 50, 10,
        67, 67, 71, 65, 65, 44, 50, 10, 67, 67, 71, 65, 71, 44, 50, 10, 67, 67, 71, 67, 67, 44, 50,
        10, 67, 67, 71, 67, 84, 44, 50, 10, 67, 67, 71, 84, 67, 44, 50, 10, 67, 67, 71, 84, 84, 44,
        50, 10, 67, 67, 71, 71, 65, 44, 50, 10, 67, 67, 71, 71, 71, 44, 50, 10, 67, 84, 65, 65, 65,
        44, 50, 10, 67, 84, 65, 65, 71, 44, 50, 10, 67, 84, 65, 67, 67, 44, 50, 10, 67, 84, 65, 67,
        84, 44, 50, 10, 67, 84, 65, 84, 67, 44, 50, 10, 67, 84, 65, 84, 84, 44, 50, 10, 67, 84, 65,
        71, 65, 44, 50, 10, 67, 84, 65, 71, 71, 44, 50, 10, 67, 84, 67, 65, 67, 44, 50, 10, 67, 84,
        67, 65, 84, 44, 50, 10, 67, 84, 67, 67, 65, 44, 50, 10, 67, 84, 67, 67, 71, 44, 50, 10, 67,
        84, 67, 84, 65, 44, 50, 10, 67, 84, 67, 84, 71, 44, 50, 10, 67, 84, 67, 71, 67, 44, 50, 10,
        67, 84, 67, 71, 84, 44, 50, 10, 67, 84, 84, 65, 67, 44, 50, 10, 67, 84, 84, 65, 84, 44, 50,
        10, 67, 84, 84, 67, 65, 44, 50, 10, 67, 84, 84, 67, 71, 44, 50, 10, 67, 84, 84, 84, 65, 44,
        50, 10, 67, 84, 84, 84, 71, 44, 50, 10, 67, 84, 84, 71, 67, 44, 50, 10, 67, 84, 84, 71, 84,
        44, 50, 10, 67, 84, 71, 65, 65, 44, 50, 10, 67, 84, 71, 65, 71, 44, 50, 10, 67, 84, 71, 67,
        67, 44, 50, 10, 67, 84, 71, 67, 84, 44, 50, 10, 67, 84, 71, 84, 67, 44, 50, 10, 67, 84, 71,
        84, 84, 44, 50, 10, 67, 84, 71, 71, 65, 44, 50, 10, 67, 84, 71, 71, 71, 44, 50, 10, 67, 71,
        65, 65, 67, 44, 50, 10, 67, 71, 65, 65, 84, 44, 50, 10, 67, 71, 65, 67, 65, 44, 50, 10, 67,
        71, 65, 67, 71, 44, 50, 10, 67, 71, 65, 84, 65, 44, 50, 10, 67, 71, 65, 84, 71, 44, 50, 10,
        67, 71, 65, 71, 67, 44, 50, 10, 67, 71, 65, 71, 84, 44, 50, 10, 67, 71, 67, 65, 65, 44, 50,
        10, 67, 71, 67, 65, 71, 44, 50, 10, 67, 71, 67, 67, 67, 44, 50, 10, 67, 71, 67, 67, 84, 44,
        50, 10, 67, 71, 67, 84, 67, 44, 50, 10, 67, 71, 67, 84, 84, 44, 50, 10, 67, 71, 67, 71, 65,
        44, 50, 10, 67, 71, 67, 71, 71, 44, 50, 10, 67, 71, 84, 65, 65, 44, 50, 10, 67, 71, 84, 65,
        71, 44, 50, 10, 67, 71, 84, 67, 67, 44, 50, 10, 67, 71, 84, 67, 84, 44, 50, 10, 67, 71, 84,
        84, 67, 44, 50, 10, 67, 71, 84, 84, 84, 44, 50, 10, 67, 71, 84, 71, 65, 44, 50, 10, 67, 71,
        84, 71, 71, 44, 50, 10, 67, 71, 71, 65, 67, 44, 50, 10, 67, 71, 71, 65, 84, 44, 50, 10, 67,
        71, 71, 67, 65, 44, 50, 10, 67, 71, 71, 67, 71, 44, 50, 10, 67, 71, 71, 84, 65, 44, 50, 10,
        67, 71, 71, 84, 71, 44, 50, 10, 67, 71, 71, 71, 67, 44, 50, 10, 67, 71, 71, 71, 84, 44, 50,
        10, 84, 65, 65, 65, 67, 44, 50, 10, 84, 65, 65, 65, 84, 44, 50, 10, 84, 65, 65, 67, 65, 44,
        50, 10, 84, 65, 65, 67, 71, 44, 50, 10, 84, 65, 65, 84, 65, 44, 50, 10, 84, 65, 65, 84, 71,
        44, 50, 10, 84, 65, 65, 71, 67, 44, 50, 10, 84, 65, 65, 71, 84, 44, 50, 10, 84, 65, 67, 65,
        65, 44, 50, 10, 84, 65, 67, 65, 71, 44, 50, 10, 84, 65, 67, 67, 67, 44, 50, 10, 84, 65, 67,
        67, 84, 44, 50, 10, 84, 65, 67, 84, 67, 44, 50, 10, 84, 65, 67, 84, 84, 44, 50, 10, 84, 65,
        67, 71, 65, 44, 50, 10, 84, 65, 67, 71, 71, 44, 50, 10, 84, 65, 84, 65, 65, 44, 50, 10, 84,
        65, 84, 65, 71, 44, 50, 10, 84, 65, 84, 67, 67, 44, 50, 10, 84, 65, 84, 67, 84, 44, 50, 10,
        84, 65, 84, 84, 67, 44, 50, 10, 84, 65, 84, 84, 84, 44, 50, 10, 84, 65, 84, 71, 65, 44, 50,
        10, 84, 65, 84, 71, 71, 44, 50, 10, 84, 65, 71, 65, 67, 44, 50, 10, 84, 65, 71, 65, 84, 44,
        50, 10, 84, 65, 71, 67, 65, 44, 50, 10, 84, 65, 71, 67, 71, 44, 50, 10, 84, 65, 71, 84, 65,
        44, 50, 10, 84, 65, 71, 84, 71, 44, 50, 10, 84, 65, 71, 71, 67, 44, 50, 10, 84, 65, 71, 71,
        84, 44, 50, 10, 84, 67, 65, 65, 65, 44, 50, 10, 84, 67, 65, 65, 71, 44, 50, 10, 84, 67, 65,
        67, 67, 44, 50, 10, 84, 67, 65, 67, 84, 44, 50, 10, 84, 67, 65, 84, 67, 44, 50, 10, 84, 67,
        65, 84, 84, 44, 50, 10, 84, 67, 65, 71, 65, 44, 50, 10, 84, 67, 65, 71, 71, 44, 50, 10, 84,
        67, 67, 65, 67, 44, 50, 10, 84, 67, 67, 65, 84, 44, 50, 10, 84, 67, 67, 67, 65, 44, 50, 10,
        84, 67, 67, 67, 71, 44, 50, 10, 84, 67, 67, 84, 65, 44, 50, 10, 84, 67, 67, 84, 71, 44, 50,
        10, 84, 67, 67, 71, 67, 44, 50, 10, 84, 67, 67, 71, 84, 44, 50, 10, 84, 67, 84, 65, 67, 44,
        50, 10, 84, 67, 84, 65, 84, 44, 50, 10, 84, 67, 84, 67, 65, 44, 50, 10, 84, 67, 84, 67, 71,
        44, 50, 10, 84, 67, 84, 84, 65, 44, 50, 10, 84, 67, 84, 84, 71, 44, 50, 10, 84, 67, 84, 71,
        67, 44, 50, 10, 84, 67, 84, 71, 84, 44, 50, 10, 84, 67, 71, 65, 65, 44, 50, 10, 84, 67, 71,
        65, 71, 44, 50, 10, 84, 67, 71, 67, 67, 44, 50, 10, 84, 67, 71, 67, 84, 44, 50, 10, 84, 67,
        71, 84, 67, 44, 50, 10, 84, 67, 71, 84, 84, 44, 50, 10, 84, 67, 71, 71, 65, 44, 50, 10, 84,
        67, 71, 71, 71, 44, 50, 10, 84, 84, 65, 65, 65, 44, 50, 10, 84, 84, 65, 65, 71, 44, 50, 10,
        84, 84, 65, 67, 67, 44, 50, 10, 84, 84, 65, 67, 84, 44, 50, 10, 84, 84, 65, 84, 67, 44, 50,
        10, 84, 84, 65, 84, 84, 44, 50, 10, 84, 84, 65, 71, 65, 44, 50, 10, 84, 84, 65, 71, 71, 44,
        50, 10, 84, 84, 67, 65, 67, 44, 50, 10, 84, 84, 67, 65, 84, 44, 50, 10, 84, 84, 67, 67, 65,
        44, 50, 10, 84, 84, 67, 67, 71, 44, 50, 10, 84, 84, 67, 84, 65, 44, 50, 10, 84, 84, 67, 84,
        71, 44, 50, 10, 84, 84, 67, 71, 67, 44, 50, 10, 84, 84, 67, 71, 84, 44, 50, 10, 84, 84, 84,
        65, 67, 44, 50, 10, 84, 84, 84, 65, 84, 44, 50, 10, 84, 84, 84, 67, 65, 44, 50, 10, 84, 84,
        84, 67, 71, 44, 50, 10, 84, 84, 84, 84, 65, 44, 50, 10, 84, 84, 84, 84, 71, 44, 50, 10, 84,
        84, 84, 71, 67, 44, 50, 10, 84, 84, 84, 71, 84, 44, 50, 10, 84, 84, 71, 65, 65, 44, 50, 10,
        84, 84, 71, 65, 71, 44, 50, 10, 84, 84, 71, 67, 67, 44, 50, 10, 84, 84, 71, 67, 84, 44, 50,
        10, 84, 84, 71, 84, 67, 44, 50, 10, 84, 84, 71, 84, 84, 44, 50, 10, 84, 84, 71, 71, 65, 44,
        50, 10, 84, 84, 71, 71, 71, 44, 50, 10, 84, 71, 65, 65, 67, 44, 50, 10, 84, 71, 65, 65, 84,
        44, 50, 10, 84, 71, 65, 67, 65, 44, 50, 10, 84, 71, 65, 67, 71, 44, 50, 10, 84, 71, 65, 84,
        65, 44, 50, 10, 84, 71, 65, 84, 71, 44, 50, 10, 84, 71, 65, 71, 67, 44, 50, 10, 84, 71, 65,
        71, 84, 44, 50, 10, 84, 71, 67, 65, 65, 44, 50, 10, 84, 71, 67, 65, 71, 44, 50, 10, 84, 71,
        67, 67, 67, 44, 50, 10, 84, 71, 67, 67, 84, 44, 50, 10, 84, 71, 67, 84, 67, 44, 50, 10, 84,
        71, 67, 84, 84, 44, 50, 10, 84, 71, 67, 71, 65, 44, 50, 10, 84, 71, 67, 71, 71, 44, 50, 10,
        84, 71, 84, 65, 65, 44, 50, 10, 84, 71, 84, 65, 71, 44, 50, 10, 84, 71, 84, 67, 67, 44, 50,
        10, 84, 71, 84, 67, 84, 44, 50, 10, 84, 71, 84, 84, 67, 44, 50, 10, 84, 71, 84, 84, 84, 44,
        50, 10, 84, 71, 84, 71, 65, 44, 50, 10, 84, 71, 84, 71, 71, 44, 50, 10, 84, 71, 71, 65, 67,
        44, 50, 10, 84, 71, 71, 65, 84, 44, 50, 10, 84, 71, 71, 67, 65, 44, 50, 10, 84, 71, 71, 67,
        71, 44, 50, 10, 84, 71, 71, 84, 65, 44, 50, 10, 84, 71, 71, 84, 71, 44, 50, 10, 84, 71, 71,
        71, 67, 44, 50, 10, 84, 71, 71, 71, 84, 44, 50, 10, 71, 65, 65, 65, 65, 44, 50, 10, 71, 65,
        65, 65, 71, 44, 50, 10, 71, 65, 65, 67, 67, 44, 50, 10, 71, 65, 65, 67, 84, 44, 50, 10, 71,
        65, 65, 84, 67, 44, 50, 10, 71, 65, 65, 84, 84, 44, 50, 10, 71, 65, 65, 71, 65, 44, 50, 10,
        71, 65, 65, 71, 71, 44, 50, 10, 71, 65, 67, 65, 67, 44, 50, 10, 71, 65, 67, 65, 84, 44, 50,
        10, 71, 65, 67, 67, 65, 44, 50, 10, 71, 65, 67, 67, 71, 44, 50, 10, 71, 65, 67, 84, 65, 44,
        50, 10, 71, 65, 67, 84, 71, 44, 50, 10, 71, 65, 67, 71, 67, 44, 50, 10, 71, 65, 67, 71, 84,
        44, 50, 10, 71, 65, 84, 65, 67, 44, 50, 10, 71, 65, 84, 65, 84, 44, 50, 10, 71, 65, 84, 67,
        65, 44, 50, 10, 71, 65, 84, 67, 71, 44, 50, 10, 71, 65, 84, 84, 65, 44, 50, 10, 71, 65, 84,
        84, 71, 44, 50, 10, 71, 65, 84, 71, 67, 44, 50, 10, 71, 65, 84, 71, 84, 44, 50, 10, 71, 65,
        71, 65, 65, 44, 50, 10, 71, 65, 71, 65, 71, 44, 50, 10, 71, 65, 71, 67, 67, 44, 50, 10, 71,
        65, 71, 67, 84, 44, 50, 10, 71, 65, 71, 84, 67, 44, 50, 10, 71, 65, 71, 84, 84, 44, 50, 10,
        71, 65, 71, 71, 65, 44, 50, 10, 71, 65, 71, 71, 71, 44, 50, 10, 71, 67, 65, 65, 67, 44, 50,
        10, 71, 67, 65, 65, 84, 44, 50, 10, 71, 67, 65, 67, 65, 44, 50, 10, 71, 67, 65, 67, 71, 44,
        50, 10, 71, 67, 65, 84, 65, 44, 50, 10, 71, 67, 65, 84, 71, 44, 50, 10, 71, 67, 65, 71, 67,
        44, 50, 10, 71, 67, 65, 71, 84, 44, 50, 10, 71, 67, 67, 65, 65, 44, 50, 10, 71, 67, 67, 65,
        71, 44, 50, 10, 71, 67, 67, 67, 67, 44, 50, 10, 71, 67, 67, 67, 84, 44, 50, 10, 71, 67, 67,
        84, 67, 44, 50, 10, 71, 67, 67, 84, 84, 44, 50, 10, 71, 67, 67, 71, 65, 44, 50, 10, 71, 67,
        67, 71, 71, 44, 50, 10, 71, 67, 84, 65, 65, 44, 50, 10, 71, 67, 84, 65, 71, 44, 50, 10, 71,
        67, 84, 67, 67, 44, 50, 10, 71, 67, 84, 67, 84, 44, 50, 10, 71, 67, 84, 84, 67, 44, 50, 10,
        71, 67, 84, 84, 84, 44, 50, 10, 71, 67, 84, 71, 65, 44, 50, 10, 71, 67, 84, 71, 71, 44, 50,
        10, 71, 67, 71, 65, 67, 44, 50, 10, 71, 67, 71, 65, 84, 44, 50, 10, 71, 67, 71, 67, 65, 44,
        50, 10, 71, 67, 71, 67, 71, 44, 50, 10, 71, 67, 71, 84, 65, 44, 50, 10, 71, 67, 71, 84, 71,
        44, 50, 10, 71, 67, 71, 71, 67, 44, 50, 10, 71, 67, 71, 71, 84, 44, 50, 10, 71, 84, 65, 65,
        67, 44, 50, 10, 71, 84, 65, 65, 84, 44, 50, 10, 71, 84, 65, 67, 65, 44, 50, 10, 71, 84, 65,
        67, 71, 44, 50, 10, 71, 84, 65, 84, 65, 44, 50, 10, 71, 84, 65, 84, 71, 44, 50, 10, 71, 84,
        65, 71, 67, 44, 50, 10, 71, 84, 65, 71, 84, 44, 50, 10, 71, 84, 67, 65, 65, 44, 50, 10, 71,
        84, 67, 65, 71, 44, 50, 10, 71, 84, 67, 67, 67, 44, 50, 10, 71, 84, 67, 67, 84, 44, 50, 10,
        71, 84, 67, 84, 67, 44, 50, 10, 71, 84, 67, 84, 84, 44, 50, 10, 71, 84, 67, 71, 65, 44, 50,
        10, 71, 84, 67, 71, 71, 44, 50, 10, 71, 84, 84, 65, 65, 44, 50, 10, 71, 84, 84, 65, 71, 44,
        50, 10, 71, 84, 84, 67, 67, 44, 50, 10, 71, 84, 84, 67, 84, 44, 50, 10, 71, 84, 84, 84, 67,
        44, 50, 10, 71, 84, 84, 84, 84, 44, 50, 10, 71, 84, 84, 71, 65, 44, 50, 10, 71, 84, 84, 71,
        71, 44, 50, 10, 71, 84, 71, 65, 67, 44, 50, 10, 71, 84, 71, 65, 84, 44, 50, 10, 71, 84, 71,
        67, 65, 44, 50, 10, 71, 84, 71, 67, 71, 44, 50, 10, 71, 84, 71, 84, 65, 44, 50, 10, 71, 84,
        71, 84, 71, 44, 50, 10, 71, 84, 71, 71, 67, 44, 50, 10, 71, 84, 71, 71, 84, 44, 50, 10, 71,
        71, 65, 65, 65, 44, 50, 10, 71, 71, 65, 65, 71, 44, 50, 10, 71, 71, 65, 67, 67, 44, 50, 10,
        71, 71, 65, 67, 84, 44, 50, 10, 71, 71, 65, 84, 67, 44, 50, 10, 71, 71, 65, 84, 84, 44, 50,
        10, 71, 71, 65, 71, 65, 44, 50, 10, 71, 71, 65, 71, 71, 44, 50, 10, 71, 71, 67, 65, 67, 44,
        50, 10, 71, 71, 67, 65, 84, 44, 50, 10, 71, 71, 67, 67, 65, 44, 50, 10, 71, 71, 67, 67, 71,
        44, 50, 10, 71, 71, 67, 84, 65, 44, 50, 10, 71, 71, 67, 84, 71, 44, 50, 10, 71, 71, 67, 71,
        67, 44, 50, 10, 71, 71, 67, 71, 84, 44, 50, 10, 71, 71, 84, 65, 67, 44, 50, 10, 71, 71, 84,
        65, 84, 44, 50, 10, 71, 71, 84, 67, 65, 44, 50, 10, 71, 71, 84, 67, 71, 44, 50, 10, 71, 71,
        84, 84, 65, 44, 50, 10, 71, 71, 84, 84, 71, 44, 50, 10, 71, 71, 84, 71, 67, 44, 50, 10, 71,
        71, 84, 71, 84, 44, 50, 10, 71, 71, 71, 65, 65, 44, 50, 10, 71, 71, 71, 65, 71, 44, 50, 10,
        71, 71, 71, 67, 67, 44, 50, 10, 71, 71, 71, 67, 84, 44, 50, 10, 71, 71, 71, 84, 67, 44, 50,
        10, 71, 71, 71, 84, 84, 44, 50, 10, 71, 71, 71, 71, 65, 44, 50, 10, 71, 71, 71, 71, 71, 44,
        50, 10,
    ];

    const CSV_ABUNDANCE_MIN_2: &[u8] = &[65, 65, 65, 65, 65, 44, 51, 10];

    #[test]
    fn csv() {
        let mut outfile = Vec::new();
        let counter = &*COUNTER;

        crate::dump::csv(&mut outfile, counter, 1).unwrap();
        assert_eq!(&outfile[..], &CSV_ABUNDANCE_MIN_1[..]);

        outfile.clear();

        crate::dump::csv(&mut outfile, counter, 2).unwrap();
        assert_eq!(&outfile[..], &CSV_ABUNDANCE_MIN_2[..]);
    }

    const SOLID_ABUNDANCE_MIN_1: &[u8] = &[
        5, 0, 0, 2, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255,
    ];

    const SOLID_ABUNDANCE_MIN_2: &[u8] = &[
        5, 0, 0, 2, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ];

    #[test]
    fn solid() {
        let mut outfile = Vec::new();
        let counter = &*COUNTER;

        crate::dump::solid(&mut outfile, counter, 1).unwrap();
        assert_eq!(&outfile[..], &SOLID_ABUNDANCE_MIN_1[..]);

        outfile.clear();

        crate::dump::solid(&mut outfile, counter, 2).unwrap();
        assert_eq!(&outfile[..], &SOLID_ABUNDANCE_MIN_2[..]);
    }

    const SPECTRUM_ABUNDANCE_MIN_1: &[u8] = &[
        48, 44, 48, 10, 49, 44, 48, 10, 50, 44, 53, 49, 49, 10, 51, 44, 49, 10, 52, 44, 48, 10, 53,
        44, 48, 10, 54, 44, 48, 10, 55, 44, 48, 10, 56, 44, 48, 10, 57, 44, 48, 10, 49, 48, 44, 48,
        10, 49, 49, 44, 48, 10, 49, 50, 44, 48, 10, 49, 51, 44, 48, 10, 49, 52, 44, 48, 10, 49, 53,
        44, 48, 10, 49, 54, 44, 48, 10, 49, 55, 44, 48, 10, 49, 56, 44, 48, 10, 49, 57, 44, 48, 10,
        50, 48, 44, 48, 10, 50, 49, 44, 48, 10, 50, 50, 44, 48, 10, 50, 51, 44, 48, 10, 50, 52, 44,
        48, 10, 50, 53, 44, 48, 10, 50, 54, 44, 48, 10, 50, 55, 44, 48, 10, 50, 56, 44, 48, 10, 50,
        57, 44, 48, 10, 51, 48, 44, 48, 10, 51, 49, 44, 48, 10, 51, 50, 44, 48, 10, 51, 51, 44, 48,
        10, 51, 52, 44, 48, 10, 51, 53, 44, 48, 10, 51, 54, 44, 48, 10, 51, 55, 44, 48, 10, 51, 56,
        44, 48, 10, 51, 57, 44, 48, 10, 52, 48, 44, 48, 10, 52, 49, 44, 48, 10, 52, 50, 44, 48, 10,
        52, 51, 44, 48, 10, 52, 52, 44, 48, 10, 52, 53, 44, 48, 10, 52, 54, 44, 48, 10, 52, 55, 44,
        48, 10, 52, 56, 44, 48, 10, 52, 57, 44, 48, 10, 53, 48, 44, 48, 10, 53, 49, 44, 48, 10, 53,
        50, 44, 48, 10, 53, 51, 44, 48, 10, 53, 52, 44, 48, 10, 53, 53, 44, 48, 10, 53, 54, 44, 48,
        10, 53, 55, 44, 48, 10, 53, 56, 44, 48, 10, 53, 57, 44, 48, 10, 54, 48, 44, 48, 10, 54, 49,
        44, 48, 10, 54, 50, 44, 48, 10, 54, 51, 44, 48, 10, 54, 52, 44, 48, 10, 54, 53, 44, 48, 10,
        54, 54, 44, 48, 10, 54, 55, 44, 48, 10, 54, 56, 44, 48, 10, 54, 57, 44, 48, 10, 55, 48, 44,
        48, 10, 55, 49, 44, 48, 10, 55, 50, 44, 48, 10, 55, 51, 44, 48, 10, 55, 52, 44, 48, 10, 55,
        53, 44, 48, 10, 55, 54, 44, 48, 10, 55, 55, 44, 48, 10, 55, 56, 44, 48, 10, 55, 57, 44, 48,
        10, 56, 48, 44, 48, 10, 56, 49, 44, 48, 10, 56, 50, 44, 48, 10, 56, 51, 44, 48, 10, 56, 52,
        44, 48, 10, 56, 53, 44, 48, 10, 56, 54, 44, 48, 10, 56, 55, 44, 48, 10, 56, 56, 44, 48, 10,
        56, 57, 44, 48, 10, 57, 48, 44, 48, 10, 57, 49, 44, 48, 10, 57, 50, 44, 48, 10, 57, 51, 44,
        48, 10, 57, 52, 44, 48, 10, 57, 53, 44, 48, 10, 57, 54, 44, 48, 10, 57, 55, 44, 48, 10, 57,
        56, 44, 48, 10, 57, 57, 44, 48, 10, 49, 48, 48, 44, 48, 10, 49, 48, 49, 44, 48, 10, 49, 48,
        50, 44, 48, 10, 49, 48, 51, 44, 48, 10, 49, 48, 52, 44, 48, 10, 49, 48, 53, 44, 48, 10, 49,
        48, 54, 44, 48, 10, 49, 48, 55, 44, 48, 10, 49, 48, 56, 44, 48, 10, 49, 48, 57, 44, 48, 10,
        49, 49, 48, 44, 48, 10, 49, 49, 49, 44, 48, 10, 49, 49, 50, 44, 48, 10, 49, 49, 51, 44, 48,
        10, 49, 49, 52, 44, 48, 10, 49, 49, 53, 44, 48, 10, 49, 49, 54, 44, 48, 10, 49, 49, 55, 44,
        48, 10, 49, 49, 56, 44, 48, 10, 49, 49, 57, 44, 48, 10, 49, 50, 48, 44, 48, 10, 49, 50, 49,
        44, 48, 10, 49, 50, 50, 44, 48, 10, 49, 50, 51, 44, 48, 10, 49, 50, 52, 44, 48, 10, 49, 50,
        53, 44, 48, 10, 49, 50, 54, 44, 48, 10, 49, 50, 55, 44, 48, 10, 49, 50, 56, 44, 48, 10, 49,
        50, 57, 44, 48, 10, 49, 51, 48, 44, 48, 10, 49, 51, 49, 44, 48, 10, 49, 51, 50, 44, 48, 10,
        49, 51, 51, 44, 48, 10, 49, 51, 52, 44, 48, 10, 49, 51, 53, 44, 48, 10, 49, 51, 54, 44, 48,
        10, 49, 51, 55, 44, 48, 10, 49, 51, 56, 44, 48, 10, 49, 51, 57, 44, 48, 10, 49, 52, 48, 44,
        48, 10, 49, 52, 49, 44, 48, 10, 49, 52, 50, 44, 48, 10, 49, 52, 51, 44, 48, 10, 49, 52, 52,
        44, 48, 10, 49, 52, 53, 44, 48, 10, 49, 52, 54, 44, 48, 10, 49, 52, 55, 44, 48, 10, 49, 52,
        56, 44, 48, 10, 49, 52, 57, 44, 48, 10, 49, 53, 48, 44, 48, 10, 49, 53, 49, 44, 48, 10, 49,
        53, 50, 44, 48, 10, 49, 53, 51, 44, 48, 10, 49, 53, 52, 44, 48, 10, 49, 53, 53, 44, 48, 10,
        49, 53, 54, 44, 48, 10, 49, 53, 55, 44, 48, 10, 49, 53, 56, 44, 48, 10, 49, 53, 57, 44, 48,
        10, 49, 54, 48, 44, 48, 10, 49, 54, 49, 44, 48, 10, 49, 54, 50, 44, 48, 10, 49, 54, 51, 44,
        48, 10, 49, 54, 52, 44, 48, 10, 49, 54, 53, 44, 48, 10, 49, 54, 54, 44, 48, 10, 49, 54, 55,
        44, 48, 10, 49, 54, 56, 44, 48, 10, 49, 54, 57, 44, 48, 10, 49, 55, 48, 44, 48, 10, 49, 55,
        49, 44, 48, 10, 49, 55, 50, 44, 48, 10, 49, 55, 51, 44, 48, 10, 49, 55, 52, 44, 48, 10, 49,
        55, 53, 44, 48, 10, 49, 55, 54, 44, 48, 10, 49, 55, 55, 44, 48, 10, 49, 55, 56, 44, 48, 10,
        49, 55, 57, 44, 48, 10, 49, 56, 48, 44, 48, 10, 49, 56, 49, 44, 48, 10, 49, 56, 50, 44, 48,
        10, 49, 56, 51, 44, 48, 10, 49, 56, 52, 44, 48, 10, 49, 56, 53, 44, 48, 10, 49, 56, 54, 44,
        48, 10, 49, 56, 55, 44, 48, 10, 49, 56, 56, 44, 48, 10, 49, 56, 57, 44, 48, 10, 49, 57, 48,
        44, 48, 10, 49, 57, 49, 44, 48, 10, 49, 57, 50, 44, 48, 10, 49, 57, 51, 44, 48, 10, 49, 57,
        52, 44, 48, 10, 49, 57, 53, 44, 48, 10, 49, 57, 54, 44, 48, 10, 49, 57, 55, 44, 48, 10, 49,
        57, 56, 44, 48, 10, 49, 57, 57, 44, 48, 10, 50, 48, 48, 44, 48, 10, 50, 48, 49, 44, 48, 10,
        50, 48, 50, 44, 48, 10, 50, 48, 51, 44, 48, 10, 50, 48, 52, 44, 48, 10, 50, 48, 53, 44, 48,
        10, 50, 48, 54, 44, 48, 10, 50, 48, 55, 44, 48, 10, 50, 48, 56, 44, 48, 10, 50, 48, 57, 44,
        48, 10, 50, 49, 48, 44, 48, 10, 50, 49, 49, 44, 48, 10, 50, 49, 50, 44, 48, 10, 50, 49, 51,
        44, 48, 10, 50, 49, 52, 44, 48, 10, 50, 49, 53, 44, 48, 10, 50, 49, 54, 44, 48, 10, 50, 49,
        55, 44, 48, 10, 50, 49, 56, 44, 48, 10, 50, 49, 57, 44, 48, 10, 50, 50, 48, 44, 48, 10, 50,
        50, 49, 44, 48, 10, 50, 50, 50, 44, 48, 10, 50, 50, 51, 44, 48, 10, 50, 50, 52, 44, 48, 10,
        50, 50, 53, 44, 48, 10, 50, 50, 54, 44, 48, 10, 50, 50, 55, 44, 48, 10, 50, 50, 56, 44, 48,
        10, 50, 50, 57, 44, 48, 10, 50, 51, 48, 44, 48, 10, 50, 51, 49, 44, 48, 10, 50, 51, 50, 44,
        48, 10, 50, 51, 51, 44, 48, 10, 50, 51, 52, 44, 48, 10, 50, 51, 53, 44, 48, 10, 50, 51, 54,
        44, 48, 10, 50, 51, 55, 44, 48, 10, 50, 51, 56, 44, 48, 10, 50, 51, 57, 44, 48, 10, 50, 52,
        48, 44, 48, 10, 50, 52, 49, 44, 48, 10, 50, 52, 50, 44, 48, 10, 50, 52, 51, 44, 48, 10, 50,
        52, 52, 44, 48, 10, 50, 52, 53, 44, 48, 10, 50, 52, 54, 44, 48, 10, 50, 52, 55, 44, 48, 10,
        50, 52, 56, 44, 48, 10, 50, 52, 57, 44, 48, 10, 50, 53, 48, 44, 48, 10, 50, 53, 49, 44, 48,
        10, 50, 53, 50, 44, 48, 10, 50, 53, 51, 44, 48, 10, 50, 53, 52, 44, 48, 10, 50, 53, 53, 44,
        48, 10,
    ];

    #[test]
    fn spectrum() {
        let mut outfile = Vec::new();
        let counter = &*COUNTER;

        crate::dump::spectrum(&mut outfile, counter).unwrap();

        assert_eq!(&outfile[..], &SPECTRUM_ABUNDANCE_MIN_1[..]);
    }
}
