//! Run count command

/* std use */

/* crate use */

/* project use */
use crate::cli;
use crate::counter;
use crate::error;

/// Run count
pub fn count(params: cli::Count) -> error::Result<()> {
    log::info!("Start init counter");
    let mut counter = counter::Counter::<crate::CountType>::new(params.kmer_size());
    log::info!("End init counter");

    log::info!("Start count kmer");
    match params.format() {
        cli::Format::Fasta => counter.count_fasta(params.inputs()?, params.record_buffer()),
        #[cfg(feature = "fastq")]
        cli::Format::Fastq => counter.count_fastq(params.inputs()?, params.record_buffer()),
    }
    log::info!("End count kmer");

    let serialize = counter.serialize();

    for (out_type, output) in params.outputs().into_iter() {
        match out_type {
            cli::DumpType::Pcon => {
                log::info!("Start write count in pcon format");
                serialize.pcon(output?)?;
                log::info!("End write count in pcon format");
            }
            cli::DumpType::Csv => {
                log::info!("Start write count in csv format");
                serialize.csv(params.abundance(), output?)?;
                log::info!("End write count in csv format");
            }
            cli::DumpType::Solid => {
                log::info!("Start write count in solid format");
                serialize.solid(params.abundance(), output?)?;
                log::info!("End write count in solid format");
            }
        }
    }

    Ok(())
}
