//! Run minicount command

/* std use */

/* crate use */

/* project use */
use crate::cli;
use crate::error;
use crate::minicounter;

/// Run count
pub fn minicount(params: cli::MiniCount) -> error::Result<()> {
    log::info!("Start init counter");
    let mut counter = minicounter::MiniCounter::<crate::CountType>::new(
        params.kmer_size(),
        params.minimizer_size(),
        params.abundance(),
    );
    log::info!("End init counter");

    log::info!("Start count kmer");
    match params.format() {
        cli::Format::Fasta => counter.count_fasta(params.inputs()?, params.record_buffer()),
        #[cfg(feature = "fastq")]
        cli::Format::Fastq => counter.count_fastq(params.inputs()?, params.record_buffer()),
    }
    log::info!("End count kmer");

    for (out_type, output) in params.outputs().into_iter() {
        match out_type {
            cli::DumpType::Csv => {
                log::info!("Start write count in csv format");
                counter.serialize(params.abundance(), output?)?;
                log::info!("End write count in csv format");
            }
            _ => log::error!("Only csv dump is available for minicount"),
        }
    }

    Ok(())
}
