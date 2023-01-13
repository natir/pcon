//! Run count command

/* std use */

/* crate use */

/* project use */
use crate::cli;
use crate::counter;
use crate::error;

cfg_if::cfg_if! {
    if #[cfg(feature = "parallel")] {
    /// Run count
    pub fn count(params: cli::Count) -> error::Result<()> {
            log::info!("Start init counter");
            let mut counter = counter::Counter::<std::sync::atomic::AtomicU8>::new(params.kmer_size());
            log::info!("End init counter");

            log::info!("Start count kmer");
            counter.count_fasta(params.inputs()?, params.record_buffer());
        log::info!("End count kmer");

    log::info!("Start write count");
    let serialize = counter.serialize();

    match params.dump() {
        cli::DumpType::Pcon => serialize.pcon(params.output()?)?,
        cli::DumpType::Csv => serialize.csv(params.abundance(), params.output()?)?,
        cli::DumpType::Solid => serialize.solid(params.abundance(), params.output()?)?,
        cli::DumpType::Spectrum => todo!(),
    }
    log::info!("End write count");

        Ok(())
    }
    } else {
    /// Run count
    pub fn count(params: cli::Count) -> error::Result<()> {
        log::info!("Start init counter");
        let mut counter = counter::Counter::<u8>::new(params.kmer_size());
        log::info!("End init counter");

        log::info!("Start count kmer");
        counter.count_fasta(params.inputs()?, params.record_buffer());
        log::info!("End count kmer");

    log::info!("Start write count");
        let serialize = counter.serialize();

    match params.dump() {
        cli::DumpType::Pcon => serialize.pcon(params.output()?)?,
        cli::DumpType::Csv => serialize.csv(params.abundance(), params.output()?)?,
        cli::DumpType::Solid => serialize.solid(params.abundance(), params.output()?)?,
        cli::DumpType::Spectrum => todo!(),
    }
    log::info!("End write count");

        Ok(())
    }
    }
}
