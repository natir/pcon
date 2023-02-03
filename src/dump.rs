//! Run dump command

/* std use */

/* crate use */

/* project use */
use crate::cli;
use crate::counter;
use crate::error;

cfg_if::cfg_if! {
    if #[cfg(feature = "parallel")] {
    /// Run dump
    pub fn dump(params: cli::Dump) -> error::Result<()> {
    log::info!("Start load count");
    let counter = counter::Counter::<std::sync::atomic::AtomicU8>::from_stream(params.input()?)?;
    log::info!("End load count");

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
    /// Run dump
    pub fn dump(params: cli::Dump) -> error::Result<()> {
    log::info!("Start load count");
    let counter = counter::Counter::<u8>::from_stream(params.input()?)?;
    log::info!("End load count");

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
