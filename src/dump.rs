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

    let serialize = counter.serialize();

    for (out_type, output) in params.outputs().into_iter() {
            match out_type {
        cli::DumpType::Pcon => {
            log::info!("Start write count in pcon format");
            serialize.pcon(output?)?;
            log::info!("End write count in pcon format");
        }
        cli::DumpType::Csv =>{
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
    } else {
    /// Run dump
    pub fn dump(params: cli::Dump) -> error::Result<()> {
    log::info!("Start load count");
    let counter = counter::Counter::<u8>::from_stream(params.input()?)?;
    log::info!("End load count");


        let serialize = counter.serialize();

    for (out_type, output) in params.outputs().into_iter() {
            match out_type {
        cli::DumpType::Pcon => {
            log::info!("Start write count in pcon format");
            serialize.pcon(output?)?;
            log::info!("End write count in pcon format");
        }
        cli::DumpType::Csv =>{
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
    }
}
