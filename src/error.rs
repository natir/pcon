//! Error struct of project pcon

/* std use */

/* crate use */
use anyhow;
use thiserror;

/* project use */

/// Enum to manage error
#[derive(std::fmt::Debug, thiserror::Error)]
pub enum Error {
    /// Error in logging system configuration
    #[error(transparent)]
    Log(#[from] log::SetLoggerError),

    /// Error in rayon thread pool build
    #[cfg(feature = "parallel")]
    #[error(transparent)]
    RayonThreadPool(#[from] rayon::ThreadPoolBuildError),

    /// Cost io error
    #[error(transparent)]
    IO(#[from] std::io::Error),

    /// Error if we can't convert a DumpTypeFromStr
    #[error("Can't convert {0} in DumpType")]
    DumpTypeFromStr(String),

    /// Error durring loading count type not match
    #[error("Type use in counter not match file count")]
    TypeNotMatch,
}

/// Alias of result
pub type Result<T> = anyhow::Result<T>;
