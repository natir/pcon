//! Error struct of project pcon

/* crate use */
use anyhow;
use thiserror;

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
}

/// Alias of result
pub type Result<T> = anyhow::Result<T>;
