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
    }

/// Alias of result
pub type Result<T> = anyhow::Result<T>;
