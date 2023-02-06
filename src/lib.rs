//! Prompt COuNter, a short kmer counter.

#![warn(missing_docs)]

/* std use */

/* crate use */

/* project use */

/* mod declaration */
pub mod cli;
pub mod count;
pub mod counter;
pub mod dump;
pub mod error;
pub mod serialize;
pub mod solid;

#[cfg(not(tarpaulin_include))]
/// Define a const
type ByteOrder = byteorder::LittleEndian; // WARNING IF YOU CHANGE THIS CHECK AND CHANGE SERIALIZE.RS
