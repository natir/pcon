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
pub mod spectrum;

/// Define a const
type ByteOrder = byteorder::LittleEndian; // WARNING IF YOU CHANGE THIS CHECK AND CHANGE SERIALIZE.RS

cfg_if::cfg_if! {
    if #[cfg(all(feature = "count_u16", feature = "parallel"))] {
    /// Define count type
    pub type CountType = std::sync::atomic::AtomicU16;
    /// Define count type for all never atomic thing
    pub type CountTypeNoAtomic = u16;
    } else if #[cfg(all(feature = "count_u16", not(feature = "parallel")))] {
    /// Define count type
    pub type CountType = u16;
    /// Define count type for all never atomic thing
    pub type CountTypeNoAtomic = u16;
    } else if #[cfg(all(feature = "count_u32", feature = "parallel"))] {
    /// Define count type
    pub type CountType = std::sync::atomic::AtomicU32;
    /// Define count type for all never atomic thing
    pub type CountTypeNoAtomic = u32;
    } else if #[cfg(all(feature = "count_u32", not(feature = "parallel")))] {
    /// Define count type
    pub type CountType = u32;
    /// Define count type for all never atomic thing
    pub type CountTypeNoAtomic = u32;
    } else if #[cfg(all(feature = "count_u64", feature = "parallel"))] {
    /// Define count type
    pub type CountType = std::sync::atomic::AtomicU64;
    /// Define count type for all never atomic thing
    pub type CountTypeNoAtomic = u64;
    } else if #[cfg(all(feature = "count_u64", not(feature = "parallel")))] {
    /// Define count type
    pub type CountType = u64;
    /// Define count type for all never atomic thing
    pub type CountTypeNoAtomic = u64;
    } else if #[cfg(feature = "parallel")] {
    /// Define count type
    pub type CountType = std::sync::atomic::AtomicU8;
    /// Define count type for all never atomic thing
    pub type CountTypeNoAtomic = u8;
    } else {
    /// Define count type
    pub type CountType = u8;
    /// Define count type for all never atomic thing
    pub type CountTypeNoAtomic = u8;
    }
}
