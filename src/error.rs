/*
Copyright (c) 2020 Pierre Marijon <pmarijon@mpi-inf.mpg.de>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

/* crate use */
use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error(transparent)]
    Cli(#[from] Cli),

    #[error(transparent)]
    IO(#[from] IO),

    #[error("If you get this error please contact the author with this message and command line you use: {name:?}")]
    NotReachableCode { name: String },
}

#[derive(Debug, Error)]
pub enum Cli {
    #[error("Kmer size must be odd")]
    KMustBeOdd,

    #[error("Minimizer size must be lower than kmer size")]
    MinimizerLowerThanKmer,
}

#[derive(Debug, Error)]
pub enum IO {
    #[error("We can't create file at path {path:}")]
    CantCreateFile { path: String },

    #[error("We can't open file at path {path:}")]
    CantOpenFile { path: String },

    #[error("Error durring write: {context:}")]
    ErrorDurringWrite { context: String },

    #[error("Error durring write in {path:}")]
    ErrorDurringWriteIn { path: String },

    #[error("Error durring read: {context:}")]
    ErrorDurringRead { context: String },

    #[error("Error durring read of {path:}")]
    ErrorDurringReadIn { path: String },
}
