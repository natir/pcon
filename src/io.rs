/*
Copyright (c) 2019 Pierre Marijon <pierre.marijon@inria.fr>

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

#[derive(Clone, PartialEq)]
pub enum Mode {
    AllCounts,
    Counts,
    KmerCounts,
    Numpy,
}

impl From<&str> for Mode {
    fn from(mode: &str) -> Self {
        match mode {
            "all_counts" => return Mode::AllCounts,
            "counts" => return Mode::Counts,
            "kmer_counts" => return Mode::KmerCounts,
            "numpy" => return Mode::Numpy,
            _ => return Mode::AllCounts,
        }
    }
}

impl From<u8> for Mode {
    fn from(mode: u8) -> Self {
        match mode {
            0 => Mode::AllCounts,
            1 => Mode::Counts,
            2 => Mode::KmerCounts,
            3 => Mode::Numpy,
            _ => Mode::AllCounts,
        }
    }
}

impl Into<u8> for Mode {
    fn into(self) -> u8 {
        match self {
            Mode::AllCounts => 0,
            Mode::Counts => 1,
            Mode::KmerCounts => 2,
            Mode::Numpy => 3,
        }
    }
}

impl Into<u8> for &Mode {
    fn into(self) -> u8 {
        match self {
            Mode::AllCounts => 0,
            Mode::Counts => 1,
            Mode::KmerCounts => 2,
            Mode::Numpy => 3,
        }
    }
}
