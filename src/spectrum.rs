/*
Copyright (c) 2020 Pierre Marijon <pierre.marijon@hhu.de>

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
use anyhow::{Context, Result};
use rayon::prelude::*;

/* local use */
use crate::error::IO::*;
use crate::error::*;
use crate::*;

/// Based on Kmergenie we assume kmer spectrum is a mixture of Pareto law and some Gaussians law
/// Erroneous kmer follow Pareto law, Gaussians law represente true and repetitive kmer
/// We use this property to found the threshold to remove many Erroneous kmer and keep Many True kmer
pub enum ThresholdMethod {
    /// The first local minimum match with the intersection of Pareto and Gaussians
    FirstMinimum,

    /// More we remove kmer less we remove Erroneous kmer when remove less than n percent view before    
    Rarefaction,

    /// Remove at most n percent of total kmer
    PercentAtMost,

    /// Remove at least n percent of total kmer
    PercentAtLeast,
}

/// A struct to represent kmer spectrum and usefull corresponding function
pub struct Spectrum {
    data: Box<[u64]>,
}

impl Spectrum {
    pub fn from_counter(counter: &counter::Counter) -> Self {
        let counts = unsafe {
            &(*(counter.get_raw_count() as *const [counter::AtoCount] as *const [counter::Count]))
        };

        let data = counts
            .par_chunks(counts.len() / rayon::current_num_threads())
            .map(|chunk| {
                let mut d: Box<[u64]> = vec![0u64; 256].into_boxed_slice();

                for count in chunk.iter() {
                    d[*count as usize] = d[*count as usize].saturating_add(1);
                }

                d
            })
            .reduce(
                || vec![0u64; 256].into_boxed_slice(),
                |a, b| {
                    let mut d = Vec::with_capacity(256);

                    for x in a.iter().zip(b.iter()) {
                        d.push(x.0.saturating_add(*x.1));
                    }

                    d.into_boxed_slice()
                },
            );

        Self { data }
    }

    pub fn get_threshold(&self, method: ThresholdMethod, params: f64) -> Option<u8> {
        match method {
            ThresholdMethod::FirstMinimum => self.first_minimum(),
            ThresholdMethod::Rarefaction => self.rarefaction(params),
            ThresholdMethod::PercentAtMost => self.percent_at_most(params),
            ThresholdMethod::PercentAtLeast => self.percent_at_least(params),
        }
    }

    fn first_minimum(&self) -> Option<u8> {
        for (i, d) in self.data.windows(2).enumerate() {
            if d[1] > d[0] {
                return Some(i as u8);
            }
        }

        None
    }

    fn rarefaction(&self, limit: f64) -> Option<u8> {
        let mut cumulative_sum = 0;

        for (index, value) in self.data.iter().enumerate() {
            cumulative_sum += index as u64 * value;

            if (*value as f64 / cumulative_sum as f64) < limit {
                return Some(index as u8);
            }
        }

        None
    }

    fn percent_at_most(&self, percent: f64) -> Option<u8> {
        self.percent_at_least(percent).map(|x| x - 1)
    }

    fn percent_at_least(&self, percent: f64) -> Option<u8> {
        let total: u64 = self
            .data
            .iter()
            .enumerate()
            .map(|(index, value)| index as u64 * value)
            .sum();

        let mut cumulative_sum = 0;
        for (index, value) in self.data.iter().enumerate() {
            cumulative_sum += index as u64 * value;

            if (cumulative_sum as f64 / total as f64) > percent {
                return Some(index as u8);
            }
        }

        None
    }

    #[allow(dead_code)]
    pub(crate) fn get_raw_histogram(&self) -> &[u64] {
        &self.data
    }

    pub fn write_csv<W>(&self, mut writer: W) -> Result<()>
    where
        W: std::io::Write,
    {
        for (i, nb) in self.data.iter().enumerate() {
            writeln!(writer, "{},{}", i, nb).with_context(|| Error::IO(ErrorDurringWrite))?;
        }

        Ok(())
    }

    pub fn write_histogram<W>(&self, mut out: W, point: Option<u8>) -> std::io::Result<()>
    where
        W: std::io::Write,
    {
        // Draw kmer spectrum in console
        let shape = get_shape();

        let factor = (*self.data.iter().max().unwrap() as f64).log(10.0) / shape.1 as f64;

        let normalized: Box<[f64]> = self
            .data
            .iter()
            .map(|y| {
                if *y == 0 {
                    0.0
                } else {
                    (*y as f64).log(10.0) / factor
                }
            })
            .collect();

        for h in (1..=shape.1).rev() {
            for w in 0..shape.0 {
                if normalized[w] >= h as f64 {
                    let delta = normalized[w] - h as f64;
                    if delta > 1.0 {
                        write!(out, "\u{258c}")?;
                    } else if delta > 0.5 {
                        write!(out, "\u{2596}")?;
                    } else {
                        write!(out, " ")?;
                    }
                } else {
                    write!(out, " ")?;
                }
            }

            writeln!(out)?;
        }

        let mut last_line = vec![b' '; shape.0];
        for x in (0..shape.0).step_by(5) {
            last_line[x] = b'5'
        }
        last_line[0] = b'0';
        if let Some(pos) = point {
            if (pos as usize) < last_line.len() {
                last_line[pos as usize] = b'*';
            }
        }

        writeln!(out, "{}", std::str::from_utf8(&last_line).unwrap())?;

        Ok(())
    }
}

#[allow(dead_code)]
fn term_size() -> (usize, usize) {
    match term_size::dimensions() {
        Some((w, h)) => (w, h),
        None => (80, 24),
    }
}

#[cfg(test)]
fn get_shape() -> (usize, usize) {
    (256, 48)
}

#[cfg(not(test))]
fn get_shape() -> (usize, usize) {
    let term_size = term_size();

    (
        core::cmp::min(256, term_size.0),
        core::cmp::min(48, term_size.1),
    )
}

#[cfg(test)]
mod tests {

    use super::*;

    lazy_static::lazy_static! {
    static ref COUNTER: crate::counter::Counter = {
            let mut counter = crate::counter::Counter::new(5);

        for i in 0..cocktail::kmer::get_kmer_space_size(5) {
        counter.inc(i);
            }

            counter.inc(0);

            counter
    };
    }

    #[test]
    fn from_counter() {
        let spectrum = Spectrum::from_counter(&COUNTER);

        assert_eq!(
            spectrum.get_raw_histogram(),
            &[
                0, 0, 511, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0
            ]
        );
    }

    static SPECTRUM: [u64; 256] = [
        992273316, 64106898, 6792586, 1065818, 220444, 62400, 36748, 54062, 100806, 178868, 287058,
        424184, 568742, 705680, 805332, 871544, 874546, 827252, 744428, 636722, 523488, 418036,
        320506, 237956, 170642, 118046, 77290, 48320, 30500, 21096, 15632, 12758, 11838, 10888,
        10402, 9872, 9018, 7960, 7236, 6304, 5276, 4524, 3714, 3056, 2628, 2018, 1578, 1256, 1036,
        906, 708, 716, 592, 476, 540, 520, 446, 388, 316, 264, 258, 200, 230, 172, 164, 184, 154,
        162, 126, 124, 126, 156, 152, 98, 116, 108, 134, 116, 88, 124, 96, 94, 96, 72, 52, 56, 68,
        50, 54, 66, 54, 28, 44, 48, 30, 42, 48, 32, 38, 34, 44, 30, 32, 28, 18, 34, 20, 28, 26, 28,
        28, 32, 22, 16, 10, 26, 8, 26, 14, 14, 30, 6, 32, 38, 26, 26, 16, 30, 20, 38, 20, 22, 22,
        28, 14, 16, 20, 20, 20, 10, 12, 14, 12, 10, 18, 16, 16, 12, 18, 2, 14, 6, 12, 8, 0, 6, 2,
        4, 2, 0, 0, 2, 4, 2, 2, 6, 6, 0, 0, 2, 0, 2, 4, 0, 2, 2, 6, 2, 0, 0, 0, 2, 2, 2, 2, 2, 0,
        2, 2, 0, 2, 2, 0, 2, 2, 2, 0, 0, 2, 4, 2, 0, 2, 0, 2, 2, 2, 0, 2, 2, 2, 2, 2, 0, 0, 0, 2,
        0, 0, 2, 2, 2, 2, 4, 0, 2, 4, 4, 0, 2, 0, 0, 2, 2, 0, 0, 0, 0, 0, 4, 2, 0, 2, 0, 0, 0, 2,
        0, 4, 2, 0, 4, 2, 0, 0, 284,
    ];

    #[test]
    fn first_local_min() {
        let spectrum = Spectrum {
            data: Box::new(SPECTRUM),
        };

        assert_eq!(
            spectrum.get_threshold(ThresholdMethod::FirstMinimum, 0.1),
            Some(6)
        );
    }

    #[test]
    fn failled_first_local_min() {
        let tmp = (0..256).map(|_| 1).collect::<Box<[u64]>>();

        let spectrum = Spectrum { data: tmp };

        assert_eq!(
            spectrum.get_threshold(ThresholdMethod::FirstMinimum, 0.1),
            None
        );
    }

    #[test]
    fn rarefaction() {
        let spectrum = Spectrum {
            data: Box::new(SPECTRUM),
        };

        assert_eq!(
            spectrum.get_threshold(ThresholdMethod::Rarefaction, 0.1),
            Some(2)
        );
    }

    #[test]
    fn failled_rarefaction() {
        let tmp = (0..256).map(|_| 1).collect::<Box<[u64]>>();

        let spectrum = Spectrum { data: tmp };

        assert_eq!(
            spectrum.get_threshold(ThresholdMethod::FirstMinimum, 0.00001),
            None
        );
    }

    #[test]
    fn test_draw_hist() {
        let spectrum = Spectrum {
            data: Box::new(SPECTRUM),
        };

        let mut output = vec![];
        spectrum.write_histogram(&mut output, Some(6)).unwrap();

        let good_output = "                                                                                                                                                                                                                                                                
▖                                                                                                                                                                                                                                                               
▌                                                                                                                                                                                                                                                               
▌                                                                                                                                                                                                                                                               
▌                                                                                                                                                                                                                                                               
▌                                                                                                                                                                                                                                                               
▌                                                                                                                                                                                                                                                               
▌▖                                                                                                                                                                                                                                                              
▌▌                                                                                                                                                                                                                                                              
▌▌                                                                                                                                                                                                                                                              
▌▌                                                                                                                                                                                                                                                              
▌▌                                                                                                                                                                                                                                                              
▌▌                                                                                                                                                                                                                                                              
▌▌▌                                                                                                                                                                                                                                                             
▌▌▌                                                                                                                                                                                                                                                             
▌▌▌                                                                                                                                                                                                                                                             
▌▌▌                                                                                                                                                                                                                                                             
▌▌▌▌          ▖▖▖▖                                                                                                                                                                                                                                              
▌▌▌▌        ▖▌▌▌▌▌▌▖▖                                                                                                                                                                                                                                           
▌▌▌▌       ▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                                                          
▌▌▌▌▖     ▌▌▌▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                                                        
▌▌▌▌▌    ▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                                                       
▌▌▌▌▌   ▖▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌                                                                                                                                                                                                                                      
▌▌▌▌▌▖  ▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌                                                                                                                                                                                                                                     
▌▌▌▌▌▌ ▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                                                    
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                                                   
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌                                                                                                                                                                                                                                  
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖▖                                                                                                                                                                                                                              
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖                                                                                                                                                                                                                         
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖                                                                                                                                                                                                                      
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                                    
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                                  
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                                
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                              
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖ ▖                                                                                                                                                                                                         
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                      
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖ ▖                                                                                                                                                                                                ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖▌▖▖   ▖▖                                                                                                                                                                                      ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▌▖▌▌ ▌▖▖▖                                                                                                                                                                            ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖  ▖  ▖                                                                                                                                                                     ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌ ▖▖ ▖▖   ▖                                                                                                                                                          ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▌▌▖▌▌▌▌▌▌▖▌▖ ▌ ▖▖▖▖▌   ▖ ▖  ▖ ▌▌▖▖ ▖ ▌   ▖                                                                                                                         ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▌▖▌▌▌▌▌▌  ▌ ▌  ▌ ▌▌▌▌ ▌▖▌▖▌▌▌  ▖▖▖     ▖   ▖                                                                                                          ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌ ▌ ▌▌▌▌ ▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌ ▖▌▖ ▌▌▌▖▌ ▌ ▖                                                                                                      ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▌▌▌▌ ▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌ ▌ ▌▖                                                                                                     ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌ ▌▌▌▌ ▌         ▌▌         ▌                                                                              ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌ ▌▌▌▌ ▌ ▌    ▌  ▌▌     ▌   ▌                      ▌                       ▌  ▌▌           ▌        ▌  ▌   ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▌▌▌▌ ▌▖▌▖  ▖▌▖▖▌▌  ▖ ▖▌ ▖▖▌▖   ▖▖▖▖▖ ▖▖ ▖▖ ▖▖▖  ▖▌▖ ▖ ▖▖▖ ▖▖▖▖▖   ▖  ▖▖▖▖▌ ▖▌▌ ▖  ▖▖     ▌▖ ▖   ▖ ▌▖ ▌▖  ▌
0    5*   5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5
";

        assert_eq!(good_output, std::str::from_utf8(&output).unwrap());
    }
}
