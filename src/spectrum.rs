//! Define Spectrum struct

/* std use */

/* crate use */

/* local use */

/// Based on Kmergenie we assume kmer spectrum is a mixture of Pareto law and some Gaussians law
/// Erroneous kmer follow Pareto law, Gaussians law represente true and repetitive kmer
/// We use this property to found the threshold to remove most Erroneous kmer and keep Many True kmer
#[derive(Debug, PartialEq, Clone)]
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
    /// Create a new Spectrum with count in `counter`
    pub fn from_count<T>(counts: &[T]) -> Self
    where
        T: std::convert::Into<usize> + std::marker::Copy,
    {
        let mut data =
            vec![0u64; 2_usize.pow(8 * std::mem::size_of::<T>() as u32)].into_boxed_slice();

        for count in counts {
            data[<T as std::convert::Into<usize>>::into(*count)] =
                data[<T as std::convert::Into<usize>>::into(*count)].saturating_add(1);
        }

        Self { data }
    }

    /// Found threshold matching with method
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

    #[cfg(test)]
    pub(crate) fn get_raw_histogram(&self) -> &[u64] {
        &self.data
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::counter;

    fn generate_counter() -> counter::Counter<u8> {
        let mut counter: counter::Counter<u8> = counter::Counter::<u8>::new(5);

        for i in 0..cocktail::kmer::get_kmer_space_size(5) {
            counter::Counter::<u8>::inc(
                counter.raw_mut(),
                (cocktail::kmer::canonical(i, 5) >> 1) as usize,
            );
        }

        counter::Counter::<u8>::inc(counter.raw_mut(), 0);

        counter
    }

    #[test]
    fn from_counter() {
        let counter = generate_counter();
        let spectrum = Spectrum::from_count(counter.raw());

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
            spectrum.get_threshold(ThresholdMethod::Rarefaction, 0.00001),
            None
        );
    }

    #[test]
    fn percent_at_most() {
        let spectrum = Spectrum {
            data: Box::new(SPECTRUM),
        };

        assert_eq!(
            spectrum.get_threshold(ThresholdMethod::PercentAtMost, 0.3),
            Some(1)
        );
    }

    #[test]
    fn failled_percent_at_most() {
        let tmp = (0..256).map(|_| 1).collect::<Box<[u64]>>();

        let spectrum = Spectrum { data: tmp };

        assert_eq!(
            spectrum.get_threshold(ThresholdMethod::PercentAtMost, 1.2),
            None
        );
    }

    #[test]
    fn percent_at_least() {
        let spectrum = Spectrum {
            data: Box::new(SPECTRUM),
        };

        assert_eq!(
            spectrum.get_threshold(ThresholdMethod::PercentAtLeast, 0.3),
            Some(2)
        );
    }

    #[test]
    fn failled_percent_at_least() {
        let tmp = (0..256).map(|_| 1).collect::<Box<[u64]>>();

        let spectrum = Spectrum { data: tmp };

        assert_eq!(
            spectrum.get_threshold(ThresholdMethod::PercentAtLeast, 1.2),
            None
        );
    }
}
