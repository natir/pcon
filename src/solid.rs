//! Define Solid struct

/* std use */

/* crate use */
use bitvec::prelude::*;
use byteorder::ReadBytesExt as _;

/* local use */
use crate::error;

/// A struct to store if a kmer is Solid or not. Only kmer with abundance upper than a threshold is solid
pub struct Solid {
    k: u8,
    solid: BitBox<u8, Lsb0>,
}

impl Solid {
    /// Create a new Solid for kmer size equal to `k`
    pub fn new(k: u8) -> Self {
        Self {
            k,
            solid: bitbox![u8, Lsb0; 0; cocktail::kmer::get_hash_space_size(k) as usize],
        }
    }

    /// Create a new Solid with count in `counter` only kmer upper than `abundance` are solid
    pub fn from_count<T>(k: u8, count: &[T], abundance: T) -> Self
    where
        T: std::cmp::PartialOrd,
    {
        let mut solid = bitbox![u8, Lsb0; 0; count.len()];

        for (index, count) in count.iter().enumerate() {
            if *count > abundance {
                solid.set(index, true);
            }
        }

        Self { k, solid }
    }

    /// Create a new Solid by read
    pub fn from_stream<R>(mut input: R) -> error::Result<Self>
    where
        R: std::io::Read,
    {
        let k = input.read_u8()?;

        let mut solid = bitbox![u8, Lsb0; 0; cocktail::kmer::get_hash_space_size(k) as usize];

        input.read_exact(solid.as_raw_mut_slice())?;

        Ok(Self { k, solid })
    }

    /// Create a new Solid from path
    pub fn from_path<P>(path: P) -> error::Result<Self>
    where
        P: std::convert::AsRef<std::path::Path>,
    {
        let (readable, _compression) = niffler::get_reader(
            std::fs::File::open(path)
                .map(std::io::BufReader::new)
                .map(Box::new)?,
        )?;

        Self::from_stream(readable)
    }

    /// Get value of k
    pub fn k(&self) -> u8 {
        self.k
    }

    /// Solidity status of `kmer` is set to `value`
    pub fn set(&mut self, kmer: u64, value: bool) {
        self.set_canonic(cocktail::kmer::canonical(kmer, self.k), value);
    }

    /// Solidity status of a canonical`kmer` is set to `value`
    pub fn set_canonic(&mut self, canonical: u64, value: bool) {
        let hash = (canonical >> 1) as usize;

        if let Some(mut v) = self.solid.get_mut(hash) {
            *v = value;
        }
    }

    /// Get the solidity status of `kmer`
    pub fn get(&self, kmer: u64) -> bool {
        self.get_canonic(cocktail::kmer::canonical(kmer, self.k))
    }

    /// Get the solidity status of a canonical `kmer`
    pub fn get_canonic(&self, canonical: u64) -> bool {
        let hash = (canonical >> 1) as usize;

        self.solid[hash]
    }

    /// Extend
    pub fn extend(&mut self, rhs: Solid) {
        self.solid |= rhs.get_raw_solid()
    }

    pub(crate) fn get_raw_solid(&self) -> &BitBox<u8, Lsb0> {
        &self.solid
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const FASTA_FILE: &[u8] = b">random_seq 0
GTTCTGCAAATTAGAACAGACAATACACTGGCAGGCGTTGCGTTGGGGGAGATCTTCCGTAACGAGCCGGCATTTGTAAGAAAGAGATTTCGAGTAAATG
>random_seq 1
AGGATAGAAGCTTAAGTACAAGATAATTCCCATAGAGGAAGGGTGGTATTACAGTGCCGCCTGTTGAAAGCCCCAATCCCGCTTCAATTGTTGAGCTCAG
";

    fn get_counter() -> crate::counter::Counter<u8> {
        let mut counter = crate::counter::Counter::<u8>::new(5);

        counter.count_fasta(Box::new(FASTA_FILE), 1);

        counter
    }

    fn get_solid() -> Solid {
        let counter = get_counter();

        Solid::from_count(counter.k(), counter.raw(), 0)
    }

    const SOLID: &[u8] = &[
        112, 64, 113, 143, 130, 8, 128, 4, 6, 60, 214, 0, 243, 8, 193, 1, 30, 4, 34, 97, 4, 70,
        192, 12, 16, 144, 133, 38, 192, 41, 1, 4, 218, 179, 140, 0, 0, 140, 242, 35, 90, 56, 205,
        179, 64, 3, 25, 20, 226, 0, 32, 76, 1, 134, 48, 64, 7, 0, 200, 144, 98, 131, 2, 203,
    ];

    #[test]
    fn new_solid() {
        let mut solid = Solid::new(5);

        assert_eq!(solid.get(43), false);

        solid.set(42, true);

        assert_eq!(solid.get(42), true);
    }

    #[test]
    fn presence() {
        let solid = get_solid();

        assert_eq!(solid.get(1), false);
        assert_eq!(solid.get(724), true);

        assert_eq!(solid.get_raw_solid().as_raw_slice(), SOLID);
    }

    const SOLID_SET: &[u8] = &[
        112, 64, 113, 143, 130, 8, 128, 4, 6, 52, 214, 0, 243, 8, 193, 1, 30, 4, 2, 97, 4, 70, 192,
        12, 16, 144, 133, 36, 192, 41, 1, 4, 218, 179, 140, 0, 0, 140, 242, 35, 90, 56, 205, 179,
        64, 3, 25, 20, 226, 0, 32, 76, 1, 134, 48, 64, 7, 0, 192, 144, 98, 131, 2, 203,
    ];

    #[test]
    fn set_value() {
        let mut solid = get_solid();

        solid.set(cocktail::kmer::seq2bit(b"GTTCT"), false);
        solid.set(cocktail::kmer::seq2bit(b"AAATG"), false);
        solid.set(cocktail::kmer::seq2bit(b"AGGAT"), false);
        solid.set(cocktail::kmer::seq2bit(b"CTCAG"), false);

        assert_eq!(solid.get_raw_solid().as_raw_slice(), SOLID_SET);
    }

    #[test]
    fn extend() {
        let mut solid = get_solid();
        let mut other = Solid::new(5);

        assert_eq!(solid.get(42), true);
        assert_eq!(solid.get(44), false);
        assert_eq!(other.get(44), false);

        other.set(44, true);
        assert_eq!(other.get(44), true);

        assert_eq!(solid.get_raw_solid().count_ones(), 158);
        assert_eq!(solid.get_raw_solid().count_zeros(), 354);
        solid.extend(other);
        assert_eq!(solid.get_raw_solid().count_ones(), 159);
        assert_eq!(solid.get_raw_solid().count_zeros(), 353);
        assert_eq!(solid.get(44), true);
    }

    #[test]
    fn deserilize() -> error::Result<()> {
        let counter = get_counter();
        assert_eq!(counter.k(), 5);
        let serializer = counter.serialize();

        let temp = tempfile::NamedTempFile::new()?;
        serializer.solid(0, &temp)?;

        let path = temp.into_temp_path();

        let solid = Solid::from_path(path)?;

        assert_eq!(solid.k(), 5);
        assert_eq!(solid.get_raw_solid().as_raw_slice(), SOLID);

        Ok(())
    }
}
