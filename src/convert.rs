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

pub fn seq2bit(subseq: &[u8]) -> u64 {
    let mut kmer: u64 = 0;

    for n in subseq {
        kmer <<= 2;
        kmer |= nuc2bit(*n);
    }

    return kmer;
}

fn nuc2bit(nuc: u8) -> u64 {
    return (nuc as u64 >> 1) & 0b11;
}

pub fn bit2seq(mut kmer: u64, k: u8) -> String {
    let mut result = vec![0; k as usize];

    for i in (0..k).rev() {
        let val = kmer & 0b11;

        if val == 0b00 {
            result[i as usize] = b'A';
        } else if val == 0b01 {
            result[i as usize] = b'C';
        } else if val == 0b10 {
            result[i as usize] = b'T';
        } else {
            result[i as usize] = b'G';
        }

        kmer >>= 2;
    }

    return String::from_utf8(result).unwrap();
}

pub fn cannonical(kmer: u64, k: u8) -> u64 {
    if parity_even(kmer) {
        return kmer;
    } else {
        return revcomp(kmer, k);
    }
}

pub fn parity_even(kmer: u64) -> bool {
    return kmer.count_ones() % 2 == 0;
}

pub fn revcomp(kmer: u64, k: u8) -> u64 {
    return rev(comp(kmer), k);
}

pub fn comp(kmer: u64) -> u64 {
    return kmer ^ 0b1010101010101010101010101010101010101010101010101010101010101010;
}

pub fn rev(kmer: u64, k: u8) -> u64 {
    let clean_move = 64 - k * 2;

    let mut reverse = reverse_2(kmer, k);
    reverse <<= clean_move;
    reverse >>= clean_move;

    return reverse;
}

pub fn reverse_2(mut kmer: u64, k: u8) -> u64 {
    let mut reversed: u64 = 0;

    for _ in 0..(k - 1) {
        reversed = (reversed ^ (kmer & 0b11)) << 2;
        kmer >>= 2;
    }

    return reversed ^ (kmer & 0b11);
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn seq2bit_() {
        // TAGGC -> 1000111101
        assert_eq!(seq2bit(b"TAGGC"), 0b1000111101);

        // GCCTA -> 110101100
        assert_eq!(seq2bit(b"GCCTA"), 0b1101011000);
    }

    #[test]
    fn bit2seq_() {
        // 1000111101 -> TAGGC
        assert_eq!(bit2seq(0b1000111101, 5), "TAGGC");

        // 110101100 -> GCCTA
        assert_eq!(bit2seq(0b1101011000, 5), "GCCTA");
    }

    #[test]
    fn cannonical_() {
        // TAGGC -> 1000111101 cannonical TAGGC -> 1000111101
        assert_eq!(cannonical(0b1000111101, 5), 0b1000111101);

        // GCCTA -> 1101011000 cannonical TAGGC -> 1000111101
        assert_eq!(cannonical(0b1101011000, 5), 0b1000111101);
    }

    #[test]
    fn parity_even_() {
        assert_eq!(parity_even(0b1111), true);
        assert_eq!(parity_even(0b1110), false);
    }

    #[test]
    fn revcomp_() {
        // TAGGC -> 1000111101 revcomp GCCTA -> 1101011000
        assert_eq!(0b1000111101, revcomp(0b1101011000, 5))
    }

    #[test]
    fn comp_() {
        // TAGGC -> 1000111101 comp 0001001011
        assert_eq!(
            comp(0b1000111101),
            0b1010101010101010101010101010101010101010101010101010100010010111
        );
    }

    #[test]
    fn rev_() {
        // TAGGC -> 1000111101 rev CGGAT -> 0111110010
        let var = 0b1010101010101010101010101010101010101010101010101010101000111101;

        assert_eq!(498, rev(var, 5));
    }

    #[test]
    fn reverse_2_() {
        // TAGGC -> 1000111101 rev CGGAT -> 0111110010
        assert_eq!(498, reverse_2(573, 5));
    }

}
