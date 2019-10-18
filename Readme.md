# Prompt COuNter, is a short kmer counter

[![Build Status](https://travis-ci.org/natir/pcon.svg?branch=master)](https://travis-ci.org/natir/pcon)

pcon is a fast kmer counter but with some important limitations:

- only fasta file
- k can't be upper than 31 (k is silently reassign to 31)
- if k is even k is reduced to the nearest lower odd number
- ssik allocate 2^(k * 2 - 1) / 2 bytes (for k 19 ssik required 69 go)
- max abundance is 15
- if data contains something other than A C T or G is consider like A C T or G (check the 2nd and the 3rd bit of lettre, N was consider as G for exemple)

## Instalation

### With source

If you haven't a rust environment you can use [rustup](https://rustup.rs/) or your package manager.

When you have a rust environment setup you can run this command:

```
git clone https://github.com/natir/ssik.git
cd ssik

cargo build
cargo test
cargo install
```

## Usage

ssik are contains 5 subcommand:

- count:
  count kmer in fasta input and write result in binary format
- dump:
  convert ssik binary format in csv format, or as bit field of exist or not exist kmer

