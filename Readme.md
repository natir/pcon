# Prompt COuNter, is a short kmer counter

[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/natir/pcon/blob/master/LICENSE)
![CI](https://github.com/natir/pcon/workflows/CI/badge.svg)
[![Documentation](https://github.com/natir/pcon/workflows/Documentation/badge.svg)](https://natir.github.io/pcon/pcon)
[![CodeCov](https://codecov.io/gh/natir/pcon/branch/master/graph/badge.svg)](https://codecov.io/gh/natir/pcon)

pcon is a fast kmer counter but with some important limitations:

- only fasta file
- k can't be upper than 31 (k is silently reassign to 31)
- if k is even k is reduced to the nearest lower odd number
- pcon allocate 2^(k * 2 - 1) / 2 bytes (for k 19 pcon required 69 go)
- max abundance is 15
- if data contains something other than A C T or G is consider like A C T or G (check the 2nd and the 3rd bit of lettre, N was consider as G for exemple)

- [Instalation](#instalation)
- [Usage](#usage)
- [Minimum supported Rust version](#minimum-supported-rust-version)
- [Citation](#citation)


## Instalation

If you haven't a rust environment you can use [rustup](https://rustup.rs/) or your package manager.

### With cargo

Recommended solution.

```
cargo install --git https://github.com/natir/br.git
```

### With source

```
git clone https://github.com/natir/pcon.git
cd pcon
cargo install --path .
```

### Python binding

Give this to pip:
```
git+https://github.com/natir/pcon.git#subdirectory=dist/python
```

### C binding

```
git clone https://github.com/natir/pcon.git
cd cocktail
cargo build --release

cbindgen --config cbindgen.toml --crate pcon --output dist/c/pcon.h
cd dist/c/
make
./test
```

Dynamic and static library is avaible her `target/release/libcocktail.{a|so}` header is her `dist/c/cocktail.h`. To build a C programe you need to add `-lpthread -lm -ldl` durring linking phase.

## Usage

pcon contains 2 subcommand:

- count: count kmer in fasta input and write result in binary format
- dump: convert pcon binary format in csv format, presence absence bitfield or kmer spectrum

## Minimum supported Rust version

Currently the minimum supported Rust version is 1.45.0.

## Citation

WIP
