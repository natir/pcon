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

## Instalation

### With source

If you haven't a rust environment you can use [rustup](https://rustup.rs/) or your package manager.

When you have a rust environment setup you can run this command:

```
git clone https://github.com/natir/pcon.git
cd pcon

cargo build
cargo test
cargo install
```

## Usage

pcon contains 2 subcommand:

- count: count kmer in fasta input and write result in binary format
- dump: convert pcon binary format in csv format, presence absence bitfield or kmer spectrum

## Binding

### Python

You need [rust toolchain setup on your system](https://rustup.rs/)

Give this to pip:
```
git+https://github.com/natir/pcon.git#egg=cocktail&subdirectory=dist/python
```

### C

You need [rust toolchain setup on your system](https://rustup.rs/)

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

## Minimum supported Rust version

Currently the minimum supported Rust version is 1.45.0.
