<h1 style="text-align: center;">pcon</h1>

![Test](https://github.com/natir/pcon/workflows/Test/badge.svg)
![Lints](https://github.com/natir/pcon/workflows/Lints/badge.svg)
![MSRV](https://github.com/natir/pcon/workflows/MSRV/badge.svg)
[![CodeCov](https://codecov.io/gh/natir/pcon/branch/master/graph/badge.svg)](https://codecov.io/gh/natir/pcon)
[![Documentation](https://github.com/natir/pcon/workflows/Documentation/badge.svg)](https://natir.github.io/pcon/pcon)
[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/natir/pcon/blob/master/LICENSE)


Prompt COuNter, a short kmer counter.

- only fasta file
- if k is even k is reduced to the nearest lower odd number
- pcon allocate 2^(k * 2 - 1) times number bytes used by counter value (for k 19 pcon and one bytes counts required 69 go)

If data contains something other than A C T or G is consider like A C T or G (check the 2nd and the 3rd bit of lettre, N was consider as G for exemple)

## Installation

### Features

Pcon use rust feature to control some behavior at build time.

#### Maximum of count

Size of variable used to store kmer count is control by *count_\** features.

1. *count\_u16*: count on 2 bytes max value 65535
2. *count\_u32*: count on 4 bytes max value 4294967295
3. *count\_u64*: count on 8 bytes max value 18446744073709551615
4. *count\_u8*:  count on 1 bytes max value 255

If you set multiple count\_\* feature priority of feature follow previous list order.

#### Parallel

If you set feature parallel you activate pcon parallel feature, based on [rayon](https://docs.rs/rayon/latest/rayon/) and [rust atomic](https://doc.rust-lang.org/core/sync/atomic/index.html) type.

#### Default

*count\_u8* is the only default features.

### From source

```bash
git clone https://github.com/natir/pcon.git
cd pcon
cargo install --path . --features {features1},{features2}
```

### With cargo

```bash
cargo install https://github.com/natir/pcon --features {features1},{features2}
```

## Usage

### Count

By default `pcon count` read input fasta file from stdin and write count in stdout in pcon internal format.

```
-k, --kmer-size <KMER_SIZE>          Size of kmer
-i, --inputs <INPUTS>                Path to inputs, default read stdin
-p, --pcon <PCON>                    Path where count are store, default write in stdout
-c, --csv <CSV>                      Path where count are store
-s, --solid <SOLID>                  Path where count are store
-a, --abundance <ABUNDANCE>          Minimal abundance, default value 0
-b, --record_buffer <RECORD_BUFFER>  Number of sequence record load in buffer, default 8192
```

Count 7-mer in `example.fasta` file and write result in pcon format in `example.pcon` file:
```bash
pcon count -k 7 -i example.fasta -p example.pcon
```

### Dump

By default `pcon dump` read input pcon file from stdin and write count in csv format in stdout.

```
-i, --inputs <INPUT>         Path to inputs, default read stdin
-c, --csv <CSV>              Path where count are store, default write in stdout
-p, --pcon <PCON>            Path where count are store
-s, --solid <SOLID>          Path where count are store
-a, --abundance <ABUNDANCE>  Minimal abundance, default value 0
```

Convert 7-mer count in `example.pcon` in csv file `example.csv`:
```bash
pcon dump -i example.pcon -c example.csv
```

### Not subcommand parameter

```
-q, --quiet           Silence all output
-v, --verbosity...    Verbose mode (-v, -vv, -vvv, etc)
-T, --timestamp <TS>  Timestamp (sec, ms, ns, none)
-h, --help            Print help
-V, --version         Print version
```

## Minimum supported Rust version

Currently the minimum supported Rust version is 1.65.
