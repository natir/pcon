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

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

use pcon;

fn build_counter() -> pcon::counter::Counter {
    let mut counter = pcon::counter::Counter::new(15);

    let reader = std::fs::File::open("reads.fasta").unwrap();

    counter.count_fasta(reader, 1024);

    counter
}

fn dump_bin(c: &mut Criterion) {
    let counter = build_counter();

    let mut g = c.benchmark_group("dump_bin");

    g.sample_size(10);
    g.warm_up_time(std::time::Duration::from_secs(1));

    g.bench_function("raw", |b| {
        b.iter(|| {
            let writer = std::fs::File::create("raw.pcon").unwrap();

            counter.serialize(writer, 0).unwrap();
        })
    });

    g.bench_function("bincode", |b| {
        b.iter(|| {
            let writer = std::fs::File::create("bincode.pcon").unwrap();

            bincode::serialize_into(writer, counter).unwrap();
        })
    });

    for power in 10..17 {
        let buffer_len = 1 << power;

        g.bench_with_input(
            BenchmarkId::new("buffer", buffer_len),
            &buffer_len,
            |b, buffer_len| {
                b.iter(|| {
                    let writer = std::io::BufWriter::with_capacity(
                        *buffer_len,
                        std::fs::File::create(format!("buffer_{}.pcon", buffer_len)).unwrap(),
                    );

                    counter.serialize(writer, 0).unwrap();
                })
            },
        );
    }
}

fn dump_csv(c: &mut Criterion) {
    let counter = build_counter();

    let mut g = c.benchmark_group("dump_csv");

    g.sample_size(10);
    g.warm_up_time(std::time::Duration::from_secs(1));

    g.bench_function("raw", |b| {
        b.iter(|| {
            let writer = std::fs::File::create("raw.csv").unwrap();

            counter.serialize(writer, 0).unwrap();
        })
    });

    for power in 10..17 {
        let buffer_len = 1 << power;

        g.bench_with_input(
            BenchmarkId::new("buffer", buffer_len),
            &buffer_len,
            |b, buffer_len| {
                b.iter(|| {
                    let writer = std::io::BufWriter::with_capacity(
                        *buffer_len,
                        std::fs::File::create(format!("buffer_{}.csv", buffer_len)).unwrap(),
                    );

                    pcon::dump::csv(writer, &counter, 0).unwrap();
                })
            },
        );
    }
}

fn setup(c: &mut Criterion) {
    dump_bin(c);
    dump_csv(c);
}

criterion_group!(benches, setup);

criterion_main!(benches);
