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

fn build_counter(k: u8) -> pcon::counter::Counter {
    let mut counter = pcon::counter::Counter::new(k);

    let reader = std::fs::File::open("reads.fasta").unwrap();

    counter.count_fasta(reader, 1024);

    counter
}

fn dump_bin(c: &mut Criterion) {
    let mut g = c.benchmark_group("dump_bin");

    g.sample_size(10);
    g.warm_up_time(std::time::Duration::from_secs(1));

    for k in &[15, 17, 19] {
        let counter = build_counter(*k);

        g.bench_with_input(BenchmarkId::new("all_count", k), &k, |b, k| {
            b.iter(|| {
                let writer = std::io::BufWriter::new(
                    std::fs::File::create(format!("bincode_{}.pcon", k)).unwrap(),
                );

                bincode::serialize_into(writer, &counter).unwrap();
            })
        });

        g.bench_with_input(BenchmarkId::new("index+count_gzip", k), &k, |b, k| {
            b.iter(|| {
                let writer =
                    std::fs::File::create(format!("index+count_gzip_k{}.pcon", k)).unwrap();

                counter.serialize(writer, 0).unwrap();
            })
        });

        g.bench_with_input(BenchmarkId::new("count_rle_gzip", k), &k, |b, k| {
            b.iter(|| {
                let writer = std::fs::File::create(format!("count_rle_gzip_k{}.pcon", k)).unwrap();

                counter.serialize(writer, 0).unwrap();
            })
        });
    }
}

fn setup(c: &mut Criterion) {
    dump_bin(c);
}

criterion_group!(benches, setup);

criterion_main!(benches);
