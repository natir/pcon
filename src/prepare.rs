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

pub fn prepare(k: u8, name: &str) {
    print!("For {} = {} ssik request", name, k);

    let base: u128 = 1;
    let l: u128 = (k as u128) * 2 - 1;
    println!("\t{}\t{} b", convert(base << l), (base << l));
}

fn convert(mut value: u128) -> String {
    let mut prefix = 0;
    while value >= 1024 && prefix < 8 {
        value /= 1024;
        prefix += 1;
    }

    let prefix_table = ["b", "kb", "mb", "gb", "tb", "pb", "eb", "zb", "yb"];
    return format!("{} {}", value, prefix_table[prefix]);
}
