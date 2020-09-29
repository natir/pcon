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

pub struct Encoder<'a> {
    input: &'a [u8],
    pos_input: usize,
    buffer_len: usize,
    last_byte: u8,
    count: u8,
}

impl<'a> Encoder<'a> {
    pub fn new(input: &'a [u8], mut buffer_len: usize) -> Self {
        if buffer_len % 2 == 1 {
            buffer_len += 1;
        }

        Self {
            input,
            pos_input: 0,
            buffer_len,
            last_byte: 0,
            count: 0,
        }
    }

    fn end_run(&mut self, local_buffer: &mut Vec<u8>, new_byte: u8) {
        if self.count != 0 {
            local_buffer.push(self.count);
            local_buffer.push(self.last_byte);

            self.last_byte = new_byte;
            self.count = 1;
        }
    }
}

impl<'a> Iterator for Encoder<'a> {
    type Item = Box<[u8]>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut local_buffer = Vec::with_capacity(self.buffer_len);

        if self.pos_input >= self.input.len() {
            return None;
        }

        while local_buffer.len() < self.buffer_len {
            println!("{} {}", self.pos_input, self.input.len());
            let byte = self.input[self.pos_input];

            if byte != self.last_byte {
                if self.count != 0 {
                    self.end_run(&mut local_buffer, byte);
                } else {
                    self.last_byte = byte;
                    self.count = 1;
                }
            } else if self.count == std::u8::MAX {
                self.end_run(&mut local_buffer, byte);
            } else if self.count != 0 {
                self.count += 1;
            } else {
                self.last_byte = byte;
                self.count = 1;
            }

            self.pos_input += 1;
            if self.pos_input >= self.input.len() {
                println!("end {} {}", self.pos_input, self.input.len());
                local_buffer.push(self.count);
                local_buffer.push(self.last_byte);
                break;
            }
        }

        if local_buffer.len() == 0 {
            None
        } else {
            Some(local_buffer.into_boxed_slice())
        }
    }
}

pub struct Decoder<'a> {
    input: &'a [u8],
    pos_input: usize,
    max_pos_input: usize,
    last_byte: u8,
    count: u8,
}

impl<'a> Decoder<'a> {
    pub fn new(input: &'a [u8]) -> Self {
        Self {
            input,
            pos_input: 0,
            max_pos_input: input.len() - input.len() % 2,
            last_byte: 0,
            count: 0,
        }
    }
}

impl<'a> Iterator for Decoder<'a> {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos_input >= self.max_pos_input && self.count == 0 {
            None
        } else if self.count != 0 {
            self.count -= 1;
            Some(self.last_byte)
        } else {
            self.count = self.input[self.pos_input] - 1;
            self.last_byte = self.input[self.pos_input + 1];
            self.pos_input += 2;

            Some(self.last_byte)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const RAW: &[u8] = &[1, 1, 1, 1, 0, 3, 3, 4, 5, 6, 6];
    const COMPRESS: &[u8] = &[4, 1, 1, 0, 2, 3, 1, 4, 1, 5, 2, 6];
    const BAD_COMPRESS: &[u8] = &[4, 1, 1, 0, 2, 3, 1, 4, 1, 5, 2, 6];

    #[test]
    fn compress() {
        let mut out = Vec::new();

        for buffer in Encoder::new(RAW, 1) {
            for byte in buffer.iter() {
                out.push(*byte);
            }
        }

        assert_eq!(out, COMPRESS);
    }

    #[test]
    fn uncompress() {
        assert_eq!(Decoder::new(COMPRESS).collect::<Vec<u8>>(), RAW);
        assert_eq!(Decoder::new(BAD_COMPRESS).collect::<Vec<u8>>(), RAW);
    }
}
