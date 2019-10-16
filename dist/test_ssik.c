/*
Copyright (c) 2019 Pierre Marijon <pmarijon@mmci.univ-saarland.de>

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

#include <stdio.h>

#include "ssik.h"

int main(void) {
  printf("A: %d\n", ssik_nuc2bit('A'));
  printf("T: %d\n", ssik_nuc2bit('T'));
  printf("C: %d\n", ssik_nuc2bit('C'));
  printf("G: %d\n", ssik_nuc2bit('G'));
  printf("\n");

  unsigned char subseq[] = "ACG";
  printf("%s: 2bit: %d hash: %d\n", subseq, ssik_seq2bit(subseq, 3), ssik_seq2bit(subseq, 3) >> 1);

  unsigned char revcomp[] = "CGT";
  printf("2bit %d revcomp %d cannonical %d\n", ssik_seq2bit(subseq, 3), ssik_seq2bit(revcomp, 3), ssik_cannonical(ssik_seq2bit(subseq, 3), 3));  
  printf("\n");

  uint8_t k = 0;
  uint8_t nb_bit = 0;
  char path[] = "./dist/small.ssik";
  ssik_get_header(path, &k, &nb_bit);
  unsigned char* data = (unsigned char*) malloc(ssik_get_data_size(k, nb_bit));
  ssik_read_count(path, data);

  for(int i = 0; i != (1 << (k * 2 - 1)) / 2; i++) {
    printf("%d, ", data[i]);
  }
  printf("\n\n");
  
  printf("k: %d\n", k);
  printf("nb_bit: %d\n", nb_bit);
  printf("\n");
  
  for(int i = 0; i < (1 << (k * 2 -1)) / 2; i++) {
    printf("%d, ", data[i]);
  }
  printf("\n\n");

  for(uint64_t hash = 0; hash != ssik_get_kmer_space_size(k); hash++) {
    printf("%s %d\n", ssik_revhash(hash, k), ssik_get_count(data, hash, nb_bit));
  }
  printf("\n\n");
}
