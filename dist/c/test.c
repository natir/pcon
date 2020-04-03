/*
Copyright (c) 2020 Pierre Marijon <pmarijon@mpi-inf.mpg.de>

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

#include "pcon.h"

int main(void) {

  IO* error = pcon_error_new();
  Counter* counter = pcon_counter_new(5);

  /* Count a fasta file */
  pcon_counter_count_fasta(counter, "../data/test.fasta", error);
  if(*error != NoError) {
    printf("Error durring count of test.fasta error code %i\n", *error);
    return -1;
  }

  /* Serialize counter */
  pcon_serialize_counter(counter, "c_counter.pcon", error);
  if(*error != NoError) {
    printf("Error serialization of counter error code %i\n", *error);
    return -1;
  }
  
  printf("Kmer ACTGA (108) is present %u\n", pcon_counter_get(counter, 108));

  /* Free counter */
  pcon_counter_free(counter);

  /* Deserialize counter */
  counter = pcon_counter_new(5);
  pcon_deserialize_counter(counter, "c_counter.pcon", error);
  if(*error != NoError) {
    printf("Error durring deserialization of counter error code %i\n", *error);
    return -1;
  }
  
  printf("Kmer ACTGA (108) is present %u\n", pcon_counter_get(counter, 108));
  pcon_counter_inc(counter, 108);
  printf("Kmer ACTGA (108) is present %u\n", pcon_counter_get(counter, 108));  

  /* Convert count in solidity */ 
  Solid* solid = pcon_solid_from_counter(counter, 20);

  /* Serialize solid */
  pcon_serialize_solid(solid, "c_solid.pcon", error);
  if(*error != NoError) {
    printf("Error durring serialization of solid error code %i\n", *error);
    return -1;
  }
  printf("Kmer ACTGA (108) is solid threshold 20 %u\n", pcon_solid_get(solid, 108));  

  pcon_solid_free(solid);
  solid = pcon_solid_new(5);
  pcon_deserialize_solid(solid, "c_solid.pcon", error);
  if(*error != NoError) {
    printf("Error durring deserialization of solid error code %i\n", *error);
    return -1;
  }
  
  printf("Kmer ACTGA (108) is solid threshold 20 %u\n", pcon_solid_get(solid, 108));
  pcon_solid_set(solid, 108, 1);
  printf("Kmer ACTGA (108) is solid threshold 20 %u\n", pcon_solid_get(solid, 108));
  
  /* test dump */
  pcon_dump_csv(counter, 0, "c_counter.csv", error);
  if(*error != NoError) {
    printf("Error durring dump counter in csv error code %i\n", *error);
    return -1;
  }
  
  pcon_dump_solid(counter, 0, "c_counter.solid", error);
  if(*error != NoError) {
    printf("Error durring dump counter in solid error code %i\n", *error);
    return -1;
  }

  pcon_dump_spectrum(counter, "c_counter.spectrum.csv", error);
  if(*error != NoError) {
    printf("Error durring dump counter in spectrum error code %i\n", *error);
    return -1;
  }

  pcon_error_free(error);
  pcon_counter_free(counter);
}
