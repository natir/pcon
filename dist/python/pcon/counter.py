
from ._lowlevel import ffi, lib

from .utils import RustObject

class Counter(RustObject):
    __dealloc_func__ = lib.pcon_counter_free

    def __init__(self, k: int):
        self._objptr = lib.pcon_counter_new(k)

    def count_fasta(self, path: str):
        self._methodcall_can_produce_error(lib.pcon_counter_count_fasta, path.encode("utf-8"))

    def count_fastq(self, path: str):
        self._methodcall_can_produce_error(lib.pcon_counter_count_fastq, path.encode("utf-8"))

    def inc(self, kmer: int):
        self._methodcall(lib.pcon_counter_inc, kmer)

    def inc_canonic(self, kmer: int):
        self._methodcall(lib.pcon_counter_inc_canonic, kmer)
        
    def get(self, kmer: int) -> int:
        return self._methodcall(lib.pcon_counter_get, kmer)

    def get_canonic(self, kmer: int) -> int:
        return self._methodcall(lib.pcon_counter_get_canonic, kmer)
    
    def serialize(self, path: str):
        self._methodcall_can_produce_error(lib.pcon_serialize_counter, path.encode("utf-8"))

    @classmethod
    def deserialize(cls, path: str):
        rv = object.__new__(cls)
        rv._objptr = lib.pcon_counter_new(1)

        rv._methodcall_can_produce_error(lib.pcon_deserialize_counter, path.encode("utf-8"))

        return rv
