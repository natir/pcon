
from ._lowlevel import ffi, lib

from .utils import RustObject

from .counter import Counter

class Solid(RustObject):
    __dealloc_func__ = lib.pcon_solid_free

    def __init__(self, k: int):
        self._objptr = lib.pcon_solid_new(k)

    @classmethod
    def from_counter(cls, counter: Counter, abundance: int):
        rv = object.__new__(cls)
        rv._objptr = lib.pcon_solid_from_counter(counter._get_objptr(), abundance)

        return rv

    def set(self, kmer: int, value: bool):
        self._methodcall(lib.pcon_solid_set, kmer, value)

    def set_canonic(self, kmer: int, value: bool):
        self._methodcall(lib.pcon_solid_set_canonic, kmer, value)

    def get(self, kmer: int) -> bool: 
        return self._methodcall(lib.pcon_solid_get, kmer)

    def get_canonic(self, kmer: int) -> bool: 
        return self._methodcall(lib.pcon_solid_get_canonic, kmer)
    
    def serialize(self, path: str):
        self._methodcall_can_produce_error(lib.pcon_serialize_solid, path.encode("utf-8"))

    @classmethod
    def deserialize(cls, path: str):
        rv = object.__new__(cls)
        rv._objptr = lib.pcon_solid_new(1)

        rv._methodcall_can_produce_error(lib.pcon_deserialize_solid, path.encode("utf-8"))

        return rv
