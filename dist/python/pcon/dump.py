from ._lowlevel import ffi, lib

from .utils import error2str

from .counter import Counter

def csv(counter: Counter, abundance: int, path: str):
    _call_can_produce_error(lib.pcon_dump_csv, counter._get_objptr(), abundance, path.encode("utf-8"))

def solid(counter: Counter, abundance: int, path: str):
    _call_can_produce_error(lib.pcon_dump_solid, counter._get_objptr(), abundance, path.encode("utf-8"))
    
def spectrum(counter: Counter, path: str):
    _call_can_produce_error(lib.pcon_dump_spectrum, counter._get_objptr(), path.encode("utf-8"))
    
def _call_can_produce_error(func, *args):
    error = ffi.new("IO *")
    error[0] = lib.NoError

    ret = func(*args, error)
    
    if error[0] != lib.NoError:
        raise Exception(error2str(error[0]))

    return ret
