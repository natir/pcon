from ._lowlevel import ffi, lib

def error2str(error) -> str:
    if error == 0:
        return "Error: can't create file"
    elif error == 1:
        return "Error: can't open file"
    elif error == 2:
        return "Error: durring write file"
    elif error == 3:
        return "Error: durring read file"
    else:
        return "No error"

class RustObject(object):
    __dealloc_func__ = None
    _objptr = None
    _shared = False

    def __init__(self):
        raise TypeError("Cannot instanciate %r objects" % self.__class__.__name__)

    @classmethod
    def _from_objptr(cls, ptr, shared=False):
        rv = object.__new__(cls)
        rv._objptr = ptr
        rv._shared = shared
        return rv

    def _get_objptr(self):
        if not self._objptr:
            raise RuntimeError("Object is closed")
        return self._objptr

    def __del__(self):
        if self._objptr is None or self._shared:
            return
        f = self.__class__.__dealloc_func__
        if f is not None:
            f(self._objptr)
            self._objptr = None
    
    def _methodcall(self, func, *args):
        return func(self._get_objptr(), *args)

    def _methodcall_can_produce_error(self, func, *args):
        error = ffi.new("IO *")
        error[0] = lib.NoError
        
        ret = self._methodcall(func, *args, error)

        if error[0] != lib.NoError:
            raise Exception(error2str(error[0]))

        return ret

