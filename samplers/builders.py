# This file was created automatically by SWIG 1.3.28.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _builders
import new
new_instancemethod = new.instancemethod
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'PySwigObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


class PySwigIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, PySwigIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, PySwigIterator, name)
    def __init__(self): raise AttributeError, "No constructor defined"
    def __repr__(self):
        try: strthis = "at 0x%x" %( self.this, ) 
        except: strthis = "" 
        return "<%s.%s; proxy of C++ swig::PySwigIterator instance %s>" % (self.__class__.__module__, self.__class__.__name__, strthis,)
    __swig_destroy__ = _builders.delete_PySwigIterator
    __del__ = lambda self : None;
    def value(*args): return _builders.PySwigIterator_value(*args)
    def incr(*args): return _builders.PySwigIterator_incr(*args)
    def decr(*args): return _builders.PySwigIterator_decr(*args)
    def distance(*args): return _builders.PySwigIterator_distance(*args)
    def equal(*args): return _builders.PySwigIterator_equal(*args)
    def copy(*args): return _builders.PySwigIterator_copy(*args)
    def next(*args): return _builders.PySwigIterator_next(*args)
    def previous(*args): return _builders.PySwigIterator_previous(*args)
    def advance(*args): return _builders.PySwigIterator_advance(*args)
    def __eq__(*args): return _builders.PySwigIterator___eq__(*args)
    def __ne__(*args): return _builders.PySwigIterator___ne__(*args)
    def __iadd__(*args): return _builders.PySwigIterator___iadd__(*args)
    def __isub__(*args): return _builders.PySwigIterator___isub__(*args)
    def __add__(*args): return _builders.PySwigIterator___add__(*args)
    def __sub__(*args): return _builders.PySwigIterator___sub__(*args)
    def __iter__(self):
      return self

_builders.PySwigIterator_swigregister(PySwigIterator)

class PhipsiSampler(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, PhipsiSampler, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, PhipsiSampler, name)
    def __repr__(self):
        try: strthis = "at 0x%x" %( self.this, ) 
        except: strthis = "" 
        return "<%s.%s; proxy of C++ PhipsiSampler instance %s>" % (self.__class__.__module__, self.__class__.__name__, strthis,)
    def __init__(self, *args):
        this = _builders.new_PhipsiSampler(*args)
        try: self.this.append(this)
        except: self.this = this
    def sample(*args): return _builders.PhipsiSampler_sample(*args)
    def readPhipsiMap(*args): return _builders.PhipsiSampler_readPhipsiMap(*args)
    __swig_destroy__ = _builders.delete_PhipsiSampler
    __del__ = lambda self : None;
_builders.PhipsiSampler_swigregister(PhipsiSampler)

class OmegaSampler(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, OmegaSampler, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, OmegaSampler, name)
    def __repr__(self):
        try: strthis = "at 0x%x" %( self.this, ) 
        except: strthis = "" 
        return "<%s.%s; proxy of C++ OmegaSampler instance %s>" % (self.__class__.__module__, self.__class__.__name__, strthis,)
    def __init__(self, *args):
        this = _builders.new_OmegaSampler(*args)
        try: self.this.append(this)
        except: self.this = this
    def sample(*args): return _builders.OmegaSampler_sample(*args)
    __swig_destroy__ = _builders.delete_OmegaSampler
    __del__ = lambda self : None;
_builders.OmegaSampler_swigregister(OmegaSampler)



