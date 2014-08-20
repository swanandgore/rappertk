# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.2
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.
# This file is compatible with both classic and new-style classes.

from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_geometry', [dirname(__file__)])
        except ImportError:
            import _geometry
            return _geometry
        if fp is not None:
            try:
                _mod = imp.load_module('_geometry', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _geometry = swig_import_helper()
    del swig_import_helper
else:
    import _geometry
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
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
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class SwigPyIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SwigPyIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SwigPyIterator, name)
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _geometry.delete_SwigPyIterator
    __del__ = lambda self : None;
    def value(self): return _geometry.SwigPyIterator_value(self)
    def incr(self, n = 1): return _geometry.SwigPyIterator_incr(self, n)
    def decr(self, n = 1): return _geometry.SwigPyIterator_decr(self, n)
    def distance(self, *args): return _geometry.SwigPyIterator_distance(self, *args)
    def equal(self, *args): return _geometry.SwigPyIterator_equal(self, *args)
    def copy(self): return _geometry.SwigPyIterator_copy(self)
    def next(self): return _geometry.SwigPyIterator_next(self)
    def __next__(self): return _geometry.SwigPyIterator___next__(self)
    def previous(self): return _geometry.SwigPyIterator_previous(self)
    def advance(self, *args): return _geometry.SwigPyIterator_advance(self, *args)
    def __eq__(self, *args): return _geometry.SwigPyIterator___eq__(self, *args)
    def __ne__(self, *args): return _geometry.SwigPyIterator___ne__(self, *args)
    def __iadd__(self, *args): return _geometry.SwigPyIterator___iadd__(self, *args)
    def __isub__(self, *args): return _geometry.SwigPyIterator___isub__(self, *args)
    def __add__(self, *args): return _geometry.SwigPyIterator___add__(self, *args)
    def __sub__(self, *args): return _geometry.SwigPyIterator___sub__(self, *args)
    def __iter__(self): return self
SwigPyIterator_swigregister = _geometry.SwigPyIterator_swigregister
SwigPyIterator_swigregister(SwigPyIterator)

class VecFloat(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, VecFloat, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, VecFloat, name)
    __repr__ = _swig_repr
    def iterator(self): return _geometry.VecFloat_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _geometry.VecFloat___nonzero__(self)
    def __bool__(self): return _geometry.VecFloat___bool__(self)
    def __len__(self): return _geometry.VecFloat___len__(self)
    def pop(self): return _geometry.VecFloat_pop(self)
    def __getslice__(self, *args): return _geometry.VecFloat___getslice__(self, *args)
    def __setslice__(self, *args): return _geometry.VecFloat___setslice__(self, *args)
    def __delslice__(self, *args): return _geometry.VecFloat___delslice__(self, *args)
    def __delitem__(self, *args): return _geometry.VecFloat___delitem__(self, *args)
    def __getitem__(self, *args): return _geometry.VecFloat___getitem__(self, *args)
    def __setitem__(self, *args): return _geometry.VecFloat___setitem__(self, *args)
    def append(self, *args): return _geometry.VecFloat_append(self, *args)
    def empty(self): return _geometry.VecFloat_empty(self)
    def size(self): return _geometry.VecFloat_size(self)
    def clear(self): return _geometry.VecFloat_clear(self)
    def swap(self, *args): return _geometry.VecFloat_swap(self, *args)
    def get_allocator(self): return _geometry.VecFloat_get_allocator(self)
    def begin(self): return _geometry.VecFloat_begin(self)
    def end(self): return _geometry.VecFloat_end(self)
    def rbegin(self): return _geometry.VecFloat_rbegin(self)
    def rend(self): return _geometry.VecFloat_rend(self)
    def pop_back(self): return _geometry.VecFloat_pop_back(self)
    def erase(self, *args): return _geometry.VecFloat_erase(self, *args)
    def __init__(self, *args): 
        this = _geometry.new_VecFloat(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _geometry.VecFloat_push_back(self, *args)
    def front(self): return _geometry.VecFloat_front(self)
    def back(self): return _geometry.VecFloat_back(self)
    def assign(self, *args): return _geometry.VecFloat_assign(self, *args)
    def resize(self, *args): return _geometry.VecFloat_resize(self, *args)
    def insert(self, *args): return _geometry.VecFloat_insert(self, *args)
    def reserve(self, *args): return _geometry.VecFloat_reserve(self, *args)
    def capacity(self): return _geometry.VecFloat_capacity(self)
    __swig_destroy__ = _geometry.delete_VecFloat
    __del__ = lambda self : None;
VecFloat_swigregister = _geometry.VecFloat_swigregister
VecFloat_swigregister(VecFloat)

class VecLong(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, VecLong, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, VecLong, name)
    __repr__ = _swig_repr
    def iterator(self): return _geometry.VecLong_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _geometry.VecLong___nonzero__(self)
    def __bool__(self): return _geometry.VecLong___bool__(self)
    def __len__(self): return _geometry.VecLong___len__(self)
    def pop(self): return _geometry.VecLong_pop(self)
    def __getslice__(self, *args): return _geometry.VecLong___getslice__(self, *args)
    def __setslice__(self, *args): return _geometry.VecLong___setslice__(self, *args)
    def __delslice__(self, *args): return _geometry.VecLong___delslice__(self, *args)
    def __delitem__(self, *args): return _geometry.VecLong___delitem__(self, *args)
    def __getitem__(self, *args): return _geometry.VecLong___getitem__(self, *args)
    def __setitem__(self, *args): return _geometry.VecLong___setitem__(self, *args)
    def append(self, *args): return _geometry.VecLong_append(self, *args)
    def empty(self): return _geometry.VecLong_empty(self)
    def size(self): return _geometry.VecLong_size(self)
    def clear(self): return _geometry.VecLong_clear(self)
    def swap(self, *args): return _geometry.VecLong_swap(self, *args)
    def get_allocator(self): return _geometry.VecLong_get_allocator(self)
    def begin(self): return _geometry.VecLong_begin(self)
    def end(self): return _geometry.VecLong_end(self)
    def rbegin(self): return _geometry.VecLong_rbegin(self)
    def rend(self): return _geometry.VecLong_rend(self)
    def pop_back(self): return _geometry.VecLong_pop_back(self)
    def erase(self, *args): return _geometry.VecLong_erase(self, *args)
    def __init__(self, *args): 
        this = _geometry.new_VecLong(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _geometry.VecLong_push_back(self, *args)
    def front(self): return _geometry.VecLong_front(self)
    def back(self): return _geometry.VecLong_back(self)
    def assign(self, *args): return _geometry.VecLong_assign(self, *args)
    def resize(self, *args): return _geometry.VecLong_resize(self, *args)
    def insert(self, *args): return _geometry.VecLong_insert(self, *args)
    def reserve(self, *args): return _geometry.VecLong_reserve(self, *args)
    def capacity(self): return _geometry.VecLong_capacity(self)
    __swig_destroy__ = _geometry.delete_VecLong
    __del__ = lambda self : None;
VecLong_swigregister = _geometry.VecLong_swigregister
VecLong_swigregister(VecLong)

class VecInt(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, VecInt, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, VecInt, name)
    __repr__ = _swig_repr
    def iterator(self): return _geometry.VecInt_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _geometry.VecInt___nonzero__(self)
    def __bool__(self): return _geometry.VecInt___bool__(self)
    def __len__(self): return _geometry.VecInt___len__(self)
    def pop(self): return _geometry.VecInt_pop(self)
    def __getslice__(self, *args): return _geometry.VecInt___getslice__(self, *args)
    def __setslice__(self, *args): return _geometry.VecInt___setslice__(self, *args)
    def __delslice__(self, *args): return _geometry.VecInt___delslice__(self, *args)
    def __delitem__(self, *args): return _geometry.VecInt___delitem__(self, *args)
    def __getitem__(self, *args): return _geometry.VecInt___getitem__(self, *args)
    def __setitem__(self, *args): return _geometry.VecInt___setitem__(self, *args)
    def append(self, *args): return _geometry.VecInt_append(self, *args)
    def empty(self): return _geometry.VecInt_empty(self)
    def size(self): return _geometry.VecInt_size(self)
    def clear(self): return _geometry.VecInt_clear(self)
    def swap(self, *args): return _geometry.VecInt_swap(self, *args)
    def get_allocator(self): return _geometry.VecInt_get_allocator(self)
    def begin(self): return _geometry.VecInt_begin(self)
    def end(self): return _geometry.VecInt_end(self)
    def rbegin(self): return _geometry.VecInt_rbegin(self)
    def rend(self): return _geometry.VecInt_rend(self)
    def pop_back(self): return _geometry.VecInt_pop_back(self)
    def erase(self, *args): return _geometry.VecInt_erase(self, *args)
    def __init__(self, *args): 
        this = _geometry.new_VecInt(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _geometry.VecInt_push_back(self, *args)
    def front(self): return _geometry.VecInt_front(self)
    def back(self): return _geometry.VecInt_back(self)
    def assign(self, *args): return _geometry.VecInt_assign(self, *args)
    def resize(self, *args): return _geometry.VecInt_resize(self, *args)
    def insert(self, *args): return _geometry.VecInt_insert(self, *args)
    def reserve(self, *args): return _geometry.VecInt_reserve(self, *args)
    def capacity(self): return _geometry.VecInt_capacity(self)
    __swig_destroy__ = _geometry.delete_VecInt
    __del__ = lambda self : None;
VecInt_swigregister = _geometry.VecInt_swigregister
VecInt_swigregister(VecInt)

class VecVecInt(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, VecVecInt, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, VecVecInt, name)
    __repr__ = _swig_repr
    def iterator(self): return _geometry.VecVecInt_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _geometry.VecVecInt___nonzero__(self)
    def __bool__(self): return _geometry.VecVecInt___bool__(self)
    def __len__(self): return _geometry.VecVecInt___len__(self)
    def pop(self): return _geometry.VecVecInt_pop(self)
    def __getslice__(self, *args): return _geometry.VecVecInt___getslice__(self, *args)
    def __setslice__(self, *args): return _geometry.VecVecInt___setslice__(self, *args)
    def __delslice__(self, *args): return _geometry.VecVecInt___delslice__(self, *args)
    def __delitem__(self, *args): return _geometry.VecVecInt___delitem__(self, *args)
    def __getitem__(self, *args): return _geometry.VecVecInt___getitem__(self, *args)
    def __setitem__(self, *args): return _geometry.VecVecInt___setitem__(self, *args)
    def append(self, *args): return _geometry.VecVecInt_append(self, *args)
    def empty(self): return _geometry.VecVecInt_empty(self)
    def size(self): return _geometry.VecVecInt_size(self)
    def clear(self): return _geometry.VecVecInt_clear(self)
    def swap(self, *args): return _geometry.VecVecInt_swap(self, *args)
    def get_allocator(self): return _geometry.VecVecInt_get_allocator(self)
    def begin(self): return _geometry.VecVecInt_begin(self)
    def end(self): return _geometry.VecVecInt_end(self)
    def rbegin(self): return _geometry.VecVecInt_rbegin(self)
    def rend(self): return _geometry.VecVecInt_rend(self)
    def pop_back(self): return _geometry.VecVecInt_pop_back(self)
    def erase(self, *args): return _geometry.VecVecInt_erase(self, *args)
    def __init__(self, *args): 
        this = _geometry.new_VecVecInt(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _geometry.VecVecInt_push_back(self, *args)
    def front(self): return _geometry.VecVecInt_front(self)
    def back(self): return _geometry.VecVecInt_back(self)
    def assign(self, *args): return _geometry.VecVecInt_assign(self, *args)
    def resize(self, *args): return _geometry.VecVecInt_resize(self, *args)
    def insert(self, *args): return _geometry.VecVecInt_insert(self, *args)
    def reserve(self, *args): return _geometry.VecVecInt_reserve(self, *args)
    def capacity(self): return _geometry.VecVecInt_capacity(self)
    __swig_destroy__ = _geometry.delete_VecVecInt
    __del__ = lambda self : None;
VecVecInt_swigregister = _geometry.VecVecInt_swigregister
VecVecInt_swigregister(VecVecInt)

class VoxelIndex(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, VoxelIndex, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, VoxelIndex, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _geometry.new_VoxelIndex(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_setmethods__["i"] = _geometry.VoxelIndex_i_set
    __swig_getmethods__["i"] = _geometry.VoxelIndex_i_get
    if _newclass:i = _swig_property(_geometry.VoxelIndex_i_get, _geometry.VoxelIndex_i_set)
    __swig_setmethods__["j"] = _geometry.VoxelIndex_j_set
    __swig_getmethods__["j"] = _geometry.VoxelIndex_j_get
    if _newclass:j = _swig_property(_geometry.VoxelIndex_j_get, _geometry.VoxelIndex_j_set)
    __swig_setmethods__["k"] = _geometry.VoxelIndex_k_set
    __swig_getmethods__["k"] = _geometry.VoxelIndex_k_get
    if _newclass:k = _swig_property(_geometry.VoxelIndex_k_get, _geometry.VoxelIndex_k_set)
    __swig_destroy__ = _geometry.delete_VoxelIndex
    __del__ = lambda self : None;
VoxelIndex_swigregister = _geometry.VoxelIndex_swigregister
VoxelIndex_swigregister(VoxelIndex)

class compareVoxelIndices(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, compareVoxelIndices, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, compareVoxelIndices, name)
    __repr__ = _swig_repr
    def __call__(self, *args): return _geometry.compareVoxelIndices___call__(self, *args)
    def __init__(self): 
        this = _geometry.new_compareVoxelIndices()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _geometry.delete_compareVoxelIndices
    __del__ = lambda self : None;
compareVoxelIndices_swigregister = _geometry.compareVoxelIndices_swigregister
compareVoxelIndices_swigregister(compareVoxelIndices)

class MapIntIntFloat(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, MapIntIntFloat, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, MapIntIntFloat, name)
    __repr__ = _swig_repr
    def __init__(self): 
        this = _geometry.new_MapIntIntFloat()
        try: self.this.append(this)
        except: self.this = this
    def put(self, *args): return _geometry.MapIntIntFloat_put(self, *args)
    def get(self, *args): return _geometry.MapIntIntFloat_get(self, *args)
    def contains(self, *args): return _geometry.MapIntIntFloat_contains(self, *args)
    __swig_destroy__ = _geometry.delete_MapIntIntFloat
    __del__ = lambda self : None;
MapIntIntFloat_swigregister = _geometry.MapIntIntFloat_swigregister
MapIntIntFloat_swigregister(MapIntIntFloat)

class Grid(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Grid, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Grid, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _geometry.new_Grid(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _geometry.delete_Grid
    __del__ = lambda self : None;
    def add(self, *args): return _geometry.Grid_add(self, *args)
    def remove(self, *args): return _geometry.Grid_remove(self, *args)
    def removes(self, *args): return _geometry.Grid_removes(self, *args)
    __swig_getmethods__["xyz2ijk"] = lambda x: _geometry.Grid_xyz2ijk
    if _newclass:xyz2ijk = staticmethod(_geometry.Grid_xyz2ijk)
    def justAdd(self, *args): return _geometry.Grid_justAdd(self, *args)
    def printIJK(self, *args): return _geometry.Grid_printIJK(self, *args)
    def withinEnvelope(self, *args): return _geometry.Grid_withinEnvelope(self, *args)
    def popCopyGridPoints(self, *args): return _geometry.Grid_popCopyGridPoints(self, *args)
    def describe(self): return _geometry.Grid_describe(self)
Grid_swigregister = _geometry.Grid_swigregister
Grid_swigregister(Grid)

def Grid_xyz2ijk(*args):
  return _geometry.Grid_xyz2ijk(*args)
Grid_xyz2ijk = _geometry.Grid_xyz2ijk

class GridHelper(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, GridHelper, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, GridHelper, name)
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    def clash(self, *args): return _geometry.GridHelper_clash(self, *args)
    def cellsym(self, *args): return _geometry.GridHelper_cellsym(self, *args)
    def setCellsym(self, *args): return _geometry.GridHelper_setCellsym(self, *args)
    __swig_destroy__ = _geometry.delete_GridHelper
    __del__ = lambda self : None;
GridHelper_swigregister = _geometry.GridHelper_swigregister
GridHelper_swigregister(GridHelper)

class CAtraceGH(GridHelper):
    __swig_setmethods__ = {}
    for _s in [GridHelper]: __swig_setmethods__.update(getattr(_s,'__swig_setmethods__',{}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, CAtraceGH, name, value)
    __swig_getmethods__ = {}
    for _s in [GridHelper]: __swig_getmethods__.update(getattr(_s,'__swig_getmethods__',{}))
    __getattr__ = lambda self, name: _swig_getattr(self, CAtraceGH, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _geometry.new_CAtraceGH(*args)
        try: self.this.append(this)
        except: self.this = this
    def clash(self, *args): return _geometry.CAtraceGH_clash(self, *args)
    def cellsym(self, *args): return _geometry.CAtraceGH_cellsym(self, *args)
    def setCellsym(self, *args): return _geometry.CAtraceGH_setCellsym(self, *args)
    __swig_destroy__ = _geometry.delete_CAtraceGH
    __del__ = lambda self : None;
CAtraceGH_swigregister = _geometry.CAtraceGH_swigregister
CAtraceGH_swigregister(CAtraceGH)


def EigenSort(*args):
  return _geometry.EigenSort(*args)
EigenSort = _geometry.EigenSort

def Jacobi(*args):
  return _geometry.Jacobi(*args)
Jacobi = _geometry.Jacobi

def TRTtransform(*args):
  return _geometry.TRTtransform(*args)
TRTtransform = _geometry.TRTtransform

def TRTtransform1(*args):
  return _geometry.TRTtransform1(*args)
TRTtransform1 = _geometry.TRTtransform1

def findSuperpositionTransform(*args):
  return _geometry.findSuperpositionTransform(*args)
findSuperpositionTransform = _geometry.findSuperpositionTransform

def sphereVolSample(*args):
  return _geometry.sphereVolSample(*args)
sphereVolSample = _geometry.sphereVolSample

def dihedDiff(*args):
  return _geometry.dihedDiff(*args)
dihedDiff = _geometry.dihedDiff

def withinPlusMinus180(*args):
  return _geometry.withinPlusMinus180(*args)
withinPlusMinus180 = _geometry.withinPlusMinus180

def findTriangleBaseHeight(*args):
  return _geometry.findTriangleBaseHeight(*args)
findTriangleBaseHeight = _geometry.findTriangleBaseHeight

def calcDist(*args):
  return _geometry.calcDist(*args)
calcDist = _geometry.calcDist

def calcDihed(*args):
  return _geometry.calcDihed(*args)
calcDihed = _geometry.calcDihed

def find4thPoint(*args):
  return _geometry.find4thPoint(*args)
find4thPoint = _geometry.find4thPoint

def magnitude(*args):
  return _geometry.magnitude(*args)
magnitude = _geometry.magnitude

def normalize_point(*args):
  return _geometry.normalize_point(*args)
normalize_point = _geometry.normalize_point

def linear_combination(*args):
  return _geometry.linear_combination(*args)
linear_combination = _geometry.linear_combination

def findYZgivenX(*args):
  return _geometry.findYZgivenX(*args)
findYZgivenX = _geometry.findYZgivenX

def findSphereSphereAngleIntx(*args):
  return _geometry.findSphereSphereAngleIntx(*args)
findSphereSphereAngleIntx = _geometry.findSphereSphereAngleIntx

def randomNormalVector(*args):
  return _geometry.randomNormalVector(*args)
randomNormalVector = _geometry.randomNormalVector

def rotate(*args):
  return _geometry.rotate(*args)
rotate = _geometry.rotate

def findRotnOperator(*args):
  return _geometry.findRotnOperator(*args)
findRotnOperator = _geometry.findRotnOperator

def cross_product(*args):
  return _geometry.cross_product(*args)
cross_product = _geometry.cross_product


def calcAngle(*args):
  return _geometry.calcAngle(*args)
calcAngle = _geometry.calcAngle
