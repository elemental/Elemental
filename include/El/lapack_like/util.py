#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from ..core import *
import ctypes

# Median
# ======
lib.ElMedian_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElMedian_i.restype = c_uint
lib.ElMedian_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElMedian_s.restype = c_uint
lib.ElMedian_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElMedian_d.restype = c_uint
lib.ElMedianDist_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElMedianDist_i.restype = c_uint
lib.ElMedianDist_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElMedianDist_s.restype = c_uint
lib.ElMedianDist_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElMedianDist_d.restype = c_uint
def Median(x):
  median = TagToType(x.tag)
  if type(x) is Matrix:
    if   x.tag == iTag: lib.ElMedian_i(x.obj,pointer(median))
    elif x.tag == sTag: lib.ElMedian_s(x.obj,pointer(median))
    elif x.tag == dTag: lib.ElMedian_d(x.obj,pointer(median))
    else: raise Exception('Unsupported datatype')
  elif type(x) is DistMatrix:
    if   x.tag == iTag: lib.ElMedianDist_i(x.obj,pointer(median))
    elif x.tag == sTag: lib.ElMedianDist_s(x.obj,pointer(median))
    elif x.tag == dTag: lib.ElMedianDist_d(x.obj,pointer(median))
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')
  return median

# Sort
# ====
lib.ElSort_i.argtypes = [c_void_p,c_uint]
lib.ElSort_i.restype = c_uint
lib.ElSort_s.argtypes = [c_void_p,c_uint]
lib.ElSort_s.restype = c_uint
lib.ElSort_d.argtypes = [c_void_p,c_uint]
lib.ElSort_d.restype = c_uint
lib.ElSortDist_i.argtypes = [c_void_p,c_uint]
lib.ElSortDist_i.restype = c_uint
lib.ElSortDist_s.argtypes = [c_void_p,c_uint]
lib.ElSortDist_s.restype = c_uint
lib.ElSortDist_d.argtypes = [c_void_p,c_uint]
lib.ElSortDist_d.restype = c_uint
def Sort(X,sort=ASCENDING):
  if type(X) is Matrix:
    if   X.tag == iTag: lib.ElSort_i(X.obj,sort)
    elif X.tag == sTag: lib.ElSort_s(X.obj,sort)
    elif X.tag == dTag: lib.ElSort_d(X.obj,sort)
    else: raise Exception('Unsupported datatype')
  elif type(X) is DistMatrix:
    if   X.tag == iTag: lib.ElSortDist_i(X.obj,sort)
    elif X.tag == sTag: lib.ElSortDist_s(X.obj,sort)
    elif X.tag == dTag: lib.ElSortDist_d(X.obj,sort)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElTaggedSort_i.argtypes = [c_void_p,c_uint,POINTER(iType)]
lib.ElTaggedSort_i.restype = c_uint
lib.ElTaggedSort_s.argtypes = [c_void_p,c_uint,POINTER(sType)]
lib.ElTaggedSort_s.restype = c_uint
lib.ElTaggedSort_d.argtypes = [c_void_p,c_uint,POINTER(dType)]
lib.ElTaggedSort_d.restype = c_uint
lib.ElTaggedSortDist_i.argtypes = [c_void_p,c_uint,POINTER(iType)]
lib.ElTaggedSortDist_i.restype = c_uint
lib.ElTaggedSortDist_s.argtypes = [c_void_p,c_uint,POINTER(sType)]
lib.ElTaggedSortDist_s.restype = c_uint
lib.ElTaggedSortDist_d.argtypes = [c_void_p,c_uint,POINTER(dType)]
lib.ElTaggedSortDist_d.restype = c_uint
def TaggedSort(x,sort):
  taggedOrder = (TagToType(x.tag)*x.Height())()
  if type(x) is Matrix:
    if   x.tag == iTag: lib.ElTaggedSort_i(x.obj,sort,taggedOrder)
    elif x.tag == sTag: lib.ElTaggedSort_s(x.obj,sort,taggedOrder)
    elif x.tag == dTag: lib.ElTaggedSort_d(x.obj,sort,taggedOrder)
    else: raise Exception('Unsupported datatype')
  elif type(x) is DistMatrix:
    if   x.tag == iTag: lib.ElTaggedSortDist_i(x.obj,sort,taggedOrder)
    elif x.tag == sTag: lib.ElTaggedSortDist_s(x.obj,sort,taggedOrder)
    elif x.tag == dTag: lib.ElTaggedSortDist_d(x.obj,sort,taggedOrder)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')
  return taggedOrder
