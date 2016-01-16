#
#  Copyright (c) 2009-2016, Jack Poulson
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
lib.ElMedian_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElMedian_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElMedianDist_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElMedianDist_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElMedianDist_d.argtypes = [c_void_p,POINTER(dType)]
def Median(x):
  median = TagToType(x.tag)
  args = [x.obj,pointer(median)]
  if type(x) is Matrix:
    if   x.tag == iTag: lib.ElMedian_i(*args)
    elif x.tag == sTag: lib.ElMedian_s(*args)
    elif x.tag == dTag: lib.ElMedian_d(*args)
    else: DataExcept()
  elif type(x) is DistMatrix:
    if   x.tag == iTag: lib.ElMedianDist_i(*args)
    elif x.tag == sTag: lib.ElMedianDist_s(*args)
    elif x.tag == dTag: lib.ElMedianDist_d(*args)
    else: DataExcept()
  else: TypeExcept()
  return median

# Sort
# ====
lib.ElSort_i.argtypes = [c_void_p,c_uint]
lib.ElSort_s.argtypes = [c_void_p,c_uint]
lib.ElSort_d.argtypes = [c_void_p,c_uint]
lib.ElSortDist_i.argtypes = [c_void_p,c_uint]
lib.ElSortDist_s.argtypes = [c_void_p,c_uint]
lib.ElSortDist_d.argtypes = [c_void_p,c_uint]
def Sort(X,sort=ASCENDING):
  args = [X.obj,sort]
  if type(X) is Matrix:
    if   X.tag == iTag: lib.ElSort_i(*args)
    elif X.tag == sTag: lib.ElSort_s(*args)
    elif X.tag == dTag: lib.ElSort_d(*args)
    else: DataExcept()
  elif type(X) is DistMatrix:
    if   X.tag == iTag: lib.ElSortDist_i(*args)
    elif X.tag == sTag: lib.ElSortDist_s(*args)
    elif X.tag == dTag: lib.ElSortDist_d(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElTaggedSort_i.argtypes = [c_void_p,c_uint,POINTER(iType)]
lib.ElTaggedSort_s.argtypes = [c_void_p,c_uint,POINTER(sType)]
lib.ElTaggedSort_d.argtypes = [c_void_p,c_uint,POINTER(dType)]
lib.ElTaggedSortDist_i.argtypes = [c_void_p,c_uint,POINTER(iType)]
lib.ElTaggedSortDist_s.argtypes = [c_void_p,c_uint,POINTER(sType)]
lib.ElTaggedSortDist_d.argtypes = [c_void_p,c_uint,POINTER(dType)]
def TaggedSort(x,sort):
  taggedOrder = (TagToType(x.tag)*x.Height())()
  args = [x.obj,sort,taggedOrder]
  if type(x) is Matrix:
    if   x.tag == iTag: lib.ElTaggedSort_i(*args)
    elif x.tag == sTag: lib.ElTaggedSort_s(*args)
    elif x.tag == dTag: lib.ElTaggedSort_d(*args)
    else: DataExcept()
  elif type(x) is DistMatrix:
    if   x.tag == iTag: lib.ElTaggedSortDist_i(*args)
    elif x.tag == sTag: lib.ElTaggedSortDist_s(*args)
    elif x.tag == dTag: lib.ElTaggedSortDist_d(*args)
    else: DataExcept()
  else: TypeExcept()
  return taggedOrder
