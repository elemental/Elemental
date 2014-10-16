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

# Hermitian tridiagonal eigensolvers
# ==================================
class HermitianEigSubset_s(ctypes.Structure):
  _fields_ = [("indexSubset",bType),
              ("lowerIndex",iType),("upperIndex",iType),
              ("rangeSubset",bType),
              ("lowerBound",sType),("upperBound",sType)]
class HermitianEigSubset_d(ctypes.Structure):
  _fields_ = [("indexSubset",bType),
              ("lowerIndex",iType),("upperIndex",iType),
              ("rangeSubset",bType),
              ("lowerBound",dType),("upperBound",dType)]

lib.ElHermitianTridiagEig_s.argtypes = [c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianTridiagEig_s.restype = c_uint
lib.ElHermitianTridiagEig_d.argtypes = [c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianTridiagEig_d.restype = c_uint
lib.ElHermitianTridiagEig_c.argtypes = [c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianTridiagEig_c.restype = c_uint
lib.ElHermitianTridiagEig_z.argtypes = [c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianTridiagEig_z.restype = c_uint
lib.ElHermitianTridiagEigDist_s.argtypes = [c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianTridiagEigDist_s.restype = c_uint
lib.ElHermitianTridiagEigDist_d.argtypes = [c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianTridiagEigDist_d.restype = c_uint
lib.ElHermitianTridiagEigDist_c.argtypes = [c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianTridiagEigDist_c.restype = c_uint
lib.ElHermitianTridiagEigDist_z.argtypes = [c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianTridiagEigDist_z.restype = c_uint
lib.ElHermitianTridiagEigPair_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianTridiagEigPair_s.restype = c_uint
lib.ElHermitianTridiagEigPair_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianTridiagEigPair_d.restype = c_uint
lib.ElHermitianTridiagEigPair_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianTridiagEigPair_c.restype = c_uint
lib.ElHermitianTridiagEigPair_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianTridiagEigPair_z.restype = c_uint
lib.ElHermitianTridiagEigPairDist_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianTridiagEigPairDist_s.restype = c_uint
lib.ElHermitianTridiagEigPairDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianTridiagEigPairDist_d.restype = c_uint
lib.ElHermitianTridiagEigPairDist_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianTridiagEigPairDist_c.restype = c_uint
lib.ElHermitianTridiagEigPairDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianTridiagEigPairDist_z.restype = c_uint

lib.ElHermitianTridiagEigPartial_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]
lib.ElHermitianTridiagEigPartial_s.restype = c_uint
lib.ElHermitianTridiagEigPartial_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]
lib.ElHermitianTridiagEigPartial_d.restype = c_uint
lib.ElHermitianTridiagEigPartial_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]
lib.ElHermitianTridiagEigPartial_c.restype = c_uint
lib.ElHermitianTridiagEigPartial_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]
lib.ElHermitianTridiagEigPartial_z.restype = c_uint
lib.ElHermitianTridiagEigPartialDist_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]
lib.ElHermitianTridiagEigPartialDist_s.restype = c_uint
lib.ElHermitianTridiagEigPartialDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]
lib.ElHermitianTridiagEigPartialDist_d.restype = c_uint
lib.ElHermitianTridiagEigPartialDist_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]
lib.ElHermitianTridiagEigPartialDist_c.restype = c_uint
lib.ElHermitianTridiagEigPartialDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]
lib.ElHermitianTridiagEigPartialDist_z.restype = c_uint
lib.ElHermitianTridiagEigPairPartial_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]
lib.ElHermitianTridiagEigPairPartial_s.restype = c_uint
lib.ElHermitianTridiagEigPairPartial_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]
lib.ElHermitianTridiagEigPairPartial_d.restype = c_uint
lib.ElHermitianTridiagEigPairPartial_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]
lib.ElHermitianTridiagEigPairPartial_c.restype = c_uint
lib.ElHermitianTridiagEigPairPartial_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]
lib.ElHermitianTridiagEigPairPartial_z.restype = c_uint
lib.ElHermitianTridiagEigPairPartialDist_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]
lib.ElHermitianTridiagEigPairPartialDist_s.restype = c_uint
lib.ElHermitianTridiagEigPairPartialDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]
lib.ElHermitianTridiagEigPairPartialDist_d.restype = c_uint
lib.ElHermitianTridiagEigPairPartialDist_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]
lib.ElHermitianTridiagEigPairPartialDist_c.restype = c_uint
lib.ElHermitianTridiagEigPairPartialDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]
lib.ElHermitianTridiagEigPairPartialDist_z.restype = c_uint
def HermitianTridiagEig(d,dSub,vectors=False,sort=ASCENDING,subset=None):
  if type(d) is Matrix:
    w = Matrix(d.tag)
    if vectors:
      X = Matrix(dSub.tag)
      if subset == None:
        if   dSub.tag == sTag: 
          lib.ElHermitianTridiagEigPair_s(d.obj,dSub.obj,w.obj,X.obj,sort)
        elif dSub.tag == dTag:
          lib.ElHermitianTridiagEigPair_d(d.obj,dSub.obj,w.obj,X.obj,sort)
        elif dSub.tag == cTag:
          lib.ElHermitianTridiagEigPair_c(d.obj,dSub.obj,w.obj,X.obj,sort)
        elif dSub.tag == zTag:
          lib.ElHermitianTridiagEigPair_z(d.obj,dSub.obj,w.obj,X.obj,sort)
        else: raise Exception('Unsupported datatype')
      else:
        if   dSub.tag == sTag: 
          lib.ElHermitianTridiagEigPairPartial_s \
          (d.obj,dSub.obj,w.obj,X.obj,sort,subset)
        elif dSub.tag == dTag:
          lib.ElHermitianTridiagEigPairPartial_d \
          (d.obj,dSub.obj,w.obj,X.obj,sort,subset)
        elif dSub.tag == cTag:
          lib.ElHermitianTridiagEigPairPartial_c \
          (d.obj,dSub.obj,w.obj,X.obj,sort,subset)
        elif dSub.tag == zTag:
          lib.ElHermitianTridiagEigPairPartial_z \
          (d.obj,dSub.obj,w.obj,X.obj,sort,subset)
        else: raise Exception('Unsupported datatype')
      return w, X
    else:
      if subset == None:
        if   dSub.tag == sTag: 
          lib.ElHermitianTridiagEigPair_s(d.obj,dSub.obj,w.obj,sort)
        elif dSub.tag == dTag:
          lib.ElHermitianTridiagEigPair_d(d.obj,dSub.obj,w.obj,sort)
        elif dSub.tag == cTag:
          lib.ElHermitianTridiagEigPair_c(d.obj,dSub.obj,w.obj,sort)
        elif dSub.tag == zTag:
          lib.ElHermitianTridiagEigPair_z(d.obj,dSub.obj,w.obj,sort)
        else: raise Exception('Unsupported datatype')
      else:
        if   dSub.tag == sTag: 
          lib.ElHermitianTridiagEigPairPartial_s \
          (d.obj,dSub.obj,w.obj,sort,subset)
        elif dSub.tag == dTag:
          lib.ElHermitianTridiagEigPairPartial_d \
          (d.obj,dSub.obj,w.obj,sort,subset)
        elif dSub.tag == cTag:
          lib.ElHermitianTridiagEigPairPartial_c \
          (d.obj,dSub.obj,w.obj,sort,subset)
        elif dSub.tag == zTag:
          lib.ElHermitianTridiagEigPairPartial_z \
          (d.obj,dSub.obj,w.obj,sort,subset)
        else: raise Exception('Unsupported datatype')
      return w
  elif type(d) is DistMatrix:
    w = DistMatrix(d.tag,STAR,STAR,d.Grid())
    if vectors:
      X = DistMatrix(dSub.tag,STAR,VR,dSub.Grid())
      if subset == None:
        if   dSub.tag == sTag: 
          lib.ElHermitianTridiagEigPairDist_s(d.obj,dSub.obj,w.obj,X.obj,sort)
        elif dSub.tag == dTag:
          lib.ElHermitianTridiagEigPairDist_d(d.obj,dSub.obj,w.obj,X.obj,sort)
        elif dSub.tag == cTag:
          lib.ElHermitianTridiagEigPairDist_c(d.obj,dSub.obj,w.obj,X.obj,sort)
        elif dSub.tag == zTag:
          lib.ElHermitianTridiagEigPairDist_z(d.obj,dSub.obj,w.obj,X.obj,sort)
        else: raise Exception('Unsupported datatype')
      else:
        if   dSub.tag == sTag: 
          lib.ElHermitianTridiagEigPairPartialDist_s \
          (d.obj,dSub.obj,w.obj,X.obj,sort,subset)
        elif dSub.tag == dTag:
          lib.ElHermitianTridiagEigPairPartialDist_d \
          (d.obj,dSub.obj,w.obj,X.obj,sort,subset)
        elif dSub.tag == cTag:
          lib.ElHermitianTridiagEigPairPartialDist_c \
          (d.obj,dSub.obj,w.obj,X.obj,sort,subset)
        elif dSub.tag == zTag:
          lib.ElHermitianTridiagEigPairPartialDist_z \
          (d.obj,dSub.obj,w.obj,X.obj,sort,subset)
        else: raise Exception('Unsupported datatype')
      return w, X
    else:
      if subset == None:
        if   dSub.tag == sTag: 
          lib.ElHermitianTridiagEigPairDist_s(d.obj,dSub.obj,w.obj,sort)
        elif dSub.tag == dTag:
          lib.ElHermitianTridiagEigPairDist_d(d.obj,dSub.obj,w.obj,sort)
        elif dSub.tag == cTag:
          lib.ElHermitianTridiagEigPairDist_c(d.obj,dSub.obj,w.obj,sort)
        elif dSub.tag == zTag:
          lib.ElHermitianTridiagEigPairDist_z(d.obj,dSub.obj,w.obj,sort)
        else: raise Exception('Unsupported datatype')
      else:
        if   dSub.tag == sTag: 
          lib.ElHermitianTridiagEigPairPartialDist_s \
          (d.obj,dSub.obj,w.obj,sort,subset)
        elif dSub.tag == dTag:
          lib.ElHermitianTridiagEigPairPartialDist_d \
          (d.obj,dSub.obj,w.obj,sort,subset)
        elif dSub.tag == cTag:
          lib.ElHermitianTridiagEigPairPartialDist_c \
          (d.obj,dSub.obj,w.obj,sort,subset)
        elif dSub.tag == zTag:
          lib.ElHermitianTridiagEigPairPartialDist_z \
          (d.obj,dSub.obj,w.obj,sort,subset)
        else: raise Exception('Unsupported datatype')
      return w
  else: raise Exception('Unsupported matrix type')

# Hermitian eigensolvers
# ======================
lib.ElHermitianEig_s.argtypes = [c_uint,c_void_p,c_void_p,c_uint]
lib.ElHermitianEig_s.restype = c_uint
lib.ElHermitianEig_d.argtypes = [c_uint,c_void_p,c_void_p,c_uint]
lib.ElHermitianEig_d.restype = c_uint
lib.ElHermitianEig_c.argtypes = [c_uint,c_void_p,c_void_p,c_uint]
lib.ElHermitianEig_c.restype = c_uint
lib.ElHermitianEig_z.argtypes = [c_uint,c_void_p,c_void_p,c_uint]
lib.ElHermitianEig_z.restype = c_uint
lib.ElHermitianEigDist_s.argtypes = [c_uint,c_void_p,c_void_p,c_uint]
lib.ElHermitianEigDist_s.restype = c_uint
lib.ElHermitianEigDist_d.argtypes = [c_uint,c_void_p,c_void_p,c_uint]
lib.ElHermitianEigDist_d.restype = c_uint
lib.ElHermitianEigDist_c.argtypes = [c_uint,c_void_p,c_void_p,c_uint]
lib.ElHermitianEigDist_c.restype = c_uint
lib.ElHermitianEigDist_z.argtypes = [c_uint,c_void_p,c_void_p,c_uint]
lib.ElHermitianEigDist_z.restype = c_uint
lib.ElHermitianEigPair_s.argtypes = [c_uint,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianEigPair_s.restype = c_uint
lib.ElHermitianEigPair_d.argtypes = [c_uint,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianEigPair_d.restype = c_uint
lib.ElHermitianEigPair_c.argtypes = [c_uint,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianEigPair_c.restype = c_uint
lib.ElHermitianEigPair_z.argtypes = [c_uint,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianEigPair_z.restype = c_uint
lib.ElHermitianEigPairDist_s.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianEigPairDist_s.restype = c_uint
lib.ElHermitianEigPairDist_d.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianEigPairDist_d.restype = c_uint
lib.ElHermitianEigPairDist_c.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianEigPairDist_c.restype = c_uint
lib.ElHermitianEigPairDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElHermitianEigPairDist_z.restype = c_uint

lib.ElHermitianEigPartial_s.argtypes = \
  [c_uint,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]
lib.ElHermitianEigPartial_s.restype = c_uint
lib.ElHermitianEigPartial_d.argtypes = \
  [c_uint,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]
lib.ElHermitianEigPartial_d.restype = c_uint
lib.ElHermitianEigPartial_c.argtypes = \
  [c_uint,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]
lib.ElHermitianEigPartial_c.restype = c_uint
lib.ElHermitianEigPartial_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]
lib.ElHermitianEigPartial_z.restype = c_uint
lib.ElHermitianEigPartialDist_s.argtypes = \
  [c_uint,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]
lib.ElHermitianEigPartialDist_s.restype = c_uint
lib.ElHermitianEigPartialDist_d.argtypes = \
  [c_uint,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]
lib.ElHermitianEigPartialDist_d.restype = c_uint
lib.ElHermitianEigPartialDist_c.argtypes = \
  [c_uint,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]
lib.ElHermitianEigPartialDist_c.restype = c_uint
lib.ElHermitianEigPartialDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]
lib.ElHermitianEigPartialDist_z.restype = c_uint
lib.ElHermitianEigPairPartial_s.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]
lib.ElHermitianEigPairPartial_s.restype = c_uint
lib.ElHermitianEigPairPartial_d.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]
lib.ElHermitianEigPairPartial_d.restype = c_uint
lib.ElHermitianEigPairPartial_c.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]
lib.ElHermitianEigPairPartial_c.restype = c_uint
lib.ElHermitianEigPairPartial_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]
lib.ElHermitianEigPairPartial_z.restype = c_uint
lib.ElHermitianEigPairPartialDist_s.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]
lib.ElHermitianEigPairPartialDist_s.restype = c_uint
lib.ElHermitianEigPairPartialDist_d.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]
lib.ElHermitianEigPairPartialDist_d.restype = c_uint
lib.ElHermitianEigPairPartialDist_c.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]
lib.ElHermitianEigPairPartialDist_c.restype = c_uint
lib.ElHermitianEigPairPartialDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]
lib.ElHermitianEigPairPartialDist_z.restype = c_uint
def HermitianEig(uplo,A,vectors=False,sort=ASCENDING,subset=None):
  if type(A) is Matrix:
    w = Matrix(A.tag)
    if vectors:
      X = Matrix(Base(A.tag))
      if subset == None:
        if   A.tag == sTag: 
          lib.ElHermitianEigPair_s(uplo,A.obj,w.obj,X.obj,sort)
        elif A.tag == dTag:
          lib.ElHermitianEigPair_d(uplo,A.obj,w.obj,X.obj,sort)
        elif A.tag == cTag:
          lib.ElHermitianEigPair_c(uplo,A.obj,w.obj,X.obj,sort)
        elif A.tag == zTag:
          lib.ElHermitianEigPair_z(uplo,A.obj,w.obj,X.obj,sort)
        else: raise Exception('Unsupported datatype')
      else:
        if   A.tag == sTag: 
          lib.ElHermitianEigPairPartial_s(uplo,A.obj,w.obj,X.obj,sort,subset)
        elif A.tag == dTag:
          lib.ElHermitianEigPairPartial_d(uplo,A.obj,w.obj,X.obj,sort,subset)
        elif A.tag == cTag:
          lib.ElHermitianEigPairPartial_c(uplo,A.obj,w.obj,X.obj,sort,subset)
        elif A.tag == zTag:
          lib.ElHermitianEigPairPartial_z(uplo,A.obj,w.obj,X.obj,sort,subset)
        else: raise Exception('Unsupported datatype')
      return w, X
    else:
      if subset == None:
        if   A.tag == sTag: lib.ElHermitianEigPair_s(uplo,A.obj,w.obj,sort)
        elif A.tag == dTag: lib.ElHermitianEigPair_d(uplo,A.obj,w.obj,sort)
        elif A.tag == cTag: lib.ElHermitianEigPair_c(uplo,A.obj,w.obj,sort)
        elif A.tag == zTag: lib.ElHermitianEigPair_z(uplo,A.obj,w.obj,sort)
        else: raise Exception('Unsupported datatype')
      else:
        if   A.tag == sTag: 
          lib.ElHermitianEigPairPartial_s(uplo,A.obj,w.obj,sort,subset)
        elif A.tag == dTag:
          lib.ElHermitianEigPairPartial_d(uplo,A.obj,w.obj,sort,subset)
        elif A.tag == cTag:
          lib.ElHermitianEigPairPartial_c(uplo,A.obj,w.obj,sort,subset)
        elif A.tag == zTag:
          lib.ElHermitianEigPairPartial_z(uplo,A.obj,w.obj,sort,subset)
        else: raise Exception('Unsupported datatype')
      return w
  elif type(A) is DistMatrix:
    w = DistMatrix(Base(A.tag),STAR,STAR,A.Grid())
    if vectors:
      X = DistMatrix(A.tag,MC,MR,dSub.Grid())
      if subset == None:
        if   A.tag == sTag: 
          lib.ElHermitianEigPairDist_s(uplo,A.obj,w.obj,X.obj,sort)
        elif A.tag == dTag:
          lib.ElHermitianEigPairDist_d(uplo,A.obj,w.obj,X.obj,sort)
        elif A.tag == cTag:
          lib.ElHermitianEigPairDist_c(uplo,A.obj,w.obj,X.obj,sort)
        elif A.tag == zTag:
          lib.ElHermitianEigPairDist_z(uplo,A.obj,w.obj,X.obj,sort)
        else: raise Exception('Unsupported datatype')
      else:
        if   A.tag == sTag: 
          lib.ElHermitianEigPairPartialDist_s \
          (uplo,A.obj,w.obj,X.obj,sort,subset)
        elif A.tag == dTag:
          lib.ElHermitianEigPairPartialDist_d \
          (uplo,A.obj,w.obj,X.obj,sort,subset)
        elif A.tag == cTag:
          lib.ElHermitianEigPairPartialDist_c \
          (uplo,A.obj,w.obj,X.obj,sort,subset)
        elif A.tag == zTag:
          lib.ElHermitianEigPairPartialDist_z \
          (uplo,A.obj,w.obj,X.obj,sort,subset)
        else: raise Exception('Unsupported datatype')
      return w, X
    else:
      if subset == None:
        if   A.tag == sTag: lib.ElHermitianEigPairDist_s(uplo,A.obj,w.obj,sort)
        elif A.tag == dTag: lib.ElHermitianEigPairDist_d(uplo,A.obj,w.obj,sort)
        elif A.tag == cTag: lib.ElHermitianEigPairDist_c(uplo,A.obj,w.obj,sort)
        elif A.tag == zTag: lib.ElHermitianEigPairDist_z(uplo,A.obj,w.obj,sort)
        else: raise Exception('Unsupported datatype')
      else:
        if   A.tag == sTag: 
          lib.ElHermitianEigPairPartialDist_s(uplo,A.obj,w.obj,sort,subset)
        elif A.tag == dTag:
          lib.ElHermitianEigPairPartialDist_d(uplo,A.obj,w.obj,sort,subset)
        elif A.tag == cTag:
          lib.ElHermitianEigPairPartialDist_c(uplo,A.obj,w.obj,sort,subset)
        elif A.tag == zTag:
          lib.ElHermitianEigPairPartialDist_z(uplo,A.obj,w.obj,sort,subset)
        else: raise Exception('Unsupported datatype')
      return w
  else: raise Exception('Unsupported matrix type')
