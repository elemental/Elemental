#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from ..core import *
from ..blas_like import Copy, EntrywiseMap
from ..io import ProcessEvents
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

lib.ElHermitianTridiagEig_s.argtypes = \
lib.ElHermitianTridiagEig_d.argtypes = \
lib.ElHermitianTridiagEig_c.argtypes = \
lib.ElHermitianTridiagEig_z.argtypes = \
lib.ElHermitianTridiagEigDist_s.argtypes = \
lib.ElHermitianTridiagEigDist_d.argtypes = \
lib.ElHermitianTridiagEigDist_c.argtypes = \
lib.ElHermitianTridiagEigDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_uint]

lib.ElHermitianTridiagEig_s.restype = \
lib.ElHermitianTridiagEig_d.restype = \
lib.ElHermitianTridiagEig_c.restype = \
lib.ElHermitianTridiagEig_z.restype = \
lib.ElHermitianTridiagEigDist_s.restype = \
lib.ElHermitianTridiagEigDist_d.restype = \
lib.ElHermitianTridiagEigDist_c.restype = \
lib.ElHermitianTridiagEigDist_z.restype = \
  c_uint

lib.ElHermitianTridiagEigPair_s.argtypes = \
lib.ElHermitianTridiagEigPair_d.argtypes = \
lib.ElHermitianTridiagEigPair_c.argtypes = \
lib.ElHermitianTridiagEigPair_z.argtypes = \
lib.ElHermitianTridiagEigPairDist_s.argtypes = \
lib.ElHermitianTridiagEigPairDist_d.argtypes = \
lib.ElHermitianTridiagEigPairDist_c.argtypes = \
lib.ElHermitianTridiagEigPairDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]

lib.ElHermitianTridiagEigPair_s.restype = \
lib.ElHermitianTridiagEigPair_d.restype = \
lib.ElHermitianTridiagEigPair_c.restype = \
lib.ElHermitianTridiagEigPair_z.restype = \
lib.ElHermitianTridiagEigPairDist_s.restype = \
lib.ElHermitianTridiagEigPairDist_d.restype = \
lib.ElHermitianTridiagEigPairDist_c.restype = \
lib.ElHermitianTridiagEigPairDist_z.restype = \
  c_uint

lib.ElHermitianTridiagEigPartial_s.argtypes = \
lib.ElHermitianTridiagEigPartial_c.argtypes = \
lib.ElHermitianTridiagEigPartialDist_s.argtypes = \
lib.ElHermitianTridiagEigPartialDist_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]

lib.ElHermitianTridiagEigPartial_d.argtypes = \
lib.ElHermitianTridiagEigPartial_z.argtypes = \
lib.ElHermitianTridiagEigPartialDist_d.argtypes = \
lib.ElHermitianTridiagEigPartialDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]

lib.ElHermitianTridiagEigPartial_s.restype = \
lib.ElHermitianTridiagEigPartial_d.restype = \
lib.ElHermitianTridiagEigPartial_c.restype = \
lib.ElHermitianTridiagEigPartial_z.restype = \
lib.ElHermitianTridiagEigPartialDist_s.restype = \
lib.ElHermitianTridiagEigPartialDist_d.restype = \
lib.ElHermitianTridiagEigPartialDist_c.restype = \
lib.ElHermitianTridiagEigPartialDist_z.restype = \
  c_uint

lib.ElHermitianTridiagEigPairPartial_s.argtypes = \
lib.ElHermitianTridiagEigPairPartial_c.argtypes = \
lib.ElHermitianTridiagEigPairPartialDist_s.argtypes = \
lib.ElHermitianTridiagEigPairPartialDist_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]

lib.ElHermitianTridiagEigPairPartial_d.argtypes = \
lib.ElHermitianTridiagEigPairPartial_z.argtypes = \
lib.ElHermitianTridiagEigPairPartialDist_d.argtypes = \
lib.ElHermitianTridiagEigPairPartialDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]

lib.ElHermitianTridiagEigPairPartial_s.restype = \
lib.ElHermitianTridiagEigPairPartial_d.restype = \
lib.ElHermitianTridiagEigPairPartial_c.restype = \
lib.ElHermitianTridiagEigPairPartial_z.restype = \
lib.ElHermitianTridiagEigPairPartialDist_s.restype = \
lib.ElHermitianTridiagEigPairPartialDist_d.restype = \
lib.ElHermitianTridiagEigPairPartialDist_c.restype = \
lib.ElHermitianTridiagEigPairPartialDist_z.restype = \
  c_uint

def HermitianTridiagEig(d,dSub,vectors=False,sort=ASCENDING,subset=None):
  if type(d) is Matrix:
    w = Matrix(d.tag)
    if vectors:
      X = Matrix(dSub.tag)
      if subset == None:
        args = [d.obj,dSub.obj,w.obj,X.obj,sort]
        if   dSub.tag == sTag: lib.ElHermitianTridiagEigPair_s(*args)
        elif dSub.tag == dTag: lib.ElHermitianTridiagEigPair_d(*args)
        elif dSub.tag == cTag: lib.ElHermitianTridiagEigPair_c(*args)
        elif dSub.tag == zTag: lib.ElHermitianTridiagEigPair_z(*args)
        else: DataExcept()
      else:
        args = [d.obj,dSub.obj,w.obj,X.obj,sort,subset]
        if   dSub.tag == sTag: lib.ElHermitianTridiagEigPairPartial_s(*args)
        elif dSub.tag == dTag: lib.ElHermitianTridiagEigPairPartial_d(*args)
        elif dSub.tag == cTag: lib.ElHermitianTridiagEigPairPartial_c(*args)
        elif dSub.tag == zTag: lib.ElHermitianTridiagEigPairPartial_z(*args)
        else: DataExcept()
      return w, X
    else:
      if subset == None:
        args = [d.obj,dSub.obj,w.obj,sort]
        if   dSub.tag == sTag: lib.ElHermitianTridiagEigPair_s(*args)
        elif dSub.tag == dTag: lib.ElHermitianTridiagEigPair_d(*args)
        elif dSub.tag == cTag: lib.ElHermitianTridiagEigPair_c(*args)
        elif dSub.tag == zTag: lib.ElHermitianTridiagEigPair_z(*args)
        else: DataExcept()
      else:
        args = [d.obj,dSub.obj,w.obj,sort,subset]
        if   dSub.tag == sTag: lib.ElHermitianTridiagEigPairPartial_s(*args)
        elif dSub.tag == dTag: lib.ElHermitianTridiagEigPairPartial_d(*args)
        elif dSub.tag == cTag: lib.ElHermitianTridiagEigPairPartial_c(*args)
        elif dSub.tag == zTag: lib.ElHermitianTridiagEigPairPartial_z(*args)
        else: DataExcept()
      return w
  elif type(d) is DistMatrix:
    w = DistMatrix(d.tag,STAR,STAR,d.Grid())
    if vectors:
      X = DistMatrix(dSub.tag,STAR,VR,dSub.Grid())
      if subset == None:
        args = [d.obj,dSub.obj,w.obj,X.obj,sort]
        if   dSub.tag == sTag: lib.ElHermitianTridiagEigPairDist_s(*args)
        elif dSub.tag == dTag: lib.ElHermitianTridiagEigPairDist_d(*args)
        elif dSub.tag == cTag: lib.ElHermitianTridiagEigPairDist_c(*args)
        elif dSub.tag == zTag: lib.ElHermitianTridiagEigPairDist_z(*args)
        else: DataExcept()
      else:
        args = [d.obj,dSub.obj,w.obj,X.obj,sort,subset]
        if   dSub.tag == sTag: lib.ElHermitianTridiagEigPairPartialDist_s(*args)
        elif dSub.tag == dTag: lib.ElHermitianTridiagEigPairPartialDist_d(*args)
        elif dSub.tag == cTag: lib.ElHermitianTridiagEigPairPartialDist_c(*args)
        elif dSub.tag == zTag: lib.ElHermitianTridiagEigPairPartialDist_z(*args)
        else: DataExcept()
      return w, X
    else:
      if subset == None:
        args = [d.obj,dSub.obj,w.obj,sort]
        if   dSub.tag == sTag: lib.ElHermitianTridiagEigPairDist_s(*args)
        elif dSub.tag == dTag: lib.ElHermitianTridiagEigPairDist_d(*args)
        elif dSub.tag == cTag: lib.ElHermitianTridiagEigPairDist_c(*args)
        elif dSub.tag == zTag: lib.ElHermitianTridiagEigPairDist_z(*args)
        else: DataExcept()
      else:
        args = [d.obj,dSub.obj,w.obj,sort,subset]
        if   dSub.tag == sTag: lib.ElHermitianTridiagEigPairPartialDist_s(*args)
        elif dSub.tag == dTag: lib.ElHermitianTridiagEigPairPartialDist_d(*args)
        elif dSub.tag == cTag: lib.ElHermitianTridiagEigPairPartialDist_c(*args)
        elif dSub.tag == zTag: lib.ElHermitianTridiagEigPairPartialDist_z(*args)
        else: DataExcept()
      return w
  else: TypeExcept()

# Hermitian eigensolvers
# ======================
lib.ElHermitianEig_s.argtypes = \
lib.ElHermitianEig_d.argtypes = \
lib.ElHermitianEig_c.argtypes = \
lib.ElHermitianEig_z.argtypes = \
lib.ElHermitianEigDist_s.argtypes = \
lib.ElHermitianEigDist_d.argtypes = \
lib.ElHermitianEigDist_c.argtypes = \
lib.ElHermitianEigDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_uint]

lib.ElHermitianEig_s.restype = \
lib.ElHermitianEig_d.restype = \
lib.ElHermitianEig_c.restype = \
lib.ElHermitianEig_z.restype = \
lib.ElHermitianEigDist_s.restype = \
lib.ElHermitianEigDist_d.restype = \
lib.ElHermitianEigDist_c.restype = \
lib.ElHermitianEigDist_z.restype = \
  c_uint

lib.ElHermitianEigPair_s.argtypes = \
lib.ElHermitianEigPair_d.argtypes = \
lib.ElHermitianEigPair_c.argtypes = \
lib.ElHermitianEigPair_z.argtypes = \
lib.ElHermitianEigPairDist_s.argtypes = \
lib.ElHermitianEigPairDist_d.argtypes = \
lib.ElHermitianEigPairDist_c.argtypes = \
lib.ElHermitianEigPairDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint]

lib.ElHermitianEigPair_s.restype = \
lib.ElHermitianEigPair_d.restype = \
lib.ElHermitianEigPair_c.restype = \
lib.ElHermitianEigPair_z.restype = \
lib.ElHermitianEigPairDist_s.restype = \
lib.ElHermitianEigPairDist_d.restype = \
lib.ElHermitianEigPairDist_c.restype = \
lib.ElHermitianEigPairDist_z.restype = \
  c_uint

lib.ElHermitianEigPartial_s.argtypes = \
lib.ElHermitianEigPartial_c.argtypes = \
lib.ElHermitianEigPartialDist_s.argtypes = \
lib.ElHermitianEigPartialDist_c.argtypes = \
  [c_uint,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]

lib.ElHermitianEigPartial_d.argtypes = \
lib.ElHermitianEigPartial_z.argtypes = \
lib.ElHermitianEigPartialDist_d.argtypes = \
lib.ElHermitianEigPartialDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]

lib.ElHermitianEigPartial_s.restype = \
lib.ElHermitianEigPartial_d.restype = \
lib.ElHermitianEigPartial_c.restype = \
lib.ElHermitianEigPartial_z.restype = \
lib.ElHermitianEigPartialDist_s.restype = \
lib.ElHermitianEigPartialDist_d.restype = \
lib.ElHermitianEigPartialDist_c.restype = \
lib.ElHermitianEigPartialDist_z.restype = \
  c_uint

lib.ElHermitianEigPairPartial_s.argtypes = \
lib.ElHermitianEigPairPartial_c.argtypes = \
lib.ElHermitianEigPairPartialDist_s.argtypes = \
lib.ElHermitianEigPairPartialDist_c.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]

lib.ElHermitianEigPairPartial_d.argtypes = \
lib.ElHermitianEigPairPartial_z.argtypes = \
lib.ElHermitianEigPairPartialDist_d.argtypes = \
lib.ElHermitianEigPairPartialDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]

lib.ElHermitianEigPairPartial_s.restype = \
lib.ElHermitianEigPairPartial_d.restype = \
lib.ElHermitianEigPairPartial_c.restype = \
lib.ElHermitianEigPairPartial_z.restype = \
lib.ElHermitianEigPairPartialDist_s.restype = \
lib.ElHermitianEigPairPartialDist_d.restype = \
lib.ElHermitianEigPairPartialDist_c.restype = \
lib.ElHermitianEigPairPartialDist_z.restype = \
  c_uint

def HermitianEig(uplo,A,vectors=False,sort=ASCENDING,subset=None):
  if type(A) is Matrix:
    w = Matrix(A.tag)
    if vectors:
      X = Matrix(Base(A.tag))
      if subset == None:
        args = [uplo,A.obj,w.obj,X.obj,sort]
        if   A.tag == sTag: lib.ElHermitianEigPair_s(*args)
        elif A.tag == dTag: lib.ElHermitianEigPair_d(*args)
        elif A.tag == cTag: lib.ElHermitianEigPair_c(*args)
        elif A.tag == zTag: lib.ElHermitianEigPair_z(*args)
        else: DataExcept()
      else:
        args = [uplo,A.obj,w.obj,X.obj,sort,subset]
        if   A.tag == sTag: lib.ElHermitianEigPairPartial_s(*args)
        elif A.tag == dTag: lib.ElHermitianEigPairPartial_d(*args)
        elif A.tag == cTag: lib.ElHermitianEigPairPartial_c(*args)
        elif A.tag == zTag: lib.ElHermitianEigPairPartial_z(*args)
        else: DataExcept()
      return w, X
    else:
      if subset == None:
        args = [uplo,A.obj,w.obj,sort]
        if   A.tag == sTag: lib.ElHermitianEigPair_s(*args)
        elif A.tag == dTag: lib.ElHermitianEigPair_d(*args)
        elif A.tag == cTag: lib.ElHermitianEigPair_c(*args)
        elif A.tag == zTag: lib.ElHermitianEigPair_z(*args)
        else: DataExcept()
      else:
        args = [uplo,A.obj,w.obj,sort,subset]
        if   A.tag == sTag: lib.ElHermitianEigPairPartial_s(*args)
        elif A.tag == dTag: lib.ElHermitianEigPairPartial_d(*args)
        elif A.tag == cTag: lib.ElHermitianEigPairPartial_c(*args)
        elif A.tag == zTag: lib.ElHermitianEigPairPartial_z(*args)
        else: DataExcept()
      return w
  elif type(A) is DistMatrix:
    w = DistMatrix(Base(A.tag),STAR,STAR,A.Grid())
    if vectors:
      X = DistMatrix(A.tag,MC,MR,dSub.Grid())
      if subset == None:
        args = [uplo,A.obj,w.obj,X.obj,sort]
        if   A.tag == sTag: lib.ElHermitianEigPairDist_s(*args)
        elif A.tag == dTag: lib.ElHermitianEigPairDist_d(*args)
        elif A.tag == cTag: lib.ElHermitianEigPairDist_c(*args)
        elif A.tag == zTag: lib.ElHermitianEigPairDist_z(*args)
        else: DataExcept()
      else:
        args = [uplo,A.obj,w.obj,X.obj,sort,subset]
        if   A.tag == sTag: lib.ElHermitianEigPairPartialDist_s(*args)
        elif A.tag == dTag: lib.ElHermitianEigPairPartialDist_d(*args)
        elif A.tag == cTag: lib.ElHermitianEigPairPartialDist_c(*args)
        elif A.tag == zTag: lib.ElHermitianEigPairPartialDist_z(*args)
        else: DataExcept()
      return w, X
    else:
      if subset == None:
        args = [uplo,A.obj,w.obj,sort]
        if   A.tag == sTag: lib.ElHermitianEigPairDist_s(*args)
        elif A.tag == dTag: lib.ElHermitianEigPairDist_d(*args)
        elif A.tag == cTag: lib.ElHermitianEigPairDist_c(*args)
        elif A.tag == zTag: lib.ElHermitianEigPairDist_z(*args)
        else: DataExcept()
      else:
        args = [uplo,A.obj,w.obj,sort,subset]
        if   A.tag == sTag: lib.ElHermitianEigPairPartialDist_s(*args)
        elif A.tag == dTag: lib.ElHermitianEigPairPartialDist_d(*args)
        elif A.tag == cTag: lib.ElHermitianEigPairPartialDist_c(*args)
        elif A.tag == zTag: lib.ElHermitianEigPairPartialDist_z(*args)
        else: DataExcept()
      return w
  else: TypeExcept()

# Skew-Hermitian eigensolvers
# ===========================

lib.ElSkewHermitianEig_s.argtypes = \
lib.ElSkewHermitianEig_d.argtypes = \
lib.ElSkewHermitianEig_c.argtypes = \
lib.ElSkewHermitianEig_z.argtypes = \
lib.ElSkewHermitianEigDist_s.argtypes = \
lib.ElSkewHermitianEigDist_d.argtypes = \
lib.ElSkewHermitianEigDist_c.argtypes = \
lib.ElSkewHermitianEigDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_uint]

lib.ElSkewHermitianEig_s.restype = \
lib.ElSkewHermitianEig_d.restype = \
lib.ElSkewHermitianEig_c.restype = \
lib.ElSkewHermitianEig_z.restype = \
lib.ElSkewHermitianEigDist_s.restype = \
lib.ElSkewHermitianEigDist_d.restype = \
lib.ElSkewHermitianEigDist_c.restype = \
lib.ElSkewHermitianEigDist_z.restype = \
  c_uint

lib.ElSkewHermitianEigPair_s.argtypes = \
lib.ElSkewHermitianEigPair_d.argtypes = \
lib.ElSkewHermitianEigPair_c.argtypes = \
lib.ElSkewHermitianEigPair_z.argtypes = \
lib.ElSkewHermitianEigPairDist_s.argtypes = \
lib.ElSkewHermitianEigPairDist_d.argtypes = \
lib.ElSkewHermitianEigPairDist_c.argtypes = \
lib.ElSkewHermitianEigPairDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint]

lib.ElSkewHermitianEigPair_s.restype = \
lib.ElSkewHermitianEigPair_d.restype = \
lib.ElSkewHermitianEigPair_c.restype = \
lib.ElSkewHermitianEigPair_z.restype = \
lib.ElSkewHermitianEigPairDist_s.restype = \
lib.ElSkewHermitianEigPairDist_d.restype = \
lib.ElSkewHermitianEigPairDist_c.restype = \
lib.ElSkewHermitianEigPairDist_z.restype = \
  c_uint

lib.ElSkewHermitianEigPartial_s.argtypes = \
lib.ElSkewHermitianEigPartial_c.argtypes = \
lib.ElSkewHermitianEigPartialDist_s.argtypes = \
lib.ElSkewHermitianEigPartialDist_c.argtypes = \
  [c_uint,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]

lib.ElSkewHermitianEigPartial_d.argtypes = \
lib.ElSkewHermitianEigPartial_z.argtypes = \
lib.ElSkewHermitianEigPartialDist_d.argtypes = \
lib.ElSkewHermitianEigPartialDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]

lib.ElSkewHermitianEigPartial_s.restype = \
lib.ElSkewHermitianEigPartial_d.restype = \
lib.ElSkewHermitianEigPartial_c.restype = \
lib.ElSkewHermitianEigPartial_z.restype = \
lib.ElSkewHermitianEigPartialDist_s.restype = \
lib.ElSkewHermitianEigPartialDist_d.restype = \
lib.ElSkewHermitianEigPartialDist_c.restype = \
lib.ElSkewHermitianEigPartialDist_z.restype = \
  c_uint

lib.ElSkewHermitianEigPairPartial_s.argtypes = \
lib.ElSkewHermitianEigPairPartial_c.argtypes = \
lib.ElSkewHermitianEigPairPartialDist_s.argtypes = \
lib.ElSkewHermitianEigPairPartialDist_c.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]

lib.ElSkewHermitianEigPairPartial_d.argtypes = \
lib.ElSkewHermitianEigPairPartial_z.argtypes = \
lib.ElSkewHermitianEigPairPartialDist_d.argtypes = \
lib.ElSkewHermitianEigPairPartialDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]

lib.ElSkewHermitianEigPairPartial_s.restype = \
lib.ElSkewHermitianEigPairPartial_d.restype = \
lib.ElSkewHermitianEigPairPartial_c.restype = \
lib.ElSkewHermitianEigPairPartial_z.restype = \
lib.ElSkewHermitianEigPairPartialDist_s.restype = \
lib.ElSkewHermitianEigPairPartialDist_d.restype = \
lib.ElSkewHermitianEigPairPartialDist_c.restype = \
lib.ElSkewHermitianEigPairPartialDist_z.restype = \
  c_uint

def SkewHermitianEig(uplo,A,vectors=False,sort=ASCENDING,subset=None):
  if type(A) is Matrix:
    w = Matrix(A.tag)
    if vectors:
      X = Matrix(Base(A.tag))
      if subset == None:
        args = [uplo,A.obj,w.obj,X.obj,sort]
        if   A.tag == sTag: lib.ElSkewHermitianEigPair_s(*args)
        elif A.tag == dTag: lib.ElSkewHermitianEigPair_d(*args)
        elif A.tag == cTag: lib.ElSkewHermitianEigPair_c(*args)
        elif A.tag == zTag: lib.ElSkewHermitianEigPair_z(*args)
        else: DataExcept()
      else:
        args = [uplo,A.obj,w.obj,X.obj,sort,subset]
        if   A.tag == sTag: lib.ElSkewHermitianEigPairPartial_s(*args)
        elif A.tag == dTag: lib.ElSkewHermitianEigPairPartial_d(*args)
        elif A.tag == cTag: lib.ElSkewHermitianEigPairPartial_c(*args)
        elif A.tag == zTag: lib.ElSkewHermitianEigPairPartial_z(*args)
        else: DataExcept()
      return w, X
    else:
      if subset == None:
        args = [uplo,A.obj,w.obj,sort]
        if   A.tag == sTag: lib.ElSkewHermitianEigPair_s(*args)
        elif A.tag == dTag: lib.ElSkewHermitianEigPair_d(*args)
        elif A.tag == cTag: lib.ElSkewHermitianEigPair_c(*args)
        elif A.tag == zTag: lib.ElSkewHermitianEigPair_z(*args)
        else: DataExcept()
      else:
        args = [uplo,A.obj,w.obj,sort,subset]
        if   A.tag == sTag: lib.ElSkewHermitianEigPairPartial_s(*args)
        elif A.tag == dTag: lib.ElSkewHermitianEigPairPartial_d(*args)
        elif A.tag == cTag: lib.ElSkewHermitianEigPairPartial_c(*args)
        elif A.tag == zTag: lib.ElSkewHermitianEigPairPartial_z(*args)
        else: DataExcept()
      return w
  elif type(A) is DistMatrix:
    w = DistMatrix(Base(A.tag),STAR,STAR,A.Grid())
    if vectors:
      X = DistMatrix(A.tag,MC,MR,dSub.Grid())
      if subset == None:
        args = [uplo,A.obj,w.obj,X.obj,sort]
        if   A.tag == sTag: lib.ElSkewHermitianEigPairDist_s(*args)
        elif A.tag == dTag: lib.ElSkewHermitianEigPairDist_d(*args)
        elif A.tag == cTag: lib.ElSkewHermitianEigPairDist_c(*args)
        elif A.tag == zTag: lib.ElSkewHermitianEigPairDist_z(*args)
        else: DataExcept()
      else:
        args = [uplo,A.obj,w.obj,X.obj,sort,subset]
        if   A.tag == sTag: lib.ElSkewHermitianEigPairPartialDist_s(*args)
        elif A.tag == dTag: lib.ElSkewHermitianEigPairPartialDist_d(*args)
        elif A.tag == cTag: lib.ElSkewHermitianEigPairPartialDist_c(*args)
        elif A.tag == zTag: lib.ElSkewHermitianEigPairPartialDist_z(*args)
        else: DataExcept()
      return w, X
    else:
      if subset == None:
        args = [uplo,A.obj,w.obj,sort]
        if   A.tag == sTag: lib.ElSkewHermitianEigPairDist_s(*args)
        elif A.tag == dTag: lib.ElSkewHermitianEigPairDist_d(*args)
        elif A.tag == cTag: lib.ElSkewHermitianEigPairDist_c(*args)
        elif A.tag == zTag: lib.ElSkewHermitianEigPairDist_z(*args)
        else: DataExcept()
      else:
        args = [uplo,A.obj,w.obj,sort,subset]
        if   A.tag == sTag: lib.ElSkewHermitianEigPairPartialDist_s(*args)
        elif A.tag == dTag: lib.ElSkewHermitianEigPairPartialDist_d(*args)
        elif A.tag == cTag: lib.ElSkewHermitianEigPairPartialDist_c(*args)
        elif A.tag == zTag: lib.ElSkewHermitianEigPairPartialDist_z(*args)
        else: DataExcept()
      return w
  else: TypeExcept()

# Hermitian generalized-definite eigensolvers
# ===========================================
lib.ElHermitianGenDefEig_s.argtypes = \
lib.ElHermitianGenDefEig_d.argtypes = \
lib.ElHermitianGenDefEig_c.argtypes = \
lib.ElHermitianGenDefEig_z.argtypes = \
lib.ElHermitianGenDefEigDist_s.argtypes = \
lib.ElHermitianGenDefEigDist_d.argtypes = \
lib.ElHermitianGenDefEigDist_c.argtypes = \
lib.ElHermitianGenDefEigDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_uint]

lib.ElHermitianGenDefEig_s.restype = \
lib.ElHermitianGenDefEig_d.restype = \
lib.ElHermitianGenDefEig_c.restype = \
lib.ElHermitianGenDefEig_z.restype = \
lib.ElHermitianGenDefEigDist_s.restype = \
lib.ElHermitianGenDefEigDist_d.restype = \
lib.ElHermitianGenDefEigDist_c.restype = \
lib.ElHermitianGenDefEigDist_z.restype = \
  c_uint

lib.ElHermitianGenDefEigPair_s.argtypes = \
lib.ElHermitianGenDefEigPair_d.argtypes = \
lib.ElHermitianGenDefEigPair_c.argtypes = \
lib.ElHermitianGenDefEigPair_z.argtypes = \
lib.ElHermitianGenDefEigPairDist_s.argtypes = \
lib.ElHermitianGenDefEigPairDist_d.argtypes = \
lib.ElHermitianGenDefEigPairDist_c.argtypes = \
lib.ElHermitianGenDefEigPairDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_uint]

lib.ElHermitianGenDefEigPair_s.restype = \
lib.ElHermitianGenDefEigPair_d.restype = \
lib.ElHermitianGenDefEigPair_c.restype = \
lib.ElHermitianGenDefEigPair_z.restype = \
lib.ElHermitianGenDefEigPairDist_s.restype = \
lib.ElHermitianGenDefEigPairDist_d.restype = \
lib.ElHermitianGenDefEigPairDist_c.restype = \
lib.ElHermitianGenDefEigPairDist_z.restype = \
  c_uint

lib.ElHermitianGenDefEigPartial_s.argtypes = \
lib.ElHermitianGenDefEigPartial_c.argtypes = \
lib.ElHermitianGenDefEigPartialDist_s.argtypes = \
lib.ElHermitianGenDefEigPartialDist_c.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_s]

lib.ElHermitianGenDefEigPartial_d.argtypes = \
lib.ElHermitianGenDefEigPartial_z.argtypes = \
lib.ElHermitianGenDefEigPartialDist_d.argtypes = \
lib.ElHermitianGenDefEigPartialDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_uint,HermitianEigSubset_d]

lib.ElHermitianGenDefEigPartial_s.restype = \
lib.ElHermitianGenDefEigPartial_d.restype = \
lib.ElHermitianGenDefEigPartial_c.restype = \
lib.ElHermitianGenDefEigPartial_z.restype = \
lib.ElHermitianGenDefEigPartialDist_s.restype = \
lib.ElHermitianGenDefEigPartialDist_d.restype = \
lib.ElHermitianGenDefEigPartialDist_c.restype = \
lib.ElHermitianGenDefEigPartialDist_z.restype = \
  c_uint

lib.ElHermitianGenDefEigPairPartial_s.argtypes = \
lib.ElHermitianGenDefEigPairPartial_c.argtypes = \
lib.ElHermitianGenDefEigPairPartialDist_s.argtypes = \
lib.ElHermitianGenDefEigPairPartialDist_c.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_uint,
   HermitianEigSubset_s]

lib.ElHermitianGenDefEigPairPartial_d.argtypes = \
lib.ElHermitianGenDefEigPairPartial_z.argtypes = \
lib.ElHermitianGenDefEigPairPartialDist_d.argtypes = \
lib.ElHermitianGenDefEigPairPartialDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_uint,
   HermitianEigSubset_d]

lib.ElHermitianGenDefEigPairPartial_s.restype = \
lib.ElHermitianGenDefEigPairPartial_d.restype = \
lib.ElHermitianGenDefEigPairPartial_c.restype = \
lib.ElHermitianGenDefEigPairPartial_z.restype = \
lib.ElHermitianGenDefEigPairPartialDist_s.restype = \
lib.ElHermitianGenDefEigPairPartialDist_d.restype = \
lib.ElHermitianGenDefEigPairPartialDist_c.restype = \
lib.ElHermitianGenDefEigPairPartialDist_z.restype = \
  c_uint

def HermitianGenDefEig(uplo,A,vectors=False,sort=ASCENDING,subset=None):
  if type(A) is Matrix:
    w = Matrix(A.tag)
    if vectors:
      X = Matrix(Base(A.tag))
      if subset == None:
        args = [pencil,uplo,A.obj,B.obj,w.obj,X.obj,sort]
        if   A.tag == sTag: lib.ElHermitianGenDefEigPair_s(*args)
        elif A.tag == dTag: lib.ElHermitianGenDefEigPair_d(*args)
        elif A.tag == cTag: lib.ElHermitianGenDefEigPair_c(*args)
        elif A.tag == zTag: lib.ElHermitianGenDefEigPair_z(*args)
        else: DataExcept()
      else:
        args = [pencil,uplo,A.obj,B.obj,w.obj,X.obj,sort,subset]
        if   A.tag == sTag: lib.ElHermitianGenDefEigPairPartial_s(*args)
        elif A.tag == dTag: lib.ElHermitianGenDefEigPairPartial_d(*args)
        elif A.tag == cTag: lib.ElHermitianGenDefEigPairPartial_c(*args)
        elif A.tag == zTag: lib.ElHermitianGenDefEigPairPartial_z(*args)
        else: DataExcept()
      return w, X
    else:
      if subset == None:
        args = [pencil,uplo,A.obj,B.obj,w.obj,sort]
        if   A.tag == sTag: lib.ElHermitianGenDefEigPair_s(*args)
        elif A.tag == dTag: lib.ElHermitianGenDefEigPair_d(*args)
        elif A.tag == cTag: lib.ElHermitianGenDefEigPair_c(*args)
        elif A.tag == zTag: lib.ElHermitianGenDefEigPair_z(*args)
        else: DataExcept()
      else:
        args = [pencil,uplo,A.obj,B.obj,w.obj,sort,subset]
        if   A.tag == sTag: lib.ElHermitianGenDefEigPairPartial_s(*args)
        elif A.tag == dTag: lib.ElHermitianGenDefEigPairPartial_d(*args)
        elif A.tag == cTag: lib.ElHermitianGenDefEigPairPartial_c(*args)
        elif A.tag == zTag: lib.ElHermitianGenDefEigPairPartial_z(*args)
        else: DataExcept()
      return w
  elif type(A) is DistMatrix:
    w = DistMatrix(Base(A.tag),STAR,STAR,A.Grid())
    if vectors:
      X = DistMatrix(A.tag,MC,MR,dSub.Grid())
      if subset == None:
        args = [pencil,uplo,A.obj,B.obj,w.obj,X.obj,sort]
        if   A.tag == sTag: lib.ElHermitianGenDefEigPairDist_s(*args)
        elif A.tag == dTag: lib.ElHermitianGenDefEigPairDist_d(*args)
        elif A.tag == cTag: lib.ElHermitianGenDefEigPairDist_c(*args)
        elif A.tag == zTag: lib.ElHermitianGenDefEigPairDist_z(*args)
        else: DataExcept()
      else:
        args = [pencil,uplo,A.obj,B.obj,w.obj,X.obj,sort,subset]
        if   A.tag == sTag: lib.ElHermitianGenDefEigPairPartialDist_s(*args)
        elif A.tag == dTag: lib.ElHermitianGenDefEigPairPartialDist_d(*args)
        elif A.tag == cTag: lib.ElHermitianGenDefEigPairPartialDist_c(*args)
        elif A.tag == zTag: lib.ElHermitianGenDefEigPairPartialDist_z(*args)
        else: DataExcept()
      return w, X
    else:
      if subset == None:
        args = [pencil,uplo,A.obj,B.obj,w.obj,sort]
        if   A.tag == sTag: lib.ElHermitianGenDefEigPairDist_s(*args)
        elif A.tag == dTag: lib.ElHermitianGenDefEigPairDist_d(*args)
        elif A.tag == cTag: lib.ElHermitianGenDefEigPairDist_c(*args)
        elif A.tag == zTag: lib.ElHermitianGenDefEigPairDist_z(*args)
        else: DataExcept()
      else:
        args = [pencil,uplo,A.obj,B.obj,w.obj,sort,subset]
        if   A.tag == sTag: lib.ElHermitianGenDefEigPairPartialDist_s(*args)
        elif A.tag == dTag: lib.ElHermitianGenDefEigPairPartialDist_d(*args)
        elif A.tag == cTag: lib.ElHermitianGenDefEigPairPartialDist_c(*args)
        elif A.tag == zTag: lib.ElHermitianGenDefEigPairPartialDist_z(*args)
        else: DataExcept()
      return w
  else: TypeExcept()

# Hermitian SVD
# =============
lib.ElHermitianSingularValues_s.argtypes = \
lib.ElHermitianSingularValues_d.argtypes = \
lib.ElHermitianSingularValues_c.argtypes = \
lib.ElHermitianSingularValues_z.argtypes = \
lib.ElHermitianSingularValuesDist_s.argtypes = \
lib.ElHermitianSingularValuesDist_d.argtypes = \
lib.ElHermitianSingularValuesDist_c.argtypes = \
lib.ElHermitianSingularValuesDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p]

lib.ElHermitianSingularValues_s.restype = \
lib.ElHermitianSingularValues_d.restype = \
lib.ElHermitianSingularValues_c.restype = \
lib.ElHermitianSingularValues_z.restype = \
lib.ElHermitianSingularValuesDist_s.restype = \
lib.ElHermitianSingularValuesDist_d.restype = \
lib.ElHermitianSingularValuesDist_c.restype = \
lib.ElHermitianSingularValuesDist_z.restype = \
  c_uint

lib.ElHermitianSVD_s.argtypes = \
lib.ElHermitianSVD_d.argtypes = \
lib.ElHermitianSVD_c.argtypes = \
lib.ElHermitianSVD_z.argtypes = \
lib.ElHermitianSVDDist_s.argtypes = \
lib.ElHermitianSVDDist_d.argtypes = \
lib.ElHermitianSVDDist_c.argtypes = \
lib.ElHermitianSVDDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]

lib.ElHermitianSVD_s.restype = \
lib.ElHermitianSVD_d.restype = \
lib.ElHermitianSVD_c.restype = \
lib.ElHermitianSVD_z.restype = \
lib.ElHermitianSVDDist_s.restype = \
lib.ElHermitianSVDDist_d.restype = \
lib.ElHermitianSVDDist_c.restype = \
lib.ElHermitianSVDDist_z.restype = \
  c_uint

def HermitianSVD(uplo,A,vectors=True):
  if type(A) is Matrix:
    s = Matrix(Base(A.tag))
    if vectors:
      U = Matrix(A.tag)
      V = Matrix(A.tag)
      args = [uplo,A.obj,s.obj,U.obj,V.obj]
      if   A.tag == sTag: lib.ElHermitianSVD_s(*args)
      elif A.tag == dTag: lib.ElHermitianSVD_d(*args)
      elif A.tag == cTag: lib.ElHermitianSVD_c(*args)
      elif A.tag == zTag: lib.ElHermitianSVD_z(*args)
      else: DataExcept()
      return U, s, V 
    else:
      args = [uplo,A.obj,s.obj]
      if   A.tag == sTag: lib.ElHermitianSingularValues_s(*args)
      elif A.tag == dTag: lib.ElHermitianSingularValues_d(*args)
      elif A.tag == cTag: lib.ElHermitianSingularValues_c(*args)
      elif A.tag == zTag: lib.ElHermitianSingularValues_z(*args)
      else: DataExcept()
    return s
  elif type(A) is DistMatrix:
    s = DistMatrix(Base(A.tag),STAR,STAR,A.Grid())
    if vectors:
      U = DistMatrix(A.tag,MC,MR,A.Grid())
      V = DistMatrix(A.tag,MC,MR,A.Grid())
      args = [uplo,A.obj,s.obj,U.obj,V.obj]
      if   A.tag == sTag: lib.ElHermitianSVDDist_s(*args)
      elif A.tag == dTag: lib.ElHermitianSVDDist_d(*args)
      elif A.tag == cTag: lib.ElHermitianSVDDist_c(*args)
      elif A.tag == zTag: lib.ElHermitianSVDDist_z(*args)
      else: DataExcept()
      return U, s, V 
    else:
      args = [uplo,A.obj,s.obj]
      if   A.tag == sTag: lib.ElHermitianSingularValuesDist_s(*args)
      elif A.tag == dTag: lib.ElHermitianSingularValuesDist_d(*args)
      elif A.tag == cTag: lib.ElHermitianSingularValuesDist_c(*args)
      elif A.tag == zTag: lib.ElHermitianSingularValuesDist_z(*args)
      else: DataExcept()
    return s
  else: TypeExcept()

# Polar decomposition
# ===================

lib.ElPolar_s.argtypes = \
lib.ElPolar_d.argtypes = \
lib.ElPolar_c.argtypes = \
lib.ElPolar_z.argtypes = \
lib.ElPolarDist_s.argtypes = \
lib.ElPolarDist_d.argtypes = \
lib.ElPolarDist_c.argtypes = \
lib.ElPolarDist_z.argtypes = \
  [c_void_p]

lib.ElPolar_s.restype = \
lib.ElPolar_d.restype = \
lib.ElPolar_c.restype = \
lib.ElPolar_z.restype = \
lib.ElPolarDist_s.restype = \
lib.ElPolarDist_d.restype = \
lib.ElPolarDist_c.restype = \
lib.ElPolarDist_z.restype = \
  c_uint

lib.ElPolarDecomp_s.argtypes = \
lib.ElPolarDecomp_d.argtypes = \
lib.ElPolarDecomp_c.argtypes = \
lib.ElPolarDecomp_z.argtypes = \
lib.ElPolarDecompDist_s.argtypes = \
lib.ElPolarDecompDist_d.argtypes = \
lib.ElPolarDecompDist_c.argtypes = \
lib.ElPolarDecompDist_z.argtypes = \
  [c_void_p,c_void_p]

lib.ElPolarDecomp_s.restype = \
lib.ElPolarDecomp_d.restype = \
lib.ElPolarDecomp_c.restype = \
lib.ElPolarDecomp_z.restype = \
lib.ElPolarDecompDist_s.restype = \
lib.ElPolarDecompDist_d.restype = \
lib.ElPolarDecompDist_c.restype = \
lib.ElPolarDecompDist_z.restype = \
  c_uint

def Polar(A,fullDecomp=False):
  if type(A) is Matrix:
    if fullDecomp:
      P = Matrix(A.tag)
      args = [A.obj,P.obj]
      if   A.tag == sTag: lib.ElPolarDecomp_s(*args)
      elif A.tag == dTag: lib.ElPolarDecomp_d(*args)
      elif A.tag == cTag: lib.ElPolarDecomp_c(*args)
      elif A.tag == zTag: lib.ElPolarDecomp_z(*args)
      else: DataExcept()
      return P
    else:
      args = [A.obj]
      if   A.tag == sTag: lib.ElPolar_s(*args)
      elif A.tag == dTag: lib.ElPolar_d(*args)
      elif A.tag == cTag: lib.ElPolar_c(*args)
      elif A.tag == zTag: lib.ElPolar_z(*args)
      else: DataExcept()
  elif type(A) is DistMatrix:
    if fullDecomp:
      P = DistMatrix(A.tag)
      args = [A.obj,P.obj]
      if   A.tag == sTag: lib.ElPolarDecompDist_s(*args)
      elif A.tag == dTag: lib.ElPolarDecompDist_d(*args)
      elif A.tag == cTag: lib.ElPolarDecompDist_c(*args)
      elif A.tag == zTag: lib.ElPolarDecompDist_z(*args)
      else: DataExcept()
      return P
    else:
      args = [A.obj]
      if   A.tag == sTag: lib.ElPolarDist_s(*args)
      elif A.tag == dTag: lib.ElPolarDist_d(*args)
      elif A.tag == cTag: lib.ElPolarDist_c(*args)
      elif A.tag == zTag: lib.ElPolarDist_z(*args)
      else: DataExcept()
  else: TypeExcept()

lib.ElHermitianPolar_s.argtypes = \
lib.ElHermitianPolar_d.argtypes = \
lib.ElHermitianPolar_c.argtypes = \
lib.ElHermitianPolar_z.argtypes = \
lib.ElHermitianPolarDist_s.argtypes = \
lib.ElHermitianPolarDist_d.argtypes = \
lib.ElHermitianPolarDist_c.argtypes = \
lib.ElHermitianPolarDist_z.argtypes = \
  [c_uint,c_void_p]

lib.ElHermitianPolar_s.restype = \
lib.ElHermitianPolar_d.restype = \
lib.ElHermitianPolar_c.restype = \
lib.ElHermitianPolar_z.restype = \
lib.ElHermitianPolarDist_s.restype = \
lib.ElHermitianPolarDist_d.restype = \
lib.ElHermitianPolarDist_c.restype = \
lib.ElHermitianPolarDist_z.restype = \
  c_uint

lib.ElHermitianPolarDecomp_s.argtypes = \
lib.ElHermitianPolarDecomp_d.argtypes = \
lib.ElHermitianPolarDecomp_c.argtypes = \
lib.ElHermitianPolarDecomp_z.argtypes = \
lib.ElHermitianPolarDecompDist_s.argtypes = \
lib.ElHermitianPolarDecompDist_d.argtypes = \
lib.ElHermitianPolarDecompDist_c.argtypes = \
lib.ElHermitianPolarDecompDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p]

lib.ElHermitianPolarDecomp_s.restype = \
lib.ElHermitianPolarDecomp_d.restype = \
lib.ElHermitianPolarDecomp_c.restype = \
lib.ElHermitianPolarDecomp_z.restype = \
lib.ElHermitianPolarDecompDist_s.restype = \
lib.ElHermitianPolarDecompDist_d.restype = \
lib.ElHermitianPolarDecompDist_c.restype = \
lib.ElHermitianPolarDecompDist_z.restype = \
  c_uint

def HermitianPolar(uplo,A,fullDecomp=False):
  if type(A) is Matrix:
    if fullDecomp:
      P = Matrix(A.tag)
      args = [uplo,A.obj,P.obj]
      if   A.tag == sTag: lib.ElHermitianPolarDecomp_s(*args)
      elif A.tag == dTag: lib.ElHermitianPolarDecomp_d(*args)
      elif A.tag == cTag: lib.ElHermitianPolarDecomp_c(*args)
      elif A.tag == zTag: lib.ElHermitianPolarDecomp_z(*args)
      else: DataExcept()
      return P
    else:
      args = [uplo,A.obj]
      if   A.tag == sTag: lib.ElHermitianPolar_s(*args)
      elif A.tag == dTag: lib.ElHermitianPolar_d(*args)
      elif A.tag == cTag: lib.ElHermitianPolar_c(*args)
      elif A.tag == zTag: lib.ElHermitianPolar_z(*args)
      else: DataExcept()
  elif type(A) is DistMatrix:
    if fullDecomp:
      P = Matrix(A.tag)
      args = [uplo,A.obj,P.obj]
      if   A.tag == sTag: lib.ElHermitianPolarDecompDist_s(*args)
      elif A.tag == dTag: lib.ElHermitianPolarDecompDist_d(*args)
      elif A.tag == cTag: lib.ElHermitianPolarDecompDist_c(*args)
      elif A.tag == zTag: lib.ElHermitianPolarDecompDist_z(*args)
      else: DataExcept()
      return P
    else:
      args = [uplo,A.obj]
      if   A.tag == sTag: lib.ElHermitianPolarDist_s(*args)
      elif A.tag == dTag: lib.ElHermitianPolarDist_d(*args)
      elif A.tag == cTag: lib.ElHermitianPolarDist_c(*args)
      elif A.tag == zTag: lib.ElHermitianPolarDist_z(*args)
      else: DataExcept()
  else: TypeExcept()

# Schur decomposition
# ===================
# Emulate an enum for the sign scaling
(SIGN_SCALE_NONE,SIGN_SCALE_DET,SIGN_SCALE_FROB)=(0,1,2)

lib.ElSignCtrlDefault_s.argtypes = [c_void_p]
lib.ElSignCtrlDefault_s.restype = c_uint
class SignCtrl_s(ctypes.Structure):
  _fields_ = [("maxIts",iType),
              ("tol",sType),
              ("power",sType),
              ("scaling",c_uint)]
  def __init__(self):
    lib.ElSignCtrlDefault_s(pointer(self))

lib.ElSignCtrlDefault_d.argtypes = [c_void_p]
lib.ElSignCtrlDefault_d.restype = c_uint
class SignCtrl_d(ctypes.Structure):
  _fields_ = [("maxIts",iType),
              ("tol",dType),
              ("power",dType),
              ("scaling",c_uint)]
  def __init__(self):
    lib.ElSignCtrlDefault_d(pointer(self))

lib.ElHessQRCtrlDefault.argtypes = [c_void_p]
lib.ElHessQRCtrlDefault.restype = c_uint
class HessQRCtrl(ctypes.Structure):
  _fields_ = [("distAED",bType),
              ("blockHeight",iType),("blockWidth",iType)]
  def __init__(self):
    lib.ElHessQRCtrlDefault(pointer(self))

lib.ElSDCCtrlDefault_s.argtypes = [c_void_p]
lib.ElSDCCtrlDefault_s.restype = c_uint
class SDCCtrl_s(ctypes.Structure):
  _fields_ = [("cutoff",iType),
              ("maxInnerIts",iType),("maxOuterIts",iType),
              ("tol",sType),
              ("spreadFactor",sType),
              ("random",bType),
              ("progress",bType),
              ("signCtrl",SignCtrl_s)]
  def __init__(self):
    lib.ElSDCCtrlDefault_s(pointer(self))

lib.ElSDCCtrlDefault_d.argtypes = [c_void_p]
lib.ElSDCCtrlDefault_d.restype = c_uint
class SDCCtrl_d(ctypes.Structure):
  _fields_ = [("cutoff",iType),
              ("maxInnerIts",iType),("maxOuterIts",iType),
              ("tol",dType),
              ("spreadFactor",dType),
              ("random",bType),
              ("progress",bType),
              ("signCtrl",SignCtrl_d)]
  def __init__(self):
    lib.ElSDCCtrlDefault_d(pointer(self))

lib.ElSchurCtrlDefault_s.argtypes = [c_void_p]
lib.ElSchurCtrlDefault_s.restype = c_uint
class SchurCtrl_s(ctypes.Structure):
  _fields_ = [("useSDC",bType),
              ("qrCtrl",HessQRCtrl),
              ("sdcCtrl",SDCCtrl_s)]
  def __init__(self):
    lib.ElSchurCtrlDefault_s(pointer(self))

lib.ElSchurCtrlDefault_d.argtypes = [c_void_p]
lib.ElSchurCtrlDefault_d.restype = c_uint
class SchurCtrl_d(ctypes.Structure):
  _fields_ = [("useSDC",bType),
              ("qrCtrl",HessQRCtrl),
              ("sdcCtrl",SDCCtrl_d)]
  def __init__(self):
    lib.ElSchurCtrlDefault_d(pointer(self))

lib.ElSchur_s.argtypes = \
lib.ElSchur_d.argtypes = \
lib.ElSchur_c.argtypes = \
lib.ElSchur_z.argtypes = \
lib.ElSchurDist_s.argtypes = \
lib.ElSchurDist_d.argtypes = \
lib.ElSchurDist_c.argtypes = \
lib.ElSchurDist_z.argtypes = \
  [c_void_p,c_void_p,bType]

lib.ElSchur_s.restype = \
lib.ElSchur_d.restype = \
lib.ElSchur_c.restype = \
lib.ElSchur_z.restype = \
lib.ElSchurDist_s.restype = \
lib.ElSchurDist_d.restype = \
lib.ElSchurDist_c.restype = \
lib.ElSchurDist_z.restype = \
  c_uint

lib.ElSchurDecomp_s.argtypes = \
lib.ElSchurDecomp_d.argtypes = \
lib.ElSchurDecomp_c.argtypes = \
lib.ElSchurDecomp_z.argtypes = \
lib.ElSchurDecompDist_s.argtypes = \
lib.ElSchurDecompDist_d.argtypes = \
lib.ElSchurDecompDist_c.argtypes = \
lib.ElSchurDecompDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,bType]

lib.ElSchurDecomp_s.restype = \
lib.ElSchurDecomp_d.restype = \
lib.ElSchurDecomp_c.restype = \
lib.ElSchurDecomp_z.restype = \
lib.ElSchurDecompDist_s.restype = \
lib.ElSchurDecompDist_d.restype = \
lib.ElSchurDecompDist_c.restype = \
lib.ElSchurDecompDist_z.restype = \
  c_uint

def Schur(A,fullTriangle=True,vectors=False):
  if type(A) is Matrix:
    w = Matrix(Complexify(A.tag))
    if vectors:
      Q = Matrix(A.tag)
      args = [A.obj,w.obj,Q.obj,fullTriangle]
      if   A.tag == sTag: lib.ElSchurDecomp_s(*args)
      elif A.tag == dTag: lib.ElSchurDecomp_d(*args)
      elif A.tag == cTag: lib.ElSchurDecomp_c(*args)
      elif A.tag == zTag: lib.ElSchurDecomp_z(*args)
      else: DataExcept()
      return w, Q
    else:
      args = [A.obj,w.obj,fullTriangle]
      if   A.tag == sTag: lib.ElSchur_s(*args)
      elif A.tag == dTag: lib.ElSchur_d(*args)
      elif A.tag == cTag: lib.ElSchur_c(*args)
      elif A.tag == zTag: lib.ElSchur_z(*args)
      else: DataExcept()
      return w
  elif type(A) is DistMatrix:
    w = DistMatrix(Complexify(A.tag),STAR,STAR,A.Grid())
    if vectors:
      Q = DistMatrix(A.tag,MC,MR,A.Grid())
      args = [A.obj,w.obj,Q.obj,fullTriangle]
      if   A.tag == sTag: lib.ElSchurDecompDist_s(*args)
      elif A.tag == dTag: lib.ElSchurDecompDist_d(*args)
      elif A.tag == cTag: lib.ElSchurDecompDist_c(*args)
      elif A.tag == zTag: lib.ElSchurDecompDist_z(*args)
      else: DataExcept()
      return w, Q
    else:
      args = [A.obj,w.obj,fullTriangle]
      if   A.tag == sTag: lib.ElSchurDist_s(*args)
      elif A.tag == dTag: lib.ElSchurDist_d(*args)
      elif A.tag == cTag: lib.ElSchurDist_c(*args)
      elif A.tag == zTag: lib.ElSchurDist_z(*args)
      else: DataExcept()
      return w
  else: TypeExcept()

# Singular value decomposition
# ============================

lib.ElSVD_s.argtypes = \
lib.ElSVD_d.argtypes = \
lib.ElSVD_c.argtypes = \
lib.ElSVD_z.argtypes = \
lib.ElSVDDist_s.argtypes = \
lib.ElSVDDist_d.argtypes = \
lib.ElSVDDist_c.argtypes = \
lib.ElSVDDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p]

lib.ElSVD_s.restype = \
lib.ElSVD_d.restype = \
lib.ElSVD_c.restype = \
lib.ElSVD_z.restype = \
lib.ElSVDDist_s.restype = \
lib.ElSVDDist_d.restype = \
lib.ElSVDDist_c.restype = \
lib.ElSVDDist_z.restype = \
  c_uint

lib.ElSingularValues_s.argtypes = \
lib.ElSingularValues_d.argtypes = \
lib.ElSingularValues_c.argtypes = \
lib.ElSingularValues_z.argtypes = \
lib.ElSingularValuesDist_s.argtypes = \
lib.ElSingularValuesDist_d.argtypes = \
lib.ElSingularValuesDist_c.argtypes = \
lib.ElSingularValuesDist_z.argtypes = \
  [c_void_p,c_void_p]

lib.ElSingularValues_s.restype = \
lib.ElSingularValues_d.restype = \
lib.ElSingularValues_c.restype = \
lib.ElSingularValues_z.restype = \
lib.ElSingularValuesDist_s.restype = \
lib.ElSingularValuesDist_d.restype = \
lib.ElSingularValuesDist_c.restype = \
lib.ElSingularValuesDist_z.restype = \
  c_uint

def SVD(A,vectors=False):
  if type(A) is Matrix:
    s = Matrix(Base(A.tag))
    if vectors:
      V = Matrix(A.tag)
      args = [A.obj,s.obj,V.obj]
      if   A.tag == sTag: lib.ElSVD_s(*args)
      elif A.tag == dTag: lib.ElSVD_d(*args)
      elif A.tag == cTag: lib.ElSVD_c(*args)
      elif A.tag == zTag: lib.ElSVD_z(*args)
      else: DataExcept()
      return s, V
    else:
      args = [A.obj,s.obj]
      if   A.tag == sTag: lib.ElSingularValues_s(*args)
      elif A.tag == dTag: lib.ElSingularValues_d(*args)
      elif A.tag == cTag: lib.ElSingularValues_c(*args)
      elif A.tag == zTag: lib.ElSingularValues_z(*args)
      else: DataExcept()
      return s
  elif type(A) is DistMatrix:
    s = DistMatrix(Base(A.tag),STAR,STAR,A.Grid())
    if vectors:
      V = DistMatrix(A.tag,MC,MR,A.Grid())
      args = [A.obj,s.obj,V.obj]
      if   A.tag == sTag: lib.ElSVDDist_s(*args)
      elif A.tag == dTag: lib.ElSVDDist_d(*args)
      elif A.tag == cTag: lib.ElSVDDist_c(*args)
      elif A.tag == zTag: lib.ElSVDDist_z(*args)
      else: DataExcept()
      return s, V
    else:
      args = [A.obj,s.obj]
      if   A.tag == sTag: lib.ElSingularValuesDist_s(*args)
      elif A.tag == dTag: lib.ElSingularValuesDist_d(*args)
      elif A.tag == cTag: lib.ElSingularValuesDist_c(*args)
      elif A.tag == zTag: lib.ElSingularValuesDist_z(*args)
      else: DataExcept()
      return s
  else: TypeExcept()

# Pseudospectra
# =============
lib.ElSnapshotCtrlDefault.argtypes = [c_void_p]
lib.ElSnapshotCtrlDefault.restype = c_uint
lib.ElSnapshotCtrlDestroy.argtypes = [c_void_p]
lib.ElSnapshotCtrlDestroy.restype = c_uint
class SnapshotCtrl(ctypes.Structure):
  _fields_ = [("realSize",iType),("imagSize",iType),
              ("imgSaveFreq",iType),("numSaveFreq",iType),
              ("imgDispFreq",iType),
              ("imgSaveCount",iType),("numSaveCount",iType),
              ("imgDispCount",iType),
              ("imgBase",c_char_p),("numBase",c_char_p),
              ("imgFormat",c_uint),("numFormat",c_uint),
              ("itCounts",bType)]
  def __init__(self):
    lib.ElSnaphsotCtrlDefault(pointer(self))
  def Destroy(self):
    lib.ElSnapshotCtrlDestroy(pointer(self))

# Emulate an enum for the pseudospectral norm
(PS_TWO_NORM,PS_ONE_NORM)=(0,1)

lib.ElPseudospecCtrlDefault_s.argtypes = [c_void_p]
lib.ElPseudospecCtrlDefault_s.restype = c_uint
lib.ElPseudospecCtrlDestroy_s.argtypes = [c_void_p]
lib.ElPseudospecCtrlDestroy_s.restype = c_uint
class PseudospecCtrl_s(ctypes.Structure):
  _fields_ = [("norm",c_uint),
              ("blockWidth",iType),
              ("schur",bType),
              ("forceComplexSchur",bType),
              ("forceComplexPs",bType),
              ("schurCtrl",SchurCtrl_s),
              ("maxIts",iType),
              ("tol",sType),
              ("deflate",bType),
              ("arnoldi",bType),
              ("basisSize",iType),
              ("reorthog",bType),
              ("progress",bType),
              ("snapCtrl",SnapshotCtrl),
              ("center",cType),
              ("realWidth",sType),
              ("imagWidth",sType)]
  def __init__(self):
    lib.ElPseudospecCtrlDefault_s(pointer(self))
  def Destroy(self):
    lib.ElPseudospecCtrlDestroy_s(pointer(self))

lib.ElPseudospecCtrlDefault_d.argtypes = [c_void_p]
lib.ElPseudospecCtrlDefault_d.restype = c_uint
lib.ElPseudospecCtrlDestroy_d.argtypes = [c_void_p]
lib.ElPseudospecCtrlDestroy_d.restype = c_uint
class PseudospecCtrl_d(ctypes.Structure):
  _fields_ = [("norm",c_uint),
              ("blockWidth",iType),
              ("schur",bType),
              ("forceComplexSchur",bType),
              ("forceComplexPs",bType),
              ("schurCtrl",SchurCtrl_d),
              ("maxIts",iType),
              ("tol",dType),
              ("deflate",bType),
              ("arnoldi",bType),
              ("basisSize",iType),
              ("reorthog",bType),
              ("progress",bType),
              ("snapCtrl",SnapshotCtrl),
              ("center",zType),
              ("realWidth",dType),
              ("imagWidth",dType)]
  def __init__(self):
    lib.ElPseudospecCtrlDefault_d(pointer(self))
  def Destroy(self):
    lib.ElPseudospecCtrlDestroy_d(pointer(self))

class SpectralBox_s(ctypes.Structure):
  _fields_ = [("center",cType),
              ("realWidth",sType),
              ("imagWidth",sType)]
class SpectralBox_d(ctypes.Structure):
  _fields_ = [("center",zType),
              ("realWidth",dType),
              ("imagWidth",dType)]

def DisplayPortrait(portrait,box,title='',tryPython=True):
  import math
  if tryPython:
    if type(portrait) is Matrix:
      EntrywiseMap(portrait,math.log10)
      try:
        import numpy as np
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        isInline = 'inline' in mpl.get_backend()
        isVec = min(portrait.Height(),portrait.Width()) == 1
        fig = plt.figure()
        axis = fig.add_axes([0.1,0.1,0.8,0.8])
        if isVec:
          axis.plot(np.squeeze(portrait.ToNumPy()),'bo-')
        else:
          lBound = box.center.real - box.realWidth/2
          rBound = box.center.real + box.realWidth/2 
          bBound = box.center.imag - box.imagWidth/2
          tBound = box.center.imag + box.imagWidth/2
          im = axis.imshow(portrait.ToNumPy(),
                           extent=[lBound,rBound,bBound,tBound])
          fig.colorbar(im,ax=axis)
        plt.title(title)
        plt.draw()
        if not isInline:
            plt.show(block=False)
        return
      except:
        print 'Could not import matplotlib.pyplot'
    elif type(portrait) is DistMatrix:
      portrait_CIRC_CIRC = DistMatrix(portrait.tag,CIRC,CIRC,portrait.Grid())
      Copy(portrait,portrait_CIRC_CIRC)
      if portrait_CIRC_CIRC.CrossRank() == portrait_CIRC_CIRC.Root():
        DisplayPortrait(portrait_CIRC_CIRC.Matrix(),box,title,True)
      return

  # Fall back to the built-in Display if we have not succeeded
  if not tryPython or type(portrait) is not Matrix:
    EntrywiseMap(portrait,math.log10)
  args = [portrait.obj,title]
  numMsExtra = 200
  if type(portrait) is Matrix:
    if   portrait.tag == sTag: lib.ElDisplay_s(*args)
    elif portrait.tag == dTag: lib.ElDisplay_d(*args)
    else: DataExcept()
    ProcessEvents(numMsExtra)
  elif type(portrait) is DistMatrix:
    if   portrait.tag == sTag: lib.ElDisplayDist_s(*args)
    elif portrait.tag == dTag: lib.ElDisplayDist_d(*args)
    else: DataExcept()
    ProcessEvents(numMsExtra)
  else: TypeExcept()

# (Pseudo-)Spectral portrait
# --------------------------
# The choice is based upon a few different norms of the Schur factor, as simply
# using the spectral radius would be insufficient for highly non-normal 
# matrices, e.g., a Jordan block with eigenvalue zero

lib.ElSpectralPortrait_s.argtypes = \
lib.ElSpectralPortrait_c.argtypes = \
lib.ElSpectralPortraitDist_s.argtypes = \
lib.ElSpectralPortraitDist_c.argtypes = \
  [c_void_p,c_void_p,iType,iType,POINTER(SpectralBox_s)]

lib.ElSpectralPortrait_d.argtypes = \
lib.ElSpectralPortrait_z.argtypes = \
lib.ElSpectralPortraitDist_d.argtypes = \
lib.ElSpectralPortraitDist_z.argtypes = \
  [c_void_p,c_void_p,iType,iType,POINTER(SpectralBox_d)]

lib.ElSpectralPortrait_s.restype = \
lib.ElSpectralPortrait_d.restype = \
lib.ElSpectralPortrait_c.restype = \
lib.ElSpectralPortrait_z.restype = \
lib.ElSpectralPortraitDist_s.restype = \
lib.ElSpectralPortraitDist_d.restype = \
lib.ElSpectralPortraitDist_c.restype = \
lib.ElSpectralPortraitDist_z.restype = \
  c_uint

lib.ElSpectralPortraitX_s.argtypes = \
lib.ElSpectralPortraitX_c.argtypes = \
lib.ElSpectralPortraitXDist_s.argtypes = \
lib.ElSpectralPortraitXDist_c.argtypes = \
  [c_void_p,c_void_p,iType,iType,POINTER(SpectralBox_s),PseudospecCtrl_s]

lib.ElSpectralPortraitX_d.argtypes = \
lib.ElSpectralPortraitX_z.argtypes = \
lib.ElSpectralPortraitXDist_d.argtypes = \
lib.ElSpectralPortraitXDist_z.argtypes = \
  [c_void_p,c_void_p,iType,iType,POINTER(SpectralBox_d),PseudospecCtrl_d]

lib.ElSpectralPortraitX_s.restype = \
lib.ElSpectralPortraitX_d.restype = \
lib.ElSpectralPortraitX_c.restype = \
lib.ElSpectralPortraitX_z.restype = \
lib.ElSpectralPortraitXDist_s.restype = \
lib.ElSpectralPortraitXDist_d.restype = \
lib.ElSpectralPortraitXDist_c.restype = \
lib.ElSpectralPortraitXDist_z.restype = \
  c_uint

def SpectralPortrait(A,realSize=200,imagSize=200,ctrl=None):
  if type(A) is Matrix:
    invNormMap = Matrix(Base(A.tag))
    if   A.tag == sTag:
      box = SpectralBox_s()
      args = [A.obj,invNormMap.obj,realSize,imagSize,pointer(box)]
      argsCtrl = [A.obj,invNormMap.obj,realSize,imagSize,pointer(box),ctrl]
      if ctrl == None: lib.ElSpectralPortrait_s(*args)
      else:            lib.ElSpectralPortraitX_s(*argsCtrl)
    elif A.tag == dTag:
      box = SpectralBox_d()
      args = [A.obj,invNormMap.obj,realSize,imagSize,pointer(box)]
      argsCtrl = [A.obj,invNormMap.obj,realSize,imagSize,pointer(box),ctrl]
      if ctrl == None: lib.ElSpectralPortrait_d(*args)
      else:            lib.ElSpectralPortraitX_d(*argsCtrl)
    elif A.tag == cTag:
      box = SpectralBox_s()
      args = [A.obj,invNormMap.obj,realSize,imagSize,pointer(box)]
      argsCtrl = [A.obj,invNormMap.obj,realSize,imagSize,pointer(box),ctrl]
      if ctrl == None: lib.ElSpectralPortrait_c(*args)
      else:            lib.ElSpectralPortraitX_c(*argsCtrl)
    elif A.tag == zTag:
      box = SpectralBox_d()
      args = [A.obj,invNormMap.obj,realSize,imagSize,pointer(box)]
      argsCtrl = [A.obj,invNormMap.obj,realSize,imagSize,pointer(box),ctrl]
      if ctrl == None: lib.ElSpectralPortrait_z(*args)
      else:            lib.ElSpectralPortraitX_z(*argsCtrl)
    else: DataExcept()
    return invNormMap, box
  elif type(A) is DistMatrix:
    invNormMap = DistMatrix(Base(A.tag),MC,MR,A.Grid())
    if   A.tag == sTag:
      box = SpectralBox_s()
      args = [A.obj,invNormMap.obj,realSize,imagSize,pointer(box)]
      argsCtrl = [A.obj,invNormMap.obj,realSize,imagSize,pointer(box),ctrl]
      if ctrl == None: lib.ElSpectralPortraitDist_s(*args)
      else:            lib.ElSpectralPortraitXDist_s(*argsCtrl)
    elif A.tag == dTag:
      box = SpectralBox_d()
      args = [A.obj,invNormMap.obj,realSize,imagSize,pointer(box)]
      argsCtrl = [A.obj,invNormMap.obj,realSize,imagSize,pointer(box),ctrl]
      if ctrl == None: lib.ElSpectralPortraitDist_d(*args)
      else:            lib.ElSpectralPortraitXDist_d(*argsCtrl)
    elif A.tag == cTag:
      box = SpectralBox_s()
      args = [A.obj,invNormMap.obj,realSize,imagSize,pointer(box)]
      argsCtrl = [A.obj,invNormMap.obj,realSize,imagSize,pointer(box),ctrl]
      if ctrl == None: lib.ElSpectralPortraitDist_c(*args)
      else:            lib.ElSpectralPortraitXDist_c(*argsCtrl)
    elif A.tag == zTag:
      box = SpectralBox_d()
      args = [A.obj,invNormMap.obj,realSize,imagSize,pointer(box)]
      argsCtrl = [A.obj,invNormMap.obj,realSize,imagSize,pointer(box),ctrl]
      if ctrl == None: lib.ElSpectralPortraitDist_z(*args)
      else:            lib.ElSpectralPortraitXDist_z(*argsCtrl)
    else: DataExcept()
    return invNormMap, box
  else: TypeExcept()

# (Pseudo-)Spectral window
# ------------------------
lib.ElSpectralWindow_s.argtypes = \
lib.ElSpectralWindowDist_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,sType,iType,iType]

lib.ElSpectralWindow_d.argtypes = \
lib.ElSpectralWindowDist_d.argtypes = \
  [c_void_p,c_void_p,dType,dType,dType,iType,iType]

lib.ElSpectralWindow_c.argtypes = \
lib.ElSpectralWindowDist_c.argtypes = \
  [c_void_p,c_void_p,cType,sType,sType,iType,iType]

lib.ElSpectralWindow_z.argtypes = \
lib.ElSpectralWindowDist_z.argtypes = \
  [c_void_p,c_void_p,zType,dType,dType,iType,iType]

lib.ElSpectralWindow_s.restype = \
lib.ElSpectralWindow_d.restype = \
lib.ElSpectralWindow_c.restype = \
lib.ElSpectralWindow_z.restype = \
lib.ElSpectralWindowDist_s.restype = \
lib.ElSpectralWindowDist_d.restype = \
lib.ElSpectralWindowDist_c.restype = \
lib.ElSpectralWindowDist_z.restype = \
  c_uint

lib.ElSpectralWindowX_s.argtypes = \
lib.ElSpectralWindowXDist_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,sType,iType,iType,PseudospecCtrl_s]

lib.ElSpectralWindowX_d.argtypes = \
lib.ElSpectralWindowXDist_d.argtypes = \
  [c_void_p,c_void_p,dType,dType,dType,iType,iType,PseudospecCtrl_d]

lib.ElSpectralWindowX_c.argtypes = \
lib.ElSpectralWindowXDist_c.argtypes = \
  [c_void_p,c_void_p,cType,sType,sType,iType,iType,PseudospecCtrl_s]

lib.ElSpectralWindowX_z.argtypes = \
lib.ElSpectralWindowXDist_z.argtypes = \
  [c_void_p,c_void_p,zType,dType,dType,iType,iType,PseudospecCtrl_d]

lib.ElSpectralWindowX_s.restype = \
lib.ElSpectralWindowX_d.restype = \
lib.ElSpectralWindowX_c.restype = \
lib.ElSpectralWindowX_z.restype = \
lib.ElSpectralWindowXDist_s.restype = \
lib.ElSpectralWindowXDist_d.restype = \
lib.ElSpectralWindowXDist_c.restype = \
lib.ElSpectralWindowXDist_z.restype = \
  c_uint

def SpectralWindow \
    (A,centerPre,realWidth,imagWidth,realSize=200,imagSize=200,ctrl=None):
  center = TagToType(A.tag)(centerPre)
  if type(A) is Matrix:
    invNormMap = Matrix(Base(A.tag))
    args = [A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize]
    argsCtrl = [A.obj,invNormMap.obj,center,realWidth,imagWidth,
                realSize,imagSize,ctrl]
    if   A.tag == sTag:
      if ctrl == None: lib.ElSpectralWindow_s(*args)
      else:            lib.ElSpectralWindowX_s(*argsCtrl)
    elif A.tag == dTag:
      if ctrl == None: lib.ElSpectralWindow_d(*args)
      else:            lib.ElSpectralWindowX_d(*argsCtrl)
    elif A.tag == cTag:
      if ctrl == None: lib.ElSpectralWindow_c(*args)
      else:            lib.ElSpectralWindowX_c(*argsCtrl)
    elif A.tag == zTag:
      if ctrl == None: lib.ElSpectralWindow_z(*args)
      else:            lib.ElSpectralWindowX_z(*argsCtrl)
    else: DataExcept()
    return invNormMap
  elif type(A) is DistMatrix:
    invNormMap = DistMatrix(Base(A.tag),MC,MR,A.Grid())
    args = [A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize]
    argsCtrl = [A.obj,invNormMap.obj,center,realWidth,imagWidth,
                realSize,imagSize,ctrl]
    if   A.tag == sTag:
      if ctrl == None: lib.ElSpectralWindowDist_s(*args)
      else:            lib.ElSpectralWindowXDist_s(*argsCtrl)
    elif A.tag == dTag:
      if ctrl == None: lib.ElSpectralWindowDist_d(*args)
      else:            lib.ElSpectralWindowXDist_d(*argsCtrl)
    elif A.tag == cTag:
      if ctrl == None: lib.ElSpectralWindowDist_c(*args)
      else:            lib.ElSpectralWindowXDist_c(*argsCtrl)
    elif A.tag == zTag:
      if ctrl == None: lib.ElSpectralWindowDist_z(*args)
      else:            lib.ElSpectralWindowXDist_z(*argsCtrl)
    else: DataExcept()
    return invNormMap
  else: TypeExcept()

# (Pseudo-)Spectral cloud
# -----------------------
lib.ElSpectralCloud_s.argtypes = \
lib.ElSpectralCloud_d.argtypes = \
lib.ElSpectralCloud_c.argtypes = \
lib.ElSpectralCloud_z.argtypes = \
lib.ElSpectralCloudDist_s.argtypes = \
lib.ElSpectralCloudDist_d.argtypes = \
lib.ElSpectralCloudDist_c.argtypes = \
lib.ElSpectralCloudDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p]

lib.ElSpectralCloud_s.restype = \
lib.ElSpectralCloud_d.restype = \
lib.ElSpectralCloud_c.restype = \
lib.ElSpectralCloud_z.restype = \
lib.ElSpectralCloudDist_s.restype = \
lib.ElSpectralCloudDist_d.restype = \
lib.ElSpectralCloudDist_c.restype = \
lib.ElSpectralCloudDist_z.restype = \
  c_uint

lib.ElSpectralCloudX_s.argtypes = \
lib.ElSpectralCloudX_c.argtypes = \
lib.ElSpectralCloudXDist_s.argtypes = \
lib.ElSpectralCloudXDist_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_s]

lib.ElSpectralCloudX_d.argtypes = \
lib.ElSpectralCloudX_z.argtypes = \
lib.ElSpectralCloudXDist_d.argtypes = \
lib.ElSpectralCloudXDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_d]

lib.ElSpectralCloudX_s.restype = \
lib.ElSpectralCloudX_d.restype = \
lib.ElSpectralCloudX_c.restype = \
lib.ElSpectralCloudX_z.restype = \
lib.ElSpectralCloudXDist_s.restype = \
lib.ElSpectralCloudXDist_d.restype = \
lib.ElSpectralCloudXDist_c.restype = \
lib.ElSpectralCloudXDist_z.restype = \
  c_uint

def SpectralCloud(A,shifts,ctrl=None):
  if type(A) is Matrix:
    invNorms = Matrix(Base(A.tag))
    args = [A.obj,shifts.obj,invNorms.obj]
    argsCtrl = [A.obj,shifts.obj,invNorms.obj,ctrl]
    if   A.tag == sTag:
      if ctrl == None: lib.ElSpectralCloud_s(*args)
      else:            lib.ElSpectralCloudX_s(*argsCtrl)
    elif A.tag == dTag:
      if ctrl == None: lib.ElSpectralCloud_d(*args)
      else:            lib.ElSpectralCloudX_d(*argsCtrl)
    elif A.tag == cTag:
      if ctrl == None: lib.ElSpectralCloud_c(*args)
      else:            lib.ElSpectralCloudX_c(*argsCtrl)
    elif A.tag == zTag:
      if ctrl == None: lib.ElSpectralCloud_z(*args)
      else:            lib.ElSpectralCloudX_z(*argsCtrl)
    else: DataExcept()
    return invNorms
  elif type(A) is DistMatrix:
    invNorms = DistMatrix(Base(A.tag),VR,STAR,A.Grid())
    args = [A.obj,shifts.obj,invNorms.obj]
    argsCtrl = [A.obj,shifts.obj,invNorms.obj,ctrl]
    if   A.tag == sTag:
      if ctrl == None: lib.ElSpectralCloudDist_s(*args)
      else:            lib.ElSpectralCloudXDist_s(*argsCtrl)
    elif A.tag == dTag:
      if ctrl == None: lib.ElSpectralCloudDist_d(*args)
      else:            lib.ElSpectralCloudXDist_d(*argsCtrl)
    elif A.tag == cTag:
      if ctrl == None: lib.ElSpectralCloudDist_c(*args)
      else:            lib.ElSpectralCloudXDist_c(*argsCtrl)
    elif A.tag == zTag:
      if ctrl == None: lib.ElSpectralCloudDist_z(*args)
      else:            lib.ElSpectralCloudXDist_z(*argsCtrl)
    else: DataExcept()
    return invNorms
  else: TypeExcept()
