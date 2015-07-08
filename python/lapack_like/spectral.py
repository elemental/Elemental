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

lib.ElHermitianTridiagEigPair_s.argtypes = \
lib.ElHermitianTridiagEigPair_d.argtypes = \
lib.ElHermitianTridiagEigPair_c.argtypes = \
lib.ElHermitianTridiagEigPair_z.argtypes = \
lib.ElHermitianTridiagEigPairDist_s.argtypes = \
lib.ElHermitianTridiagEigPairDist_d.argtypes = \
lib.ElHermitianTridiagEigPairDist_c.argtypes = \
lib.ElHermitianTridiagEigPairDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]

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

lib.ElHermitianEigPair_s.argtypes = \
lib.ElHermitianEigPair_d.argtypes = \
lib.ElHermitianEigPair_c.argtypes = \
lib.ElHermitianEigPair_z.argtypes = \
lib.ElHermitianEigPairDist_s.argtypes = \
lib.ElHermitianEigPairDist_d.argtypes = \
lib.ElHermitianEigPairDist_c.argtypes = \
lib.ElHermitianEigPairDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint]

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
      X = DistMatrix(A.tag,MC,MR,A.Grid())
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

lib.ElSkewHermitianEigPair_s.argtypes = \
lib.ElSkewHermitianEigPair_d.argtypes = \
lib.ElSkewHermitianEigPair_c.argtypes = \
lib.ElSkewHermitianEigPair_z.argtypes = \
lib.ElSkewHermitianEigPairDist_s.argtypes = \
lib.ElSkewHermitianEigPairDist_d.argtypes = \
lib.ElSkewHermitianEigPairDist_c.argtypes = \
lib.ElSkewHermitianEigPairDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_uint]

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
      X = DistMatrix(A.tag,MC,MR,A.Grid())
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

lib.ElHermitianGenDefEigPair_s.argtypes = \
lib.ElHermitianGenDefEigPair_d.argtypes = \
lib.ElHermitianGenDefEigPair_c.argtypes = \
lib.ElHermitianGenDefEigPair_z.argtypes = \
lib.ElHermitianGenDefEigPairDist_s.argtypes = \
lib.ElHermitianGenDefEigPairDist_d.argtypes = \
lib.ElHermitianGenDefEigPairDist_c.argtypes = \
lib.ElHermitianGenDefEigPairDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_uint]

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
      X = DistMatrix(A.tag,MC,MR,A.Grid())
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

lib.ElHermitianSVD_s.argtypes = \
lib.ElHermitianSVD_d.argtypes = \
lib.ElHermitianSVD_c.argtypes = \
lib.ElHermitianSVD_z.argtypes = \
lib.ElHermitianSVDDist_s.argtypes = \
lib.ElHermitianSVDDist_d.argtypes = \
lib.ElHermitianSVDDist_c.argtypes = \
lib.ElHermitianSVDDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]

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

lib.ElPolarDecomp_s.argtypes = \
lib.ElPolarDecomp_d.argtypes = \
lib.ElPolarDecomp_c.argtypes = \
lib.ElPolarDecomp_z.argtypes = \
lib.ElPolarDecompDist_s.argtypes = \
lib.ElPolarDecompDist_d.argtypes = \
lib.ElPolarDecompDist_c.argtypes = \
lib.ElPolarDecompDist_z.argtypes = \
  [c_void_p,c_void_p]

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

lib.ElHermitianPolarDecomp_s.argtypes = \
lib.ElHermitianPolarDecomp_d.argtypes = \
lib.ElHermitianPolarDecomp_c.argtypes = \
lib.ElHermitianPolarDecomp_z.argtypes = \
lib.ElHermitianPolarDecompDist_s.argtypes = \
lib.ElHermitianPolarDecompDist_d.argtypes = \
lib.ElHermitianPolarDecompDist_c.argtypes = \
lib.ElHermitianPolarDecompDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p]

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
class SignCtrl_s(ctypes.Structure):
  _fields_ = [("maxIts",iType),
              ("tol",sType),
              ("power",sType),
              ("scaling",c_uint)]
  def __init__(self):
    lib.ElSignCtrlDefault_s(pointer(self))

lib.ElSignCtrlDefault_d.argtypes = [c_void_p]
class SignCtrl_d(ctypes.Structure):
  _fields_ = [("maxIts",iType),
              ("tol",dType),
              ("power",dType),
              ("scaling",c_uint)]
  def __init__(self):
    lib.ElSignCtrlDefault_d(pointer(self))

lib.ElHessQRCtrlDefault.argtypes = [c_void_p]
class HessQRCtrl(ctypes.Structure):
  _fields_ = [("distAED",bType),
              ("blockHeight",iType),("blockWidth",iType)]
  def __init__(self):
    lib.ElHessQRCtrlDefault(pointer(self))

lib.ElSDCCtrlDefault_s.argtypes = [c_void_p]
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
class SchurCtrl_s(ctypes.Structure):
  _fields_ = [("useSDC",bType),
              ("qrCtrl",HessQRCtrl),
              ("sdcCtrl",SDCCtrl_s)]
  def __init__(self):
    lib.ElSchurCtrlDefault_s(pointer(self))

lib.ElSchurCtrlDefault_d.argtypes = [c_void_p]
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

lib.ElSchurDecomp_s.argtypes = \
lib.ElSchurDecomp_d.argtypes = \
lib.ElSchurDecomp_c.argtypes = \
lib.ElSchurDecomp_z.argtypes = \
lib.ElSchurDecompDist_s.argtypes = \
lib.ElSchurDecompDist_d.argtypes = \
lib.ElSchurDecompDist_c.argtypes = \
lib.ElSchurDecompDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,bType]

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

lib.ElSingularValues_s.argtypes = \
lib.ElSingularValues_d.argtypes = \
lib.ElSingularValues_c.argtypes = \
lib.ElSingularValues_z.argtypes = \
lib.ElSingularValuesDist_s.argtypes = \
lib.ElSingularValuesDist_d.argtypes = \
lib.ElSingularValuesDist_c.argtypes = \
lib.ElSingularValuesDist_z.argtypes = \
  [c_void_p,c_void_p]

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

# Product Lanczos
# ===============
lib.ElProductLanczosSparse_s.argtypes = \
lib.ElProductLanczosSparse_d.argtypes = \
lib.ElProductLanczosSparse_c.argtypes = \
lib.ElProductLanczosSparse_z.argtypes = \
lib.ElProductLanczosDistSparse_s.argtypes = \
lib.ElProductLanczosDistSparse_d.argtypes = \
lib.ElProductLanczosDistSparse_c.argtypes = \
lib.ElProductLanczosDistSparse_z.argtypes = \
  [c_void_p,c_void_p,iType]

def ProductLanczos(A,basisSize=20):
  T = Matrix(Base(A.tag))
  args = [A.obj,T.obj,basisSize]
  if type(A) is SparseMatrix:
    if   A.tag == sTag: lib.ElProductLanczosSparse_s(*args)
    elif A.tag == dTag: lib.ElProductLanczosSparse_d(*args)
    elif A.tag == cTag: lib.ElProductLanczosSparse_c(*args)
    elif A.tag == zTag: lib.ElProductLanczosSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == sTag: lib.ElProductLanczosDistSparse_s(*args)
    elif A.tag == dTag: lib.ElProductLanczosDistSparse_d(*args)
    elif A.tag == cTag: lib.ElProductLanczosDistSparse_c(*args)
    elif A.tag == zTag: lib.ElProductLanczosDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return T

lib.ElProductLanczosDecompSparse_s.argtypes = \
lib.ElProductLanczosDecompSparse_c.argtypes = \
lib.ElProductLanczosDecompDistSparse_s.argtypes = \
lib.ElProductLanczosDecompDistSparse_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,POINTER(sType),iType]
lib.ElProductLanczosDecompSparse_d.argtypes = \
lib.ElProductLanczosDecompSparse_z.argtypes = \
lib.ElProductLanczosDecompDistSparse_d.argtypes = \
lib.ElProductLanczosDecompDistSparse_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,POINTER(dType),iType]

def ProductLanczosDecomp(A,basisSize=20):
  T = Matrix(Base(A.tag))
  beta = TagToType(Base(A.tag))()
  if type(A) is SparseMatrix:
    V = Matrix(A.tag)
    v = Matrix(A.tag)
    args = [A.obj,V.obj,T.obj,v.obj,pointer(beta),basisSize]
    if   A.tag == sTag: lib.ElProductLanczosDecompSparse_s(*args)
    elif A.tag == dTag: lib.ElProductLanczosDecompSparse_d(*args)
    elif A.tag == cTag: lib.ElProductLanczosDecompSparse_c(*args)
    elif A.tag == zTag: lib.ElProductLanczosDecompSparse_z(*args)
    else: DataExcept()
    return V, T, v, beta.value
  elif type(A) is DistSparseMatrix:
    V = DistMultiVec(A.tag,A.Comm())
    v = DistMultiVec(A.tag,A.Comm())
    args = [A.obj,V.obj,T.obj,v.obj,pointer(beta),basisSize]
    if   A.tag == sTag: lib.ElProductLanczosDecompDistSparse_s(*args)
    elif A.tag == dTag: lib.ElProductLanczosDecompDistSparse_d(*args)
    elif A.tag == cTag: lib.ElProductLanczosDecompDistSparse_c(*args)
    elif A.tag == zTag: lib.ElProductLanczosDecompDistSparse_z(*args)
    else: DataExcept()
    return V, T, v, beta.value
  else: TypeExcept()

# Extremal singular value estimation
# ==================================
lib.ElExtremalSingValEstSparse_s.argtypes = \
lib.ElExtremalSingValEstSparse_c.argtypes = \
lib.ElExtremalSingValEstDistSparse_s.argtypes = \
lib.ElExtremalSingValEstDistSparse_c.argtypes = \
  [c_void_p,iType,POINTER(sType),POINTER(sType)]
lib.ElExtremalSingValEstSparse_d.argtypes = \
lib.ElExtremalSingValEstSparse_z.argtypes = \
lib.ElExtremalSingValEstDistSparse_d.argtypes = \
lib.ElExtremalSingValEstDistSparse_z.argtypes = \
  [c_void_p,iType,POINTER(dType),POINTER(dType)]

def ExtremalSingValEst(A,basisSize=20):
  sigMin = TagToType(Base(A.tag))()
  sigMax = TagToType(Base(A.tag))()
  args = [A.obj,basisSize,pointer(sigMin),pointer(sigMax)]
  if type(A) is SparseMatrix:
    if   A.tag == sTag: lib.ElExtremalSingValEstSparse_s(*args)
    elif A.tag == dTag: lib.ElExtremalSingValEstSparse_d(*args)
    elif A.tag == cTag: lib.ElExtremalSingValEstSparse_c(*args)
    elif A.tag == zTag: lib.ElExtremalSingValEstSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == sTag: lib.ElExtremalSingValEstDistSparse_s(*args)
    elif A.tag == dTag: lib.ElExtremalSingValEstDistSparse_d(*args)
    elif A.tag == cTag: lib.ElExtremalSingValEstDistSparse_c(*args)
    elif A.tag == zTag: lib.ElExtremalSingValEstDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return sigMin, sigMax

# Pseudospectra
# =============
lib.ElSnapshotCtrlDefault.argtypes = \
lib.ElSnapshotCtrlDestroy.argtypes = \
  [c_void_p]
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

lib.ElPseudospecCtrlDefault_s.argtypes = \
lib.ElPseudospecCtrlDestroy_s.argtypes = \
  [c_void_p]
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

lib.ElPseudospecCtrlDefault_d.argtypes = \
lib.ElPseudospecCtrlDestroy_d.argtypes = \
  [c_void_p]
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
  [c_void_p,c_void_p,cType,sType,sType,iType,iType]

lib.ElSpectralWindow_d.argtypes = \
lib.ElSpectralWindowDist_d.argtypes = \
  [c_void_p,c_void_p,zType,dType,dType,iType,iType]

lib.ElSpectralWindow_c.argtypes = \
lib.ElSpectralWindowDist_c.argtypes = \
  [c_void_p,c_void_p,cType,sType,sType,iType,iType]

lib.ElSpectralWindow_z.argtypes = \
lib.ElSpectralWindowDist_z.argtypes = \
  [c_void_p,c_void_p,zType,dType,dType,iType,iType]

lib.ElSpectralWindowX_s.argtypes = \
lib.ElSpectralWindowXDist_s.argtypes = \
  [c_void_p,c_void_p,cType,sType,sType,iType,iType,PseudospecCtrl_s]

lib.ElSpectralWindowX_d.argtypes = \
lib.ElSpectralWindowXDist_d.argtypes = \
  [c_void_p,c_void_p,zType,dType,dType,iType,iType,PseudospecCtrl_d]

lib.ElSpectralWindowX_c.argtypes = \
lib.ElSpectralWindowXDist_c.argtypes = \
  [c_void_p,c_void_p,cType,sType,sType,iType,iType,PseudospecCtrl_s]

lib.ElSpectralWindowX_z.argtypes = \
lib.ElSpectralWindowXDist_z.argtypes = \
  [c_void_p,c_void_p,zType,dType,dType,iType,iType,PseudospecCtrl_d]

def SpectralWindow \
    (A,centerPre,realWidth,imagWidth,realSize=200,imagSize=200,ctrl=None):
  center = TagToType(Complexify(A.tag))(centerPre)
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
