#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *

# Coherence
# =========
lib.ElCoherence_s.argtypes = \
lib.ElCoherence_c.argtypes = \
lib.ElCoherenceDist_s.argtypes = \
lib.ElCoherenceDist_c.argtypes = \
  [c_void_p,POINTER(sType)]
lib.ElCoherence_d.argtypes = \
lib.ElCoherence_z.argtypes = \
lib.ElCoherenceDist_d.argtypes = \
lib.ElCoherenceDist_z.argtypes = \
  [c_void_p,POINTER(dType)]
lib.ElCoherence_s.restype = \
lib.ElCoherence_d.restype = \
lib.ElCoherence_c.restype = \
lib.ElCoherence_z.restype = \
lib.ElCoherenceDist_s.restype = \
lib.ElCoherenceDist_d.restype = \
lib.ElCoherenceDist_c.restype = \
lib.ElCoherenceDist_z.restype = \
  c_uint

def Coherence(A):
  value = TagToType(Base(A.tag))()
  args = [A.obj,pointer(value)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElCoherence_s(*args)
    elif A.tag == dTag: lib.ElCoherence_d(*args)
    elif A.tag == cTag: lib.ElCoherence_c(*args)
    elif A.tag == zTag: lib.ElCoherence_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElCoherenceDist_s(*args)
    elif A.tag == dTag: lib.ElCoherenceDist_d(*args)
    elif A.tag == cTag: lib.ElCoherenceDist_c(*args)
    elif A.tag == zTag: lib.ElCoherenceDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return value

# Covariance
# ==========
lib.ElCovariance_s.argtypes = \
lib.ElCovariance_d.argtypes = \
lib.ElCovariance_c.argtypes = \
lib.ElCovariance_z.argtypes = \
lib.ElCovarianceDist_s.argtypes = \
lib.ElCovarianceDist_d.argtypes = \
lib.ElCovarianceDist_c.argtypes = \
lib.ElCovarianceDist_z.argtypes = \
  [c_void_p,c_void_p]

lib.ElCovariance_s.restype = \
lib.ElCovariance_d.restype = \
lib.ElCovariance_c.restype = \
lib.ElCovariance_z.restype = \
lib.ElCovarianceDist_s.restype = \
lib.ElCovarianceDist_d.restype = \
lib.ElCovarianceDist_c.restype = \
lib.ElCovarianceDist_z.restype = \
  c_uint

def Covariance(D):
  if type(D) is Matrix:
    S = Matrix(D.tag)
    args = [D.obj,S.obj]
    if   D.tag == sTag: lib.ElCovariance_s(*args)
    elif D.tag == dTag: lib.ElCovariance_d(*args)
    elif D.tag == cTag: lib.ElCovariance_c(*args)
    elif D.tag == zTag: lib.ElCovariance_z(*args)
    else: DataExcept()
    return S
  elif type(D) is DistMatrix:
    S = DistMatrix(D.tag,MC,MR,D.Grid())
    args = [D.obj,S.obj]
    if   D.tag == sTag: lib.ElCovarianceDist_s(*args)
    elif D.tag == dTag: lib.ElCovarianceDist_d(*args)
    elif D.tag == cTag: lib.ElCovarianceDist_c(*args)
    elif D.tag == zTag: lib.ElCovarianceDist_z(*args)
    else: DataExcept()
    return S
  else: TypeExcept()

# Log barrier
# ===========
lib.ElLogBarrier_s.argtypes = \
lib.ElLogBarrier_c.argtypes = \
lib.ElLogBarrierDist_s.argtypes = \
lib.ElLogBarrierDist_c.argtypes = \
  [c_uint,c_void_p,POINTER(sType)]
lib.ElLogBarrier_d.argtypes = \
lib.ElLogBarrier_z.argtypes = \
lib.ElLogBarrierDist_d.argtypes = \
lib.ElLogBarrierDist_z.argtypes = \
  [c_uint,c_void_p,POINTER(dType)]

lib.ElLogBarrier_s.restype = \
lib.ElLogBarrier_d.restype = \
lib.ElLogBarrier_c.restype = \
lib.ElLogBarrier_z.restype = \
lib.ElLogBarrierDist_s.restype = \
lib.ElLogBarrierDist_d.restype = \
lib.ElLogBarrierDist_c.restype = \
lib.ElLogBarrierDist_z.restype = \
  c_uint

def LogBarrier(uplo,A):
  barrier = TagToType(Base(A.tag))()
  args = [uplo,A.obj,pointer(barrier)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElLogBarrier_s(*args)
    elif A.tag == dTag: lib.ElLogBarrier_d(*args)
    elif A.tag == cTag: lib.ElLogBarrier_c(*args)
    elif A.tag == zTag: lib.ElLogBarrier_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElLogBarrierDist_s(*args)
    elif A.tag == dTag: lib.ElLogBarrierDist_d(*args)
    elif A.tag == cTag: lib.ElLogBarrierDist_c(*args)
    elif A.tag == zTag: lib.ElLogBarrierDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return barrier

# Log-det divergence
# ==================
lib.ElLogDetDiv_s.argtypes = \
lib.ElLogDetDiv_c.argtypes = \
lib.ElLogDetDivDist_s.argtypes = \
lib.ElLogDetDivDist_c.argtypes = \
  [c_uint,c_void_p,c_void_p,POINTER(sType)]
lib.ElLogDetDiv_d.argtypes = \
lib.ElLogDetDiv_z.argtypes = \
lib.ElLogDetDivDist_d.argtypes = \
lib.ElLogDetDivDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,POINTER(dType)]

lib.ElLogDetDiv_s.restype = \
lib.ElLogDetDiv_d.restype = \
lib.ElLogDetDiv_c.restype = \
lib.ElLogDetDiv_z.restype = \
lib.ElLogDetDivDist_s.restype = \
lib.ElLogDetDivDist_d.restype = \
lib.ElLogDetDivDist_c.restype = \
lib.ElLogDetDivDist_z.restype = \
  c_uint

def LogDetDiv(uplo,A,B):
  div = TagToType(Base(A.tag))()
  args = [uplo,A.obj,B.obj,pointer(div)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElLogDetDiv_s(*args)
    elif A.tag == dTag: lib.ElLogDetDiv_d(*args)
    elif A.tag == cTag: lib.ElLogDetDiv_c(*args)
    elif A.tag == zTag: lib.ElLogDetDiv_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElLogDetDivDist_s(*args)
    elif A.tag == dTag: lib.ElLogDetDivDist_d(*args)
    elif A.tag == cTag: lib.ElLogDetDivDist_c(*args)
    elif A.tag == zTag: lib.ElLogDetDivDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return div
