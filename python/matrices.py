#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *

from ctypes import CFUNCTYPE

# Special matrices
# ****************

# Deterministic
# =============

# Bull's head
# -----------
lib.ElBullsHead_c.argtypes = \
lib.ElBullsHead_z.argtypes = \
lib.ElBullsHeadDist_c.argtypes = \
lib.ElBullsHeadDist_z.argtypes = \
  [c_void_p,iType]

def BullsHead(A,n):
  args = [A.obj,n]
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElBullsHead_c(*args)
    elif A.tag == zTag: lib.ElBullsHead_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElBullsHeadDist_c(*args)
    elif A.tag == zTag: lib.ElBullsHeadDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Cauchy
# ------
lib.ElCauchy_s.argtypes = \
lib.ElCauchyDist_s.argtypes = \
  [c_void_p,iType,POINTER(sType),iType,POINTER(sType)]
lib.ElCauchy_d.argtypes = \
lib.ElCauchyDist_d.argtypes = \
  [c_void_p,iType,POINTER(dType),iType,POINTER(dType)]
lib.ElCauchy_c.argtypes = \
lib.ElCauchyDist_c.argtypes = \
  [c_void_p,iType,POINTER(cType),iType,POINTER(cType)]
lib.ElCauchy_z.argtypes = \
lib.ElCauchyDist_z.argtypes = \
  [c_void_p,iType,POINTER(zType),iType,POINTER(zType)]

def Cauchy(A,x,y):
  xLen = len(x)
  yLen = len(y)
  xBuf = (TagToType(A.tag)*xLen)(*x)
  yBuf = (TagToType(A.tag)*yLen)(*y)
  args = [A.obj,xLen,xBuf.yLen,yBuf]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElCauchy_s(*args)
    elif A.tag == dTag: lib.ElCauchy_d(*args)
    elif A.tag == cTag: lib.ElCauchy_c(*args)
    elif A.tag == zTag: lib.ElCauchy_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElCauchyDist_s(*args)
    elif A.tag == dTag: lib.ElCauchyDist_d(*args)
    elif A.tag == cTag: lib.ElCauchyDist_c(*args)
    elif A.tag == zTag: lib.ElCauchyDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Cauchy-like
# -----------
lib.ElCauchyLike_s.argtypes = \
lib.ElCauchyLikeDist_s.argtypes = \
  [c_void_p,iType,POINTER(sType),iType,POINTER(sType),
            iType,POINTER(sType),iType,POINTER(sType)]
lib.ElCauchyLike_d.argtypes = \
lib.ElCauchyLikeDist_d.argtypes = \
  [c_void_p,iType,POINTER(dType),iType,POINTER(dType),
            iType,POINTER(dType),iType,POINTER(dType)]
lib.ElCauchyLike_c.argtypes = \
lib.ElCauchyLikeDist_c.argtypes = \
  [c_void_p,iType,POINTER(cType),iType,POINTER(cType),
            iType,POINTER(cType),iType,POINTER(cType)]
lib.ElCauchyLike_z.argtypes = \
lib.ElCauchyLikeDist_z.argtypes = \
  [c_void_p,iType,POINTER(zType),iType,POINTER(zType),
            iType,POINTER(zType),iType,POINTER(zType)]

def CauchyLike(A,r,s,x,y):
  rLen = len(r)
  sLen = len(s)
  xLen = len(x)
  yLen = len(y)
  rBuf = (TagToType(A.tag)*rLen)(*r)
  sBuf = (TagToType(A.tag)*sLen)(*s)
  xBuf = (TagToType(A.tag)*xLen)(*x)
  yBuf = (TagToType(A.tag)*yLen)(*y)
  args = [A.obj,rLen,rBuf,sLen,sBuf,xLen,xBuf,yLen,yBuf]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElCauchyLike_s(*args)
    elif A.tag == dTag: lib.ElCauchyLike_d(*args)
    elif A.tag == cTag: lib.ElCauchyLike_c(*args)
    elif A.tag == zTag: lib.ElCauchyLike_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElCauchyLikeDist_s(*args)
    elif A.tag == dTag: lib.ElCauchyLikeDist_d(*args)
    elif A.tag == cTag: lib.ElCauchyLikeDist_c(*args)
    elif A.tag == zTag: lib.ElCauchyLikeDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Circulant
# ---------
lib.ElCirculant_i.argtypes = \
lib.ElCirculantDist_i.argtypes = \
  [c_void_p,iType,POINTER(iType)]

lib.ElCirculant_s.argtypes = \
lib.ElCirculantDist_s.argtypes = \
  [c_void_p,iType,POINTER(sType)]

lib.ElCirculant_d.argtypes = \
lib.ElCirculantDist_d.argtypes = \
  [c_void_p,iType,POINTER(dType)]

lib.ElCirculant_c.argtypes = \
lib.ElCirculantDist_c.argtypes = \
  [c_void_p,iType,POINTER(cType)]

lib.ElCirculant_z.argtypes = \
lib.ElCirculantDist_z.argtypes = \
  [c_void_p,iType,POINTER(zType)]

def Circulant(A,a):
  aLen = len(a)
  aBuf = (TagToType(A.tag)*aLen)(*a)
  args = [A.obj,aLen,aBuf]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElCirculant_i(*args)
    elif A.tag == sTag: lib.ElCirculant_s(*args)
    elif A.tag == dTag: lib.ElCirculant_d(*args)
    elif A.tag == cTag: lib.ElCirculant_c(*args)
    elif A.tag == zTag: lib.ElCirculant_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElCirculantDist_i(*args)
    elif A.tag == sTag: lib.ElCirculantDist_s(*args)
    elif A.tag == dTag: lib.ElCirculantDist_d(*args)
    elif A.tag == cTag: lib.ElCirculantDist_c(*args)
    elif A.tag == zTag: lib.ElCirculantDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Demmel
# ------
lib.ElDemmel_s.argtypes = \
lib.ElDemmel_d.argtypes = \
lib.ElDemmel_c.argtypes = \
lib.ElDemmel_z.argtypes = \
lib.ElDemmelDist_s.argtypes = \
lib.ElDemmelDist_d.argtypes = \
lib.ElDemmelDist_c.argtypes = \
lib.ElDemmelDist_z.argtypes = \
  [c_void_p,iType]

def Demmel(A,n):
  args = [A.obj,n]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElDemmel_s(*args)
    elif A.tag == dTag: lib.ElDemmel_d(*args)
    elif A.tag == cTag: lib.ElDemmel_c(*args)
    elif A.tag == zTag: lib.ElDemmel_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElDemmelDist_s(*args)
    elif A.tag == dTag: lib.ElDemmelDist_d(*args)
    elif A.tag == cTag: lib.ElDemmelDist_c(*args)
    elif A.tag == zTag: lib.ElDemmelDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Diagonal
# --------
lib.ElDiagonal_i.argtypes = \
lib.ElDiagonal_s.argtypes = \
lib.ElDiagonal_d.argtypes = \
lib.ElDiagonal_c.argtypes = \
lib.ElDiagonal_z.argtypes = \
lib.ElDiagonalDist_i.argtypes = \
lib.ElDiagonalDist_s.argtypes = \
lib.ElDiagonalDist_d.argtypes = \
lib.ElDiagonalDist_c.argtypes = \
lib.ElDiagonalDist_z.argtypes = \
lib.ElDiagonalSparse_i.argtypes = \
lib.ElDiagonalSparse_s.argtypes = \
lib.ElDiagonalSparse_d.argtypes = \
lib.ElDiagonalSparse_c.argtypes = \
lib.ElDiagonalSparse_z.argtypes = \
lib.ElDiagonalDistSparse_i.argtypes = \
lib.ElDiagonalDistSparse_s.argtypes = \
lib.ElDiagonalDistSparse_d.argtypes = \
lib.ElDiagonalDistSparse_c.argtypes = \
lib.ElDiagonalDistSparse_z.argtypes = \
  [c_void_p,c_void_p]

def Diagonal(A,d):
  args = [A.obj,d.obj]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElDiagonal_i(*args)
    elif A.tag == sTag: lib.ElDiagonal_s(*args)
    elif A.tag == dTag: lib.ElDiagonal_d(*args)
    elif A.tag == cTag: lib.ElDiagonal_c(*args)
    elif A.tag == zTag: lib.ElDiagonal_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElDiagonalDist_i(*args)
    elif A.tag == sTag: lib.ElDiagonalDist_s(*args)
    elif A.tag == dTag: lib.ElDiagonalDist_d(*args)
    elif A.tag == cTag: lib.ElDiagonalDist_c(*args)
    elif A.tag == zTag: lib.ElDiagonalDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == iTag: lib.ElDiagonalSparse_i(*args)
    elif A.tag == sTag: lib.ElDiagonalSparse_s(*args)
    elif A.tag == dTag: lib.ElDiagonalSparse_d(*args)
    elif A.tag == cTag: lib.ElDiagonalSparse_c(*args)
    elif A.tag == zTag: lib.ElDiagonalSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == iTag: lib.ElDiagonalDistSparse_i(*args)
    elif A.tag == sTag: lib.ElDiagonalDistSparse_s(*args)
    elif A.tag == dTag: lib.ElDiagonalDistSparse_d(*args)
    elif A.tag == cTag: lib.ElDiagonalDistSparse_c(*args)
    elif A.tag == zTag: lib.ElDiagonalDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

# DruinskyToledo
# --------------
lib.ElDruinskyToledo_s.argtypes = \
lib.ElDruinskyToledo_d.argtypes = \
lib.ElDruinskyToledo_c.argtypes = \
lib.ElDruinskyToledo_z.argtypes = \
lib.ElDruinskyToledoDist_s.argtypes = \
lib.ElDruinskyToledoDist_d.argtypes = \
lib.ElDruinskyToledoDist_c.argtypes = \
lib.ElDruinskyToledoDist_z.argtypes = \
  [c_void_p,iType]

def DruinskyToledo(A,k):
  args = [A.obj,k]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElDruinskyToledo_s(*args)
    elif A.tag == dTag: lib.ElDruinskyToledo_d(*args)
    elif A.tag == cTag: lib.ElDruinskyToledo_c(*args)
    elif A.tag == zTag: lib.ElDruinskyToledo_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElDruinskyToledoDist_s(*args)
    elif A.tag == dTag: lib.ElDruinskyToledoDist_d(*args)
    elif A.tag == cTag: lib.ElDruinskyToledoDist_c(*args)
    elif A.tag == zTag: lib.ElDruinskyToledoDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Dynamic regularization counter-example
# --------------------------------------
lib.ElDynamicRegCounter_s.argtypes = \
lib.ElDynamicRegCounter_d.argtypes = \
lib.ElDynamicRegCounter_c.argtypes = \
lib.ElDynamicRegCounter_z.argtypes = \
lib.ElDynamicRegCounterDist_s.argtypes = \
lib.ElDynamicRegCounterDist_d.argtypes = \
lib.ElDynamicRegCounterDist_c.argtypes = \
lib.ElDynamicRegCounterDist_z.argtypes = \
lib.ElDynamicRegCounterSparse_s.argtypes = \
lib.ElDynamicRegCounterSparse_d.argtypes = \
lib.ElDynamicRegCounterSparse_c.argtypes = \
lib.ElDynamicRegCounterSparse_z.argtypes = \
lib.ElDynamicRegCounterDistSparse_s.argtypes = \
lib.ElDynamicRegCounterDistSparse_d.argtypes = \
lib.ElDynamicRegCounterDistSparse_c.argtypes = \
lib.ElDynamicRegCounterDistSparse_z.argtypes = \
  [c_void_p,iType]

def DynamicRegCounter(A,n):
  args = [A.obj,n]
  if type(A) is Matrix: 
    if   A.tag == sTag: lib.ElDynamicRegCounter_s(*args)
    elif A.tag == dTag: lib.ElDynamicRegCounter_d(*args)
    elif A.tag == cTag: lib.ElDynamicRegCounter_c(*args)
    elif A.tag == zTag: lib.ElDynamicRegCounter_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix: 
    if   A.tag == sTag: lib.ElDynamicRegCounterDist_s(*args)
    elif A.tag == dTag: lib.ElDynamicRegCounterDist_d(*args)
    elif A.tag == cTag: lib.ElDynamicRegCounterDist_c(*args)
    elif A.tag == zTag: lib.ElDynamicRegCounterDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix: 
    if   A.tag == sTag: lib.ElDynamicRegCounterSparse_s(*args)
    elif A.tag == dTag: lib.ElDynamicRegCounterSparse_d(*args)
    elif A.tag == cTag: lib.ElDynamicRegCounterSparse_c(*args)
    elif A.tag == zTag: lib.ElDynamicRegCounterSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix: 
    if   A.tag == sTag: lib.ElDynamicRegCounterDistSparse_s(*args)
    elif A.tag == dTag: lib.ElDynamicRegCounterDistSparse_d(*args)
    elif A.tag == cTag: lib.ElDynamicRegCounterDistSparse_c(*args)
    elif A.tag == zTag: lib.ElDynamicRegCounterDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Egorov
# ------
lib.ElEgorov_c.argtypes = \
lib.ElEgorovDist_c.argtypes = \
  [c_void_p,CFUNCTYPE(sType,iType,iType),iType]

lib.ElEgorov_z.argtypes = \
lib.ElEgorovDist_z.argtypes = \
  [c_void_p,CFUNCTYPE(dType,iType,iType),iType]

def Egorov(A,phase,n):
  cPhase = CFUNCTYPE(TagToType(Base(A.tag)),iType,iType)(phase)
  args = [A.obj,cPhase,n]
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElEgorov_c(*args)
    elif A.tag == zTag: lib.ElEgorov_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElEgorovDist_c(*args)
    elif A.tag == zTag: lib.ElEgorovDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Ehrenfest
# ---------
lib.ElEhrenfest_s.argtypes = \
lib.ElEhrenfest_d.argtypes = \
lib.ElEhrenfest_c.argtypes = \
lib.ElEhrenfest_z.argtypes = \
lib.ElEhrenfestDist_s.argtypes = \
lib.ElEhrenfestDist_d.argtypes = \
lib.ElEhrenfestDist_c.argtypes = \
lib.ElEhrenfestDist_z.argtypes = \
  [c_void_p,iType]

def Ehrenfest(P,n):
  args = [P.obj,n]
  if type(P) is Matrix:
    if   P.tag == sTag: lib.ElEhrenfest_s(*args)
    elif P.tag == dTag: lib.ElEhrenfest_d(*args)
    elif P.tag == cTag: lib.ElEhrenfest_c(*args)
    elif P.tag == zTag: lib.ElEhrenfest_z(*args)
    else: DataExcept()
  elif type(P) is DistMatrix:
    if   P.tag == sTag: lib.ElEhrenfestDist_s(*args)
    elif P.tag == dTag: lib.ElEhrenfestDist_d(*args)
    elif P.tag == cTag: lib.ElEhrenfestDist_c(*args)
    elif P.tag == zTag: lib.ElEhrenfestDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElEhrenfestStationary_s.argtypes = \
lib.ElEhrenfestStationary_d.argtypes = \
lib.ElEhrenfestStationary_c.argtypes = \
lib.ElEhrenfestStationary_z.argtypes = \
lib.ElEhrenfestStationaryDist_s.argtypes = \
lib.ElEhrenfestStationaryDist_d.argtypes = \
lib.ElEhrenfestStationaryDist_c.argtypes = \
lib.ElEhrenfestStationaryDist_z.argtypes = \
  [c_void_p,iType]

def EhrenfestStationary(PInf,n):
  args = [PInf.obj,n]
  if type(PInf) is Matrix:
    if   PInf.tag == sTag: lib.ElEhrenfestStationary_s(*args)
    elif PInf.tag == dTag: lib.ElEhrenfestStationary_d(*args)
    elif PInf.tag == cTag: lib.ElEhrenfestStationary_c(*args)
    elif PInf.tag == zTag: lib.ElEhrenfestStationary_z(*args)
    else: DataExcept()
  elif type(PInf) is DistMatrix:
    if   PInf.tag == sTag: lib.ElEhrenfestStationaryDist_s(*args)
    elif PInf.tag == dTag: lib.ElEhrenfestStationaryDist_d(*args)
    elif PInf.tag == cTag: lib.ElEhrenfestStationaryDist_c(*args)
    elif PInf.tag == zTag: lib.ElEhrenfestStationaryDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElEhrenfestDecay_s.argtypes = \
lib.ElEhrenfestDecay_d.argtypes = \
lib.ElEhrenfestDecay_c.argtypes = \
lib.ElEhrenfestDecay_z.argtypes = \
lib.ElEhrenfestDecayDist_s.argtypes = \
lib.ElEhrenfestDecayDist_d.argtypes = \
lib.ElEhrenfestDecayDist_c.argtypes = \
lib.ElEhrenfestDecayDist_z.argtypes = \
  [c_void_p,iType]

def EhrenfestDecay(PInf,n):
  args = [PInf.obj,n]
  if type(PInf) is Matrix:
    if   PInf.tag == sTag: lib.ElEhrenfestDecay_s(*args)
    elif PInf.tag == dTag: lib.ElEhrenfestDecay_d(*args)
    elif PInf.tag == cTag: lib.ElEhrenfestDecay_c(*args)
    elif PInf.tag == zTag: lib.ElEhrenfestDecay_z(*args)
    else: DataExcept()
  elif type(PInf) is DistMatrix:
    if   PInf.tag == sTag: lib.ElEhrenfestDecayDist_s(*args)
    elif PInf.tag == dTag: lib.ElEhrenfestDecayDist_d(*args)
    elif PInf.tag == cTag: lib.ElEhrenfestDecayDist_c(*args)
    elif PInf.tag == zTag: lib.ElEhrenfestDecayDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Extended Kahan
# --------------
lib.ElExtendedKahan_s.argtypes = \
lib.ElExtendedKahan_c.argtypes = \
lib.ElExtendedKahanDist_s.argtypes = \
lib.ElExtendedKahanDist_c.argtypes = \
  [c_void_p,iType,sType,sType]

lib.ElExtendedKahan_d.argtypes = \
lib.ElExtendedKahan_z.argtypes = \
lib.ElExtendedKahanDist_d.argtypes = \
lib.ElExtendedKahanDist_z.argtypes = \
  [c_void_p,iType,dType,dType]

def ExtendedKahan(A,k,phi,mu):
  args = [A.obj,k,phi,mu]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElExtendedKahan_s(*args)
    elif A.tag == dTag: lib.ElExtendedKahan_d(*args)
    elif A.tag == cTag: lib.ElExtendedKahan_c(*args)
    elif A.tag == zTag: lib.ElExtendedKahan_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElExtendedKahanDist_s(*args)
    elif A.tag == dTag: lib.ElExtendedKahanDist_d(*args)
    elif A.tag == cTag: lib.ElExtendedKahanDist_c(*args)
    elif A.tag == zTag: lib.ElExtendedKahanDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Fiedler
# -------
lib.ElFiedler_s.argtypes = \
lib.ElFiedlerDist_s.argtypes = \
  [c_void_p,iType,POINTER(sType)]

lib.ElFiedler_d.argtypes = \
lib.ElFiedlerDist_d.argtypes = \
  [c_void_p,iType,POINTER(dType)]

lib.ElFiedler_c.argtypes = \
lib.ElFiedlerDist_c.argtypes = \
  [c_void_p,iType,POINTER(cType)]

lib.ElFiedler_z.argtypes = \
lib.ElFiedlerDist_z.argtypes = \
  [c_void_p,iType,POINTER(zType)]

def Fiedler(A,c):
  cLen = len(c)
  cBuf = (TagToType(A.tag)*cLen)(*c)
  args = [A.obj,cLen,cBuf]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElFiedler_s(*args)
    elif A.tag == dTag: lib.ElFiedler_d(*args)
    elif A.tag == cTag: lib.ElFiedler_c(*args)
    elif A.tag == zTag: lib.ElFiedler_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElFiedlerDist_s(*args)
    elif A.tag == dTag: lib.ElFiedlerDist_d(*args)
    elif A.tag == cTag: lib.ElFiedlerDist_c(*args)
    elif A.tag == zTag: lib.ElFiedlerDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Forsythe
# --------
lib.ElForsythe_i.argtypes = \
lib.ElForsytheDist_i.argtypes = \
  [c_void_p,iType,iType,iType]

lib.ElForsythe_s.argtypes = \
lib.ElForsytheDist_s.argtypes = \
  [c_void_p,iType,sType,sType]

lib.ElForsythe_d.argtypes = \
lib.ElForsytheDist_d.argtypes = \
  [c_void_p,iType,dType,dType]

lib.ElForsythe_c.argtypes = \
lib.ElForsytheDist_c.argtypes = \
  [c_void_p,iType,cType,cType]

lib.ElForsythe_z.argtypes = \
lib.ElForsytheDist_z.argtypes = \
  [c_void_p,iType,zType,zType]

def Forsythe(J,n,alpha,lamb):
  args = [A.obj,n,alpha,lamb]
  if type(J) is Matrix:
    if   J.tag == iTag: lib.ElForsythe_i(*args)
    elif J.tag == sTag: lib.ElForsythe_s(*args)
    elif J.tag == dTag: lib.ElForsythe_d(*args)
    elif J.tag == cTag: lib.ElForsythe_c(*args)
    elif J.tag == zTag: lib.ElForsythe_z(*args)
    else: DataExcept()
  elif type(J) is DistMatrix:
    if   J.tag == iTag: lib.ElForsytheDist_i(*args)
    elif J.tag == sTag: lib.ElForsytheDist_s(*args)
    elif J.tag == dTag: lib.ElForsytheDist_d(*args)
    elif J.tag == cTag: lib.ElForsytheDist_c(*args)
    elif J.tag == zTag: lib.ElForsytheDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Fox-Li
# ------
lib.ElFoxLi_c.argtypes = \
lib.ElFoxLiDist_c.argtypes = \
  [c_void_p,iType,sType]

lib.ElFoxLi_z.argtypes = \
lib.ElFoxLiDist_z.argtypes = \
  [c_void_p,iType,dType]

def FoxLi(A,n,omega=48.):
  args = [A.obj,n,omega]
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElFoxLi_c(*args)
    elif A.tag == zTag: lib.ElFoxLi_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElFoxLiDist_c(*args)
    elif A.tag == zTag: lib.ElFoxLiDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Fourier
# -------
lib.ElFourier_c.argtypes = \
lib.ElFourier_z.argtypes = \
lib.ElFourierDist_c.argtypes = \
lib.ElFourierDist_z.argtypes = \
  [c_void_p,iType]

def Fourier(A,n):
  args = [A.obj,n]
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElFourier_c(*args)
    elif A.tag == zTag: lib.ElFourier_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElFourierDist_c(*args)
    elif A.tag == zTag: lib.ElFourierDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Fourier-Identity
# ----------------
lib.ElFourierIdentity_c.argtypes = \
lib.ElFourierIdentity_z.argtypes = \
lib.ElFourierIdentityDist_c.argtypes = \
lib.ElFourierIdentityDist_z.argtypes = \
  [c_void_p,iType]

def FourierIdentity(A,n):
  args = [A.obj,n]
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElFourierIdentity_c(*args)
    elif A.tag == zTag: lib.ElFourierIdentity_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElFourierIdentityDist_c(*args)
    elif A.tag == zTag: lib.ElFourierIdentityDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# GCD matrix
# ----------
lib.ElGCDMatrix_i.argtypes = \
lib.ElGCDMatrix_s.argtypes = \
lib.ElGCDMatrix_d.argtypes = \
lib.ElGCDMatrix_c.argtypes = \
lib.ElGCDMatrix_z.argtypes = \
lib.ElGCDMatrixDist_i.argtypes = \
lib.ElGCDMatrixDist_s.argtypes = \
lib.ElGCDMatrixDist_d.argtypes = \
lib.ElGCDMatrixDist_c.argtypes = \
lib.ElGCDMatrixDist_z.argtypes = \
  [c_void_p,iType,iType]

def GCDMatrix(G,m,n):
  args = [G.obj,m,n]
  if type(G) is Matrix:
    if   G.tag == iTag: lib.ElGCDMatrix_i(*args)
    elif G.tag == sTag: lib.ElGCDMatrix_s(*args)
    elif G.tag == dTag: lib.ElGCDMatrix_d(*args)
    elif G.tag == cTag: lib.ElGCDMatrix_c(*args)
    elif G.tag == zTag: lib.ElGCDMatrix_z(*args)
    else: DataExcept()
  elif type(G) is DistMatrix:
    if   G.tag == iTag: lib.ElGCDMatrixDist_i(*args)
    elif G.tag == sTag: lib.ElGCDMatrixDist_s(*args)
    elif G.tag == dTag: lib.ElGCDMatrixDist_d(*args)
    elif G.tag == cTag: lib.ElGCDMatrixDist_c(*args)
    elif G.tag == zTag: lib.ElGCDMatrixDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Gear matrix
# -----------
lib.ElGear_i.argtypes = \
lib.ElGear_s.argtypes = \
lib.ElGear_d.argtypes = \
lib.ElGear_c.argtypes = \
lib.ElGear_z.argtypes = \
lib.ElGearDist_i.argtypes = \
lib.ElGearDist_s.argtypes = \
lib.ElGearDist_d.argtypes = \
lib.ElGearDist_c.argtypes = \
lib.ElGearDist_z.argtypes = \
  [c_void_p,iType,iType,iType]

def Gear(G,n,s,t):
  args = [G.obj,n,s,t]
  if type(G) is Matrix:
    if   G.tag == iTag: lib.ElGear_i(*args)
    elif G.tag == sTag: lib.ElGear_s(*args)
    elif G.tag == dTag: lib.ElGear_d(*args)
    elif G.tag == cTag: lib.ElGear_c(*args)
    elif G.tag == zTag: lib.ElGear_z(*args)
    else: DataExcept()
  elif type(G) is DistMatrix:
    if   G.tag == iTag: lib.ElGearDist_i(*args)
    elif G.tag == sTag: lib.ElGearDist_s(*args)
    elif G.tag == dTag: lib.ElGearDist_d(*args)
    elif G.tag == cTag: lib.ElGearDist_c(*args)
    elif G.tag == zTag: lib.ElGearDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# GEPP Growth
# -----------
lib.ElGEPPGrowth_s.argtypes = \
lib.ElGEPPGrowth_d.argtypes = \
lib.ElGEPPGrowth_c.argtypes = \
lib.ElGEPPGrowth_z.argtypes = \
lib.ElGEPPGrowthDist_s.argtypes = \
lib.ElGEPPGrowthDist_d.argtypes = \
lib.ElGEPPGrowthDist_c.argtypes = \
lib.ElGEPPGrowthDist_z.argtypes = \
  [c_void_p,iType]

def GEPPGrowth(A,n):
  args = [A.obj,n]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElGEPPGrowth_s(*args)
    elif A.tag == dTag: lib.ElGEPPGrowth_d(*args)
    elif A.tag == cTag: lib.ElGEPPGrowth_c(*args)
    elif A.tag == zTag: lib.ElGEPPGrowth_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElGEPPGrowthDist_s(*args)
    elif A.tag == dTag: lib.ElGEPPGrowthDist_d(*args)
    elif A.tag == cTag: lib.ElGEPPGrowthDist_c(*args)
    elif A.tag == zTag: lib.ElGEPPGrowthDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Golub/Klema/Stewart
# -------------------
lib.ElGKS_s.argtypes = \
lib.ElGKS_d.argtypes = \
lib.ElGKS_c.argtypes = \
lib.ElGKS_z.argtypes = \
lib.ElGKSDist_s.argtypes = \
lib.ElGKSDist_d.argtypes = \
lib.ElGKSDist_c.argtypes = \
lib.ElGKSDist_z.argtypes = \
  [c_void_p,iType]

def GKS(A,n):
  args = [A.obj,n]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElGKS_s(*args)
    elif A.tag == dTag: lib.ElGKS_d(*args)
    elif A.tag == cTag: lib.ElGKS_c(*args)
    elif A.tag == zTag: lib.ElGKS_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElGKSDist_s(*args)
    elif A.tag == dTag: lib.ElGKSDist_d(*args)
    elif A.tag == cTag: lib.ElGKSDist_c(*args)
    elif A.tag == zTag: lib.ElGKSDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Grcar
# -----
lib.ElGrcar_i.argtypes = \
lib.ElGrcar_s.argtypes = \
lib.ElGrcar_d.argtypes = \
lib.ElGrcar_c.argtypes = \
lib.ElGrcar_z.argtypes = \
lib.ElGrcarDist_i.argtypes = \
lib.ElGrcarDist_s.argtypes = \
lib.ElGrcarDist_d.argtypes = \
lib.ElGrcarDist_c.argtypes = \
lib.ElGrcarDist_z.argtypes = \
  [c_void_p,iType,iType]

def Grcar(A,n,k=3):
  args = [A.obj,n,k]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElGrcar_i(*args)
    elif A.tag == sTag: lib.ElGrcar_s(*args)
    elif A.tag == dTag: lib.ElGrcar_d(*args)
    elif A.tag == cTag: lib.ElGrcar_c(*args)
    elif A.tag == zTag: lib.ElGrcar_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElGrcarDist_i(*args)
    elif A.tag == sTag: lib.ElGrcarDist_s(*args)
    elif A.tag == dTag: lib.ElGrcarDist_d(*args)
    elif A.tag == cTag: lib.ElGrcarDist_c(*args)
    elif A.tag == zTag: lib.ElGrcarDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Haar
# ----
lib.ElHaar_s.argtypes = \
lib.ElHaar_d.argtypes = \
lib.ElHaar_c.argtypes = \
lib.ElHaar_z.argtypes = \
lib.ElHaarDist_s.argtypes = \
lib.ElHaarDist_d.argtypes = \
lib.ElHaarDist_c.argtypes = \
lib.ElHaarDist_z.argtypes = \
  [c_void_p,iType]

def Haar(A,n):
  args = [A.obj,n]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHaar_s(*args)
    elif A.tag == dTag: lib.ElHaar_d(*args)
    elif A.tag == cTag: lib.ElHaar_c(*args)
    elif A.tag == zTag: lib.ElHaar_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHaarDist_s(*args)
    elif A.tag == dTag: lib.ElHaarDist_d(*args)
    elif A.tag == cTag: lib.ElHaarDist_c(*args)
    elif A.tag == zTag: lib.ElHaarDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElImplicitHaar_s.argtypes = \
lib.ElImplicitHaar_d.argtypes = \
lib.ElImplicitHaar_c.argtypes = \
lib.ElImplicitHaar_z.argtypes = \
lib.ElImplicitHaarDist_s.argtypes = \
lib.ElImplicitHaarDist_d.argtypes = \
lib.ElImplicitHaarDist_c.argtypes = \
lib.ElImplicitHaarDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,iType]

def ImplicitHaar(A,n):
  if type(A) is Matrix:
    t = Matrix(A.tag)
    d = Matrix(Base(A.tag))
    args = [A.obj,t.obj,d.obj,n]
    if   A.tag == sTag: lib.ElImplicitHaar_s(*args)
    elif A.tag == dTag: lib.ElImplicitHaar_d(*args)
    elif A.tag == cTag: lib.ElImplicitHaar_c(*args)
    elif A.tag == zTag: lib.ElImplicitHaar_z(*args)
    else: DataExcept()
    return t, d
  elif type(A) is DistMatrix:
    t = DistMatrix(A.tag,MC,STAR,A.Grid())
    d = DistMatrix(Base(A.tag),MC,STAR,A.Grid())
    args = [A.obj,t.obj,d.obj,n]
    if   A.tag == sTag: lib.ElImplicitHaarDist_s(*args)
    elif A.tag == dTag: lib.ElImplicitHaarDist_d(*args)
    elif A.tag == cTag: lib.ElImplicitHaarDist_c(*args)
    elif A.tag == zTag: lib.ElImplicitHaarDist_z(*args)
    else: DataExcept()
    return t, d
  else: TypeExcept()

# Hankel
# ------
lib.ElHankel_i.argtypes = \
lib.ElHankelDist_i.argtypes = \
  [c_void_p,iType,iType,iType,POINTER(iType)]

lib.ElHankel_s.argtypes = \
lib.ElHankelDist_s.argtypes = \
  [c_void_p,iType,iType,iType,POINTER(sType)]

lib.ElHankel_d.argtypes = \
lib.ElHankelDist_d.argtypes = \
  [c_void_p,iType,iType,iType,POINTER(dType)]

lib.ElHankel_c.argtypes = \
lib.ElHankelDist_c.argtypes = \
  [c_void_p,iType,iType,iType,POINTER(cType)]

lib.ElHankel_z.argtypes = \
lib.ElHankelDist_z.argtypes = \
  [c_void_p,iType,iType,iType,POINTER(zType)]

def Hankel(A,m,n,a):
  aLen = len(a)
  aBuf = (TagToType(A.tag)*aLen)(*a)
  args = [A.obj,m,n,aLen,aBuf]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElHankel_i(*args)
    elif A.tag == sTag: lib.ElHankel_s(*args)
    elif A.tag == dTag: lib.ElHankel_d(*args)
    elif A.tag == cTag: lib.ElHankel_c(*args)
    elif A.tag == zTag: lib.ElHankel_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElHankelDist_i(*args)
    elif A.tag == sTag: lib.ElHankelDist_s(*args)
    elif A.tag == dTag: lib.ElHankelDist_d(*args)
    elif A.tag == cTag: lib.ElHankelDist_c(*args)
    elif A.tag == zTag: lib.ElHankelDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Hanowa
# ------
lib.ElHanowa_i.argtypes = \
lib.ElHanowaDist_i.argtypes = \
  [c_void_p,iType,iType]

lib.ElHanowa_s.argtypes = \
lib.ElHanowaDist_s.argtypes = \
  [c_void_p,iType,sType]

lib.ElHanowa_d.argtypes = \
lib.ElHanowaDist_d.argtypes = \
  [c_void_p,iType,dType]

lib.ElHanowa_c.argtypes = \
lib.ElHanowaDist_c.argtypes = \
  [c_void_p,iType,cType]

lib.ElHanowa_z.argtypes = \
lib.ElHanowaDist_z.argtypes = \
  [c_void_p,iType,zType]

def Hanowa(A,n,mu):
  args = [A.obj,n,mu]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElHanowa_i(*args)
    elif A.tag == sTag: lib.ElHanowa_s(*args)
    elif A.tag == dTag: lib.ElHanowa_d(*args)
    elif A.tag == cTag: lib.ElHanowa_c(*args)
    elif A.tag == zTag: lib.ElHanowa_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElHanowaDist_i(*args)
    elif A.tag == sTag: lib.ElHanowaDist_s(*args)
    elif A.tag == dTag: lib.ElHanowaDist_d(*args)
    elif A.tag == cTag: lib.ElHanowaDist_c(*args)
    elif A.tag == zTag: lib.ElHanowaDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Hatano-Nelson
# -------------
lib.ElHatanoNelson_s.argtypes = \
lib.ElHatanoNelsonDist_s.argtypes = \
  [c_void_p,iType,sType,sType,sType,bType]

lib.ElHatanoNelson_d.argtypes = \
lib.ElHatanoNelsonDist_d.argtypes = \
  [c_void_p,iType,dType,dType,dType,bType]

lib.ElHatanoNelson_c.argtypes = \
lib.ElHatanoNelsonDist_c.argtypes = \
  [c_void_p,iType,cType,sType,cType,bType]

lib.ElHatanoNelson_z.argtypes = \
lib.ElHatanoNelsonDist_z.argtypes = \
  [c_void_p,iType,zType,dType,zType,bType]

def HatanoNelson(A,n,center,radius,g,periodic=True):
  args = [A.obj,n,center,radius,g,periodic]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHatanoNelson_s(*args)
    elif A.tag == dTag: lib.ElHatanoNelson_d(*args)
    elif A.tag == cTag: lib.ElHatanoNelson_c(*args)
    elif A.tag == zTag: lib.ElHatanoNelson_z(*args)
    else: DataExcept()
  elif tyep(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHatanoNelsonDist_s(*args)
    elif A.tag == dTag: lib.ElHatanoNelsonDist_d(*args)
    elif A.tag == cTag: lib.ElHatanoNelsonDist_c(*args)
    elif A.tag == zTag: lib.ElHatanoNelsonDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Helmholtz
# ---------
lib.ElHelmholtz1D_s.argtypes = \
lib.ElHelmholtz1DDist_s.argtypes = \
lib.ElHelmholtz1DSparse_s.argtypes = \
lib.ElHelmholtz1DDistSparse_s.argtypes = \
  [c_void_p,iType,sType]

lib.ElHelmholtz1D_d.argtypes = \
lib.ElHelmholtz1DDist_d.argtypes = \
lib.ElHelmholtz1DSparse_d.argtypes = \
lib.ElHelmholtz1DDistSparse_d.argtypes = \
  [c_void_p,iType,dType]

lib.ElHelmholtz1D_c.argtypes = \
lib.ElHelmholtz1DDist_c.argtypes = \
lib.ElHelmholtz1DSparse_c.argtypes = \
lib.ElHelmholtz1DDistSparse_c.argtypes = \
  [c_void_p,iType,cType]

lib.ElHelmholtz1D_z.argtypes = \
lib.ElHelmholtz1DDist_z.argtypes = \
lib.ElHelmholtz1DSparse_z.argtypes = \
lib.ElHelmholtz1DDistSparse_z.argtypes = \
  [c_void_p,iType,zType]

def Helmholtz1D(H,nx,shift):
  args = [H.obj,nx,shift]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHelmholtz1D_s(*args)
    elif A.tag == dTag: lib.ElHelmholtz1D_d(*args)
    elif A.tag == cTag: lib.ElHelmholtz1D_c(*args)
    elif A.tag == zTag: lib.ElHelmholtz1D_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHelmholtz1DDist_s(*args)
    elif A.tag == dTag: lib.ElHelmholtz1DDist_d(*args)
    elif A.tag == cTag: lib.ElHelmholtz1DDist_c(*args)
    elif A.tag == zTag: lib.ElHelmholtz1DDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == sTag: lib.ElHelmholtz1DSparse_s(*args)
    elif A.tag == dTag: lib.ElHelmholtz1DSparse_d(*args)
    elif A.tag == cTag: lib.ElHelmholtz1DSparse_c(*args)
    elif A.tag == zTag: lib.ElHelmholtz1DSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == sTag: lib.ElHelmholtz1DDistSparse_s(*args)
    elif A.tag == dTag: lib.ElHelmholtz1DDistSparse_d(*args)
    elif A.tag == cTag: lib.ElHelmholtz1DDistSparse_c(*args)
    elif A.tag == zTag: lib.ElHelmholtz1DDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElHelmholtz2D_s.argtypes = \
lib.ElHelmholtz2DDist_s.argtypes = \
lib.ElHelmholtz2DSparse_s.argtypes = \
lib.ElHelmholtz2DDistSparse_s.argtypes = \
  [c_void_p,iType,iType,sType]

lib.ElHelmholtz2D_d.argtypes = \
lib.ElHelmholtz2DDist_d.argtypes = \
lib.ElHelmholtz2DSparse_d.argtypes = \
lib.ElHelmholtz2DDistSparse_d.argtypes = \
  [c_void_p,iType,iType,dType]

lib.ElHelmholtz2D_c.argtypes = \
lib.ElHelmholtz2DDist_c.argtypes = \
lib.ElHelmholtz2DSparse_c.argtypes = \
lib.ElHelmholtz2DDistSparse_c.argtypes = \
  [c_void_p,iType,iType,cType]

lib.ElHelmholtz2D_z.argtypes = \
lib.ElHelmholtz2DDist_z.argtypes = \
lib.ElHelmholtz2DSparse_z.argtypes = \
lib.ElHelmholtz2DDistSparse_z.argtypes = \
  [c_void_p,iType,iType,zType]

def Helmholtz2D(H,nx,ny,shift):
  args = [H.obj,nx,ny,shift]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHelmholtz2D_s(*args)
    elif A.tag == dTag: lib.ElHelmholtz2D_d(*args)
    elif A.tag == cTag: lib.ElHelmholtz2D_c(*args)
    elif A.tag == zTag: lib.ElHelmholtz2D_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHelmholtz2DDist_s(*args)
    elif A.tag == dTag: lib.ElHelmholtz2DDist_d(*args)
    elif A.tag == cTag: lib.ElHelmholtz2DDist_c(*args)
    elif A.tag == zTag: lib.ElHelmholtz2DDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == sTag: lib.ElHelmholtz2DSparse_s(*args)
    elif A.tag == dTag: lib.ElHelmholtz2DSparse_d(*args)
    elif A.tag == cTag: lib.ElHelmholtz2DSparse_c(*args)
    elif A.tag == zTag: lib.ElHelmholtz2DSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == sTag: lib.ElHelmholtz2DSparse_s(*args)
    elif A.tag == dTag: lib.ElHelmholtz2DSparse_d(*args)
    elif A.tag == cTag: lib.ElHelmholtz2DSparse_c(*args)
    elif A.tag == zTag: lib.ElHelmholtz2DSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElHelmholtz3D_s.argtypes = \
lib.ElHelmholtz3DDist_s.argtypes = \
  [c_void_p,iType,iType,iType,sType]

lib.ElHelmholtz3D_d.argtypes = \
lib.ElHelmholtz3DDist_d.argtypes = \
  [c_void_p,iType,iType,iType,dType]

lib.ElHelmholtz3D_c.argtypes = \
lib.ElHelmholtz3DDist_c.argtypes = \
  [c_void_p,iType,iType,iType,cType]

lib.ElHelmholtz3D_z.argtypes = \
lib.ElHelmholtz3DDist_z.argtypes = \
  [c_void_p,iType,iType,iType,zType]

def Helmholtz3D(H,nx,ny,nz,shift):
  args = [H.obj,nx,ny,nz,shift]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHelmholtz3D_s(*args)
    elif A.tag == dTag: lib.ElHelmholtz3D_d(*args)
    elif A.tag == cTag: lib.ElHelmholtz3D_c(*args)
    elif A.tag == zTag: lib.ElHelmholtz3D_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHelmholtz3DDist_s(*args)
    elif A.tag == dTag: lib.ElHelmholtz3DDist_d(*args)
    elif A.tag == cTag: lib.ElHelmholtz3DDist_c(*args)
    elif A.tag == zTag: lib.ElHelmholtz3DDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == sTag: lib.ElHelmholtz3DSparse_s(*args)
    elif A.tag == dTag: lib.ElHelmholtz3DSparse_d(*args)
    elif A.tag == cTag: lib.ElHelmholtz3DSparse_c(*args)
    elif A.tag == zTag: lib.ElHelmholtz3DSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == sTag: lib.ElHelmholtz3DDistSparse_s(*args)
    elif A.tag == dTag: lib.ElHelmholtz3DDistSparse_d(*args)
    elif A.tag == cTag: lib.ElHelmholtz3DDistSparse_c(*args)
    elif A.tag == zTag: lib.ElHelmholtz3DDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Helmholtz with PML
# ------------------
lib.ElHelmholtzPML1D_c.argtypes = \
lib.ElHelmholtzPML1DDist_c.argtypes = \
lib.ElHelmholtzPML1DSparse_c.argtypes = \
lib.ElHelmholtzPML1DDistSparse_c.argtypes = \
  [c_void_p,iType,cType,iType,sType,sType]

lib.ElHelmholtzPML1D_z.argtypes = \
lib.ElHelmholtzPML1DDist_z.argtypes = \
lib.ElHelmholtzPML1DSparse_z.argtypes = \
lib.ElHelmholtzPML1DDistSparse_z.argtypes = \
  [c_void_p,iType,zType,iType,dType,dType]

def HelmholtzPML1D(H,nx,omega,numPml,sigma,pmlExp):
  args = [H.obj,nx,omega,numPml,sigma,pmlExp]
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElHelmholtzPML1D_c(*args)
    elif A.tag == zTag: lib.ElHelmholtzPML1D_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElHelmholtzPML1DDist_c(*args)
    elif A.tag == zTag: lib.ElHelmholtzPML1DDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == cTag: lib.ElHelmholtzPML1DSparse_c(*args)
    elif A.tag == zTag: lib.ElHelmholtzPML1DSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == cTag: lib.ElHelmholtzPML1DDistSparse_c(*args)
    elif A.tag == zTag: lib.ElHelmholtzPML1DDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElHelmholtzPML2D_c.argtypes = \
lib.ElHelmholtzPML2DDist_c.argtypes = \
lib.ElHelmholtzPML2DSparse_c.argtypes = \
lib.ElHelmholtzPML2DDistSparse_c.argtypes = \
  [c_void_p,iType,iType,cType,iType,sType,sType]

lib.ElHelmholtzPML2D_z.argtypes = \
lib.ElHelmholtzPML2DDist_z.argtypes = \
lib.ElHelmholtzPML2DSparse_z.argtypes = \
lib.ElHelmholtzPML2DDistSparse_z.argtypes = \
  [c_void_p,iType,iType,zType,iType,dType,dType]

def HelmholtzPML2D(H,nx,ny,omega,numPml,sigma,pmlExp):
  args = [H.obj,nx,ny,omega,numPml,sigma,pmlExp]
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElHelmholtzPML2D_c(*args)
    elif A.tag == zTag: lib.ElHelmholtzPML2D_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElHelmholtzPML2DDist_c(*args)
    elif A.tag == zTag: lib.ElHelmholtzPML2DDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == cTag: lib.ElHelmholtzPML2DSparse_c(*args)
    elif A.tag == zTag: lib.ElHelmholtzPML2DSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == cTag: lib.ElHelmholtzPML2DDistSparse_c(*args)
    elif A.tag == zTag: lib.ElHelmholtzPML2DDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElHelmholtzPML3D_c.argtypes = \
lib.ElHelmholtzPML3DDist_c.argtypes = \
lib.ElHelmholtzPML3DSparse_c.argtypes = \
lib.ElHelmholtzPML3DDistSparse_c.argtypes = \
  [c_void_p,iType,iType,iType,cType,iType,sType,sType]

lib.ElHelmholtzPML3D_z.argtypes = \
lib.ElHelmholtzPML3DDist_z.argtypes = \
lib.ElHelmholtzPML3DSparse_z.argtypes = \
lib.ElHelmholtzPML3DDistSparse_z.argtypes = \
  [c_void_p,iType,iType,iType,zType,iType,dType,dType]

def HelmholtzPML3D(H,nx,ny,nz,omega,numPml,sigma,pmlExp):
  args = [H.obj,nx,ny,nz,omega,numPml,sigma,pmlExp]
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElHelmholtzPML3D_c(*args)
    elif A.tag == zTag: lib.ElHelmholtzPML3D_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElHelmholtzPML3DDist_c(*args)
    elif A.tag == zTag: lib.ElHelmholtzPML3DDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == cTag: lib.ElHelmholtzPML3DSparse_c(*args)
    elif A.tag == zTag: lib.ElHelmholtzPML3DSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == cTag: lib.ElHelmholtzPML3DDistSparse_c(*args)
    elif A.tag == zTag: lib.ElHelmholtzPML3DDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Hermitian from EVD
# ------------------
lib.ElHermitianFromEVD_s.argtypes = \
lib.ElHermitianFromEVD_d.argtypes = \
lib.ElHermitianFromEVD_c.argtypes = \
lib.ElHermitianFromEVD_z.argtypes = \
lib.ElHermitianFromEVDDist_s.argtypes = \
lib.ElHermitianFromEVDDist_d.argtypes = \
lib.ElHermitianFromEVDDist_c.argtypes = \
lib.ElHermitianFromEVDDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p]

def HermitianFromEVD(uplo,A,w,Z):
  if type(A) is not type(w) or type(w) is not type(Z):
    raise Exception('Types of {A,w,Z} must match')
  if A.tag != Z.tag:
    raise Exception('Datatypes of A and Z must match')
  if w.tag != Base(Z.tag):
    raise Exception('w must be of the base datatype of Z')
  args = [uplo,A.obj,w.obj,Z.obj]
  if type(Z) is Matrix:
    if   Z.tag == sTag: lib.ElHermitianFromEVD_s(*args)
    elif Z.tag == dTag: lib.ElHermitianFromEVD_d(*args)
    elif Z.tag == cTag: lib.ElHermitianFromEVD_c(*args)
    elif Z.tag == zTag: lib.ElHermitianFromEVD_z(*args)
    else: DataExcept()
  elif type(Z) is DistMatrix:
    if   Z.tag == sTag: lib.ElHermitianFromEVDDist_s(*args)
    elif Z.tag == dTag: lib.ElHermitianFromEVDDist_d(*args)
    elif Z.tag == cTag: lib.ElHermitianFromEVDDist_c(*args)
    elif Z.tag == zTag: lib.ElHermitianFromEVDDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Hermitian uniform spectrum
# --------------------------
lib.ElHermitianUniformSpectrum_s.argtypes = \
lib.ElHermitianUniformSpectrum_c.argtypes = \
lib.ElHermitianUniformSpectrumDist_s.argtypes = \
lib.ElHermitianUniformSpectrumDist_c.argtypes = \
  [c_void_p,iType,sType,sType]

lib.ElHermitianUniformSpectrum_d.argtypes = \
lib.ElHermitianUniformSpectrum_z.argtypes = \
lib.ElHermitianUniformSpectrumDist_d.argtypes = \
lib.ElHermitianUniformSpectrumDist_z.argtypes = \
  [c_void_p,iType,dType,dType]

def HermitianUniformSpectrum(A,n,lower=0,upper=1):
  args = [A.obj,n,lower,upper]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHermitianUniformSpectrum_s(*args)
    elif A.tag == dTag: lib.ElHermitianUniformSpectrum_d(*args)
    elif A.tag == cTag: lib.ElHermitianUniformSpectrum_c(*args)
    elif A.tag == zTag: lib.ElHermitianUniformSpectrum_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHermitianUniformSpectrumDist_s(*args)
    elif A.tag == dTag: lib.ElHermitianUniformSpectrumDist_d(*args)
    elif A.tag == cTag: lib.ElHermitianUniformSpectrumDist_c(*args)
    elif A.tag == zTag: lib.ElHermitianUniformSpectrumDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Hilbert
# -------
lib.ElHilbert_s.argtypes = \
lib.ElHilbert_d.argtypes = \
lib.ElHilbert_c.argtypes = \
lib.ElHilbert_z.argtypes = \
lib.ElHilbertDist_s.argtypes = \
lib.ElHilbertDist_d.argtypes = \
lib.ElHilbertDist_c.argtypes = \
lib.ElHilbertDist_z.argtypes = \
  [c_void_p,iType]

def Hilbert(A,n):
  args = [A.obj,n]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHilbert_s(*args)
    elif A.tag == dTag: lib.ElHilbert_d(*args)
    elif A.tag == cTag: lib.ElHilbert_c(*args)
    elif A.tag == zTag: lib.ElHilbert_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHilbertDist_s(*args)
    elif A.tag == dTag: lib.ElHilbertDist_d(*args)
    elif A.tag == cTag: lib.ElHilbertDist_c(*args)
    elif A.tag == zTag: lib.ElHilbertDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Identity
# --------
lib.ElIdentity_i.argtypes = \
lib.ElIdentity_s.argtypes = \
lib.ElIdentity_d.argtypes = \
lib.ElIdentity_c.argtypes = \
lib.ElIdentity_z.argtypes = \
lib.ElIdentityDist_i.argtypes = \
lib.ElIdentityDist_s.argtypes = \
lib.ElIdentityDist_d.argtypes = \
lib.ElIdentityDist_c.argtypes = \
lib.ElIdentityDist_z.argtypes = \
lib.ElIdentitySparse_i.argtypes = \
lib.ElIdentitySparse_s.argtypes = \
lib.ElIdentitySparse_d.argtypes = \
lib.ElIdentitySparse_c.argtypes = \
lib.ElIdentitySparse_z.argtypes = \
lib.ElIdentityDistSparse_i.argtypes = \
lib.ElIdentityDistSparse_s.argtypes = \
lib.ElIdentityDistSparse_d.argtypes = \
lib.ElIdentityDistSparse_c.argtypes = \
lib.ElIdentityDistSparse_z.argtypes = \
  [c_void_p,iType,iType]

def Identity(A,m,n):
  args = [A.obj,m,n]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElIdentity_i(*args)
    elif A.tag == sTag: lib.ElIdentity_s(*args)
    elif A.tag == dTag: lib.ElIdentity_d(*args)
    elif A.tag == cTag: lib.ElIdentity_c(*args)
    elif A.tag == zTag: lib.ElIdentity_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElIdentityDist_i(*args)
    elif A.tag == sTag: lib.ElIdentityDist_s(*args)
    elif A.tag == dTag: lib.ElIdentityDist_d(*args)
    elif A.tag == cTag: lib.ElIdentityDist_c(*args)
    elif A.tag == zTag: lib.ElIdentityDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == iTag: lib.ElIdentitySparse_i(*args)
    elif A.tag == sTag: lib.ElIdentitySparse_s(*args)
    elif A.tag == dTag: lib.ElIdentitySparse_d(*args)
    elif A.tag == cTag: lib.ElIdentitySparse_c(*args)
    elif A.tag == zTag: lib.ElIdentitySparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == iTag: lib.ElIdentityDistSparse_i(*args)
    elif A.tag == sTag: lib.ElIdentityDistSparse_s(*args)
    elif A.tag == dTag: lib.ElIdentityDistSparse_d(*args)
    elif A.tag == cTag: lib.ElIdentityDistSparse_c(*args)
    elif A.tag == zTag: lib.ElIdentityDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Jordan
# ------
lib.ElJordan_i.argtypes = \
lib.ElJordanDist_i.argtypes = \
  [c_void_p,iType,iType]

lib.ElJordan_s.argtypes = \
lib.ElJordanDist_s.argtypes = \
  [c_void_p,iType,sType]

lib.ElJordan_d.argtypes = \
lib.ElJordanDist_d.argtypes = \
  [c_void_p,iType,dType]

lib.ElJordan_c.argtypes = \
lib.ElJordanDist_c.argtypes = \
  [c_void_p,iType,cType]

lib.ElJordan_z.argtypes = \
lib.ElJordanDist_z.argtypes = \
  [c_void_p,iType,zType]

def Jordan(J,n,lambPre):
  lamb = TagToType(J.tag)(lambPre)
  args = [J.obj,n,lamb]
  if type(J) is Matrix: 
    if   J.tag == iTag: lib.ElJordan_i(*args)
    elif J.tag == sTag: lib.ElJordan_s(*args)
    elif J.tag == dTag: lib.ElJordan_d(*args)
    elif J.tag == cTag: lib.ElJordan_c(*args)
    elif J.tag == zTag: lib.ElJordan_z(*args)
    else: DataExcept()
  elif type(J) is DistMatrix:
    if   J.tag == iTag: lib.ElJordanDist_i(*args)
    elif J.tag == sTag: lib.ElJordanDist_s(*args)
    elif J.tag == dTag: lib.ElJordanDist_d(*args)
    elif J.tag == cTag: lib.ElJordanDist_c(*args)
    elif J.tag == zTag: lib.ElJordanDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Jordan-Cholesky
# ---------------
lib.ElJordanCholesky_s.argtypes = \
lib.ElJordanCholesky_d.argtypes = \
lib.ElJordanCholesky_c.argtypes = \
lib.ElJordanCholesky_z.argtypes = \
lib.ElJordanCholeskyDist_s.argtypes = \
lib.ElJordanCholeskyDist_d.argtypes = \
lib.ElJordanCholeskyDist_c.argtypes = \
lib.ElJordanCholeskyDist_z.argtypes = \
lib.ElJordanCholeskySparse_s.argtypes = \
lib.ElJordanCholeskySparse_d.argtypes = \
lib.ElJordanCholeskySparse_c.argtypes = \
lib.ElJordanCholeskySparse_z.argtypes = \
lib.ElJordanCholeskyDistSparse_s.argtypes = \
lib.ElJordanCholeskyDistSparse_d.argtypes = \
lib.ElJordanCholeskyDistSparse_c.argtypes = \
lib.ElJordanCholeskyDistSparse_z.argtypes = \
  [c_void_p,iType]

def JordanCholesky(A,n):
  args = [A.obj,n]
  if type(A) is Matrix: 
    if   A.tag == sTag: lib.ElJordanCholesky_s(*args)
    elif A.tag == dTag: lib.ElJordanCholesky_d(*args)
    elif A.tag == cTag: lib.ElJordanCholesky_c(*args)
    elif A.tag == zTag: lib.ElJordanCholesky_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix: 
    if   A.tag == sTag: lib.ElJordanCholeskyDist_s(*args)
    elif A.tag == dTag: lib.ElJordanCholeskyDist_d(*args)
    elif A.tag == cTag: lib.ElJordanCholeskyDist_c(*args)
    elif A.tag == zTag: lib.ElJordanCholeskyDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix: 
    if   A.tag == sTag: lib.ElJordanCholeskySparse_s(*args)
    elif A.tag == dTag: lib.ElJordanCholeskySparse_d(*args)
    elif A.tag == cTag: lib.ElJordanCholeskySparse_c(*args)
    elif A.tag == zTag: lib.ElJordanCholeskySparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix: 
    if   A.tag == sTag: lib.ElJordanCholeskyDistSparse_s(*args)
    elif A.tag == dTag: lib.ElJordanCholeskyDistSparse_d(*args)
    elif A.tag == cTag: lib.ElJordanCholeskyDistSparse_c(*args)
    elif A.tag == zTag: lib.ElJordanCholeskyDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Kahan
# -----
lib.ElKahan_s.argtypes = \
lib.ElKahanDist_s.argtypes = \
  [c_void_p,iType,sType]

lib.ElKahan_d.argtypes = \
lib.ElKahanDist_d.argtypes = \
  [c_void_p,iType,dType]

lib.ElKahan_c.argtypes = \
lib.ElKahanDist_c.argtypes = \
  [c_void_p,iType,cType]

lib.ElKahan_z.argtypes = \
lib.ElKahanDist_z.argtypes = \
  [c_void_p,iType,zType]

def Kahan(A,n,phi):
  args = [A.obj,n,phi]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElKahan_s(*args)
    elif A.tag == dTag: lib.ElKahan_d(*args)
    elif A.tag == cTag: lib.ElKahan_c(*args)
    elif A.tag == zTag: lib.ElKahan_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElKahanDist_s(*args)
    elif A.tag == dTag: lib.ElKahanDist_d(*args)
    elif A.tag == cTag: lib.ElKahanDist_c(*args)
    elif A.tag == zTag: lib.ElKahanDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# KMS
# ---
lib.ElKMS_i.argtypes = \
lib.ElKMSDist_i.argtypes = \
  [c_void_p,iType,iType]

lib.ElKMS_s.argtypes = \
lib.ElKMSDist_s.argtypes = \
  [c_void_p,iType,sType]

lib.ElKMS_d.argtypes = \
lib.ElKMSDist_d.argtypes = \
  [c_void_p,iType,dType]

lib.ElKMS_c.argtypes = \
lib.ElKMSDist_c.argtypes = \
  [c_void_p,iType,cType]

lib.ElKMS_z.argtypes = \
lib.ElKMSDist_z.argtypes = \
  [c_void_p,iType,zType]

def KMS(K,n,rho):
  args = [K.obj,n,rho]
  if type(K) is Matrix:
    if   K.tag == iTag: lib.ElKMS_i(*args)
    elif K.tag == sTag: lib.ElKMS_s(*args)
    elif K.tag == dTag: lib.ElKMS_d(*args)
    elif K.tag == cTag: lib.ElKMS_c(*args)
    elif K.tag == zTag: lib.ElKMS_z(*args)
    else: DataExcept()
  elif type(K) is DistMatrix:
    if   K.tag == iTag: lib.ElKMSDist_i(*args)
    elif K.tag == sTag: lib.ElKMSDist_s(*args)
    elif K.tag == dTag: lib.ElKMSDist_d(*args)
    elif K.tag == cTag: lib.ElKMSDist_c(*args)
    elif K.tag == zTag: lib.ElKMSDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Laplacian
# ---------
lib.ElLaplacian1D_s.argtypes = \
lib.ElLaplacian1D_d.argtypes = \
lib.ElLaplacian1D_c.argtypes = \
lib.ElLaplacian1D_z.argtypes = \
lib.ElLaplacian1DDist_s.argtypes = \
lib.ElLaplacian1DDist_d.argtypes = \
lib.ElLaplacian1DDist_c.argtypes = \
lib.ElLaplacian1DDist_z.argtypes = \
lib.ElLaplacian1DSparse_s.argtypes = \
lib.ElLaplacian1DSparse_d.argtypes = \
lib.ElLaplacian1DSparse_c.argtypes = \
lib.ElLaplacian1DSparse_z.argtypes = \
lib.ElLaplacian1DDistSparse_s.argtypes = \
lib.ElLaplacian1DDistSparse_d.argtypes = \
lib.ElLaplacian1DDistSparse_c.argtypes = \
lib.ElLaplacian1DDistSparse_z.argtypes = \
  [c_void_p,iType]

def Laplacian1D(L,nx):
  args = [L.obj,nx]
  if type(L) is Matrix:
    if   L.tag == sTag: lib.ElLaplacian1D_s(*args)
    elif L.tag == dTag: lib.ElLaplacian1D_d(*args)
    elif L.tag == cTag: lib.ElLaplacian1D_c(*args)
    elif L.tag == zTag: lib.ElLaplacian1D_z(*args)
    else: DataExcept()
  elif type(L) is DistMatrix:
    if   L.tag == sTag: lib.ElLaplacian1DDist_s(*args)
    elif L.tag == dTag: lib.ElLaplacian1DDist_d(*args)
    elif L.tag == cTag: lib.ElLaplacian1DDist_c(*args)
    elif L.tag == zTag: lib.ElLaplacian1DDist_z(*args)
    else: DataExcept()
  elif type(L) is SparseMatrix:
    if   L.tag == sTag: lib.ElLaplacian1DSparse_s(*args)
    elif L.tag == dTag: lib.ElLaplacian1DSparse_d(*args)
    elif L.tag == cTag: lib.ElLaplacian1DSparse_c(*args)
    elif L.tag == zTag: lib.ElLaplacian1DSparse_z(*args)
    else: DataExcept()
  elif type(L) is DistSparseMatrix:
    if   L.tag == sTag: lib.ElLaplacian1DDistSparse_s(*args)
    elif L.tag == dTag: lib.ElLaplacian1DDistSparse_d(*args)
    elif L.tag == cTag: lib.ElLaplacian1DDistSparse_c(*args)
    elif L.tag == zTag: lib.ElLaplacian1DDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

# LEFT OFF HERE (TODO: Add sparse wrappers)
lib.ElLaplacian2D_s.argtypes = \
lib.ElLaplacian2D_d.argtypes = \
lib.ElLaplacian2D_c.argtypes = \
lib.ElLaplacian2D_z.argtypes = \
lib.ElLaplacian2DDist_s.argtypes = \
lib.ElLaplacian2DDist_d.argtypes = \
lib.ElLaplacian2DDist_c.argtypes = \
lib.ElLaplacian2DDist_z.argtypes = \
  [c_void_p,iType,iType]

def Laplacian2D(L,nx,ny):
  args = [L.obj,nx,ny]
  if type(L) is Matrix:
    if   L.tag == sTag: lib.ElLaplacian2D_s(*args)
    elif L.tag == dTag: lib.ElLaplacian2D_d(*args)
    elif L.tag == cTag: lib.ElLaplacian2D_c(*args)
    elif L.tag == zTag: lib.ElLaplacian2D_z(*args)
    else: DataExcept()
  elif type(L) is DistMatrix:
    if   L.tag == sTag: lib.ElLaplacian2DDist_s(*args)
    elif L.tag == dTag: lib.ElLaplacian2DDist_d(*args)
    elif L.tag == cTag: lib.ElLaplacian2DDist_c(*args)
    elif L.tag == zTag: lib.ElLaplacian2DDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElLaplacian3D_s.argtypes = \
lib.ElLaplacian3D_d.argtypes = \
lib.ElLaplacian3D_c.argtypes = \
lib.ElLaplacian3D_z.argtypes = \
lib.ElLaplacian3DDist_s.argtypes = \
lib.ElLaplacian3DDist_d.argtypes = \
lib.ElLaplacian3DDist_c.argtypes = \
lib.ElLaplacian3DDist_z.argtypes = \
  [c_void_p,iType,iType,iType]

def Laplacian3D(L,nx,ny,nz):
  args = [L.obj,nx,ny,nz]
  if type(L) is Matrix:
    if   L.tag == sTag: lib.ElLaplacian3D_s(*args)
    elif L.tag == dTag: lib.ElLaplacian3D_d(*args)
    elif L.tag == cTag: lib.ElLaplacian3D_c(*args)
    elif L.tag == zTag: lib.ElLaplacian3D_z(*args)
    else: DataExcept()
  elif type(L) is DistMatrix:
    if   L.tag == sTag: lib.ElLaplacian3DDist_s(*args)
    elif L.tag == dTag: lib.ElLaplacian3DDist_d(*args)
    elif L.tag == cTag: lib.ElLaplacian3DDist_c(*args)
    elif L.tag == zTag: lib.ElLaplacian3DDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Lauchli
# -------
lib.ElLauchli_i.argtypes = \
lib.ElLauchliDist_i.argtypes = \
  [c_void_p,iType,iType]

lib.ElLauchli_s.argtypes = \
lib.ElLauchliDist_s.argtypes = \
  [c_void_p,iType,sType]

lib.ElLauchli_d.argtypes = \
lib.ElLauchliDist_d.argtypes = \
  [c_void_p,iType,dType]

lib.ElLauchli_c.argtypes = \
lib.ElLauchliDist_c.argtypes = \
  [c_void_p,iType,cType]

lib.ElLauchli_z.argtypes = \
lib.ElLauchliDist_z.argtypes = \
  [c_void_p,iType,zType]

def Lauchli(A,n,mu):
  args = [A.obj,n,mu]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElLauchli_i(*args)
    elif A.tag == sTag: lib.ElLauchli_s(*args)
    elif A.tag == dTag: lib.ElLauchli_d(*args)
    elif A.tag == cTag: lib.ElLauchli_c(*args)
    elif A.tag == zTag: lib.ElLauchli_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElLauchliDist_i(*args)
    elif A.tag == sTag: lib.ElLauchliDist_s(*args)
    elif A.tag == dTag: lib.ElLauchliDist_d(*args)
    elif A.tag == cTag: lib.ElLauchliDist_c(*args)
    elif A.tag == zTag: lib.ElLauchliDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Legendre
# --------
lib.ElLegendre_s.argtypes = \
lib.ElLegendre_d.argtypes = \
lib.ElLegendre_c.argtypes = \
lib.ElLegendre_z.argtypes = \
lib.ElLegendreDist_s.argtypes = \
lib.ElLegendreDist_d.argtypes = \
lib.ElLegendreDist_c.argtypes = \
lib.ElLegendreDist_z.argtypes = \
  [c_void_p,iType]

def Legendre(A,n):
  args = [A.obj,n]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElLegendre_s(*args)
    elif A.tag == dTag: lib.ElLegendre_d(*args)
    elif A.tag == cTag: lib.ElLegendre_c(*args)
    elif A.tag == zTag: lib.ElLegendre_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElLegendreDist_s(*args)
    elif A.tag == dTag: lib.ElLegendreDist_d(*args)
    elif A.tag == cTag: lib.ElLegendreDist_c(*args)
    elif A.tag == zTag: lib.ElLegendreDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Lehmer
# ------
lib.ElLehmer_s.argtypes = \
lib.ElLehmer_d.argtypes = \
lib.ElLehmer_c.argtypes = \
lib.ElLehmer_z.argtypes = \
lib.ElLehmerDist_s.argtypes = \
lib.ElLehmerDist_d.argtypes = \
lib.ElLehmerDist_c.argtypes = \
lib.ElLehmerDist_z.argtypes = \
  [c_void_p,iType]

def Lehmer(L,n):
  args = [L.obj,n]
  if type(L) is Matrix:
    if   L.tag == sTag: lib.ElLehmer_s(*args)
    elif L.tag == dTag: lib.ElLehmer_d(*args)
    elif L.tag == cTag: lib.ElLehmer_c(*args)
    elif L.tag == zTag: lib.ElLehmer_z(*args)
    else: DataExcept()
  elif type(L) is DistMatrix:
    if   L.tag == sTag: lib.ElLehmerDist_s(*args)
    elif L.tag == dTag: lib.ElLehmerDist_d(*args)
    elif L.tag == cTag: lib.ElLehmerDist_c(*args)
    elif L.tag == zTag: lib.ElLehmerDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Lotkin
# ------
lib.ElLotkin_s.argtypes = \
lib.ElLotkin_d.argtypes = \
lib.ElLotkin_c.argtypes = \
lib.ElLotkin_z.argtypes = \
lib.ElLotkinDist_s.argtypes = \
lib.ElLotkinDist_d.argtypes = \
lib.ElLotkinDist_c.argtypes = \
lib.ElLotkinDist_z.argtypes = \
  [c_void_p,iType]

def Lotkin(A,n):
  args = [A.obj,n]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElLotkin_s(*args)
    elif A.tag == dTag: lib.ElLotkin_d(*args)
    elif A.tag == cTag: lib.ElLotkin_c(*args)
    elif A.tag == zTag: lib.ElLotkin_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElLotkinDist_s(*args)
    elif A.tag == dTag: lib.ElLotkinDist_d(*args)
    elif A.tag == cTag: lib.ElLotkinDist_c(*args)
    elif A.tag == zTag: lib.ElLotkinDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# MinIJ
# -----
lib.ElMinIJ_i.argtypes = \
lib.ElMinIJ_s.argtypes = \
lib.ElMinIJ_d.argtypes = \
lib.ElMinIJ_c.argtypes = \
lib.ElMinIJ_z.argtypes = \
lib.ElMinIJDist_i.argtypes = \
lib.ElMinIJDist_s.argtypes = \
lib.ElMinIJDist_d.argtypes = \
lib.ElMinIJDist_c.argtypes = \
lib.ElMinIJDist_z.argtypes = \
  [c_void_p,iType]

def MinIJ(A,n):
  args = [A.obj,n]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElMinIJ_i(*args)
    elif A.tag == sTag: lib.ElMinIJ_s(*args)
    elif A.tag == dTag: lib.ElMinIJ_d(*args)
    elif A.tag == cTag: lib.ElMinIJ_c(*args)
    elif A.tag == zTag: lib.ElMinIJ_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElMinIJDist_i(*args)
    elif A.tag == sTag: lib.ElMinIJDist_s(*args)
    elif A.tag == dTag: lib.ElMinIJDist_d(*args)
    elif A.tag == cTag: lib.ElMinIJDist_c(*args)
    elif A.tag == zTag: lib.ElMinIJDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Normal from EVD
# ---------------
lib.ElNormalFromEVD_c.argtypes = \
lib.ElNormalFromEVD_z.argtypes = \
lib.ElNormalFromEVDDist_c.argtypes = \
lib.ElNormalFromEVDDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p]

def NormalFromEVD(A,w,Z):
  if type(A) is not type(w): raise Exception('Types of A and w must match')
  if type(A) is not type(Z): raise Exception('Types of A and Z must match')
  if Z.tag != A.tag: raise Exception('Datatypes of A and Z must match')
  if w.tag != Base(A.tag): raise Exception('Base datatype of A must match w')
  args = [A.obj,w.obj,Z.obj]
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElNormalFromEVD_c(*args)
    elif A.tag == zTag: lib.ElNormalFromEVD_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElNormalFromEVDDist_c(*args)
    elif A.tag == zTag: lib.ElNormalFromEVDDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Ones
# ----
lib.ElOnes_i.argtypes = \
lib.ElOnes_s.argtypes = \
lib.ElOnes_d.argtypes = \
lib.ElOnes_c.argtypes = \
lib.ElOnes_z.argtypes = \
lib.ElOnesDist_i.argtypes = \
lib.ElOnesDist_s.argtypes = \
lib.ElOnesDist_d.argtypes = \
lib.ElOnesDist_c.argtypes = \
lib.ElOnesDist_z.argtypes = \
lib.ElOnesDistMultiVec_i.argtypes = \
lib.ElOnesDistMultiVec_s.argtypes = \
lib.ElOnesDistMultiVec_d.argtypes = \
lib.ElOnesDistMultiVec_c.argtypes = \
lib.ElOnesDistMultiVec_z.argtypes = \
lib.ElOnesSparse_i.argtypes = \
lib.ElOnesSparse_s.argtypes = \
lib.ElOnesSparse_d.argtypes = \
lib.ElOnesSparse_c.argtypes = \
lib.ElOnesSparse_z.argtypes = \
lib.ElOnesDistSparse_i.argtypes = \
lib.ElOnesDistSparse_s.argtypes = \
lib.ElOnesDistSparse_d.argtypes = \
lib.ElOnesDistSparse_c.argtypes = \
lib.ElOnesDistSparse_z.argtypes = \
  [c_void_p,iType,iType]

def Ones(A,m,n):
  args = [A.obj,m,n]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElOnes_i(*args)
    elif A.tag == sTag: lib.ElOnes_s(*args)
    elif A.tag == dTag: lib.ElOnes_d(*args)
    elif A.tag == cTag: lib.ElOnes_c(*args)
    elif A.tag == zTag: lib.ElOnes_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElOnesDist_i(*args)
    elif A.tag == sTag: lib.ElOnesDist_s(*args)
    elif A.tag == dTag: lib.ElOnesDist_d(*args)
    elif A.tag == cTag: lib.ElOnesDist_c(*args)
    elif A.tag == zTag: lib.ElOnesDist_z(*args)
    else: DataExcept()
  elif type(A) is DistMultiVec:
    if   A.tag == iTag: lib.ElOnesDistMultiVec_i(*args)
    elif A.tag == sTag: lib.ElOnesDistMultiVec_s(*args)
    elif A.tag == dTag: lib.ElOnesDistMultiVec_d(*args)
    elif A.tag == cTag: lib.ElOnesDistMultiVec_c(*args)
    elif A.tag == zTag: lib.ElOnesDistMultiVec_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == iTag: lib.ElOnesSparse_i(*args)
    elif A.tag == sTag: lib.ElOnesSparse_s(*args)
    elif A.tag == dTag: lib.ElOnesSparse_d(*args)
    elif A.tag == cTag: lib.ElOnesSparse_c(*args)
    elif A.tag == zTag: lib.ElOnesSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == iTag: lib.ElOnesDistSparse_i(*args)
    elif A.tag == sTag: lib.ElOnesDistSparse_s(*args)
    elif A.tag == dTag: lib.ElOnesDistSparse_d(*args)
    elif A.tag == cTag: lib.ElOnesDistSparse_c(*args)
    elif A.tag == zTag: lib.ElOnesDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

# 1-2-1 matrix
# ------------
lib.ElOneTwoOne_i.argtypes = \
lib.ElOneTwoOne_s.argtypes = \
lib.ElOneTwoOne_d.argtypes = \
lib.ElOneTwoOne_c.argtypes = \
lib.ElOneTwoOne_z.argtypes = \
lib.ElOneTwoOneDist_i.argtypes = \
lib.ElOneTwoOneDist_s.argtypes = \
lib.ElOneTwoOneDist_d.argtypes = \
lib.ElOneTwoOneDist_c.argtypes = \
lib.ElOneTwoOneDist_z.argtypes = \
  [c_void_p,iType]

def OneTwoOne(A,n):
  args = [A.obj,n]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElOneTwoOne_i(*args)
    elif A.tag == sTag: lib.ElOneTwoOne_s(*args)
    elif A.tag == dTag: lib.ElOneTwoOne_d(*args)
    elif A.tag == cTag: lib.ElOneTwoOne_c(*args)
    elif A.tag == zTag: lib.ElOneTwoOne_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElOneTwoOneDist_i(*args)
    elif A.tag == sTag: lib.ElOneTwoOneDist_s(*args)
    elif A.tag == dTag: lib.ElOneTwoOneDist_d(*args)
    elif A.tag == cTag: lib.ElOneTwoOneDist_c(*args)
    elif A.tag == zTag: lib.ElOneTwoOneDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Parter
# ------
lib.ElParter_s.argtypes = \
lib.ElParter_d.argtypes = \
lib.ElParter_c.argtypes = \
lib.ElParter_z.argtypes = \
lib.ElParterDist_s.argtypes = \
lib.ElParterDist_d.argtypes = \
lib.ElParterDist_c.argtypes = \
lib.ElParterDist_z.argtypes = \
  [c_void_p,iType]

def Parter(A,n):
  args = [A.obj,n]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElParter_s(*args)
    elif A.tag == dTag: lib.ElParter_d(*args)
    elif A.tag == cTag: lib.ElParter_c(*args)
    elif A.tag == zTag: lib.ElParter_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElParterDist_s(*args)
    elif A.tag == dTag: lib.ElParterDist_d(*args)
    elif A.tag == cTag: lib.ElParterDist_c(*args)
    elif A.tag == zTag: lib.ElParterDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Pei
# ---
lib.ElPei_s.argtypes = \
lib.ElPeiDist_s.argtypes = \
  [c_void_p,iType,sType]
lib.ElPei_d.argtypes = \
lib.ElPeiDist_d.argtypes = \
  [c_void_p,iType,dType]
lib.ElPei_c.argtypes = \
lib.ElPeiDist_c.argtypes = \
  [c_void_p,iType,cType]
lib.ElPei_z.argtypes = \
lib.ElPeiDist_z.argtypes = \
  [c_void_p,iType,zType]
def Pei(A,n,alpha):
  args = [A.obj,n,alpha]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElPei_s(*args)
    elif A.tag == dTag: lib.ElPei_d(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElPeiDist_s(*args)
    elif A.tag == dTag: lib.ElPeiDist_d(*args)
    else: DataExcept()
  else: TypeExcept()

# Redheffer
# ---------
lib.ElRedheffer_i.argtypes = \
lib.ElRedheffer_s.argtypes = \
lib.ElRedheffer_d.argtypes = \
lib.ElRedheffer_c.argtypes = \
lib.ElRedheffer_z.argtypes = \
lib.ElRedhefferDist_i.argtypes = \
lib.ElRedhefferDist_s.argtypes = \
lib.ElRedhefferDist_d.argtypes = \
lib.ElRedhefferDist_c.argtypes = \
lib.ElRedhefferDist_z.argtypes = \
  [c_void_p,iType]
def Redheffer(A,n):
  args = [A.obj,n]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElRedheffer_i(*args)
    elif A.tag == sTag: lib.ElRedheffer_s(*args)
    elif A.tag == dTag: lib.ElRedheffer_d(*args)
    elif A.tag == cTag: lib.ElRedheffer_c(*args)
    elif A.tag == zTag: lib.ElRedheffer_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElRedhefferDist_i(*args)
    elif A.tag == sTag: lib.ElRedhefferDist_s(*args)
    elif A.tag == dTag: lib.ElRedhefferDist_d(*args)
    elif A.tag == cTag: lib.ElRedhefferDist_c(*args)
    elif A.tag == zTag: lib.ElRedhefferDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Riffle
# ------
lib.ElRiffle_s.argtypes = \
lib.ElRiffle_d.argtypes = \
lib.ElRiffle_c.argtypes = \
lib.ElRiffle_z.argtypes = \
lib.ElRiffleDist_s.argtypes = \
lib.ElRiffleDist_d.argtypes = \
lib.ElRiffleDist_c.argtypes = \
lib.ElRiffleDist_z.argtypes = \
  [c_void_p,iType]
def Riffle(P,n):
  args = [P.obj,n]
  if type(P) is Matrix:
    if   P.tag == sTag: lib.ElRiffle_s(*args)
    elif P.tag == dTag: lib.ElRiffle_d(*args)
    elif P.tag == cTag: lib.ElRiffle_c(*args)
    elif P.tag == zTag: lib.ElRiffle_z(*args)
    else: DataExcept()
  elif type(P) is DistMatrix:
    if   P.tag == sTag: lib.ElRiffleDist_s(*args)
    elif P.tag == dTag: lib.ElRiffleDist_d(*args)
    elif P.tag == cTag: lib.ElRiffleDist_c(*args)
    elif P.tag == zTag: lib.ElRiffleDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElRiffleStationary_s.argtypes = \
lib.ElRiffleStationary_d.argtypes = \
lib.ElRiffleStationary_c.argtypes = \
lib.ElRiffleStationary_z.argtypes = \
lib.ElRiffleStationaryDist_s.argtypes = \
lib.ElRiffleStationaryDist_d.argtypes = \
lib.ElRiffleStationaryDist_c.argtypes = \
lib.ElRiffleStationaryDist_z.argtypes = \
  [c_void_p,iType]
def RiffleStationary(P,n):
  args = [P.obj,n]
  if type(P) is Matrix:
    if   P.tag == sTag: lib.ElRiffleStationary_s(*args)
    elif P.tag == dTag: lib.ElRiffleStationary_d(*args)
    elif P.tag == cTag: lib.ElRiffleStationary_c(*args)
    elif P.tag == zTag: lib.ElRiffleStationary_z(*args)
    else: DataExcept()
  elif type(P) is DistMatrix:
    if   P.tag == sTag: lib.ElRiffleStationaryDist_s(*args)
    elif P.tag == dTag: lib.ElRiffleStationaryDist_d(*args)
    elif P.tag == cTag: lib.ElRiffleStationaryDist_c(*args)
    elif P.tag == zTag: lib.ElRiffleStationaryDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElRiffleDecay_s.argtypes = \
lib.ElRiffleDecay_d.argtypes = \
lib.ElRiffleDecay_c.argtypes = \
lib.ElRiffleDecay_z.argtypes = \
lib.ElRiffleDecayDist_s.argtypes = \
lib.ElRiffleDecayDist_d.argtypes = \
lib.ElRiffleDecayDist_c.argtypes = \
lib.ElRiffleDecayDist_z.argtypes = \
  [c_void_p,iType]
def RiffleDecay(P,n):
  args = [P.obj,n]
  if type(P) is Matrix:
    if   P.tag == sTag: lib.ElRiffleDecay_s(*args)
    elif P.tag == dTag: lib.ElRiffleDecay_d(*args)
    elif P.tag == cTag: lib.ElRiffleDecay_c(*args)
    elif P.tag == zTag: lib.ElRiffleDecay_z(*args)
    else: DataExcept()
  elif type(P) is DistMatrix:
    if   P.tag == sTag: lib.ElRiffleDecayDist_s(*args)
    elif P.tag == dTag: lib.ElRiffleDecayDist_d(*args)
    elif P.tag == cTag: lib.ElRiffleDecayDist_c(*args)
    elif P.tag == zTag: lib.ElRiffleDecayDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Ris
# ---
lib.ElRis_s.argtypes = \
lib.ElRis_d.argtypes = \
lib.ElRis_c.argtypes = \
lib.ElRis_z.argtypes = \
lib.ElRisDist_s.argtypes = \
lib.ElRisDist_d.argtypes = \
lib.ElRisDist_c.argtypes = \
lib.ElRisDist_z.argtypes = \
  [c_void_p,iType]
def Ris(A,n):
  args = [A.obj,n]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElRis_s(*args)
    elif A.tag == dTag: lib.ElRis_d(*args)
    elif A.tag == cTag: lib.ElRis_c(*args)
    elif A.tag == zTag: lib.ElRis_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElRisDist_s(*args)
    elif A.tag == dTag: lib.ElRisDist_d(*args)
    elif A.tag == cTag: lib.ElRisDist_c(*args)
    elif A.tag == zTag: lib.ElRisDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Toeplitz
# --------
lib.ElToeplitz_i.argtypes = [c_void_p,iType,iType,iType,POINTER(iType)]
lib.ElToeplitzDist_i.argtypes = [c_void_p,iType,iType,iType,POINTER(iType)]

lib.ElToeplitz_s.argtypes = [c_void_p,iType,iType,iType,POINTER(sType)]
lib.ElToeplitzDist_s.argtypes = [c_void_p,iType,iType,iType,POINTER(sType)]

lib.ElToeplitz_d.argtypes = [c_void_p,iType,iType,iType,POINTER(dType)]
lib.ElToeplitzDist_d.argtypes = [c_void_p,iType,iType,iType,POINTER(dType)]

lib.ElToeplitz_c.argtypes = [c_void_p,iType,iType,iType,POINTER(cType)]
lib.ElToeplitzDist_c.argtypes = [c_void_p,iType,iType,iType,POINTER(cType)]

lib.ElToeplitz_z.argtypes = [c_void_p,iType,iType,iType,POINTER(zType)]
lib.ElToeplitzDist_z.argtypes = [c_void_p,iType,iType,iType,POINTER(zType)]

def Toeplitz(A,m,n,a):
  aLen = len(a)
  aBuf = (TagToType(A.tag)*aLen)(*a) 
  args = [A.obj,m,n,aLen,aBuf]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElToeplitz_i(*args)
    elif A.tag == sTag: lib.ElToeplitz_s(*args)
    elif A.tag == dTag: lib.ElToeplitz_d(*args)
    elif A.tag == cTag: lib.ElToeplitz_c(*args)
    elif A.tag == zTag: lib.ElToeplitz_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElToeplitzDist_i(*args)
    elif A.tag == sTag: lib.ElToeplitzDist_s(*args)
    elif A.tag == dTag: lib.ElToeplitzDist_d(*args)
    elif A.tag == cTag: lib.ElToeplitzDist_c(*args)
    elif A.tag == zTag: lib.ElToeplitzDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Trefethen-Embree
# ----------------
lib.ElTrefethenEmbree_c.argtypes = \
lib.ElTrefethenEmbree_z.argtypes = \
lib.ElTrefethenEmbreeDist_c.argtypes = \
lib.ElTrefethenEmbreeDist_z.argtypes = \
  [c_void_p,iType]
def TrefethenEmbree(A,n):
  args = [A.obj,n]
  if type(A) is Matrix: 
    if   A.tag == cTag: lib.ElTrefethenEmbree_c(*args)
    elif A.tag == zTag: lib.ElTrefethenEmbree_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElTrefethenEmbreeDist_c(*args)
    elif A.tag == zTag: lib.ElTrefethenEmbreeDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Triangle
# --------
lib.ElTriangle_c.argtypes = \
lib.ElTriangle_z.argtypes = \
lib.ElTriangleDist_c.argtypes = \
lib.ElTriangleDist_z.argtypes = \
  [c_void_p,iType]
def Triangle(A,n):
  args = [A.obj,n]
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElTriangle_c(*args)
    elif A.tag == zTag: lib.ElTriangle_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElTriangleDist_c(*args)
    elif A.tag == zTag: lib.ElTriangleDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# TriW
# ----
lib.ElTriW_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElTriW_s.argtypes = [c_void_p,iType,sType,iType]
lib.ElTriW_d.argtypes = [c_void_p,iType,dType,iType]
lib.ElTriW_c.argtypes = [c_void_p,iType,cType,iType]
lib.ElTriW_z.argtypes = [c_void_p,iType,zType,iType]
lib.ElTriWDist_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElTriWDist_s.argtypes = [c_void_p,iType,sType,iType]
lib.ElTriWDist_d.argtypes = [c_void_p,iType,dType,iType]
lib.ElTriWDist_c.argtypes = [c_void_p,iType,cType,iType]
lib.ElTriWDist_z.argtypes = [c_void_p,iType,zType,iType]
def TriW(A,n,alpha,k):
  args = [A.obj,n,alpha,k]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElTriW_i(*args)
    elif A.tag == sTag: lib.ElTriW_s(*args)
    elif A.tag == dTag: lib.ElTriW_d(*args)
    elif A.tag == cTag: lib.ElTriW_c(*args)
    elif A.tag == zTag: lib.ElTriW_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElTriWDist_i(*args)
    elif A.tag == sTag: lib.ElTriWDist_s(*args)
    elif A.tag == dTag: lib.ElTriWDist_d(*args)
    elif A.tag == cTag: lib.ElTriWDist_c(*args)
    elif A.tag == zTag: lib.ElTriWDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Walsh
# -----
lib.ElWalsh_i.argtypes = \
lib.ElWalsh_s.argtypes = \
lib.ElWalsh_d.argtypes = \
lib.ElWalsh_c.argtypes = \
lib.ElWalsh_z.argtypes = \
lib.ElWalshDist_i.argtypes = \
lib.ElWalshDist_s.argtypes = \
lib.ElWalshDist_d.argtypes = \
lib.ElWalshDist_c.argtypes = \
lib.ElWalshDist_z.argtypes = \
  [c_void_p,iType,bType]
def Walsh(A,k,binary=False):
  args = [A.obj,k,binary]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElWalsh_i(*args)
    elif A.tag == sTag: lib.ElWalsh_s(*args)
    elif A.tag == dTag: lib.ElWalsh_d(*args)
    elif A.tag == cTag: lib.ElWalsh_c(*args)
    elif A.tag == zTag: lib.ElWalsh_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElWalshDist_i(*args)
    elif A.tag == sTag: lib.ElWalshDist_s(*args)
    elif A.tag == dTag: lib.ElWalshDist_d(*args)
    elif A.tag == cTag: lib.ElWalshDist_c(*args)
    elif A.tag == zTag: lib.ElWalshDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Walsh-Identity
# --------------
lib.ElWalshIdentity_i.argtypes = \
lib.ElWalshIdentity_s.argtypes = \
lib.ElWalshIdentity_d.argtypes = \
lib.ElWalshIdentity_c.argtypes = \
lib.ElWalshIdentity_z.argtypes = \
lib.ElWalshIdentityDist_i.argtypes = \
lib.ElWalshIdentityDist_s.argtypes = \
lib.ElWalshIdentityDist_d.argtypes = \
lib.ElWalshIdentityDist_c.argtypes = \
lib.ElWalshIdentityDist_z.argtypes = \
  [c_void_p,iType,bType]
def WalshIdentity(A,k,binary=False):
  args = [A.obj,k,binary]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElWalshIdentity_i(*args)
    elif A.tag == sTag: lib.ElWalshIdentity_s(*args)
    elif A.tag == dTag: lib.ElWalshIdentity_d(*args)
    elif A.tag == cTag: lib.ElWalshIdentity_c(*args)
    elif A.tag == zTag: lib.ElWalshIdentity_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElWalshIdentityDist_i(*args)
    elif A.tag == sTag: lib.ElWalshIdentityDist_s(*args)
    elif A.tag == dTag: lib.ElWalshIdentityDist_d(*args)
    elif A.tag == cTag: lib.ElWalshIdentityDist_c(*args)
    elif A.tag == zTag: lib.ElWalshIdentityDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Whale
# -----
lib.ElWhale_c.argtypes = \
lib.ElWhale_z.argtypes = \
lib.ElWhaleDist_c.argtypes = \
lib.ElWhaleDist_z.argtypes = \
  [c_void_p,iType]
def Whale(A,n):
  args = [A.obj,n]
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElWhale_c(*args)
    elif A.tag == zTag: lib.ElWhale_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElWhaleDist_c(*args)
    elif A.tag == zTag: lib.ElWhaleDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Wilkinson
# ---------
lib.ElWilkinson_i.argtypes = \
lib.ElWilkinson_s.argtypes = \
lib.ElWilkinson_d.argtypes = \
lib.ElWilkinson_c.argtypes = \
lib.ElWilkinson_z.argtypes = \
lib.ElWilkinsonDist_i.argtypes = \
lib.ElWilkinsonDist_s.argtypes = \
lib.ElWilkinsonDist_d.argtypes = \
lib.ElWilkinsonDist_c.argtypes = \
lib.ElWilkinsonDist_z.argtypes = \
  [c_void_p,iType]
def Wilkinson(A,k):
  args = [A.obj,k]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElWilkinson_i(*args)
    elif A.tag == sTag: lib.ElWilkinson_s(*args)
    elif A.tag == dTag: lib.ElWilkinson_d(*args)
    elif A.tag == cTag: lib.ElWilkinson_c(*args)
    elif A.tag == zTag: lib.ElWilkinson_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElWilkinsonDist_i(*args)
    elif A.tag == sTag: lib.ElWilkinsonDist_s(*args)
    elif A.tag == dTag: lib.ElWilkinsonDist_d(*args)
    elif A.tag == cTag: lib.ElWilkinsonDist_c(*args)
    elif A.tag == zTag: lib.ElWilkinsonDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Zeros
# -----
lib.ElZeros_i.argtypes = \
lib.ElZeros_s.argtypes = \
lib.ElZeros_d.argtypes = \
lib.ElZeros_c.argtypes = \
lib.ElZeros_z.argtypes = \
lib.ElZerosDist_i.argtypes = \
lib.ElZerosDist_s.argtypes = \
lib.ElZerosDist_d.argtypes = \
lib.ElZerosDist_c.argtypes = \
lib.ElZerosDist_z.argtypes = \
lib.ElZerosSparse_i.argtypes = \
lib.ElZerosSparse_s.argtypes = \
lib.ElZerosSparse_d.argtypes = \
lib.ElZerosSparse_c.argtypes = \
lib.ElZerosSparse_z.argtypes = \
lib.ElZerosDistSparse_i.argtypes = \
lib.ElZerosDistSparse_s.argtypes = \
lib.ElZerosDistSparse_d.argtypes = \
lib.ElZerosDistSparse_c.argtypes = \
lib.ElZerosDistSparse_z.argtypes = \
lib.ElZerosDistMultiVec_i.argtypes = \
lib.ElZerosDistMultiVec_s.argtypes = \
lib.ElZerosDistMultiVec_d.argtypes = \
lib.ElZerosDistMultiVec_c.argtypes = \
lib.ElZerosDistMultiVec_z.argtypes = \
  [c_void_p,iType,iType]
def Zeros(A,m,n):
  args = [A.obj,m,n]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElZeros_i(*args)
    elif A.tag == sTag: lib.ElZeros_s(*args)
    elif A.tag == dTag: lib.ElZeros_d(*args)
    elif A.tag == cTag: lib.ElZeros_c(*args)
    elif A.tag == zTag: lib.ElZeros_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElZerosDist_i(*args)
    elif A.tag == sTag: lib.ElZerosDist_s(*args)
    elif A.tag == dTag: lib.ElZerosDist_d(*args)
    elif A.tag == cTag: lib.ElZerosDist_c(*args)
    elif A.tag == zTag: lib.ElZerosDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == iTag: lib.ElZerosSparse_i(*args)
    elif A.tag == sTag: lib.ElZerosSparse_s(*args)
    elif A.tag == dTag: lib.ElZerosSparse_d(*args)
    elif A.tag == cTag: lib.ElZerosSparse_c(*args)
    elif A.tag == zTag: lib.ElZerosSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == iTag: lib.ElZerosDistSparse_i(*args)
    elif A.tag == sTag: lib.ElZerosDistSparse_s(*args)
    elif A.tag == dTag: lib.ElZerosDistSparse_d(*args)
    elif A.tag == cTag: lib.ElZerosDistSparse_c(*args)
    elif A.tag == zTag: lib.ElZerosDistSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistMultiVec:
    if   A.tag == iTag: lib.ElZerosDistMultiVec_i(*args)
    elif A.tag == sTag: lib.ElZerosDistMultiVec_s(*args)
    elif A.tag == dTag: lib.ElZerosDistMultiVec_d(*args)
    elif A.tag == cTag: lib.ElZerosDistMultiVec_c(*args)
    elif A.tag == zTag: lib.ElZerosDistMultiVec_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Random
# ======

# Bernoulli
# ---------
lib.ElBernoulli_i.argtypes = \
lib.ElBernoulli_s.argtypes = \
lib.ElBernoulli_d.argtypes = \
lib.ElBernoulli_c.argtypes = \
lib.ElBernoulli_z.argtypes = \
lib.ElBernoulliDist_i.argtypes = \
lib.ElBernoulliDist_s.argtypes = \
lib.ElBernoulliDist_d.argtypes = \
lib.ElBernoulliDist_c.argtypes = \
lib.ElBernoulliDist_z.argtypes = \
  [c_void_p,iType,iType,dType]
def Bernoulli(A,m,n,p=0.5):
  args = [A.obj,m,n,p]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElBernoulli_i(*args)
    elif A.tag == sTag: lib.ElBernoulli_s(*args)
    elif A.tag == dTag: lib.ElBernoulli_d(*args)
    elif A.tag == cTag: lib.ElBernoulli_c(*args)
    elif A.tag == zTag: lib.ElBernoulli_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElBernoulliDist_i(*args)
    elif A.tag == sTag: lib.ElBernoulliDist_s(*args)
    elif A.tag == dTag: lib.ElBernoulliDist_d(*args)
    elif A.tag == cTag: lib.ElBernoulliDist_c(*args)
    elif A.tag == zTag: lib.ElBernoulliDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Gaussian
# --------
lib.ElGaussian_s.argtypes = [c_void_p,iType,iType,sType,sType]
lib.ElGaussian_d.argtypes = [c_void_p,iType,iType,dType,dType]
lib.ElGaussian_c.argtypes = [c_void_p,iType,iType,cType,sType]
lib.ElGaussian_z.argtypes = [c_void_p,iType,iType,zType,dType]
lib.ElGaussianDist_s.argtypes = [c_void_p,iType,iType,sType,sType]
lib.ElGaussianDist_d.argtypes = [c_void_p,iType,iType,dType,dType]
lib.ElGaussianDist_c.argtypes = [c_void_p,iType,iType,cType,sType]
lib.ElGaussianDist_z.argtypes = [c_void_p,iType,iType,zType,dType]
lib.ElGaussianDistMultiVec_s.argtypes = [c_void_p,iType,iType,sType,sType]
lib.ElGaussianDistMultiVec_d.argtypes = [c_void_p,iType,iType,dType,dType]
lib.ElGaussianDistMultiVec_c.argtypes = [c_void_p,iType,iType,cType,sType]
lib.ElGaussianDistMultiVec_z.argtypes = [c_void_p,iType,iType,zType,dType]
def Gaussian(A,m,n,meanPre=0,stddev=1):
  mean = TagToType(A.tag)(meanPre)
  args = [A.obj,m,n,mean,stddev]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElGaussian_s(*args)
    elif A.tag == dTag: lib.ElGaussian_d(*args)
    elif A.tag == cTag: lib.ElGaussian_c(*args)
    elif A.tag == zTag: lib.ElGaussian_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElGaussianDist_s(*args)
    elif A.tag == dTag: lib.ElGaussianDist_d(*args)
    elif A.tag == cTag: lib.ElGaussianDist_c(*args)
    elif A.tag == zTag: lib.ElGaussianDist_z(*args)
    else: DataExcept()
  elif type(A) is DistMultiVec:
    if   A.tag == sTag: lib.ElGaussianDistMultiVec_s(*args)
    elif A.tag == dTag: lib.ElGaussianDistMultiVec_d(*args)
    elif A.tag == cTag: lib.ElGaussianDistMultiVec_c(*args)
    elif A.tag == zTag: lib.ElGaussianDistMultiVec_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Normal uniform spectrum
# -----------------------
lib.ElNormalUniformSpectrum_c.argtypes = [c_void_p,iType,cType,sType]
lib.ElNormalUniformSpectrum_z.argtypes = [c_void_p,iType,zType,dType]
lib.ElNormalUniformSpectrumDist_c.argtypes = [c_void_p,iType,cType,sType]
lib.ElNormalUniformSpectrumDist_z.argtypes = [c_void_p,iType,zType,dType]
def NormalUniformSpectrum(A,n,centerPre=0,radius=1):
  center = TagToType(A.tag)(centerPre)
  args = [A.obj,n,center,radius]
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElNormalUniformSpectrum_c(*args)
    elif A.tag == zTag: lib.ElNormalUniformSpectrum_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElNormalUniformSpectrumDist_c(*args)
    elif A.tag == zTag: lib.ElNormalUniformSpectrumDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Rademacher
# ----------
lib.ElRademacher_i.argtypes = \
lib.ElRademacher_s.argtypes = \
lib.ElRademacher_d.argtypes = \
lib.ElRademacher_c.argtypes = \
lib.ElRademacher_z.argtypes = \
lib.ElRademacherDist_i.argtypes = \
lib.ElRademacherDist_s.argtypes = \
lib.ElRademacherDist_d.argtypes = \
lib.ElRademacherDist_c.argtypes = \
lib.ElRademacherDist_z.argtypes = \
  [c_void_p,iType,iType]
def Rademacher(A,m,n):
  args = [A.obj,m,n]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElRademacher_i(*args)
    elif A.tag == sTag: lib.ElRademacher_s(*args)
    elif A.tag == dTag: lib.ElRademacher_d(*args)
    elif A.tag == cTag: lib.ElRademacher_c(*args)
    elif A.tag == zTag: lib.ElRademacher_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElRademacherDist_i(*args)
    elif A.tag == sTag: lib.ElRademacherDist_s(*args)
    elif A.tag == dTag: lib.ElRademacherDist_d(*args)
    elif A.tag == cTag: lib.ElRademacherDist_c(*args)
    elif A.tag == zTag: lib.ElRademacherDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Three-valued
# ------------
lib.ElThreeValued_i.argtypes = \
lib.ElThreeValued_s.argtypes = \
lib.ElThreeValued_d.argtypes = \
lib.ElThreeValued_c.argtypes = \
lib.ElThreeValued_z.argtypes = \
lib.ElThreeValuedDist_i.argtypes = \
lib.ElThreeValuedDist_s.argtypes = \
lib.ElThreeValuedDist_d.argtypes = \
lib.ElThreeValuedDist_c.argtypes = \
lib.ElThreeValuedDist_z.argtypes = \
  [c_void_p,iType,iType,dType]

def ThreeValued(A,m,n,p=2./3.):
  args = [A.obj,m,n,p]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElThreeValued_i(*args)
    elif A.tag == sTag: lib.ElThreeValued_s(*args)
    elif A.tag == dTag: lib.ElThreeValued_d(*args)
    elif A.tag == cTag: lib.ElThreeValued_c(*args)
    elif A.tag == zTag: lib.ElThreeValued_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElThreeValuedDist_i(*args)
    elif A.tag == sTag: lib.ElThreeValuedDist_s(*args)
    elif A.tag == dTag: lib.ElThreeValuedDist_d(*args)
    elif A.tag == cTag: lib.ElThreeValuedDist_c(*args)
    elif A.tag == zTag: lib.ElThreeValuedDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Uniform
# -------
lib.ElUniform_i.argtypes = \
lib.ElUniformDist_i.argtypes = \
lib.ElUniformDistMultiVec_i.argtypes = \
  [c_void_p,iType,iType,iType,iType]

lib.ElUniform_s.argtypes = \
lib.ElUniformDist_s.argtypes = \
lib.ElUniformDistMultiVec_s.argtypes = \
  [c_void_p,iType,iType,sType,sType]

lib.ElUniform_d.argtypes = \
lib.ElUniformDist_d.argtypes = \
lib.ElUniformDistMultiVec_d.argtypes = \
  [c_void_p,iType,iType,dType,dType]

lib.ElUniform_c.argtypes = \
lib.ElUniformDist_c.argtypes = \
lib.ElUniformDistMultiVec_c.argtypes = \
  [c_void_p,iType,iType,cType,sType]

lib.ElUniform_z.argtypes = \
lib.ElUniformDist_z.argtypes = \
lib.ElUniformDistMultiVec_z.argtypes = \
  [c_void_p,iType,iType,zType,dType]

def Uniform(A,m,n,centerPre=0,radius=1):
  center = TagToType(A.tag)(centerPre) 
  args = [A.obj,m,n,center,radius]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElUniform_i(*args)
    elif A.tag == sTag: lib.ElUniform_s(*args)
    elif A.tag == dTag: lib.ElUniform_d(*args)
    elif A.tag == cTag: lib.ElUniform_c(*args)
    elif A.tag == zTag: lib.ElUniform_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElUniformDist_i(*args)
    elif A.tag == sTag: lib.ElUniformDist_s(*args)
    elif A.tag == dTag: lib.ElUniformDist_d(*args)
    elif A.tag == cTag: lib.ElUniformDist_c(*args)
    elif A.tag == zTag: lib.ElUniformDist_z(*args)
    else: DataExcept()
  elif type(A) is DistMultiVec:
    if   A.tag == iTag: lib.ElUniformDistMultiVec_i(*args)
    elif A.tag == sTag: lib.ElUniformDistMultiVec_s(*args)
    elif A.tag == dTag: lib.ElUniformDistMultiVec_d(*args)
    elif A.tag == cTag: lib.ElUniformDistMultiVec_c(*args)
    elif A.tag == zTag: lib.ElUniformDistMultiVec_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Uniform Helmholtz Green's
# -------------------------
lib.ElUniformHelmholtzGreens_c.argtypes = \
lib.ElUniformHelmholtzGreensDist_c.argtypes = \
  [c_void_p,iType,sType]

lib.ElUniformHelmholtzGreens_z.argtypes = \
lib.ElUniformHelmholtzGreensDist_z.argtypes = \
  [c_void_p,iType,dType]

def UniformHelmholtzGreens(A,n,lamb):
  args = [A.obj,n,lamb]
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElUniformHelmholtzGreens_c(*args)
    elif A.tag == zTag: lib.ElUniformHelmholtzGreens_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElUniformHelmholtzGreensDist_c(*args)
    elif A.tag == zTag: lib.ElUniformHelmholtzGreensDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Wigner
# ------
lib.ElWigner_s.argtypes = \
lib.ElWignerDist_s.argtypes = \
  [c_void_p,iType,sType,sType]

lib.ElWigner_d.argtypes = \
lib.ElWignerDist_d.argtypes = \
  [c_void_p,iType,dType,dType]

lib.ElWigner_c.argtypes = \
lib.ElWignerDist_c.argtypes = \
  [c_void_p,iType,cType,sType]

lib.ElWigner_z.argtypes = \
lib.ElWignerDist_z.argtypes = \
  [c_void_p,iType,zType,dType]

def Wigner(A,n,meanPre=0,stddev=1):
  mean = TagToType(A.tag)(meanPre) 
  args = [A.obj,n,mean,stddev]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElWigner_s(*args)
    elif A.tag == dTag: lib.ElWigner_d(*args)
    elif A.tag == cTag: lib.ElWigner_c(*args)
    elif A.tag == zTag: lib.ElWigner_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElWignerDist_s(*args)
    elif A.tag == dTag: lib.ElWignerDist_d(*args)
    elif A.tag == cTag: lib.ElWignerDist_c(*args)
    elif A.tag == zTag: lib.ElWignerDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
