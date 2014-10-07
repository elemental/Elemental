#
#  Copyright (c) 2009-2014, Jack Poulson
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
lib.ElBullsHead_c.argtypes = [c_void_p,iType]
lib.ElBullsHead_c.restype = c_uint
lib.ElBullsHead_z.argtypes = [c_void_p,iType]
lib.ElBullsHead_z.restype = c_uint
lib.ElBullsHeadDist_c.argtypes = [c_void_p,iType]
lib.ElBullsHeadDist_c.restype = c_uint
lib.ElBullsHeadDist_z.argtypes = [c_void_p,iType]
lib.ElBullsHeadDist_z.restype = c_uint
def BullsHead(A,n):
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElBullsHead_c(A.obj,n)
    elif A.tag == zTag: lib.ElBullsHead_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElBullsHeadDist_c(A.obj,n)
    elif A.tag == zTag: lib.ElBullsHeadDist_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Cauchy
# ------
lib.ElCauchy_s.argtypes = [c_void_p,iType,POINTER(sType),iType,POINTER(sType)]
lib.ElCauchy_s.restype = c_uint
lib.ElCauchy_d.argtypes = [c_void_p,iType,POINTER(dType),iType,POINTER(dType)]
lib.ElCauchy_d.restype = c_uint
lib.ElCauchy_c.argtypes = [c_void_p,iType,POINTER(cType),iType,POINTER(cType)]
lib.ElCauchy_c.restype = c_uint
lib.ElCauchy_z.argtypes = [c_void_p,iType,POINTER(zType),iType,POINTER(zType)]
lib.ElCauchy_z.restype = c_uint
lib.ElCauchyDist_s.argtypes = \
  [c_void_p,iType,POINTER(sType),iType,POINTER(sType)]
lib.ElCauchyDist_s.restype = c_uint
lib.ElCauchyDist_d.argtypes = \
  [c_void_p,iType,POINTER(dType),iType,POINTER(dType)]
lib.ElCauchyDist_d.restype = c_uint
lib.ElCauchyDist_c.argtypes = \
  [c_void_p,iType,POINTER(cType),iType,POINTER(cType)]
lib.ElCauchyDist_c.restype = c_uint
lib.ElCauchyDist_z.argtypes = \
  [c_void_p,iType,POINTER(zType),iType,POINTER(zType)]
lib.ElCauchyDist_z.restype = c_uint
def Cauchy(A,x,y):
  xLen = len(x)
  yLen = len(y)
  xBuf = (TagToType(A.tag)*xLen)(*x)
  yBuf = (TagToType(A.tag)*yLen)(*y)
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElCauchy_s(A.obj,xLen,xBuf,yLen,yBuf)
    elif A.tag == dTag: lib.ElCauchy_d(A.obj,xLen,xBuf,yLen,yBuf)
    elif A.tag == cTag: lib.ElCauchy_c(A.obj,xLen,xBuf,yLen,yBuf)
    elif A.tag == zTag: lib.ElCauchy_z(A.obj,xLen,xBuf,yLen,yBuf)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElCauchyDist_s(A.obj,xLen,xBuf,yLen,yBuf)
    elif A.tag == dTag: lib.ElCauchyDist_d(A.obj,xLen,xBuf,yLen,yBuf)
    elif A.tag == cTag: lib.ElCauchyDist_c(A.obj,xLen,xBuf,yLen,yBuf)
    elif A.tag == zTag: lib.ElCauchyDist_z(A.obj,xLen,xBuf,yLen,yBuf)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Cauchy-like
# -----------
lib.ElCauchyLike_s.argtypes = \
  [c_void_p,iType,POINTER(sType),iType,POINTER(sType),
            iType,POINTER(sType),iType,POINTER(sType)]
lib.ElCauchyLike_s.restype = c_uint
lib.ElCauchyLike_d.argtypes = \
  [c_void_p,iType,POINTER(dType),iType,POINTER(dType),
            iType,POINTER(dType),iType,POINTER(dType)]
lib.ElCauchyLike_d.restype = c_uint
lib.ElCauchyLike_c.argtypes = \
  [c_void_p,iType,POINTER(cType),iType,POINTER(cType),
            iType,POINTER(cType),iType,POINTER(cType)]
lib.ElCauchyLike_c.restype = c_uint
lib.ElCauchyLike_z.argtypes = \
  [c_void_p,iType,POINTER(zType),iType,POINTER(zType),
            iType,POINTER(zType),iType,POINTER(zType)]
lib.ElCauchyLike_z.restype = c_uint
lib.ElCauchyLikeDist_s.argtypes = \
  [c_void_p,iType,POINTER(sType),iType,POINTER(sType),
            iType,POINTER(sType),iType,POINTER(sType)]
lib.ElCauchyLikeDist_s.restype = c_uint
lib.ElCauchyLikeDist_d.argtypes = \
  [c_void_p,iType,POINTER(dType),iType,POINTER(dType),
            iType,POINTER(dType),iType,POINTER(dType)]
lib.ElCauchyLikeDist_d.restype = c_uint
lib.ElCauchyLikeDist_c.argtypes = \
  [c_void_p,iType,POINTER(cType),iType,POINTER(cType),
            iType,POINTER(cType),iType,POINTER(cType)]
lib.ElCauchyLikeDist_c.restype = c_uint
lib.ElCauchyLikeDist_z.argtypes = \
  [c_void_p,iType,POINTER(zType),iType,POINTER(zType),
            iType,POINTER(zType),iType,POINTER(zType)]
lib.ElCauchyLikeDist_z.restype = c_uint
def CauchyLike(A,r,s,x,y):
  rLen = len(r)
  sLen = len(s)
  xLen = len(x)
  yLen = len(y)
  rBuf = (TagToType(A.tag)*rLen)(*r)
  sBuf = (TagToType(A.tag)*sLen)(*s)
  xBuf = (TagToType(A.tag)*xLen)(*x)
  yBuf = (TagToType(A.tag)*yLen)(*y)
  if type(A) is Matrix:
    if   A.tag == sTag: 
      lib.ElCauchyLike_s(A.obj,rLen,rBuf,sLen,sBuf,xLen,xBuf,yLen,yBuf)
    elif A.tag == dTag:
      lib.ElCauchyLike_d(A.obj,rLen,rBuf,sLen,sBuf,xLen,xBuf,yLen,yBuf)
    elif A.tag == cTag:
      lib.ElCauchyLike_c(A.obj,rLen,rBuf,sLen,sBuf,xLen,xBuf,yLen,yBuf)
    elif A.tag == zTag:
      lib.ElCauchyLike_z(A.obj,rLen,rBuf,sLen,sBuf,xLen,xBuf,yLen,yBuf)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElCauchyLikeDist_s(A.obj,rLen,rBuf,sLen,sBuf,xLen,xBuf,yLen,yBuf)
    elif A.tag == dTag:
      lib.ElCauchyLikeDist_d(A.obj,rLen,rBuf,sLen,sBuf,xLen,xBuf,yLen,yBuf)
    elif A.tag == cTag:
      lib.ElCauchyLikeDist_c(A.obj,rLen,rBuf,sLen,sBuf,xLen,xBuf,yLen,yBuf)
    elif A.tag == zTag:
      lib.ElCauchyLikeDist_z(A.obj,rLen,rBuf,sLen,sBuf,xLen,xBuf,yLen,yBuf)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Circulant
# ---------
lib.ElCirculant_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElCirculant_i.restype = c_uint
lib.ElCirculant_s.argtypes = [c_void_p,iType,POINTER(sType)]
lib.ElCirculant_s.restype = c_uint
lib.ElCirculant_d.argtypes = [c_void_p,iType,POINTER(dType)]
lib.ElCirculant_d.restype = c_uint
lib.ElCirculant_c.argtypes = [c_void_p,iType,POINTER(cType)]
lib.ElCirculant_c.restype = c_uint
lib.ElCirculant_z.argtypes = [c_void_p,iType,POINTER(zType)]
lib.ElCirculant_z.restype = c_uint
lib.ElCirculantDist_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElCirculantDist_i.restype = c_uint
lib.ElCirculantDist_s.argtypes = [c_void_p,iType,POINTER(sType)]
lib.ElCirculantDist_s.restype = c_uint
lib.ElCirculantDist_d.argtypes = [c_void_p,iType,POINTER(dType)]
lib.ElCirculantDist_d.restype = c_uint
lib.ElCirculantDist_c.argtypes = [c_void_p,iType,POINTER(cType)]
lib.ElCirculantDist_c.restype = c_uint
lib.ElCirculantDist_z.argtypes = [c_void_p,iType,POINTER(zType)]
lib.ElCirculantDist_z.restype = c_uint
def Circulant(A,a):
  aLen = len(a)
  aBuf = (TagToType(A.tag)*aLen)(*a)
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElCirculant_i(A.obj,aLen,aBuf)
    elif A.tag == sTag: lib.ElCirculant_s(A.obj,aLen,aBuf)
    elif A.tag == dTag: lib.ElCirculant_d(A.obj,aLen,aBuf)
    elif A.tag == cTag: lib.ElCirculant_c(A.obj,aLen,aBuf)
    elif A.tag == zTag: lib.ElCirculant_z(A.obj,aLen,aBuf)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElCirculantDist_i(A.obj,aLen,aBuf)
    elif A.tag == sTag: lib.ElCirculantDist_s(A.obj,aLen,aBuf)
    elif A.tag == dTag: lib.ElCirculantDist_d(A.obj,aLen,aBuf)
    elif A.tag == cTag: lib.ElCirculantDist_c(A.obj,aLen,aBuf)
    elif A.tag == zTag: lib.ElCirculantDist_z(A.obj,aLen,aBuf)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Demmel
# ------
lib.ElDemmel_s.argtypes = [c_void_p,iType]
lib.ElDemmel_s.restype = c_uint
lib.ElDemmel_d.argtypes = [c_void_p,iType]
lib.ElDemmel_d.restype = c_uint
lib.ElDemmel_c.argtypes = [c_void_p,iType]
lib.ElDemmel_c.restype = c_uint
lib.ElDemmel_z.argtypes = [c_void_p,iType]
lib.ElDemmel_z.restype = c_uint
lib.ElDemmelDist_s.argtypes = [c_void_p,iType]
lib.ElDemmelDist_s.restype = c_uint
lib.ElDemmelDist_d.argtypes = [c_void_p,iType]
lib.ElDemmelDist_d.restype = c_uint
lib.ElDemmelDist_c.argtypes = [c_void_p,iType]
lib.ElDemmelDist_c.restype = c_uint
lib.ElDemmelDist_z.argtypes = [c_void_p,iType]
lib.ElDemmelDist_z.restype = c_uint
def Demmel(A,n):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElDemmel_s(A.obj,n)
    elif A.tag == dTag: lib.ElDemmel_d(A.obj,n)
    elif A.tag == cTag: lib.ElDemmel_c(A.obj,n)
    elif A.tag == zTag: lib.ElDemmel_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElDemmelDist_s(A.obj,n)
    elif A.tag == dTag: lib.ElDemmelDist_d(A.obj,n)
    elif A.tag == cTag: lib.ElDemmelDist_c(A.obj,n)
    elif A.tag == zTag: lib.ElDemmelDist_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Diagonal
# --------
lib.ElDiagonal_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDiagonal_i.restype = c_uint
lib.ElDiagonal_s.argtypes = [c_void_p,iType,POINTER(sType)]
lib.ElDiagonal_s.restype = c_uint
lib.ElDiagonal_d.argtypes = [c_void_p,iType,POINTER(dType)]
lib.ElDiagonal_d.restype = c_uint
lib.ElDiagonal_c.argtypes = [c_void_p,iType,POINTER(cType)]
lib.ElDiagonal_c.restype = c_uint
lib.ElDiagonal_z.argtypes = [c_void_p,iType,POINTER(zType)]
lib.ElDiagonal_z.restype = c_uint
lib.ElDiagonalDist_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDiagonalDist_i.restype = c_uint
lib.ElDiagonalDist_s.argtypes = [c_void_p,iType,POINTER(sType)]
lib.ElDiagonalDist_s.restype = c_uint
lib.ElDiagonalDist_d.argtypes = [c_void_p,iType,POINTER(dType)]
lib.ElDiagonalDist_d.restype = c_uint
lib.ElDiagonalDist_c.argtypes = [c_void_p,iType,POINTER(cType)]
lib.ElDiagonalDist_c.restype = c_uint
lib.ElDiagonalDist_z.argtypes = [c_void_p,iType,POINTER(zType)]
lib.ElDiagonalDist_z.restype = c_uint
def Diagonal(A,d):
  dLen = len(d)
  dBuf = (TagToType(A.tag)*dLen)(*d)
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElDiagonal_i(A.obj,dLen,dBuf)
    elif A.tag == sTag: lib.ElDiagonal_s(A.obj,dLen,dBuf)
    elif A.tag == dTag: lib.ElDiagonal_d(A.obj,dLen,dBuf)
    elif A.tag == cTag: lib.ElDiagonal_c(A.obj,dLen,dBuf)
    elif A.tag == zTag: lib.ElDiagonal_z(A.obj,dLen,dBuf)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElDiagonalDist_i(A.obj,dLen,dBuf)
    elif A.tag == sTag: lib.ElDiagonalDist_s(A.obj,dLen,dBuf)
    elif A.tag == dTag: lib.ElDiagonalDist_d(A.obj,dLen,dBuf)
    elif A.tag == cTag: lib.ElDiagonalDist_c(A.obj,dLen,dBuf)
    elif A.tag == zTag: lib.ElDiagonalDist_z(A.obj,dLen,dBuf)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Egorov
# ------
lib.ElEgorov_c.argtypes = [c_void_p,CFUNCTYPE(sType,iType,iType),iType]
lib.ElEgorov_c.restype = c_uint
lib.ElEgorov_z.argtypes = [c_void_p,CFUNCTYPE(dType,iType,iType),iType]
lib.ElEgorov_z.restype = c_uint
lib.ElEgorovDist_c.argtypes = [c_void_p,CFUNCTYPE(sType,iType,iType),iType]
lib.ElEgorovDist_c.restype = c_uint
lib.ElEgorovDist_z.argtypes = [c_void_p,CFUNCTYPE(dType,iType,iType),iType]
lib.ElEgorovDist_z.restype = c_uint
def Egorov(A,phase,n):
  cPhase = CFUNCTYPE(TagToType(Base(A.tag)),iType,iType)(phase)
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElEgorov_c(A.obj,cPhase,n)
    elif A.tag == zTag: lib.ElEgorov_z(A.obj,cPhase,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElEgorovDist_c(A.obj,cPhase,n)
    elif A.tag == zTag: lib.ElEgorovDist_z(A.obj,cPhase,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Ehrenfest
# ---------
lib.ElEhrenfest_s.argtypes = [c_void_p,iType]
lib.ElEhrenfest_s.restype = c_uint
lib.ElEhrenfest_d.argtypes = [c_void_p,iType]
lib.ElEhrenfest_d.restype = c_uint
lib.ElEhrenfest_c.argtypes = [c_void_p,iType]
lib.ElEhrenfest_c.restype = c_uint
lib.ElEhrenfest_z.argtypes = [c_void_p,iType]
lib.ElEhrenfest_z.restype = c_uint
lib.ElEhrenfestDist_s.argtypes = [c_void_p,iType]
lib.ElEhrenfestDist_s.restype = c_uint
lib.ElEhrenfestDist_d.argtypes = [c_void_p,iType]
lib.ElEhrenfestDist_d.restype = c_uint
lib.ElEhrenfestDist_c.argtypes = [c_void_p,iType]
lib.ElEhrenfestDist_c.restype = c_uint
lib.ElEhrenfestDist_z.argtypes = [c_void_p,iType]
lib.ElEhrenfestDist_z.restype = c_uint
def Ehrenfest(P,n):
  if type(P) is Matrix:
    if   P.tag == sTag: lib.ElEhrenfest_s(P.obj,n)
    elif P.tag == dTag: lib.ElEhrenfest_d(P.obj,n)
    elif P.tag == cTag: lib.ElEhrenfest_c(P.obj,n)
    elif P.tag == zTag: lib.ElEhrenfest_z(P.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(P) is DistMatrix:
    if   P.tag == sTag: lib.ElEhrenfestDist_s(P.obj,n)
    elif P.tag == dTag: lib.ElEhrenfestDist_d(P.obj,n)
    elif P.tag == cTag: lib.ElEhrenfestDist_c(P.obj,n)
    elif P.tag == zTag: lib.ElEhrenfestDist_z(P.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElEhrenfestStationary_s.argtypes = [c_void_p,iType]
lib.ElEhrenfestStationary_s.restype = c_uint
lib.ElEhrenfestStationary_d.argtypes = [c_void_p,iType]
lib.ElEhrenfestStationary_d.restype = c_uint
lib.ElEhrenfestStationary_c.argtypes = [c_void_p,iType]
lib.ElEhrenfestStationary_c.restype = c_uint
lib.ElEhrenfestStationary_z.argtypes = [c_void_p,iType]
lib.ElEhrenfestStationary_z.restype = c_uint
lib.ElEhrenfestStationaryDist_s.argtypes = [c_void_p,iType]
lib.ElEhrenfestStationaryDist_s.restype = c_uint
lib.ElEhrenfestStationaryDist_d.argtypes = [c_void_p,iType]
lib.ElEhrenfestStationaryDist_d.restype = c_uint
lib.ElEhrenfestStationaryDist_c.argtypes = [c_void_p,iType]
lib.ElEhrenfestStationaryDist_c.restype = c_uint
lib.ElEhrenfestStationaryDist_z.argtypes = [c_void_p,iType]
lib.ElEhrenfestStationaryDist_z.restype = c_uint
def EhrenfestStationary(PInf,n):
  if type(PInf) is Matrix:
    if   PInf.tag == sTag: lib.ElEhrenfestStationary_s(PInf.obj,n)
    elif PInf.tag == dTag: lib.ElEhrenfestStationary_d(PInf.obj,n)
    elif PInf.tag == cTag: lib.ElEhrenfestStationary_c(PInf.obj,n)
    elif PInf.tag == zTag: lib.ElEhrenfestStationary_z(PInf.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(PInf) is DistMatrix:
    if   PInf.tag == sTag: lib.ElEhrenfestStationaryDist_s(PInf.obj,n)
    elif PInf.tag == dTag: lib.ElEhrenfestStationaryDist_d(PInf.obj,n)
    elif PInf.tag == cTag: lib.ElEhrenfestStationaryDist_c(PInf.obj,n)
    elif PInf.tag == zTag: lib.ElEhrenfestStationaryDist_z(PInf.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElEhrenfestDecay_s.argtypes = [c_void_p,iType]
lib.ElEhrenfestDecay_s.restype = c_uint
lib.ElEhrenfestDecay_d.argtypes = [c_void_p,iType]
lib.ElEhrenfestDecay_d.restype = c_uint
lib.ElEhrenfestDecay_c.argtypes = [c_void_p,iType]
lib.ElEhrenfestDecay_c.restype = c_uint
lib.ElEhrenfestDecay_z.argtypes = [c_void_p,iType]
lib.ElEhrenfestDecay_z.restype = c_uint
lib.ElEhrenfestDecayDist_s.argtypes = [c_void_p,iType]
lib.ElEhrenfestDecayDist_s.restype = c_uint
lib.ElEhrenfestDecayDist_d.argtypes = [c_void_p,iType]
lib.ElEhrenfestDecayDist_d.restype = c_uint
lib.ElEhrenfestDecayDist_c.argtypes = [c_void_p,iType]
lib.ElEhrenfestDecayDist_c.restype = c_uint
lib.ElEhrenfestDecayDist_z.argtypes = [c_void_p,iType]
lib.ElEhrenfestDecayDist_z.restype = c_uint
def EhrenfestDecay(PInf,n):
  if type(PInf) is Matrix:
    if   PInf.tag == sTag: lib.ElEhrenfestDecay_s(PInf.obj,n)
    elif PInf.tag == dTag: lib.ElEhrenfestDecay_d(PInf.obj,n)
    elif PInf.tag == cTag: lib.ElEhrenfestDecay_c(PInf.obj,n)
    elif PInf.tag == zTag: lib.ElEhrenfestDecay_z(PInf.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(PInf) is DistMatrix:
    if   PInf.tag == sTag: lib.ElEhrenfestDecayDist_s(PInf.obj,n)
    elif PInf.tag == dTag: lib.ElEhrenfestDecayDist_d(PInf.obj,n)
    elif PInf.tag == cTag: lib.ElEhrenfestDecayDist_c(PInf.obj,n)
    elif PInf.tag == zTag: lib.ElEhrenfestDecayDist_z(PInf.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Extended Kahan
# --------------
lib.ElExtendedKahan_s.argtypes = [c_void_p,iType,sType,sType]
lib.ElExtendedKahan_s.restype = c_uint
lib.ElExtendedKahan_d.argtypes = [c_void_p,iType,dType,dType]
lib.ElExtendedKahan_d.restype = c_uint
lib.ElExtendedKahan_c.argtypes = [c_void_p,iType,sType,sType]
lib.ElExtendedKahan_c.restype = c_uint
lib.ElExtendedKahan_z.argtypes = [c_void_p,iType,dType,dType]
lib.ElExtendedKahan_z.restype = c_uint
lib.ElExtendedKahanDist_s.argtypes = [c_void_p,iType,sType,sType]
lib.ElExtendedKahanDist_s.restype = c_uint
lib.ElExtendedKahanDist_d.argtypes = [c_void_p,iType,dType,dType]
lib.ElExtendedKahanDist_d.restype = c_uint
lib.ElExtendedKahanDist_c.argtypes = [c_void_p,iType,sType,sType]
lib.ElExtendedKahanDist_c.restype = c_uint
lib.ElExtendedKahanDist_z.argtypes = [c_void_p,iType,dType,dType]
lib.ElExtendedKahanDist_z.restype = c_uint
def ExtendedKahan(A,k,phi,mu):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElExtendedKahan_s(A.obj,k,phi,mu)
    elif A.tag == dTag: lib.ElExtendedKahan_d(A.obj,k,phi,mu)
    elif A.tag == cTag: lib.ElExtendedKahan_c(A.obj,k,phi,mu)
    elif A.tag == zTag: lib.ElExtendedKahan_z(A.obj,k,phi,mu)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElExtendedKahanDist_s(A.obj,k,phi,mu)
    elif A.tag == dTag: lib.ElExtendedKahanDist_d(A.obj,k,phi,mu)
    elif A.tag == cTag: lib.ElExtendedKahanDist_c(A.obj,k,phi,mu)
    elif A.tag == zTag: lib.ElExtendedKahanDist_z(A.obj,k,phi,mu)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Fiedler
# -------
lib.ElFiedler_s.argtypes = [c_void_p,iType,POINTER(sType)]
lib.ElFiedler_s.restype = c_uint
lib.ElFiedler_d.argtypes = [c_void_p,iType,POINTER(dType)]
lib.ElFiedler_d.restype = c_uint
lib.ElFiedler_c.argtypes = [c_void_p,iType,POINTER(cType)]
lib.ElFiedler_c.restype = c_uint
lib.ElFiedler_z.argtypes = [c_void_p,iType,POINTER(zType)]
lib.ElFiedler_z.restype = c_uint
lib.ElFiedlerDist_s.argtypes = [c_void_p,iType,POINTER(sType)]
lib.ElFiedlerDist_s.restype = c_uint
lib.ElFiedlerDist_d.argtypes = [c_void_p,iType,POINTER(dType)]
lib.ElFiedlerDist_d.restype = c_uint
lib.ElFiedlerDist_c.argtypes = [c_void_p,iType,POINTER(cType)]
lib.ElFiedlerDist_c.restype = c_uint
lib.ElFiedlerDist_z.argtypes = [c_void_p,iType,POINTER(zType)]
lib.ElFiedlerDist_z.restype = c_uint
def Fiedler(A,c):
  cLen = len(c)
  cBuf = (TagToType(A.tag)*cLen)(*c)
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElFiedler_s(A.obj,cLen,cBuf)
    elif A.tag == dTag: lib.ElFiedler_d(A.obj,cLen,cBuf)
    elif A.tag == cTag: lib.ElFiedler_c(A.obj,cLen,cBuf)
    elif A.tag == zTag: lib.ElFiedler_z(A.obj,cLen,cBuf)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElFiedlerDist_s(A.obj,cLen,cBuf)
    elif A.tag == dTag: lib.ElFiedlerDist_d(A.obj,cLen,cBuf)
    elif A.tag == cTag: lib.ElFiedlerDist_c(A.obj,cLen,cBuf)
    elif A.tag == zTag: lib.ElFiedlerDist_z(A.obj,cLen,cBuf)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Forsythe
# --------
lib.ElForsythe_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElForsythe_i.restype = c_uint
lib.ElForsythe_s.argtypes = [c_void_p,iType,sType,sType]
lib.ElForsythe_s.restype = c_uint
lib.ElForsythe_d.argtypes = [c_void_p,iType,dType,dType]
lib.ElForsythe_d.restype = c_uint
lib.ElForsythe_c.argtypes = [c_void_p,iType,cType,cType]
lib.ElForsythe_c.restype = c_uint
lib.ElForsythe_z.argtypes = [c_void_p,iType,zType,zType]
lib.ElForsythe_z.restype = c_uint
lib.ElForsytheDist_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElForsytheDist_i.restype = c_uint
lib.ElForsytheDist_s.argtypes = [c_void_p,iType,sType,sType]
lib.ElForsytheDist_s.restype = c_uint
lib.ElForsytheDist_d.argtypes = [c_void_p,iType,dType,dType]
lib.ElForsytheDist_d.restype = c_uint
lib.ElForsytheDist_c.argtypes = [c_void_p,iType,cType,cType]
lib.ElForsytheDist_c.restype = c_uint
lib.ElForsytheDist_z.argtypes = [c_void_p,iType,zType,zType]
lib.ElForsytheDist_z.restype = c_uint
def Forsythe(J,n,alpha,lamb):
  if type(J) is Matrix:
    if   J.tag == iTag: lib.ElForsythe_i(A.obj,n,alpha,lamb)
    elif J.tag == sTag: lib.ElForsythe_s(A.obj,n,alpha,lamb)
    elif J.tag == dTag: lib.ElForsythe_d(A.obj,n,alpha,lamb)
    elif J.tag == cTag: lib.ElForsythe_c(A.obj,n,alpha,lamb)
    elif J.tag == zTag: lib.ElForsythe_z(A.obj,n,alpha,lamb)
    else: raise Exception('Unsupported datatype')
  elif type(J) is DistMatrix:
    if   J.tag == iTag: lib.ElForsytheDist_i(A.obj,n,alpha,lamb)
    elif J.tag == sTag: lib.ElForsytheDist_s(A.obj,n,alpha,lamb)
    elif J.tag == dTag: lib.ElForsytheDist_d(A.obj,n,alpha,lamb)
    elif J.tag == cTag: lib.ElForsytheDist_c(A.obj,n,alpha,lamb)
    elif J.tag == zTag: lib.ElForsytheDist_z(A.obj,n,alpha,lamb)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Fox-Li
# ------
lib.ElFoxLi_c.argtypes = [c_void_p,iType,sType]
lib.ElFoxLi_c.restype = c_uint
lib.ElFoxLi_z.argtypes = [c_void_p,iType,dType]
lib.ElFoxLi_z.restype = c_uint
lib.ElFoxLiDist_c.argtypes = [c_void_p,iType,sType]
lib.ElFoxLiDist_c.restype = c_uint
lib.ElFoxLiDist_z.argtypes = [c_void_p,iType,dType]
lib.ElFoxLiDist_z.restype = c_uint
def FoxLi(A,n,omega):
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElFoxLi_c(A.obj,n,omega)
    elif A.tag == zTag: lib.ElFoxLi_z(A.obj,n,omega)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElFoxLiDist_c(A.obj,n,omega)
    elif A.tag == zTag: lib.ElFoxLiDist_z(A.obj,n,omega)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Fourier
# -------
lib.ElFourier_c.argtypes = [c_void_p,iType]
lib.ElFourier_c.restype = c_uint
lib.ElFourier_z.argtypes = [c_void_p,iType]
lib.ElFourier_z.restype = c_uint
lib.ElFourierDist_c.argtypes = [c_void_p,iType]
lib.ElFourierDist_c.restype = c_uint
lib.ElFourierDist_z.argtypes = [c_void_p,iType]
lib.ElFourierDist_z.restype = c_uint
def Fourier(A,n):
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElFourier_c(A.obj,n)
    elif A.tag == zTag: lib.ElFourier_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElFourierDist_c(A.obj,n)
    elif A.tag == zTag: lib.ElFourierDist_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Fourier-Identity
# ----------------
lib.ElFourierIdentity_c.argtypes = [c_void_p,iType]
lib.ElFourierIdentity_c.restype = c_uint
lib.ElFourierIdentity_z.argtypes = [c_void_p,iType]
lib.ElFourierIdentity_z.restype = c_uint
lib.ElFourierIdentityDist_c.argtypes = [c_void_p,iType]
lib.ElFourierIdentityDist_c.restype = c_uint
lib.ElFourierIdentityDist_z.argtypes = [c_void_p,iType]
lib.ElFourierIdentityDist_z.restype = c_uint
def FourierIdentity(A,n):
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElFourierIdentity_c(A.obj,n)
    elif A.tag == zTag: lib.ElFourierIdentity_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElFourierIdentityDist_c(A.obj,n)
    elif A.tag == zTag: lib.ElFourierIdentityDist_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# GCD matrix
# ----------
lib.ElGCDMatrix_i.argtypes = [c_void_p,iType,iType]
lib.ElGCDMatrix_i.restype = c_uint
lib.ElGCDMatrix_s.argtypes = [c_void_p,iType,iType]
lib.ElGCDMatrix_s.restype = c_uint
lib.ElGCDMatrix_d.argtypes = [c_void_p,iType,iType]
lib.ElGCDMatrix_d.restype = c_uint
lib.ElGCDMatrix_c.argtypes = [c_void_p,iType,iType]
lib.ElGCDMatrix_c.restype = c_uint
lib.ElGCDMatrix_z.argtypes = [c_void_p,iType,iType]
lib.ElGCDMatrix_z.restype = c_uint
lib.ElGCDMatrixDist_i.argtypes = [c_void_p,iType,iType]
lib.ElGCDMatrixDist_i.restype = c_uint
lib.ElGCDMatrixDist_s.argtypes = [c_void_p,iType,iType]
lib.ElGCDMatrixDist_s.restype = c_uint
lib.ElGCDMatrixDist_d.argtypes = [c_void_p,iType,iType]
lib.ElGCDMatrixDist_d.restype = c_uint
lib.ElGCDMatrixDist_c.argtypes = [c_void_p,iType,iType]
lib.ElGCDMatrixDist_c.restype = c_uint
lib.ElGCDMatrixDist_z.argtypes = [c_void_p,iType,iType]
lib.ElGCDMatrixDist_z.restype = c_uint
def GCDMatrix(G,m,n):
  if type(G) is Matrix:
    if   G.tag == iTag: lib.ElGCDMatrix_i(G.obj,m,n)
    elif G.tag == sTag: lib.ElGCDMatrix_s(G.obj,m,n)
    elif G.tag == dTag: lib.ElGCDMatrix_d(G.obj,m,n)
    elif G.tag == cTag: lib.ElGCDMatrix_c(G.obj,m,n)
    elif G.tag == zTag: lib.ElGCDMatrix_z(G.obj,m,n)
    else: raise Exception('Unsupported datatype')
  elif type(G) is DistMatrix:
    if   G.tag == iTag: lib.ElGCDMatrixDist_i(G.obj,m,n)
    elif G.tag == sTag: lib.ElGCDMatrixDist_s(G.obj,m,n)
    elif G.tag == dTag: lib.ElGCDMatrixDist_d(G.obj,m,n)
    elif G.tag == cTag: lib.ElGCDMatrixDist_c(G.obj,m,n)
    elif G.tag == zTag: lib.ElGCDMatrixDist_z(G.obj,m,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Gear matrix
# -----------
lib.ElGear_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElGear_i.restype = c_uint
lib.ElGear_s.argtypes = [c_void_p,iType,iType,iType]
lib.ElGear_s.restype = c_uint
lib.ElGear_d.argtypes = [c_void_p,iType,iType,iType]
lib.ElGear_d.restype = c_uint
lib.ElGear_c.argtypes = [c_void_p,iType,iType,iType]
lib.ElGear_c.restype = c_uint
lib.ElGear_z.argtypes = [c_void_p,iType,iType,iType]
lib.ElGear_z.restype = c_uint
lib.ElGearDist_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElGearDist_i.restype = c_uint
lib.ElGearDist_s.argtypes = [c_void_p,iType,iType,iType]
lib.ElGearDist_s.restype = c_uint
lib.ElGearDist_d.argtypes = [c_void_p,iType,iType,iType]
lib.ElGearDist_d.restype = c_uint
lib.ElGearDist_c.argtypes = [c_void_p,iType,iType,iType]
lib.ElGearDist_c.restype = c_uint
lib.ElGearDist_z.argtypes = [c_void_p,iType,iType,iType]
lib.ElGearDist_z.restype = c_uint
def Gear(G,n,s,t):
  if type(G) is Matrix:
    if   G.tag == iTag: lib.ElGear_i(G.obj,n,s,t)
    elif G.tag == sTag: lib.ElGear_s(G.obj,n,s,t)
    elif G.tag == dTag: lib.ElGear_d(G.obj,n,s,t)
    elif G.tag == cTag: lib.ElGear_c(G.obj,n,s,t)
    elif G.tag == zTag: lib.ElGear_z(G.obj,n,s,t)
    else: raise Exception('Unsupported datatype')
  elif type(G) is DistMatrix:
    if   G.tag == iTag: lib.ElGearDist_i(G.obj,n,s,t)
    elif G.tag == sTag: lib.ElGearDist_s(G.obj,n,s,t)
    elif G.tag == dTag: lib.ElGearDist_d(G.obj,n,s,t)
    elif G.tag == cTag: lib.ElGearDist_c(G.obj,n,s,t)
    elif G.tag == zTag: lib.ElGearDist_z(G.obj,n,s,t)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# GEPP Growth
# -----------
lib.ElGEPPGrowth_s.argtypes = [c_void_p,iType]
lib.ElGEPPGrowth_s.restype = c_uint
lib.ElGEPPGrowth_d.argtypes = [c_void_p,iType]
lib.ElGEPPGrowth_d.restype = c_uint
lib.ElGEPPGrowth_c.argtypes = [c_void_p,iType]
lib.ElGEPPGrowth_c.restype = c_uint
lib.ElGEPPGrowth_z.argtypes = [c_void_p,iType]
lib.ElGEPPGrowth_z.restype = c_uint
lib.ElGEPPGrowthDist_s.argtypes = [c_void_p,iType]
lib.ElGEPPGrowthDist_s.restype = c_uint
lib.ElGEPPGrowthDist_d.argtypes = [c_void_p,iType]
lib.ElGEPPGrowthDist_d.restype = c_uint
lib.ElGEPPGrowthDist_c.argtypes = [c_void_p,iType]
lib.ElGEPPGrowthDist_c.restype = c_uint
lib.ElGEPPGrowthDist_z.argtypes = [c_void_p,iType]
lib.ElGEPPGrowthDist_z.restype = c_uint
def GEPPGrowth(A,n):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElGEPPGrowth_s(A.obj,n)
    elif A.tag == dTag: lib.ElGEPPGrowth_d(A.obj,n)
    elif A.tag == cTag: lib.ElGEPPGrowth_c(A.obj,n)
    elif A.tag == zTag: lib.ElGEPPGrowth_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElGEPPGrowthDist_s(A.obj,n)
    elif A.tag == dTag: lib.ElGEPPGrowthDist_d(A.obj,n)
    elif A.tag == cTag: lib.ElGEPPGrowthDist_c(A.obj,n)
    elif A.tag == zTag: lib.ElGEPPGrowthDist_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Golub/Klema/Stewart
# -------------------
lib.ElGKS_s.argtypes = [c_void_p,iType]
lib.ElGKS_s.restype = c_uint
lib.ElGKS_d.argtypes = [c_void_p,iType]
lib.ElGKS_d.restype = c_uint
lib.ElGKS_c.argtypes = [c_void_p,iType]
lib.ElGKS_c.restype = c_uint
lib.ElGKS_z.argtypes = [c_void_p,iType]
lib.ElGKS_z.restype = c_uint
lib.ElGKSDist_s.argtypes = [c_void_p,iType]
lib.ElGKSDist_s.restype = c_uint
lib.ElGKSDist_d.argtypes = [c_void_p,iType]
lib.ElGKSDist_d.restype = c_uint
lib.ElGKSDist_c.argtypes = [c_void_p,iType]
lib.ElGKSDist_c.restype = c_uint
lib.ElGKSDist_z.argtypes = [c_void_p,iType]
lib.ElGKSDist_z.restype = c_uint
def GKS(A,n):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElGKS_s(A.obj,n)
    elif A.tag == dTag: lib.ElGKS_d(A.obj,n)
    elif A.tag == cTag: lib.ElGKS_c(A.obj,n)
    elif A.tag == zTag: lib.ElGKS_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElGKSDist_s(A.obj,n)
    elif A.tag == dTag: lib.ElGKSDist_d(A.obj,n)
    elif A.tag == cTag: lib.ElGKSDist_c(A.obj,n)
    elif A.tag == zTag: lib.ElGKSDist_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Grcar
# -----
lib.ElGrcar_i.argtypes = [c_void_p,iType,iType]
lib.ElGrcar_i.restype = c_uint
lib.ElGrcar_s.argtypes = [c_void_p,iType,iType]
lib.ElGrcar_s.restype = c_uint
lib.ElGrcar_d.argtypes = [c_void_p,iType,iType]
lib.ElGrcar_d.restype = c_uint
lib.ElGrcar_c.argtypes = [c_void_p,iType,iType]
lib.ElGrcar_c.restype = c_uint
lib.ElGrcar_z.argtypes = [c_void_p,iType,iType]
lib.ElGrcar_z.restype = c_uint
lib.ElGrcarDist_i.argtypes = [c_void_p,iType,iType]
lib.ElGrcarDist_i.restype = c_uint
lib.ElGrcarDist_s.argtypes = [c_void_p,iType,iType]
lib.ElGrcarDist_s.restype = c_uint
lib.ElGrcarDist_d.argtypes = [c_void_p,iType,iType]
lib.ElGrcarDist_d.restype = c_uint
lib.ElGrcarDist_c.argtypes = [c_void_p,iType,iType]
lib.ElGrcarDist_c.restype = c_uint
lib.ElGrcarDist_z.argtypes = [c_void_p,iType,iType]
lib.ElGrcarDist_z.restype = c_uint
def Grcar(A,n,k=3):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElGrcar_i(A.obj,n,k)
    elif A.tag == sTag: lib.ElGrcar_s(A.obj,n,k)
    elif A.tag == dTag: lib.ElGrcar_d(A.obj,n,k)
    elif A.tag == cTag: lib.ElGrcar_c(A.obj,n,k)
    elif A.tag == zTag: lib.ElGrcar_z(A.obj,n,k)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElGrcarDist_i(A.obj,n,k)
    elif A.tag == sTag: lib.ElGrcarDist_s(A.obj,n,k)
    elif A.tag == dTag: lib.ElGrcarDist_d(A.obj,n,k)
    elif A.tag == cTag: lib.ElGrcarDist_c(A.obj,n,k)
    elif A.tag == zTag: lib.ElGrcarDist_z(A.obj,n,k)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Haar
# ----
lib.ElHaar_s.argtypes = [c_void_p,iType]
lib.ElHaar_s.restype = c_uint
lib.ElHaar_d.argtypes = [c_void_p,iType]
lib.ElHaar_d.restype = c_uint
lib.ElHaar_c.argtypes = [c_void_p,iType]
lib.ElHaar_c.restype = c_uint
lib.ElHaar_z.argtypes = [c_void_p,iType]
lib.ElHaar_z.restype = c_uint
lib.ElHaarDist_s.argtypes = [c_void_p,iType]
lib.ElHaarDist_s.restype = c_uint
lib.ElHaarDist_d.argtypes = [c_void_p,iType]
lib.ElHaarDist_d.restype = c_uint
lib.ElHaarDist_c.argtypes = [c_void_p,iType]
lib.ElHaarDist_c.restype = c_uint
lib.ElHaarDist_z.argtypes = [c_void_p,iType]
lib.ElHaarDist_z.restype = c_uint
def Haar(A,n):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHaar_s(A.obj,n)
    elif A.tag == dTag: lib.ElHaar_d(A.obj,n)
    elif A.tag == cTag: lib.ElHaar_c(A.obj,n)
    elif A.tag == zTag: lib.ElHaar_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHaarDist_s(A.obj,n)
    elif A.tag == dTag: lib.ElHaarDist_d(A.obj,n)
    elif A.tag == cTag: lib.ElHaarDist_c(A.obj,n)
    elif A.tag == zTag: lib.ElHaarDist_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElImplicitHaar_s.argtypes = [c_void_p,c_void_p,c_void_p,iType]
lib.ElImplicitHaar_s.restype = c_uint
lib.ElImplicitHaar_d.argtypes = [c_void_p,c_void_p,c_void_p,iType]
lib.ElImplicitHaar_d.restype = c_uint
lib.ElImplicitHaar_c.argtypes = [c_void_p,c_void_p,c_void_p,iType]
lib.ElImplicitHaar_c.restype = c_uint
lib.ElImplicitHaar_z.argtypes = [c_void_p,c_void_p,c_void_p,iType]
lib.ElImplicitHaar_z.restype = c_uint
lib.ElImplicitHaarDist_s.argtypes = [c_void_p,c_void_p,c_void_p,iType]
lib.ElImplicitHaarDist_s.restype = c_uint
lib.ElImplicitHaarDist_d.argtypes = [c_void_p,c_void_p,c_void_p,iType]
lib.ElImplicitHaarDist_d.restype = c_uint
lib.ElImplicitHaarDist_c.argtypes = [c_void_p,c_void_p,c_void_p,iType]
lib.ElImplicitHaarDist_c.restype = c_uint
lib.ElImplicitHaarDist_z.argtypes = [c_void_p,c_void_p,c_void_p,iType]
lib.ElImplicitHaarDist_z.restype = c_uint
def ImplicitHaar(A,t,d,n):
  if type(A) is not type(d): raise Exception('A and t must be of the same type')
  if type(A) is not type(t): raise Exception('A and d must be of the same type')
  if t.tag != A.tag: raise Exception('A and t must have the same datatype')
  if d.tag != Base(A.tag): raise Exception('d must have the base datatype of A')
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElImplicitHaar_s(A.obj,t.obj,d.obj,n)
    elif A.tag == dTag: lib.ElImplicitHaar_d(A.obj,t.obj,d.obj,n)
    elif A.tag == cTag: lib.ElImplicitHaar_c(A.obj,t.obj,d.obj,n)
    elif A.tag == zTag: lib.ElImplicitHaar_z(A.obj,t.obj,d.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElImplicitHaarDist_s(A.obj,t.obj,d.obj,n)
    elif A.tag == dTag: lib.ElImplicitHaarDist_d(A.obj,t.obj,d.obj,n)
    elif A.tag == cTag: lib.ElImplicitHaarDist_c(A.obj,t.obj,d.obj,n)
    elif A.tag == zTag: lib.ElImplicitHaarDist_z(A.obj,t.obj,d.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Hankel
# ------
lib.ElHankel_i.argtypes = [c_void_p,iType,iType,iType,POINTER(iType)]
lib.ElHankel_i.restype = c_uint
lib.ElHankel_s.argtypes = [c_void_p,iType,iType,iType,POINTER(sType)]
lib.ElHankel_s.restype = c_uint
lib.ElHankel_d.argtypes = [c_void_p,iType,iType,iType,POINTER(dType)]
lib.ElHankel_d.restype = c_uint
lib.ElHankel_c.argtypes = [c_void_p,iType,iType,iType,POINTER(cType)]
lib.ElHankel_c.restype = c_uint
lib.ElHankel_z.argtypes = [c_void_p,iType,iType,iType,POINTER(zType)]
lib.ElHankel_z.restype = c_uint
lib.ElHankelDist_i.argtypes = [c_void_p,iType,iType,iType,POINTER(iType)]
lib.ElHankelDist_i.restype = c_uint
lib.ElHankelDist_s.argtypes = [c_void_p,iType,iType,iType,POINTER(sType)]
lib.ElHankelDist_s.restype = c_uint
lib.ElHankelDist_d.argtypes = [c_void_p,iType,iType,iType,POINTER(dType)]
lib.ElHankelDist_d.restype = c_uint
lib.ElHankelDist_c.argtypes = [c_void_p,iType,iType,iType,POINTER(cType)]
lib.ElHankelDist_c.restype = c_uint
lib.ElHankelDist_z.argtypes = [c_void_p,iType,iType,iType,POINTER(zType)]
lib.ElHankelDist_z.restype = c_uint
def Hankel(A,m,n,a):
  aLen = len(a)
  aBuf = (TagToType(A.tag)*aLen)(*a)
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElHankel_i(A.obj,m,n,aLen,aBuf)
    elif A.tag == sTag: lib.ElHankel_s(A.obj,m,n,aLen,aBuf)
    elif A.tag == dTag: lib.ElHankel_d(A.obj,m,n,aLen,aBuf)
    elif A.tag == cTag: lib.ElHankel_c(A.obj,m,n,aLen,aBuf)
    elif A.tag == zTag: lib.ElHankel_z(A.obj,m,n,aLen,aBuf)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElHankelDist_i(A.obj,m,n,aLen,aBuf)
    elif A.tag == sTag: lib.ElHankelDist_s(A.obj,m,n,aLen,aBuf)
    elif A.tag == dTag: lib.ElHankelDist_d(A.obj,m,n,aLen,aBuf)
    elif A.tag == cTag: lib.ElHankelDist_c(A.obj,m,n,aLen,aBuf)
    elif A.tag == zTag: lib.ElHankelDist_z(A.obj,m,n,aLen,aBuf)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Hanowa
# ------
lib.ElHanowa_i.argtypes = [c_void_p,iType,iType]
lib.ElHanowa_i.restype = c_uint
lib.ElHanowa_s.argtypes = [c_void_p,iType,sType]
lib.ElHanowa_s.restype = c_uint
lib.ElHanowa_d.argtypes = [c_void_p,iType,dType]
lib.ElHanowa_d.restype = c_uint
lib.ElHanowa_c.argtypes = [c_void_p,iType,cType]
lib.ElHanowa_c.restype = c_uint
lib.ElHanowa_z.argtypes = [c_void_p,iType,zType]
lib.ElHanowa_z.restype = c_uint
lib.ElHanowaDist_i.argtypes = [c_void_p,iType,iType]
lib.ElHanowaDist_i.restype = c_uint
lib.ElHanowaDist_s.argtypes = [c_void_p,iType,sType]
lib.ElHanowaDist_s.restype = c_uint
lib.ElHanowaDist_d.argtypes = [c_void_p,iType,dType]
lib.ElHanowaDist_d.restype = c_uint
lib.ElHanowaDist_c.argtypes = [c_void_p,iType,cType]
lib.ElHanowaDist_c.restype = c_uint
lib.ElHanowaDist_z.argtypes = [c_void_p,iType,zType]
lib.ElHanowaDist_z.restype = c_uint
def Hanowa(A,n,mu):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElHanowa_i(A.obj,n,mu)
    elif A.tag == sTag: lib.ElHanowa_s(A.obj,n,mu)
    elif A.tag == dTag: lib.ElHanowa_d(A.obj,n,mu)
    elif A.tag == cTag: lib.ElHanowa_c(A.obj,n,mu)
    elif A.tag == zTag: lib.ElHanowa_z(A.obj,n,mu)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElHanowaDist_i(A.obj,n,mu)
    elif A.tag == sTag: lib.ElHanowaDist_s(A.obj,n,mu)
    elif A.tag == dTag: lib.ElHanowaDist_d(A.obj,n,mu)
    elif A.tag == cTag: lib.ElHanowaDist_c(A.obj,n,mu)
    elif A.tag == zTag: lib.ElHanowaDist_z(A.obj,n,mu)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Hatano-Nelson
# -------------
lib.ElHatanoNelson_s.argtypes = [c_void_p,iType,sType,sType,sType,bType]
lib.ElHatanoNelson_s.restype = c_uint
lib.ElHatanoNelson_d.argtypes = [c_void_p,iType,dType,dType,dType,bType]
lib.ElHatanoNelson_d.restype = c_uint
lib.ElHatanoNelson_c.argtypes = [c_void_p,iType,cType,sType,cType,bType]
lib.ElHatanoNelson_c.restype = c_uint
lib.ElHatanoNelson_z.argtypes = [c_void_p,iType,zType,dType,zType,bType]
lib.ElHatanoNelson_z.restype = c_uint
lib.ElHatanoNelsonDist_s.argtypes = [c_void_p,iType,sType,sType,sType,bType]
lib.ElHatanoNelsonDist_s.restype = c_uint
lib.ElHatanoNelsonDist_d.argtypes = [c_void_p,iType,dType,dType,dType,bType]
lib.ElHatanoNelsonDist_d.restype = c_uint
lib.ElHatanoNelsonDist_c.argtypes = [c_void_p,iType,cType,sType,cType,bType]
lib.ElHatanoNelsonDist_c.restype = c_uint
lib.ElHatanoNelsonDist_z.argtypes = [c_void_p,iType,zType,dType,zType,bType]
lib.ElHatanoNelsonDist_z.restype = c_uint
def HatanoNelson(A,n,center,radius,g,periodic=True):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHatanoNelson_s(A.obj,n,center,radius,g,periodic)
    elif A.tag == dTag: lib.ElHatanoNelson_d(A.obj,n,center,radius,g,periodic)
    elif A.tag == cTag: lib.ElHatanoNelson_c(A.obj,n,center,radius,g,periodic)
    elif A.tag == zTag: lib.ElHatanoNelson_z(A.obj,n,center,radius,g,periodic)
    else: raise Exception('Unsupported datatype')
  elif tyep(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElHatanoNelsonDist_s(A.obj,n,center,radius,g,periodic)
    elif A.tag == dTag: 
      lib.ElHatanoNelsonDist_d(A.obj,n,center,radius,g,periodic)
    elif A.tag == cTag: 
      lib.ElHatanoNelsonDist_c(A.obj,n,center,radius,g,periodic)
    elif A.tag == zTag: 
      lib.ElHatanoNelsonDist_z(A.obj,n,center,radius,g,periodic)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Helmholtz
# ---------
lib.ElHelmholtz1D_s.argtypes = [c_void_p,iType,sType]
lib.ElHelmholtz1D_s.restype = c_uint
lib.ElHelmholtz1D_d.argtypes = [c_void_p,iType,dType]
lib.ElHelmholtz1D_d.restype = c_uint
lib.ElHelmholtz1D_c.argtypes = [c_void_p,iType,cType]
lib.ElHelmholtz1D_c.restype = c_uint
lib.ElHelmholtz1D_z.argtypes = [c_void_p,iType,zType]
lib.ElHelmholtz1D_z.restype = c_uint
lib.ElHelmholtz1DDist_s.argtypes = [c_void_p,iType,sType]
lib.ElHelmholtz1DDist_s.restype = c_uint
lib.ElHelmholtz1DDist_d.argtypes = [c_void_p,iType,dType]
lib.ElHelmholtz1DDist_d.restype = c_uint
lib.ElHelmholtz1DDist_c.argtypes = [c_void_p,iType,cType]
lib.ElHelmholtz1DDist_c.restype = c_uint
lib.ElHelmholtz1DDist_z.argtypes = [c_void_p,iType,zType]
lib.ElHelmholtz1DDist_z.restype = c_uint
def Helmholtz1D(H,nx,shift):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHelmholtz1D_s(H.obj,nx,shift)
    elif A.tag == dTag: lib.ElHelmholtz1D_d(H.obj,nx,shift)
    elif A.tag == cTag: lib.ElHelmholtz1D_c(H.obj,nx,shift)
    elif A.tag == zTag: lib.ElHelmholtz1D_z(H.obj,nx,shift)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHelmholtz1DDist_s(H.obj,nx,shift)
    elif A.tag == dTag: lib.ElHelmholtz1DDist_d(H.obj,nx,shift)
    elif A.tag == cTag: lib.ElHelmholtz1DDist_c(H.obj,nx,shift)
    elif A.tag == zTag: lib.ElHelmholtz1DDist_z(H.obj,nx,shift)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElHelmholtz2D_s.argtypes = [c_void_p,iType,iType,sType]
lib.ElHelmholtz2D_s.restype = c_uint
lib.ElHelmholtz2D_d.argtypes = [c_void_p,iType,iType,dType]
lib.ElHelmholtz2D_d.restype = c_uint
lib.ElHelmholtz2D_c.argtypes = [c_void_p,iType,iType,cType]
lib.ElHelmholtz2D_c.restype = c_uint
lib.ElHelmholtz2D_z.argtypes = [c_void_p,iType,iType,zType]
lib.ElHelmholtz2D_z.restype = c_uint
lib.ElHelmholtz2DDist_s.argtypes = [c_void_p,iType,iType,sType]
lib.ElHelmholtz2DDist_s.restype = c_uint
lib.ElHelmholtz2DDist_d.argtypes = [c_void_p,iType,iType,dType]
lib.ElHelmholtz2DDist_d.restype = c_uint
lib.ElHelmholtz2DDist_c.argtypes = [c_void_p,iType,iType,cType]
lib.ElHelmholtz2DDist_c.restype = c_uint
lib.ElHelmholtz2DDist_z.argtypes = [c_void_p,iType,iType,zType]
lib.ElHelmholtz2DDist_z.restype = c_uint
def Helmholtz2D(H,nx,ny,shift):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHelmholtz2D_s(H.obj,nx,ny,shift)
    elif A.tag == dTag: lib.ElHelmholtz2D_d(H.obj,nx,ny,shift)
    elif A.tag == cTag: lib.ElHelmholtz2D_c(H.obj,nx,ny,shift)
    elif A.tag == zTag: lib.ElHelmholtz2D_z(H.obj,nx,ny,shift)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHelmholtz2DDist_s(H.obj,nx,ny,shift)
    elif A.tag == dTag: lib.ElHelmholtz2DDist_d(H.obj,nx,ny,shift)
    elif A.tag == cTag: lib.ElHelmholtz2DDist_c(H.obj,nx,ny,shift)
    elif A.tag == zTag: lib.ElHelmholtz2DDist_z(H.obj,nx,ny,shift)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElHelmholtz3D_s.argtypes = [c_void_p,iType,iType,iType,sType]
lib.ElHelmholtz3D_s.restype = c_uint
lib.ElHelmholtz3D_d.argtypes = [c_void_p,iType,iType,iType,dType]
lib.ElHelmholtz3D_d.restype = c_uint
lib.ElHelmholtz3D_c.argtypes = [c_void_p,iType,iType,iType,cType]
lib.ElHelmholtz3D_c.restype = c_uint
lib.ElHelmholtz3D_z.argtypes = [c_void_p,iType,iType,iType,zType]
lib.ElHelmholtz3D_z.restype = c_uint
lib.ElHelmholtz3DDist_s.argtypes = [c_void_p,iType,iType,iType,sType]
lib.ElHelmholtz3DDist_s.restype = c_uint
lib.ElHelmholtz3DDist_d.argtypes = [c_void_p,iType,iType,iType,dType]
lib.ElHelmholtz3DDist_d.restype = c_uint
lib.ElHelmholtz3DDist_c.argtypes = [c_void_p,iType,iType,iType,cType]
lib.ElHelmholtz3DDist_c.restype = c_uint
lib.ElHelmholtz3DDist_z.argtypes = [c_void_p,iType,iType,iType,zType]
lib.ElHelmholtz3DDist_z.restype = c_uint
def Helmholtz3D(H,nx,ny,nz,shift):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHelmholtz3D_s(H.obj,nx,ny,nz,shift)
    elif A.tag == dTag: lib.ElHelmholtz3D_d(H.obj,nx,ny,nz,shift)
    elif A.tag == cTag: lib.ElHelmholtz3D_c(H.obj,nx,ny,nz,shift)
    elif A.tag == zTag: lib.ElHelmholtz3D_z(H.obj,nx,ny,nz,shift)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHelmholtz3DDist_s(H.obj,nx,ny,nz,shift)
    elif A.tag == dTag: lib.ElHelmholtz3DDist_d(H.obj,nx,ny,nz,shift)
    elif A.tag == cTag: lib.ElHelmholtz3DDist_c(H.obj,nx,ny,nz,shift)
    elif A.tag == zTag: lib.ElHelmholtz3DDist_z(H.obj,nx,ny,nz,shift)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Helmholtz with PML
# ------------------
lib.ElHelmholtzPML1D_c.argtypes = [c_void_p,iType,cType,iType,sType,sType]
lib.ElHelmholtzPML1D_c.restype = c_uint
lib.ElHelmholtzPML1D_z.argtypes = [c_void_p,iType,zType,iType,dType,dType]
lib.ElHelmholtzPML1D_z.restype = c_uint
lib.ElHelmholtzPML1DDist_c.argtypes = [c_void_p,iType,cType,iType,sType,sType]
lib.ElHelmholtzPML1DDist_c.restype = c_uint
lib.ElHelmholtzPML1DDist_z.argtypes = [c_void_p,iType,zType,iType,dType,dType]
lib.ElHelmholtzPML1DDist_z.restype = c_uint
def HelmholtzPML1D(H,nx,omega,numPml,sigma,pmlExp):
  if type(A) is Matrix:
    if   A.tag == cTag: 
      lib.ElHelmholtzPML1D_c(H.obj,nx,omega,numPml,sigma,pmlExp)
    elif A.tag == zTag:
      lib.ElHelmholtzPML1D_z(H.obj,nx,omega,numPml,sigma,pmlExp)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == cTag: 
      lib.ElHelmholtzPML1DDist_c(H.obj,nx,omega,numPml,sigma,pmlExp)
    elif A.tag == zTag:
      lib.ElHelmholtzPML1DDist_z(H.obj,nx,omega,numPml,sigma,pmlExp)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElHelmholtzPML2D_c.argtypes = [c_void_p,iType,iType,cType,iType,sType,sType]
lib.ElHelmholtzPML2D_c.restype = c_uint
lib.ElHelmholtzPML2D_z.argtypes = [c_void_p,iType,iType,zType,iType,dType,dType]
lib.ElHelmholtzPML2D_z.restype = c_uint
lib.ElHelmholtzPML2DDist_c.argtypes = \
  [c_void_p,iType,iType,cType,iType,sType,sType]
lib.ElHelmholtzPML2DDist_c.restype = c_uint
lib.ElHelmholtzPML2DDist_z.argtypes = \
  [c_void_p,iType,iType,zType,iType,dType,dType]
lib.ElHelmholtzPML2DDist_z.restype = c_uint
def HelmholtzPML2D(H,nx,ny,omega,numPml,sigma,pmlExp):
  if type(A) is Matrix:
    if   A.tag == cTag: 
      lib.ElHelmholtzPML2D_c(H.obj,nx,ny,omega,numPml,sigma,pmlExp)
    elif A.tag == zTag:
      lib.ElHelmholtzPML2D_z(H.obj,nx,ny,omega,numPml,sigma,pmlExp)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == cTag: 
      lib.ElHelmholtzPML2DDist_c(H.obj,nx,ny,omega,numPml,sigma,pmlExp)
    elif A.tag == zTag:
      lib.ElHelmholtzPML2DDist_z(H.obj,nx,ny,omega,numPml,sigma,pmlExp)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElHelmholtzPML3D_c.argtypes = \
  [c_void_p,iType,iType,iType,cType,iType,sType,sType]
lib.ElHelmholtzPML3D_c.restype = c_uint
lib.ElHelmholtzPML3D_z.argtypes = \
  [c_void_p,iType,iType,iType,zType,iType,dType,dType]
lib.ElHelmholtzPML3D_z.restype = c_uint
lib.ElHelmholtzPML3DDist_c.argtypes = \
  [c_void_p,iType,iType,iType,cType,iType,sType,sType]
lib.ElHelmholtzPML3DDist_c.restype = c_uint
lib.ElHelmholtzPML3DDist_z.argtypes = \
  [c_void_p,iType,iType,iType,zType,iType,dType,dType]
lib.ElHelmholtzPML3DDist_z.restype = c_uint
def HelmholtzPML3D(H,nx,ny,nz,omega,numPml,sigma,pmlExp):
  if type(A) is Matrix:
    if   A.tag == cTag: 
      lib.ElHelmholtzPML3D_c(H.obj,nx,ny,nz,omega,numPml,sigma,pmlExp)
    elif A.tag == zTag:
      lib.ElHelmholtzPML3D_z(H.obj,nx,ny,nz,omega,numPml,sigma,pmlExp)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == cTag: 
      lib.ElHelmholtzPML3DDist_c(H.obj,nx,ny,nz,omega,numPml,sigma,pmlExp)
    elif A.tag == zTag:
      lib.ElHelmholtzPML3DDist_z(H.obj,nx,ny,nz,omega,numPml,sigma,pmlExp)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Hermitian from EVD
# ------------------
lib.ElHermitianFromEVD_s.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElHermitianFromEVD_s.restype = c_uint
lib.ElHermitianFromEVD_d.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElHermitianFromEVD_d.restype = c_uint
lib.ElHermitianFromEVD_c.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElHermitianFromEVD_c.restype = c_uint
lib.ElHermitianFromEVD_z.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElHermitianFromEVD_z.restype = c_uint
lib.ElHermitianFromEVDDist_s.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElHermitianFromEVDDist_s.restype = c_uint
lib.ElHermitianFromEVDDist_d.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElHermitianFromEVDDist_d.restype = c_uint
lib.ElHermitianFromEVDDist_c.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElHermitianFromEVDDist_c.restype = c_uint
lib.ElHermitianFromEVDDist_z.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElHermitianFromEVDDist_z.restype = c_uint
def HermitianFromEVD(uplo,A,w,Z):
  if Z.tag != A.tag: raise Exception('A and Z must have the same datatype')
  if w.tag != Base(A.tag): raise Exception('w must be of the base type of A')
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHermitianFromEVD_s(uplo,A.obj,w.obj,Z.obj)
    elif A.tag == dTag: lib.ElHermitianFromEVD_d(uplo,A.obj,w.obj,Z.obj)
    elif A.tag == cTag: lib.ElHermitianFromEVD_c(uplo,A.obj,w.obj,Z.obj)
    elif A.tag == zTag: lib.ElHermitianFromEVD_z(uplo,A.obj,w.obj,Z.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHermitianFromEVDDist_s(uplo,A.obj,w.obj,Z.obj)
    elif A.tag == dTag: lib.ElHermitianFromEVDDist_d(uplo,A.obj,w.obj,Z.obj)
    elif A.tag == cTag: lib.ElHermitianFromEVDDist_c(uplo,A.obj,w.obj,Z.obj)
    elif A.tag == zTag: lib.ElHermitianFromEVDDist_z(uplo,A.obj,w.obj,Z.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Hermitian uniform spectrum
# --------------------------
lib.ElHermitianUniformSpectrum_s.argtypes = [c_void_p,iType,sType,sType]
lib.ElHermitianUniformSpectrum_s.restype = c_uint
lib.ElHermitianUniformSpectrum_d.argtypes = [c_void_p,iType,dType,dType]
lib.ElHermitianUniformSpectrum_d.restype = c_uint
lib.ElHermitianUniformSpectrum_c.argtypes = [c_void_p,iType,sType,sType]
lib.ElHermitianUniformSpectrum_c.restype = c_uint
lib.ElHermitianUniformSpectrum_z.argtypes = [c_void_p,iType,dType,dType]
lib.ElHermitianUniformSpectrum_z.restype = c_uint
lib.ElHermitianUniformSpectrumDist_s.argtypes = [c_void_p,iType,sType,sType]
lib.ElHermitianUniformSpectrumDist_s.restype = c_uint
lib.ElHermitianUniformSpectrumDist_d.argtypes = [c_void_p,iType,dType,dType]
lib.ElHermitianUniformSpectrumDist_d.restype = c_uint
lib.ElHermitianUniformSpectrumDist_c.argtypes = [c_void_p,iType,sType,sType]
lib.ElHermitianUniformSpectrumDist_c.restype = c_uint
lib.ElHermitianUniformSpectrumDist_z.argtypes = [c_void_p,iType,dType,dType]
lib.ElHermitianUniformSpectrumDist_z.restype = c_uint
def HermitianUniformSpectrum(A,n,lower=0,upper=1):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHermitianUniformSpectrum_s(A.obj,n,lower,upper)
    elif A.tag == dTag: lib.ElHermitianUniformSpectrum_d(A.obj,n,lower,upper)
    elif A.tag == cTag: lib.ElHermitianUniformSpectrum_c(A.obj,n,lower,upper)
    elif A.tag == zTag: lib.ElHermitianUniformSpectrum_z(A.obj,n,lower,upper)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElHermitianUniformSpectrumDist_s(A.obj,n,lower,upper)
    elif A.tag == dTag: 
      lib.ElHermitianUniformSpectrumDist_d(A.obj,n,lower,upper)
    elif A.tag == cTag: 
      lib.ElHermitianUniformSpectrumDist_c(A.obj,n,lower,upper)
    elif A.tag == zTag: 
      lib.ElHermitianUniformSpectrumDist_z(A.obj,n,lower,upper)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Hilbert
# -------
lib.ElHilbert_s.argtypes = [c_void_p,iType]
lib.ElHilbert_s.restype = c_uint
lib.ElHilbert_d.argtypes = [c_void_p,iType]
lib.ElHilbert_d.restype = c_uint
lib.ElHilbert_c.argtypes = [c_void_p,iType]
lib.ElHilbert_c.restype = c_uint
lib.ElHilbert_z.argtypes = [c_void_p,iType]
lib.ElHilbert_z.restype = c_uint
lib.ElHilbertDist_s.argtypes = [c_void_p,iType]
lib.ElHilbertDist_s.restype = c_uint
lib.ElHilbertDist_d.argtypes = [c_void_p,iType]
lib.ElHilbertDist_d.restype = c_uint
lib.ElHilbertDist_c.argtypes = [c_void_p,iType]
lib.ElHilbertDist_c.restype = c_uint
lib.ElHilbertDist_z.argtypes = [c_void_p,iType]
lib.ElHilbertDist_z.restype = c_uint
def Hilbert(A,n):
  if type(A) is Matrix:
    if   A.tag == sType: lib.ElHilbert_s(A.obj,n)
    elif A.tag == dType: lib.ElHilbert_d(A.obj,n)
    elif A.tag == cType: lib.ElHilbert_c(A.obj,n)
    elif A.tag == zType: lib.ElHilbert_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sType: lib.ElHilbertDist_s(A.obj,n)
    elif A.tag == dType: lib.ElHilbertDist_d(A.obj,n)
    elif A.tag == cType: lib.ElHilbertDist_c(A.obj,n)
    elif A.tag == zType: lib.ElHilbertDist_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Identity
# --------
lib.ElIdentity_i.argtypes = [c_void_p,iType,iType]
lib.ElIdentity_i.restype = c_uint
lib.ElIdentity_s.argtypes = [c_void_p,iType,iType]
lib.ElIdentity_s.restype = c_uint
lib.ElIdentity_d.argtypes = [c_void_p,iType,iType]
lib.ElIdentity_d.restype = c_uint
lib.ElIdentity_c.argtypes = [c_void_p,iType,iType]
lib.ElIdentity_c.restype = c_uint
lib.ElIdentity_z.argtypes = [c_void_p,iType,iType]
lib.ElIdentity_z.restype = c_uint
lib.ElIdentityDist_i.argtypes = [c_void_p,iType,iType]
lib.ElIdentityDist_i.restype = c_uint
lib.ElIdentityDist_s.argtypes = [c_void_p,iType,iType]
lib.ElIdentityDist_s.restype = c_uint
lib.ElIdentityDist_d.argtypes = [c_void_p,iType,iType]
lib.ElIdentityDist_d.restype = c_uint
lib.ElIdentityDist_c.argtypes = [c_void_p,iType,iType]
lib.ElIdentityDist_c.restype = c_uint
lib.ElIdentityDist_z.argtypes = [c_void_p,iType,iType]
lib.ElIdentityDist_z.restype = c_uint
def Identity(A,m,n):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElIdentity_i(A.obj,m,n)
    elif A.tag == sTag: lib.ElIdentity_s(A.obj,m,n)
    elif A.tag == dTag: lib.ElIdentity_d(A.obj,m,n)
    elif A.tag == cTag: lib.ElIdentity_c(A.obj,m,n)
    elif A.tag == zTag: lib.ElIdentity_z(A.obj,m,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElIdentityDist_i(A.obj,m,n)
    elif A.tag == sTag: lib.ElIdentityDist_s(A.obj,m,n)
    elif A.tag == dTag: lib.ElIdentityDist_d(A.obj,m,n)
    elif A.tag == cTag: lib.ElIdentityDist_c(A.obj,m,n)
    elif A.tag == zTag: lib.ElIdentityDist_z(A.obj,m,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Jordan
# ------
lib.ElJordan_i.argtypes = [c_void_p,iType,iType]
lib.ElJordan_i.restype = c_uint
lib.ElJordan_s.argtypes = [c_void_p,iType,sType]
lib.ElJordan_s.restype = c_uint
lib.ElJordan_d.argtypes = [c_void_p,iType,dType]
lib.ElJordan_d.restype = c_uint
lib.ElJordan_c.argtypes = [c_void_p,iType,cType]
lib.ElJordan_c.restype = c_uint
lib.ElJordan_z.argtypes = [c_void_p,iType,zType]
lib.ElJordan_z.restype = c_uint
lib.ElJordanDist_i.argtypes = [c_void_p,iType,iType]
lib.ElJordanDist_i.restype = c_uint
lib.ElJordanDist_s.argtypes = [c_void_p,iType,sType]
lib.ElJordanDist_s.restype = c_uint
lib.ElJordanDist_d.argtypes = [c_void_p,iType,dType]
lib.ElJordanDist_d.restype = c_uint
lib.ElJordanDist_c.argtypes = [c_void_p,iType,cType]
lib.ElJordanDist_c.restype = c_uint
lib.ElJordanDist_z.argtypes = [c_void_p,iType,zType]
lib.ElJordanDist_z.restype = c_uint
def Jordan(J,n,lamb):
  if type(J) is Matrix: 
    if   J.tag == iTag: lib.ElJordan_i(J.obj,n,lamb)
    elif J.tag == sTag: lib.ElJordan_s(J.obj,n,lamb)
    elif J.tag == dTag: lib.ElJordan_d(J.obj,n,lamb)
    elif J.tag == cTag: lib.ElJordan_c(J.obj,n,lamb)
    elif J.tag == zTag: lib.ElJordan_z(J.obj,n,lamb)
    else: raise Exception('Unsupported datatype')
  elif type(J) is DistMatrix:
    if   J.tag == iTag: lib.ElJordanDist_i(J.obj,n,lamb)
    elif J.tag == sTag: lib.ElJordanDist_s(J.obj,n,lamb)
    elif J.tag == dTag: lib.ElJordanDist_d(J.obj,n,lamb)
    elif J.tag == cTag: lib.ElJordanDist_c(J.obj,n,lamb)
    elif J.tag == zTag: lib.ElJordanDist_z(J.obj,n,lamb)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Kahan
# -----
lib.ElKahan_s.argtypes = [c_void_p,iType,sType]
lib.ElKahan_s.restype = c_uint
lib.ElKahan_d.argtypes = [c_void_p,iType,dType]
lib.ElKahan_d.restype = c_uint
lib.ElKahan_c.argtypes = [c_void_p,iType,cType]
lib.ElKahan_c.restype = c_uint
lib.ElKahan_z.argtypes = [c_void_p,iType,zType]
lib.ElKahan_z.restype = c_uint
lib.ElKahanDist_s.argtypes = [c_void_p,iType,sType]
lib.ElKahanDist_s.restype = c_uint
lib.ElKahanDist_d.argtypes = [c_void_p,iType,dType]
lib.ElKahanDist_d.restype = c_uint
lib.ElKahanDist_c.argtypes = [c_void_p,iType,cType]
lib.ElKahanDist_c.restype = c_uint
lib.ElKahanDist_z.argtypes = [c_void_p,iType,zType]
lib.ElKahanDist_z.restype = c_uint
def Kahan(A,n,phi):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElKahan_s(A.obj,n,phi)
    elif A.tag == dTag: lib.ElKahan_d(A.obj,n,phi)
    elif A.tag == cTag: lib.ElKahan_c(A.obj,n,phi)
    elif A.tag == zTag: lib.ElKahan_z(A.obj,n,phi)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElKahanDist_s(A.obj,n,phi)
    elif A.tag == dTag: lib.ElKahanDist_d(A.obj,n,phi)
    elif A.tag == cTag: lib.ElKahanDist_c(A.obj,n,phi)
    elif A.tag == zTag: lib.ElKahanDist_z(A.obj,n,phi)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# KMS
# ---
lib.ElKMS_i.argtypes = [c_void_p,iType,iType]
lib.ElKMS_i.restype = c_uint
lib.ElKMS_s.argtypes = [c_void_p,iType,sType]
lib.ElKMS_s.restype = c_uint
lib.ElKMS_d.argtypes = [c_void_p,iType,dType]
lib.ElKMS_d.restype = c_uint
lib.ElKMS_c.argtypes = [c_void_p,iType,cType]
lib.ElKMS_c.restype = c_uint
lib.ElKMS_z.argtypes = [c_void_p,iType,zType]
lib.ElKMS_z.restype = c_uint
lib.ElKMSDist_i.argtypes = [c_void_p,iType,iType]
lib.ElKMSDist_i.restype = c_uint
lib.ElKMSDist_s.argtypes = [c_void_p,iType,sType]
lib.ElKMSDist_s.restype = c_uint
lib.ElKMSDist_d.argtypes = [c_void_p,iType,dType]
lib.ElKMSDist_d.restype = c_uint
lib.ElKMSDist_c.argtypes = [c_void_p,iType,cType]
lib.ElKMSDist_c.restype = c_uint
lib.ElKMSDist_z.argtypes = [c_void_p,iType,zType]
lib.ElKMSDist_z.restype = c_uint
def KMS(K,n,rho):
  if type(K) is Matrix:
    if   K.tag == iTag: lib.ElKMS_i(K.obj,n,rho)
    elif K.tag == sTag: lib.ElKMS_s(K.obj,n,rho)
    elif K.tag == dTag: lib.ElKMS_d(K.obj,n,rho)
    elif K.tag == cTag: lib.ElKMS_c(K.obj,n,rho)
    elif K.tag == zTag: lib.ElKMS_z(K.obj,n,rho)
    else: raise Exception('Unsupported datatype')
  elif type(K) is DistMatrix:
    if   K.tag == iTag: lib.ElKMSDist_i(K.obj,n,rho)
    elif K.tag == sTag: lib.ElKMSDist_s(K.obj,n,rho)
    elif K.tag == dTag: lib.ElKMSDist_d(K.obj,n,rho)
    elif K.tag == cTag: lib.ElKMSDist_c(K.obj,n,rho)
    elif K.tag == zTag: lib.ElKMSDist_z(K.obj,n,rho)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Laplacian
# ---------
lib.ElLaplacian1D_s.argtypes = [c_void_p,iType]
lib.ElLaplacian1D_s.restype = c_uint
lib.ElLaplacian1D_d.argtypes = [c_void_p,iType]
lib.ElLaplacian1D_d.restype = c_uint
lib.ElLaplacian1D_c.argtypes = [c_void_p,iType]
lib.ElLaplacian1D_c.restype = c_uint
lib.ElLaplacian1D_z.argtypes = [c_void_p,iType]
lib.ElLaplacian1D_z.restype = c_uint
lib.ElLaplacian1DDist_s.argtypes = [c_void_p,iType]
lib.ElLaplacian1DDist_s.restype = c_uint
lib.ElLaplacian1DDist_d.argtypes = [c_void_p,iType]
lib.ElLaplacian1DDist_d.restype = c_uint
lib.ElLaplacian1DDist_c.argtypes = [c_void_p,iType]
lib.ElLaplacian1DDist_c.restype = c_uint
lib.ElLaplacian1DDist_z.argtypes = [c_void_p,iType]
lib.ElLaplacian1DDist_z.restype = c_uint
def Laplacian1D(L,nx):
  if type(L) is Matrix:
    if   L.tag == sTag: lib.ElLaplacian1D_s(L.obj,nx)
    elif L.tag == dTag: lib.ElLaplacian1D_d(L.obj,nx)
    elif L.tag == cTag: lib.ElLaplacian1D_c(L.obj,nx)
    elif L.tag == zTag: lib.ElLaplacian1D_z(L.obj,nx)
    else: raise Exception('Unsupported datatype')
  elif type(L) is DistMatrix:
    if   L.tag == sTag: lib.ElLaplacian1DDist_s(L.obj,nx)
    elif L.tag == dTag: lib.ElLaplacian1DDist_d(L.obj,nx)
    elif L.tag == cTag: lib.ElLaplacian1DDist_c(L.obj,nx)
    elif L.tag == zTag: lib.ElLaplacian1DDist_z(L.obj,nx)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElLaplacian2D_s.argtypes = [c_void_p,iType,iType]
lib.ElLaplacian2D_s.restype = c_uint
lib.ElLaplacian2D_d.argtypes = [c_void_p,iType,iType]
lib.ElLaplacian2D_d.restype = c_uint
lib.ElLaplacian2D_c.argtypes = [c_void_p,iType,iType]
lib.ElLaplacian2D_c.restype = c_uint
lib.ElLaplacian2D_z.argtypes = [c_void_p,iType,iType]
lib.ElLaplacian2D_z.restype = c_uint
lib.ElLaplacian2DDist_s.argtypes = [c_void_p,iType,iType]
lib.ElLaplacian2DDist_s.restype = c_uint
lib.ElLaplacian2DDist_d.argtypes = [c_void_p,iType,iType]
lib.ElLaplacian2DDist_d.restype = c_uint
lib.ElLaplacian2DDist_c.argtypes = [c_void_p,iType,iType]
lib.ElLaplacian2DDist_c.restype = c_uint
lib.ElLaplacian2DDist_z.argtypes = [c_void_p,iType,iType]
lib.ElLaplacian2DDist_z.restype = c_uint
def Laplacian2D(L,nx,ny):
  if type(L) is Matrix:
    if   L.tag == sTag: lib.ElLaplacian2D_s(L.obj,nx,ny)
    elif L.tag == dTag: lib.ElLaplacian2D_d(L.obj,nx,ny)
    elif L.tag == cTag: lib.ElLaplacian2D_c(L.obj,nx,ny)
    elif L.tag == zTag: lib.ElLaplacian2D_z(L.obj,nx,ny)
    else: raise Exception('Unsupported datatype')
  elif type(L) is DistMatrix:
    if   L.tag == sTag: lib.ElLaplacian2DDist_s(L.obj,nx,ny)
    elif L.tag == dTag: lib.ElLaplacian2DDist_d(L.obj,nx,ny)
    elif L.tag == cTag: lib.ElLaplacian2DDist_c(L.obj,nx,ny)
    elif L.tag == zTag: lib.ElLaplacian2DDist_z(L.obj,nx,ny)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElLaplacian3D_s.argtypes = [c_void_p,iType,iType,iType]
lib.ElLaplacian3D_s.restype = c_uint
lib.ElLaplacian3D_d.argtypes = [c_void_p,iType,iType,iType]
lib.ElLaplacian3D_d.restype = c_uint
lib.ElLaplacian3D_c.argtypes = [c_void_p,iType,iType,iType]
lib.ElLaplacian3D_c.restype = c_uint
lib.ElLaplacian3D_z.argtypes = [c_void_p,iType,iType,iType]
lib.ElLaplacian3D_z.restype = c_uint
lib.ElLaplacian3DDist_s.argtypes = [c_void_p,iType,iType,iType]
lib.ElLaplacian3DDist_s.restype = c_uint
lib.ElLaplacian3DDist_d.argtypes = [c_void_p,iType,iType,iType]
lib.ElLaplacian3DDist_d.restype = c_uint
lib.ElLaplacian3DDist_c.argtypes = [c_void_p,iType,iType,iType]
lib.ElLaplacian3DDist_c.restype = c_uint
lib.ElLaplacian3DDist_z.argtypes = [c_void_p,iType,iType,iType]
lib.ElLaplacian3DDist_z.restype = c_uint
def Laplacian3D(L,nx,ny,nz):
  if type(L) is Matrix:
    if   L.tag == sTag: lib.ElLaplacian3D_s(L.obj,nx,ny,nz)
    elif L.tag == dTag: lib.ElLaplacian3D_d(L.obj,nx,ny,nz)
    elif L.tag == cTag: lib.ElLaplacian3D_c(L.obj,nx,ny,nz)
    elif L.tag == zTag: lib.ElLaplacian3D_z(L.obj,nx,ny,nz)
    else: raise Exception('Unsupported datatype')
  elif type(L) is DistMatrix:
    if   L.tag == sTag: lib.ElLaplacian3DDist_s(L.obj,nx,ny,nz)
    elif L.tag == dTag: lib.ElLaplacian3DDist_d(L.obj,nx,ny,nz)
    elif L.tag == cTag: lib.ElLaplacian3DDist_c(L.obj,nx,ny,nz)
    elif L.tag == zTag: lib.ElLaplacian3DDist_z(L.obj,nx,ny,nz)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Lauchli
# -------
lib.ElLauchli_i.argtypes = [c_void_p,iType,iType]
lib.ElLauchli_i.restype = c_uint
lib.ElLauchli_s.argtypes = [c_void_p,iType,sType]
lib.ElLauchli_s.restype = c_uint
lib.ElLauchli_d.argtypes = [c_void_p,iType,dType]
lib.ElLauchli_d.restype = c_uint
lib.ElLauchli_c.argtypes = [c_void_p,iType,cType]
lib.ElLauchli_c.restype = c_uint
lib.ElLauchli_z.argtypes = [c_void_p,iType,zType]
lib.ElLauchli_z.restype = c_uint
lib.ElLauchliDist_i.argtypes = [c_void_p,iType,iType]
lib.ElLauchliDist_i.restype = c_uint
lib.ElLauchliDist_s.argtypes = [c_void_p,iType,sType]
lib.ElLauchliDist_s.restype = c_uint
lib.ElLauchliDist_d.argtypes = [c_void_p,iType,dType]
lib.ElLauchliDist_d.restype = c_uint
lib.ElLauchliDist_c.argtypes = [c_void_p,iType,cType]
lib.ElLauchliDist_c.restype = c_uint
lib.ElLauchliDist_z.argtypes = [c_void_p,iType,zType]
lib.ElLauchliDist_z.restype = c_uint
def Lauchli(A,n,mu):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElLauchli_i(A.obj,n,mu)
    elif A.tag == sTag: lib.ElLauchli_s(A.obj,n,mu)
    elif A.tag == dTag: lib.ElLauchli_d(A.obj,n,mu)
    elif A.tag == cTag: lib.ElLauchli_c(A.obj,n,mu)
    elif A.tag == zTag: lib.ElLauchli_z(A.obj,n,mu)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElLauchliDist_i(A.obj,n,mu)
    elif A.tag == sTag: lib.ElLauchliDist_s(A.obj,n,mu)
    elif A.tag == dTag: lib.ElLauchliDist_d(A.obj,n,mu)
    elif A.tag == cTag: lib.ElLauchliDist_c(A.obj,n,mu)
    elif A.tag == zTag: lib.ElLauchliDist_z(A.obj,n,mu)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Legendre
# --------
lib.ElLegendre_s.argtypes = [c_void_p,iType]
lib.ElLegendre_s.restype = c_uint
lib.ElLegendre_d.argtypes = [c_void_p,iType]
lib.ElLegendre_d.restype = c_uint
lib.ElLegendre_c.argtypes = [c_void_p,iType]
lib.ElLegendre_c.restype = c_uint
lib.ElLegendre_z.argtypes = [c_void_p,iType]
lib.ElLegendre_z.restype = c_uint
lib.ElLegendreDist_s.argtypes = [c_void_p,iType]
lib.ElLegendreDist_s.restype = c_uint
lib.ElLegendreDist_d.argtypes = [c_void_p,iType]
lib.ElLegendreDist_d.restype = c_uint
lib.ElLegendreDist_c.argtypes = [c_void_p,iType]
lib.ElLegendreDist_c.restype = c_uint
lib.ElLegendreDist_z.argtypes = [c_void_p,iType]
lib.ElLegendreDist_z.restype = c_uint
def Legendre(A,n):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElLegendre_s(A.obj,n)
    elif A.tag == dTag: lib.ElLegendre_d(A.obj,n)
    elif A.tag == cTag: lib.ElLegendre_c(A.obj,n)
    elif A.tag == zTag: lib.ElLegendre_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElLegendreDist_s(A.obj,n)
    elif A.tag == dTag: lib.ElLegendreDist_d(A.obj,n)
    elif A.tag == cTag: lib.ElLegendreDist_c(A.obj,n)
    elif A.tag == zTag: lib.ElLegendreDist_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Lehmer
# ------
lib.ElLehmer_s.argtypes = [c_void_p,iType]
lib.ElLehmer_s.restype = c_uint
lib.ElLehmer_d.argtypes = [c_void_p,iType]
lib.ElLehmer_d.restype = c_uint
lib.ElLehmer_c.argtypes = [c_void_p,iType]
lib.ElLehmer_c.restype = c_uint
lib.ElLehmer_z.argtypes = [c_void_p,iType]
lib.ElLehmer_z.restype = c_uint
lib.ElLehmerDist_s.argtypes = [c_void_p,iType]
lib.ElLehmerDist_s.restype = c_uint
lib.ElLehmerDist_d.argtypes = [c_void_p,iType]
lib.ElLehmerDist_d.restype = c_uint
lib.ElLehmerDist_c.argtypes = [c_void_p,iType]
lib.ElLehmerDist_c.restype = c_uint
lib.ElLehmerDist_z.argtypes = [c_void_p,iType]
lib.ElLehmerDist_z.restype = c_uint
def Lehmer(L,n):
  if type(L) is Matrix:
    if   L.tag == sTag: lib.ElLehmer_s(L.tag,n)
    elif L.tag == dTag: lib.ElLehmer_d(L.tag,n)
    elif L.tag == cTag: lib.ElLehmer_c(L.tag,n)
    elif L.tag == zTag: lib.ElLehmer_z(L.tag,n)
    else: raise Exception('Unsupported datatype')
  elif type(L) is DistMatrix:
    if   L.tag == sTag: lib.ElLehmerDist_s(L.tag,n)
    elif L.tag == dTag: lib.ElLehmerDist_d(L.tag,n)
    elif L.tag == cTag: lib.ElLehmerDist_c(L.tag,n)
    elif L.tag == zTag: lib.ElLehmerDist_z(L.tag,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Lotkin
# ------
lib.ElLotkin_s.argtypes = [c_void_p,iType]
lib.ElLotkin_s.restype = c_uint
lib.ElLotkin_d.argtypes = [c_void_p,iType]
lib.ElLotkin_d.restype = c_uint
lib.ElLotkin_c.argtypes = [c_void_p,iType]
lib.ElLotkin_c.restype = c_uint
lib.ElLotkin_z.argtypes = [c_void_p,iType]
lib.ElLotkin_z.restype = c_uint
lib.ElLotkinDist_s.argtypes = [c_void_p,iType]
lib.ElLotkinDist_s.restype = c_uint
lib.ElLotkinDist_d.argtypes = [c_void_p,iType]
lib.ElLotkinDist_d.restype = c_uint
lib.ElLotkinDist_c.argtypes = [c_void_p,iType]
lib.ElLotkinDist_c.restype = c_uint
lib.ElLotkinDist_z.argtypes = [c_void_p,iType]
lib.ElLotkinDist_z.restype = c_uint
def Lotkin(A,n):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElLotkin_s(A.obj,n)
    elif A.tag == dTag: lib.ElLotkin_d(A.obj,n)
    elif A.tag == cTag: lib.ElLotkin_c(A.obj,n)
    elif A.tag == zTag: lib.ElLotkin_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElLotkinDist_s(A.obj,n)
    elif A.tag == dTag: lib.ElLotkinDist_d(A.obj,n)
    elif A.tag == cTag: lib.ElLotkinDist_c(A.obj,n)
    elif A.tag == zTag: lib.ElLotkinDist_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# MinIJ
# -----
lib.ElMinIJ_i.argtypes = [c_void_p,iType]
lib.ElMinIJ_i.restype = c_uint
lib.ElMinIJ_s.argtypes = [c_void_p,iType]
lib.ElMinIJ_s.restype = c_uint
lib.ElMinIJ_d.argtypes = [c_void_p,iType]
lib.ElMinIJ_d.restype = c_uint
lib.ElMinIJ_c.argtypes = [c_void_p,iType]
lib.ElMinIJ_c.restype = c_uint
lib.ElMinIJ_z.argtypes = [c_void_p,iType]
lib.ElMinIJ_z.restype = c_uint
lib.ElMinIJDist_i.argtypes = [c_void_p,iType]
lib.ElMinIJDist_i.restype = c_uint
lib.ElMinIJDist_s.argtypes = [c_void_p,iType]
lib.ElMinIJDist_s.restype = c_uint
lib.ElMinIJDist_d.argtypes = [c_void_p,iType]
lib.ElMinIJDist_d.restype = c_uint
lib.ElMinIJDist_c.argtypes = [c_void_p,iType]
lib.ElMinIJDist_c.restype = c_uint
lib.ElMinIJDist_z.argtypes = [c_void_p,iType]
lib.ElMinIJDist_z.restype = c_uint
def MinIJ(A,n):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElMinIJ_i(A.obj,n)
    elif A.tag == sTag: lib.ElMinIJ_s(A.obj,n)
    elif A.tag == dTag: lib.ElMinIJ_d(A.obj,n)
    elif A.tag == cTag: lib.ElMinIJ_c(A.obj,n)
    elif A.tag == zTag: lib.ElMinIJ_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElMinIJDist_i(A.obj,n)
    elif A.tag == sTag: lib.ElMinIJDist_s(A.obj,n)
    elif A.tag == dTag: lib.ElMinIJDist_d(A.obj,n)
    elif A.tag == cTag: lib.ElMinIJDist_c(A.obj,n)
    elif A.tag == zTag: lib.ElMinIJDist_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Normal from EVD
# ---------------
lib.ElNormalFromEVD_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElNormalFromEVD_c.restype = c_uint
lib.ElNormalFromEVD_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElNormalFromEVD_z.restype = c_uint
lib.ElNormalFromEVDDist_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElNormalFromEVDDist_c.restype = c_uint
lib.ElNormalFromEVDDist_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElNormalFromEVDDist_z.restype = c_uint
def NormalFromEVD(A,w,Z):
  if type(A) is not type(w): raise Exception('Types of A and w must match')
  if type(A) is not type(Z): raise Exception('Types of A and Z must match')
  if Z.tag != A.tag: raise Exception('Datatypes of A and Z must match')
  if w.tag != Base(A.tag): raise Exception('Base datatype of A must match w')
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElNormalFromEVD_c(A.obj,w.obj,Z.obj)
    elif A.tag == zTag: lib.ElNormalFromEVD_z(A.obj,w.obj,Z.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElNormalFromEVDDist_c(A.obj,w.obj,Z.obj)
    elif A.tag == zTag: lib.ElNormalFromEVDDist_z(A.obj,w.obj,Z.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Normal uniform spectrum
# -----------------------
lib.ElNormalUniformSpectrum_c.argtypes = [c_void_p,iType,cType,sType]
lib.ElNormalUniformSpectrum_c.restype = c_uint
lib.ElNormalUniformSpectrum_z.argtypes = [c_void_p,iType,zType,dType]
lib.ElNormalUniformSpectrum_z.restype = c_uint
lib.ElNormalUniformSpectrumDist_c.argtypes = [c_void_p,iType,cType,sType]
lib.ElNormalUniformSpectrumDist_c.restype = c_uint
lib.ElNormalUniformSpectrumDist_z.argtypes = [c_void_p,iType,zType,dType]
lib.ElNormalUniformSpectrumDist_z.restype = c_uint
def NormalUniformSpectrum(A,n,centerPre=0,radius=1):
  center = TagToType(A.tag)(centerPre)
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElNormalUniformSpectrum_c(A.obj,n,center,radius)
    elif A.tag == zTag: lib.ElNormalUniformSpectrum_z(A.obj,n,center,radius)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElNormalUniformSpectrumDist_c(A.obj,n,center,radius)
    elif A.tag == zTag: lib.ElNormalUniformSpectrumDist_z(A.obj,n,center,radius)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Ones
# ----
lib.ElOnes_c.argtypes = [c_void_p,iType,iType]
lib.ElOnes_c.restype = c_uint
lib.ElOnes_z.argtypes = [c_void_p,iType,iType]
lib.ElOnes_z.restype = c_uint
lib.ElOnesDist_c.argtypes = [c_void_p,iType,iType]
lib.ElOnesDist_c.restype = c_uint
lib.ElOnesDist_z.argtypes = [c_void_p,iType,iType]
lib.ElOnesDist_z.restype = c_uint
def Ones(A,m,n):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElOnes_i(A.obj,m,n)
    elif A.tag == sTag: lib.ElOnes_s(A.obj,m,n)
    elif A.tag == dTag: lib.ElOnes_d(A.obj,m,n)
    elif A.tag == cTag: lib.ElOnes_c(A.obj,m,n)
    elif A.tag == zTag: lib.ElOnes_z(A.obj,m,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElOnesDist_i(A.obj,m,n)
    elif A.tag == sTag: lib.ElOnesDist_s(A.obj,m,n)
    elif A.tag == dTag: lib.ElOnesDist_d(A.obj,m,n)
    elif A.tag == cTag: lib.ElOnesDist_c(A.obj,m,n)
    elif A.tag == zTag: lib.ElOnesDist_z(A.obj,m,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# 1-2-1 matrix
# ------------
lib.ElOneTwoOne_i.argtypes = [c_void_p,iType]
lib.ElOneTwoOne_i.restype = c_uint
lib.ElOneTwoOne_s.argtypes = [c_void_p,iType]
lib.ElOneTwoOne_s.restype = c_uint
lib.ElOneTwoOne_d.argtypes = [c_void_p,iType]
lib.ElOneTwoOne_d.restype = c_uint
lib.ElOneTwoOne_c.argtypes = [c_void_p,iType]
lib.ElOneTwoOne_c.restype = c_uint
lib.ElOneTwoOne_z.argtypes = [c_void_p,iType]
lib.ElOneTwoOne_z.restype = c_uint
lib.ElOneTwoOneDist_i.argtypes = [c_void_p,iType]
lib.ElOneTwoOneDist_i.restype = c_uint
lib.ElOneTwoOneDist_s.argtypes = [c_void_p,iType]
lib.ElOneTwoOneDist_s.restype = c_uint
lib.ElOneTwoOneDist_d.argtypes = [c_void_p,iType]
lib.ElOneTwoOneDist_d.restype = c_uint
lib.ElOneTwoOneDist_c.argtypes = [c_void_p,iType]
lib.ElOneTwoOneDist_c.restype = c_uint
lib.ElOneTwoOneDist_z.argtypes = [c_void_p,iType]
lib.ElOneTwoOneDist_z.restype = c_uint
def OneTwoOne(A,n):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElOneTwoOne_i(A.obj,n)
    elif A.tag == sTag: lib.ElOneTwoOne_s(A.obj,n)
    elif A.tag == dTag: lib.ElOneTwoOne_d(A.obj,n)
    elif A.tag == cTag: lib.ElOneTwoOne_c(A.obj,n)
    elif A.tag == zTag: lib.ElOneTwoOne_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElOneTwoOneDist_i(A.obj,n)
    elif A.tag == sTag: lib.ElOneTwoOneDist_s(A.obj,n)
    elif A.tag == dTag: lib.ElOneTwoOneDist_d(A.obj,n)
    elif A.tag == cTag: lib.ElOneTwoOneDist_c(A.obj,n)
    elif A.tag == zTag: lib.ElOneTwoOneDist_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Parter
# ------
lib.ElParter_s.argtypes = [c_void_p,iType]
lib.ElParter_s.restype = c_uint
lib.ElParter_d.argtypes = [c_void_p,iType]
lib.ElParter_d.restype = c_uint
lib.ElParter_c.argtypes = [c_void_p,iType]
lib.ElParter_c.restype = c_uint
lib.ElParter_z.argtypes = [c_void_p,iType]
lib.ElParter_z.restype = c_uint
lib.ElParterDist_s.argtypes = [c_void_p,iType]
lib.ElParterDist_s.restype = c_uint
lib.ElParterDist_d.argtypes = [c_void_p,iType]
lib.ElParterDist_d.restype = c_uint
lib.ElParterDist_c.argtypes = [c_void_p,iType]
lib.ElParterDist_c.restype = c_uint
lib.ElParterDist_z.argtypes = [c_void_p,iType]
lib.ElParterDist_z.restype = c_uint
def Parter(A,n):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElParter_s(A.obj,n)
    elif A.tag == dTag: lib.ElParter_d(A.obj,n)
    elif A.tag == cTag: lib.ElParter_c(A.obj,n)
    elif A.tag == zTag: lib.ElParter_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElParterDist_s(A.obj,n)
    elif A.tag == dTag: lib.ElParterDist_d(A.obj,n)
    elif A.tag == cTag: lib.ElParterDist_c(A.obj,n)
    elif A.tag == zTag: lib.ElParterDist_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Pei
# ---
lib.ElPei_s.argtypes = [c_void_p,iType,sType]
lib.ElPei_s.restype = c_uint
lib.ElPei_d.argtypes = [c_void_p,iType,dType]
lib.ElPei_d.restype = c_uint
lib.ElPei_c.argtypes = [c_void_p,iType,cType]
lib.ElPei_c.restype = c_uint
lib.ElPei_z.argtypes = [c_void_p,iType,zType]
lib.ElPei_z.restype = c_uint
lib.ElPeiDist_s.argtypes = [c_void_p,iType,sType]
lib.ElPeiDist_s.restype = c_uint
lib.ElPeiDist_d.argtypes = [c_void_p,iType,dType]
lib.ElPeiDist_d.restype = c_uint
lib.ElPeiDist_c.argtypes = [c_void_p,iType,cType]
lib.ElPeiDist_c.restype = c_uint
lib.ElPeiDist_z.argtypes = [c_void_p,iType,zType]
lib.ElPeiDist_z.restype = c_uint
def Pei(A,n,alpha):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElPei_s(A.obj,n,alpha)
    elif A.tag == dTag: lib.ElPei_d(A.obj,n,alpha)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElPeiDist_s(A.obj,n,alpha)
    elif A.tag == dTag: lib.ElPeiDist_d(A.obj,n,alpha)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Redheffer
# ---------
lib.ElRedheffer_i.argtypes = [c_void_p,iType]
lib.ElRedheffer_i.restype = c_uint
lib.ElRedheffer_s.argtypes = [c_void_p,iType]
lib.ElRedheffer_s.restype = c_uint
lib.ElRedheffer_d.argtypes = [c_void_p,iType]
lib.ElRedheffer_d.restype = c_uint
lib.ElRedheffer_c.argtypes = [c_void_p,iType]
lib.ElRedheffer_c.restype = c_uint
lib.ElRedheffer_z.argtypes = [c_void_p,iType]
lib.ElRedheffer_z.restype = c_uint
lib.ElRedhefferDist_i.argtypes = [c_void_p,iType]
lib.ElRedhefferDist_i.restype = c_uint
lib.ElRedhefferDist_s.argtypes = [c_void_p,iType]
lib.ElRedhefferDist_s.restype = c_uint
lib.ElRedhefferDist_d.argtypes = [c_void_p,iType]
lib.ElRedhefferDist_d.restype = c_uint
lib.ElRedhefferDist_c.argtypes = [c_void_p,iType]
lib.ElRedhefferDist_c.restype = c_uint
lib.ElRedhefferDist_z.argtypes = [c_void_p,iType]
lib.ElRedhefferDist_z.restype = c_uint
def Redheffer(A,n):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElRedheffer_i(A.obj,n)
    elif A.tag == sTag: lib.ElRedheffer_s(A.obj,n)
    elif A.tag == dTag: lib.ElRedheffer_d(A.obj,n)
    elif A.tag == cTag: lib.ElRedheffer_c(A.obj,n)
    elif A.tag == zTag: lib.ElRedheffer_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElRedhefferDist_i(A.obj,n)
    elif A.tag == sTag: lib.ElRedhefferDist_s(A.obj,n)
    elif A.tag == dTag: lib.ElRedhefferDist_d(A.obj,n)
    elif A.tag == cTag: lib.ElRedhefferDist_c(A.obj,n)
    elif A.tag == zTag: lib.ElRedhefferDist_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Riffle
# ------
lib.ElRiffle_s.argtypes = [c_void_p,iType]
lib.ElRiffle_s.restype = c_uint
lib.ElRiffle_d.argtypes = [c_void_p,iType]
lib.ElRiffle_d.restype = c_uint
lib.ElRiffle_c.argtypes = [c_void_p,iType]
lib.ElRiffle_c.restype = c_uint
lib.ElRiffle_z.argtypes = [c_void_p,iType]
lib.ElRiffle_z.restype = c_uint
lib.ElRiffleDist_s.argtypes = [c_void_p,iType]
lib.ElRiffleDist_s.restype = c_uint
lib.ElRiffleDist_d.argtypes = [c_void_p,iType]
lib.ElRiffleDist_d.restype = c_uint
lib.ElRiffleDist_c.argtypes = [c_void_p,iType]
lib.ElRiffleDist_c.restype = c_uint
lib.ElRiffleDist_z.argtypes = [c_void_p,iType]
lib.ElRiffleDist_z.restype = c_uint
def Riffle(P,n):
  if type(P) is Matrix:
    if   P.tag == sTag: lib.ElRiffle_s(P.obj,n)
    elif P.tag == dTag: lib.ElRiffle_d(P.obj,n)
    elif P.tag == cTag: lib.ElRiffle_c(P.obj,n)
    elif P.tag == zTag: lib.ElRiffle_z(P.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(P) is DistMatrix:
    if   P.tag == sTag: lib.ElRiffleDist_s(P.obj,n)
    elif P.tag == dTag: lib.ElRiffleDist_d(P.obj,n)
    elif P.tag == cTag: lib.ElRiffleDist_c(P.obj,n)
    elif P.tag == zTag: lib.ElRiffleDist_z(P.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElRiffleStationary_s.argtypes = [c_void_p,iType]
lib.ElRiffleStationary_s.restype = c_uint
lib.ElRiffleStationary_d.argtypes = [c_void_p,iType]
lib.ElRiffleStationary_d.restype = c_uint
lib.ElRiffleStationary_c.argtypes = [c_void_p,iType]
lib.ElRiffleStationary_c.restype = c_uint
lib.ElRiffleStationary_z.argtypes = [c_void_p,iType]
lib.ElRiffleStationary_z.restype = c_uint
lib.ElRiffleStationaryDist_s.argtypes = [c_void_p,iType]
lib.ElRiffleStationaryDist_s.restype = c_uint
lib.ElRiffleStationaryDist_d.argtypes = [c_void_p,iType]
lib.ElRiffleStationaryDist_d.restype = c_uint
lib.ElRiffleStationaryDist_c.argtypes = [c_void_p,iType]
lib.ElRiffleStationaryDist_c.restype = c_uint
lib.ElRiffleStationaryDist_z.argtypes = [c_void_p,iType]
lib.ElRiffleStationaryDist_z.restype = c_uint
def RiffleStationary(P,n):
  if type(P) is Matrix:
    if   P.tag == sTag: lib.ElRiffleStationary_s(P.obj,n)
    elif P.tag == dTag: lib.ElRiffleStationary_d(P.obj,n)
    elif P.tag == cTag: lib.ElRiffleStationary_c(P.obj,n)
    elif P.tag == zTag: lib.ElRiffleStationary_z(P.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(P) is DistMatrix:
    if   P.tag == sTag: lib.ElRiffleStationaryDist_s(P.obj,n)
    elif P.tag == dTag: lib.ElRiffleStationaryDist_d(P.obj,n)
    elif P.tag == cTag: lib.ElRiffleStationaryDist_c(P.obj,n)
    elif P.tag == zTag: lib.ElRiffleStationaryDist_z(P.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElRiffleDecay_s.argtypes = [c_void_p,iType]
lib.ElRiffleDecay_s.restype = c_uint
lib.ElRiffleDecay_d.argtypes = [c_void_p,iType]
lib.ElRiffleDecay_d.restype = c_uint
lib.ElRiffleDecay_c.argtypes = [c_void_p,iType]
lib.ElRiffleDecay_c.restype = c_uint
lib.ElRiffleDecay_z.argtypes = [c_void_p,iType]
lib.ElRiffleDecay_z.restype = c_uint
lib.ElRiffleDecayDist_s.argtypes = [c_void_p,iType]
lib.ElRiffleDecayDist_s.restype = c_uint
lib.ElRiffleDecayDist_d.argtypes = [c_void_p,iType]
lib.ElRiffleDecayDist_d.restype = c_uint
lib.ElRiffleDecayDist_c.argtypes = [c_void_p,iType]
lib.ElRiffleDecayDist_c.restype = c_uint
lib.ElRiffleDecayDist_z.argtypes = [c_void_p,iType]
lib.ElRiffleDecayDist_z.restype = c_uint
def RiffleDecay(P,n):
  if type(P) is Matrix:
    if   P.tag == sTag: lib.ElRiffleDecay_s(P.obj,n)
    elif P.tag == dTag: lib.ElRiffleDecay_d(P.obj,n)
    elif P.tag == cTag: lib.ElRiffleDecay_c(P.obj,n)
    elif P.tag == zTag: lib.ElRiffleDecay_z(P.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(P) is DistMatrix:
    if   P.tag == sTag: lib.ElRiffleDecayDist_s(P.obj,n)
    elif P.tag == dTag: lib.ElRiffleDecayDist_d(P.obj,n)
    elif P.tag == cTag: lib.ElRiffleDecayDist_c(P.obj,n)
    elif P.tag == zTag: lib.ElRiffleDecayDist_z(P.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Ris
# ---
lib.ElRis_s.argtypes = [c_void_p,iType]
lib.ElRis_s.restype = c_uint
lib.ElRis_d.argtypes = [c_void_p,iType]
lib.ElRis_d.restype = c_uint
lib.ElRis_c.argtypes = [c_void_p,iType]
lib.ElRis_c.restype = c_uint
lib.ElRis_z.argtypes = [c_void_p,iType]
lib.ElRis_z.restype = c_uint
lib.ElRisDist_s.argtypes = [c_void_p,iType]
lib.ElRisDist_s.restype = c_uint
lib.ElRisDist_d.argtypes = [c_void_p,iType]
lib.ElRisDist_d.restype = c_uint
lib.ElRisDist_c.argtypes = [c_void_p,iType]
lib.ElRisDist_c.restype = c_uint
lib.ElRisDist_z.argtypes = [c_void_p,iType]
lib.ElRisDist_z.restype = c_uint
def Ris(A,n):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElRis_s(A.obj,n)
    elif A.tag == dTag: lib.ElRis_d(A.obj,n)
    elif A.tag == cTag: lib.ElRis_c(A.obj,n)
    elif A.tag == zTag: lib.ElRis_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElRisDist_s(A.obj,n)
    elif A.tag == dTag: lib.ElRisDist_d(A.obj,n)
    elif A.tag == cTag: lib.ElRisDist_c(A.obj,n)
    elif A.tag == zTag: lib.ElRisDist_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Toeplitz
# --------
lib.ElToeplitz_i.argtypes = [c_void_p,iType,iType,iType,POINTER(iType)]
lib.ElToeplitz_i.restype = c_uint
lib.ElToeplitz_s.argtypes = [c_void_p,iType,iType,iType,POINTER(sType)]
lib.ElToeplitz_s.restype = c_uint
lib.ElToeplitz_d.argtypes = [c_void_p,iType,iType,iType,POINTER(dType)]
lib.ElToeplitz_d.restype = c_uint
lib.ElToeplitz_c.argtypes = [c_void_p,iType,iType,iType,POINTER(cType)]
lib.ElToeplitz_c.restype = c_uint
lib.ElToeplitz_z.argtypes = [c_void_p,iType,iType,iType,POINTER(zType)]
lib.ElToeplitz_z.restype = c_uint
lib.ElToeplitzDist_i.argtypes = [c_void_p,iType,iType,iType,POINTER(iType)]
lib.ElToeplitzDist_i.restype = c_uint
lib.ElToeplitzDist_s.argtypes = [c_void_p,iType,iType,iType,POINTER(sType)]
lib.ElToeplitzDist_s.restype = c_uint
lib.ElToeplitzDist_d.argtypes = [c_void_p,iType,iType,iType,POINTER(dType)]
lib.ElToeplitzDist_d.restype = c_uint
lib.ElToeplitzDist_c.argtypes = [c_void_p,iType,iType,iType,POINTER(cType)]
lib.ElToeplitzDist_c.restype = c_uint
lib.ElToeplitzDist_z.argtypes = [c_void_p,iType,iType,iType,POINTER(zType)]
lib.ElToeplitzDist_z.restype = c_uint
def Toeplitz(A,m,n,a):
  aLen = len(a)
  aBuf = (TagToType(A.tag)*aLen)(*a) 
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElToeplitz_i(A.obj,m,n,aLen,aBuf)
    elif A.tag == sTag: lib.ElToeplitz_s(A.obj,m,n,aLen,aBuf)
    elif A.tag == dTag: lib.ElToeplitz_d(A.obj,m,n,aLen,aBuf)
    elif A.tag == cTag: lib.ElToeplitz_c(A.obj,m,n,aLen,aBuf)
    elif A.tag == zTag: lib.ElToeplitz_z(A.obj,m,n,aLen,aBuf)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElToeplitzDist_i(A.obj,m,n,aLen,aBuf)
    elif A.tag == sTag: lib.ElToeplitzDist_s(A.obj,m,n,aLen,aBuf)
    elif A.tag == dTag: lib.ElToeplitzDist_d(A.obj,m,n,aLen,aBuf)
    elif A.tag == cTag: lib.ElToeplitzDist_c(A.obj,m,n,aLen,aBuf)
    elif A.tag == zTag: lib.ElToeplitzDist_z(A.obj,m,n,aLen,aBuf)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Trefethen-Embree
# ----------------
lib.ElTrefethenEmbree_c.argtypes = [c_void_p,iType]
lib.ElTrefethenEmbree_c.restype = c_uint
lib.ElTrefethenEmbree_z.argtypes = [c_void_p,iType]
lib.ElTrefethenEmbree_z.restype = c_uint
lib.ElTrefethenEmbreeDist_c.argtypes = [c_void_p,iType]
lib.ElTrefethenEmbreeDist_c.restype = c_uint
lib.ElTrefethenEmbreeDist_z.argtypes = [c_void_p,iType]
lib.ElTrefethenEmbreeDist_z.restype = c_uint
def TrefethenEmbree(A,n):
  if type(A) is Matrix: 
    if   A.tag == cTag: lib.ElTrefethenEmbree_c(A.obj,n)
    elif A.tag == zTag: lib.ElTrefethenEmbree_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElTrefethenEmbreeDist_c(A.obj,n)
    elif A.tag == zTag: lib.ElTrefethenEmbreeDist_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Triangle
# --------
lib.ElTriangle_c.argtypes = [c_void_p,iType]
lib.ElTriangle_c.restype = c_uint
lib.ElTriangle_z.argtypes = [c_void_p,iType]
lib.ElTriangle_z.restype = c_uint
lib.ElTriangleDist_c.argtypes = [c_void_p,iType]
lib.ElTriangleDist_c.restype = c_uint
lib.ElTriangleDist_z.argtypes = [c_void_p,iType]
lib.ElTriangleDist_z.restype = c_uint
def Triangle(A,n):
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElTriangle_c(A.obj,n)
    elif A.tag == zTag: lib.ElTriangle_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElTriangleDist_c(A.obj,n)
    elif A.tag == zTag: lib.ElTriangleDist_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# TriW
# ----
lib.ElTriW_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElTriW_i.restype = c_uint
lib.ElTriW_s.argtypes = [c_void_p,iType,sType,iType]
lib.ElTriW_s.restype = c_uint
lib.ElTriW_d.argtypes = [c_void_p,iType,dType,iType]
lib.ElTriW_d.restype = c_uint
lib.ElTriW_c.argtypes = [c_void_p,iType,cType,iType]
lib.ElTriW_c.restype = c_uint
lib.ElTriW_z.argtypes = [c_void_p,iType,zType,iType]
lib.ElTriW_z.restype = c_uint
lib.ElTriWDist_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElTriWDist_i.restype = c_uint
lib.ElTriWDist_s.argtypes = [c_void_p,iType,sType,iType]
lib.ElTriWDist_s.restype = c_uint
lib.ElTriWDist_d.argtypes = [c_void_p,iType,dType,iType]
lib.ElTriWDist_d.restype = c_uint
lib.ElTriWDist_c.argtypes = [c_void_p,iType,cType,iType]
lib.ElTriWDist_c.restype = c_uint
lib.ElTriWDist_z.argtypes = [c_void_p,iType,zType,iType]
lib.ElTriWDist_z.restype = c_uint
def TriW(A,n,alpha,k):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElTriW_i(A.obj,n,alpha,k)
    elif A.tag == sTag: lib.ElTriW_s(A.obj,n,alpha,k)
    elif A.tag == dTag: lib.ElTriW_d(A.obj,n,alpha,k)
    elif A.tag == cTag: lib.ElTriW_c(A.obj,n,alpha,k)
    elif A.tag == zTag: lib.ElTriW_z(A.obj,n,alpha,k)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElTriWDist_i(A.obj,n,alpha,k)
    elif A.tag == sTag: lib.ElTriWDist_s(A.obj,n,alpha,k)
    elif A.tag == dTag: lib.ElTriWDist_d(A.obj,n,alpha,k)
    elif A.tag == cTag: lib.ElTriWDist_c(A.obj,n,alpha,k)
    elif A.tag == zTag: lib.ElTriWDist_z(A.obj,n,alpha,k)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Walsh
# -----
lib.ElWalsh_i.argtypes = [c_void_p,iType,bType]
lib.ElWalsh_i.restype = c_uint
lib.ElWalsh_s.argtypes = [c_void_p,iType,bType]
lib.ElWalsh_s.restype = c_uint
lib.ElWalsh_d.argtypes = [c_void_p,iType,bType]
lib.ElWalsh_d.restype = c_uint
lib.ElWalsh_c.argtypes = [c_void_p,iType,bType]
lib.ElWalsh_c.restype = c_uint
lib.ElWalsh_z.argtypes = [c_void_p,iType,bType]
lib.ElWalsh_z.restype = c_uint
lib.ElWalshDist_i.argtypes = [c_void_p,iType,bType]
lib.ElWalshDist_i.restype = c_uint
lib.ElWalshDist_s.argtypes = [c_void_p,iType,bType]
lib.ElWalshDist_s.restype = c_uint
lib.ElWalshDist_d.argtypes = [c_void_p,iType,bType]
lib.ElWalshDist_d.restype = c_uint
lib.ElWalshDist_c.argtypes = [c_void_p,iType,bType]
lib.ElWalshDist_c.restype = c_uint
lib.ElWalshDist_z.argtypes = [c_void_p,iType,bType]
lib.ElWalshDist_z.restype = c_uint
def Walsh(A,k,binary=False):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElWalsh_i(A.obj,k,binary)
    elif A.tag == sTag: lib.ElWalsh_s(A.obj,k,binary)
    elif A.tag == dTag: lib.ElWalsh_d(A.obj,k,binary)
    elif A.tag == cTag: lib.ElWalsh_c(A.obj,k,binary)
    elif A.tag == zTag: lib.ElWalsh_z(A.obj,k,binary)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElWalshDist_i(A.obj,k,binary)
    elif A.tag == sTag: lib.ElWalshDist_s(A.obj,k,binary)
    elif A.tag == dTag: lib.ElWalshDist_d(A.obj,k,binary)
    elif A.tag == cTag: lib.ElWalshDist_c(A.obj,k,binary)
    elif A.tag == zTag: lib.ElWalshDist_z(A.obj,k,binary)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Whale
# -----
lib.ElWhale_c.argtypes = [c_void_p,iType]
lib.ElWhale_c.restype = c_uint
lib.ElWhale_z.argtypes = [c_void_p,iType]
lib.ElWhale_z.restype = c_uint
lib.ElWhaleDist_c.argtypes = [c_void_p,iType]
lib.ElWhaleDist_c.restype = c_uint
lib.ElWhaleDist_z.argtypes = [c_void_p,iType]
lib.ElWhaleDist_z.restype = c_uint
def Whale(A,n):
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElWhale_c(A.obj,n)
    elif A.tag == zTag: lib.ElWhale_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElWhaleDist_c(A.obj,n)
    elif A.tag == zTag: lib.ElWhaleDist_z(A.obj,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Wilkinson
# ---------
lib.ElWilkinson_i.argtypes = [c_void_p,iType]
lib.ElWilkinson_i.restype = c_uint
lib.ElWilkinson_s.argtypes = [c_void_p,iType]
lib.ElWilkinson_s.restype = c_uint
lib.ElWilkinson_d.argtypes = [c_void_p,iType]
lib.ElWilkinson_d.restype = c_uint
lib.ElWilkinson_c.argtypes = [c_void_p,iType]
lib.ElWilkinson_c.restype = c_uint
lib.ElWilkinson_z.argtypes = [c_void_p,iType]
lib.ElWilkinson_z.restype = c_uint
lib.ElWilkinsonDist_i.argtypes = [c_void_p,iType]
lib.ElWilkinsonDist_i.restype = c_uint
lib.ElWilkinsonDist_s.argtypes = [c_void_p,iType]
lib.ElWilkinsonDist_s.restype = c_uint
lib.ElWilkinsonDist_d.argtypes = [c_void_p,iType]
lib.ElWilkinsonDist_d.restype = c_uint
lib.ElWilkinsonDist_c.argtypes = [c_void_p,iType]
lib.ElWilkinsonDist_c.restype = c_uint
lib.ElWilkinsonDist_z.argtypes = [c_void_p,iType]
lib.ElWilkinsonDist_z.restype = c_uint
def Wilkinson(A,k):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElWilkinson_i(A.obj,k)
    elif A.tag == sTag: lib.ElWilkinson_s(A.obj,k)
    elif A.tag == dTag: lib.ElWilkinson_d(A.obj,k)
    elif A.tag == cTag: lib.ElWilkinson_c(A.obj,k)
    elif A.tag == zTag: lib.ElWilkinson_z(A.obj,k)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElWilkinsonDist_i(A.obj,k)
    elif A.tag == sTag: lib.ElWilkinsonDist_s(A.obj,k)
    elif A.tag == dTag: lib.ElWilkinsonDist_d(A.obj,k)
    elif A.tag == cTag: lib.ElWilkinsonDist_c(A.obj,k)
    elif A.tag == zTag: lib.ElWilkinsonDist_z(A.obj,k)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Zeros
# -----
lib.ElZeros_i.argtypes = [c_void_p,iType,iType]
lib.ElZeros_i.restype = c_uint
lib.ElZeros_s.argtypes = [c_void_p,iType,iType]
lib.ElZeros_s.restype = c_uint
lib.ElZeros_d.argtypes = [c_void_p,iType,iType]
lib.ElZeros_d.restype = c_uint
lib.ElZeros_c.argtypes = [c_void_p,iType,iType]
lib.ElZeros_c.restype = c_uint
lib.ElZeros_z.argtypes = [c_void_p,iType,iType]
lib.ElZeros_z.restype = c_uint
lib.ElZerosDist_i.argtypes = [c_void_p,iType,iType]
lib.ElZerosDist_i.restype = c_uint
lib.ElZerosDist_s.argtypes = [c_void_p,iType,iType]
lib.ElZerosDist_s.restype = c_uint
lib.ElZerosDist_d.argtypes = [c_void_p,iType,iType]
lib.ElZerosDist_d.restype = c_uint
lib.ElZerosDist_c.argtypes = [c_void_p,iType,iType]
lib.ElZerosDist_c.restype = c_uint
lib.ElZerosDist_z.argtypes = [c_void_p,iType,iType]
lib.ElZerosDist_z.restype = c_uint
def Zeros(A,m,n):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElZeros_i(A.obj,m,n)
    elif A.tag == sTag: lib.ElZeros_s(A.obj,m,n)
    elif A.tag == dTag: lib.ElZeros_d(A.obj,m,n)
    elif A.tag == cTag: lib.ElZeros_c(A.obj,m,n)
    elif A.tag == zTag: lib.ElZeros_z(A.obj,m,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElZerosDist_i(A.obj,m,n)
    elif A.tag == sTag: lib.ElZerosDist_s(A.obj,m,n)
    elif A.tag == dTag: lib.ElZerosDist_d(A.obj,m,n)
    elif A.tag == cTag: lib.ElZerosDist_c(A.obj,m,n)
    elif A.tag == zTag: lib.ElZerosDist_z(A.obj,m,n)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Random
# ======

# Gaussian
# --------
lib.ElGaussian_s.argtypes = [c_void_p,iType,iType,sType,sType]
lib.ElGaussian_s.restype = c_uint
lib.ElGaussian_d.argtypes = [c_void_p,iType,iType,dType,dType]
lib.ElGaussian_d.restype = c_uint
lib.ElGaussian_c.argtypes = [c_void_p,iType,iType,cType,sType]
lib.ElGaussian_c.restype = c_uint
lib.ElGaussian_z.argtypes = [c_void_p,iType,iType,zType,dType]
lib.ElGaussian_z.restype = c_uint
lib.ElGaussianDist_s.argtypes = [c_void_p,iType,iType,sType,sType]
lib.ElGaussianDist_s.restype = c_uint
lib.ElGaussianDist_d.argtypes = [c_void_p,iType,iType,dType,dType]
lib.ElGaussianDist_d.restype = c_uint
lib.ElGaussianDist_c.argtypes = [c_void_p,iType,iType,cType,sType]
lib.ElGaussianDist_c.restype = c_uint
lib.ElGaussianDist_z.argtypes = [c_void_p,iType,iType,zType,dType]
lib.ElGaussian_z.restype = c_uint
def Gaussian(A,m,n,meanPre=0,stddev=1):
  mean = TagToType(A.tag)(meanPre)
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElGaussian_s(A.obj,m,n,mean,stddev)
    elif A.tag == dTag: lib.ElGaussian_d(A.obj,m,n,mean,stddev)
    elif A.tag == cTag: lib.ElGaussian_c(A.obj,m,n,mean,stddev)
    elif A.tag == zTag: lib.ElGaussian_z(A.obj,m,n,mean,stddev)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElGaussianDist_s(A.obj,m,n,mean,stddev)
    elif A.tag == dTag: lib.ElGaussianDist_d(A.obj,m,n,mean,stddev)
    elif A.tag == cTag: lib.ElGaussianDist_c(A.obj,m,n,mean,stddev)
    elif A.tag == zTag: lib.ElGaussianDist_z(A.obj,m,n,mean,stddev)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Uniform
# -------
lib.ElUniform_i.argtypes = [c_void_p,iType,iType,iType,iType]
lib.ElUniform_i.restype = c_uint
lib.ElUniform_s.argtypes = [c_void_p,iType,iType,sType,sType]
lib.ElUniform_s.restype = c_uint
lib.ElUniform_d.argtypes = [c_void_p,iType,iType,dType,dType]
lib.ElUniform_d.restype = c_uint
lib.ElUniform_c.argtypes = [c_void_p,iType,iType,cType,sType]
lib.ElUniform_c.restype = c_uint
lib.ElUniform_z.argtypes = [c_void_p,iType,iType,zType,dType]
lib.ElUniform_z.restype = c_uint
lib.ElUniformDist_i.argtypes = [c_void_p,iType,iType,iType,iType]
lib.ElUniformDist_i.restype = c_uint
lib.ElUniformDist_s.argtypes = [c_void_p,iType,iType,sType,sType]
lib.ElUniformDist_s.restype = c_uint
lib.ElUniformDist_d.argtypes = [c_void_p,iType,iType,dType,dType]
lib.ElUniformDist_d.restype = c_uint
lib.ElUniformDist_c.argtypes = [c_void_p,iType,iType,cType,sType]
lib.ElUniformDist_c.restype = c_uint
lib.ElUniformDist_z.argtypes = [c_void_p,iType,iType,zType,dType]
lib.ElUniformDist_z.restype = c_uint
def Uniform(A,m,n,centerPre=0,radius=1):
  center = TagToType(A.tag)(centerPre) 
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElUniform_i(A.obj,m,n,center,radius)
    elif A.tag == sTag: lib.ElUniform_s(A.obj,m,n,center,radius)
    elif A.tag == dTag: lib.ElUniform_d(A.obj,m,n,center,radius)
    elif A.tag == cTag: lib.ElUniform_c(A.obj,m,n,center,radius)
    elif A.tag == zTag: lib.ElUniform_z(A.obj,m,n,center,radius)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElUniformDist_i(A.obj,m,n,center,radius)
    elif A.tag == sTag: lib.ElUniformDist_s(A.obj,m,n,center,radius)
    elif A.tag == dTag: lib.ElUniformDist_d(A.obj,m,n,center,radius)
    elif A.tag == cTag: lib.ElUniformDist_c(A.obj,m,n,center,radius)
    elif A.tag == zTag: lib.ElUniformDist_z(A.obj,m,n,center,radius)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Uniform Helmholtz Green's
# -------------------------
lib.ElUniformHelmholtzGreens_c.argtypes = [c_void_p,iType,sType]
lib.ElUniformHelmholtzGreens_c.restype = c_uint
lib.ElUniformHelmholtzGreens_z.argtypes = [c_void_p,iType,dType]
lib.ElUniformHelmholtzGreens_z.restype = c_uint
lib.ElUniformHelmholtzGreensDist_c.argtypes = [c_void_p,iType,sType]
lib.ElUniformHelmholtzGreensDist_c.restype = c_uint
lib.ElUniformHelmholtzGreensDist_z.argtypes = [c_void_p,iType,dType]
lib.ElUniformHelmholtzGreensDist_z.restype = c_uint
def UniformHelmholtzGreens(A,n,lamb):
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElUniformHelmholtzGreens_c(A.obj,n,lamb)
    elif A.tag == zTag: lib.ElUniformHelmholtzGreens_z(A.obj,n,lamb)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElUniformHelmholtzGreensDist_c(A.obj,n,lamb)
    elif A.tag == zTag: lib.ElUniformHelmholtzGreensDist_z(A.obj,n,lamb)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Wigner
# ------
lib.ElWigner_s.argtypes = [c_void_p,iType,sType,sType]
lib.ElWigner_s.restype = c_uint
lib.ElWigner_d.argtypes = [c_void_p,iType,dType,dType]
lib.ElWigner_d.restype = c_uint
lib.ElWigner_c.argtypes = [c_void_p,iType,cType,sType]
lib.ElWigner_c.restype = c_uint
lib.ElWigner_z.argtypes = [c_void_p,iType,zType,dType]
lib.ElWigner_z.restype = c_uint
lib.ElWignerDist_s.argtypes = [c_void_p,iType,sType,sType]
lib.ElWignerDist_s.restype = c_uint
lib.ElWignerDist_d.argtypes = [c_void_p,iType,dType,dType]
lib.ElWignerDist_d.restype = c_uint
lib.ElWignerDist_c.argtypes = [c_void_p,iType,cType,sType]
lib.ElWignerDist_c.restype = c_uint
lib.ElWignerDist_z.argtypes = [c_void_p,iType,zType,dType]
lib.ElWignerDist_z.restype = c_uint
def Wigner(A,n,meanPre=0,stddev=1):
  mean = TagToType(A.tag)(meanPre) 
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElWigner_s(A.obj,n,mean,stddev)
    elif A.tag == dTag: lib.ElWigner_d(A.obj,n,mean,stddev)
    elif A.tag == cTag: lib.ElWigner_c(A.obj,n,mean,stddev)
    elif A.tag == zTag: lib.ElWigner_z(A.obj,n,mean,stddev)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElWignerDist_s(A.obj,n,mean,stddev)
    elif A.tag == dTag: lib.ElWignerDist_d(A.obj,n,mean,stddev)
    elif A.tag == cTag: lib.ElWignerDist_c(A.obj,n,mean,stddev)
    elif A.tag == zTag: lib.ElWignerDist_z(A.obj,n,mean,stddev)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')
