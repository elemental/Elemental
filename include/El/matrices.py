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
  if type(A) is Matrix:
    if   A.tag == sTag: 
      xBuf = (sType*xLen)(*x)
      yBuf = (sType*yLen)(*y)
      lib.ElCauchy_s(A.obj,xLen,xBuf,yLen,yBuf)
    elif A.tag == dTag:
      xBuf = (dType*xLen)(*x)
      yBuf = (dType*yLen)(*y)
      lib.ElCauchy_d(A.obj,xLen,xBuf,yLen,yBuf)
    elif A.tag == cTag:
      xBuf = (cType*xLen)(*x)
      yBuf = (cType*yLen)(*y)
      lib.ElCauchy_c(A.obj,xLen,xBuf,yLen,yBuf)
    elif A.tag == zTag:
      xBuf = (zType*xLen)(*x)
      yBuf = (zType*yLen)(*y)
      lib.ElCauchy_z(A.obj,xLen,xBuf,yLen,yBuf)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      xBuf = (sType*xLen)(*x)
      yBuf = (sType*yLen)(*y)
      lib.ElCauchyDist_s(A.obj,xLen,xBuf,yLen,yBuf)
    elif A.tag == dTag:
      xBuf = (dType*xLen)(*x)
      yBuf = (dType*yLen)(*y)
      lib.ElCauchyDist_d(A.obj,xLen,xBuf,yLen,yBuf)
    elif A.tag == cTag:
      xBuf = (cType*xLen)(*x)
      yBuf = (cType*yLen)(*y)
      lib.ElCauchyDist_c(A.obj,xLen,xBuf,yLen,yBuf)
    elif A.tag == zTag:
      xBuf = (zType*xLen)(*x)
      yBuf = (zType*yLen)(*y)
      lib.ElCauchyDist_z(A.obj,xLen,xBuf,yLen,yBuf)
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
  if type(A) is Matrix:
    if   A.tag == sTag: 
      rBuf = (sType*rLen)(*r)
      sBuf = (sType*sLen)(*s)
      xBuf = (sType*xLen)(*x)
      yBuf = (sType*yLen)(*y)
      lib.ElCauchyLike_s(A.obj,rLen,rBuf,sLen,sBuf,xLen,xBuf,yLen,yBuf)
    elif A.tag == dTag:
      rBuf = (dType*rLen)(*r)
      sBuf = (dType*sLen)(*s)
      xBuf = (dType*xLen)(*x)
      yBuf = (dType*yLen)(*y)
      lib.ElCauchyLike_d(A.obj,rLen,rBuf,sLen,sBuf,xLen,xBuf,yLen,yBuf)
    elif A.tag == cTag:
      rBuf = (cType*rLen)(*r)
      sBuf = (cType*sLen)(*s)
      xBuf = (cType*xLen)(*x)
      yBuf = (cType*yLen)(*y)
      lib.ElCauchyLike_c(A.obj,rLen,rBuf,sLen,sBuf,xLen,xBuf,yLen,yBuf)
    elif A.tag == zTag:
      rBuf = (zType*rLen)(*r)
      sBuf = (zType*sLen)(*s)
      xBuf = (zType*xLen)(*x)
      yBuf = (zType*yLen)(*y)
      lib.ElCauchyLike_z(A.obj,rLen,rBuf,sLen,sBuf,xLen,xBuf,yLen,yBuf)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      rBuf = (sType*rLen)(*r)
      sBuf = (sType*sLen)(*s)
      xBuf = (sType*xLen)(*x)
      yBuf = (sType*yLen)(*y)
      lib.ElCauchyLikeDist_s(A.obj,rLen,rBuf,sLen,sBuf,xLen,xBuf,yLen,yBuf)
    elif A.tag == dTag:
      rBuf = (dType*rLen)(*r)
      sBuf = (dType*sLen)(*s)
      xBuf = (dType*xLen)(*x)
      yBuf = (dType*yLen)(*y)
      lib.ElCauchyLikeDist_d(A.obj,rLen,rBuf,sLen,sBuf,xLen,xBuf,yLen,yBuf)
    elif A.tag == cTag:
      rBuf = (cType*rLen)(*r)
      sBuf = (cType*sLen)(*s)
      xBuf = (cType*xLen)(*x)
      yBuf = (cType*yLen)(*y)
      lib.ElCauchyLikeDist_c(A.obj,rLen,rBuf,sLen,sBuf,xLen,xBuf,yLen,yBuf)
    elif A.tag == zTag:
      rBuf = (zType*rLen)(*r)
      sBuf = (zType*sLen)(*s)
      xBuf = (zType*xLen)(*x)
      yBuf = (zType*yLen)(*y)
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
  if type(A) is Matrix:
    if   A.tag == iTag: 
      aBuf = (iType*aLen)(*a)
      lib.ElCirculant_i(A.obj,aLen,aBuf)
    elif A.tag == sTag:
      aBuf = (sType*aLen)(*a)
      lib.ElCirculant_s(A.obj,aLen,aBuf)
    elif A.tag == dTag:
      aBuf = (dType*aLen)(*a)
      lib.ElCirculant_d(A.obj,aLen,aBuf)
    elif A.tag == cTag:
      aBuf = (cType*aLen)(*a)
      lib.ElCirculant_c(A.obj,aLen,aBuf)
    elif A.tag == zTag:
      aBuf = (zType*aLen)(*a)
      lib.ElCirculant_z(A.obj,aLen,aBuf)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: 
      aBuf = (iType*aLen)(*a)
      lib.ElCirculantDist_i(A.obj,aLen,aBuf)
    elif A.tag == sTag:
      aBuf = (sType*aLen)(*a)
      lib.ElCirculantDist_s(A.obj,aLen,aBuf)
    elif A.tag == dTag:
      aBuf = (dType*aLen)(*a)
      lib.ElCirculantDist_d(A.obj,aLen,aBuf)
    elif A.tag == cTag:
      aBuf = (cType*aLen)(*a)
      lib.ElCirculantDist_c(A.obj,aLen,aBuf)
    elif A.tag == zTag:
      aBuf = (zType*aLen)(*a)
      lib.ElCirculantDist_z(A.obj,aLen,aBuf)
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
  if type(A) is Matrix:
    if   A.tag == iTag: 
      dBuf = (iType*dLen)(*d)
      lib.ElDiagonal_i(A.obj,dLen,dBuf)
    elif A.tag == sTag:
      dBuf = (sType*dLen)(*d)
      lib.ElDiagonal_s(A.obj,dLen,dBuf)
    elif A.tag == dTag:
      dBuf = (dType*dLen)(*d)
      lib.ElDiagonal_d(A.obj,dLen,dBuf)
    elif A.tag == cTag:
      dBuf = (cType*dLen)(*d)
      lib.ElDiagonal_c(A.obj,dLen,dBuf)
    elif A.tag == zTag:
      dBuf = (zType*zLen)(*d)
      lib.ElDiagonal_z(A.obj,dLen,dBuf)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: 
      dBuf = (iType*dLen)(*d)
      lib.ElDiagonalDist_i(A.obj,dLen,dBuf)
    elif A.tag == sTag:
      dBuf = (sType*dLen)(*d)
      lib.ElDiagonalDist_s(A.obj,dLen,dBuf)
    elif A.tag == dTag:
      dBuf = (dType*dLen)(*d)
      lib.ElDiagonalDist_d(A.obj,dLen,dBuf)
    elif A.tag == cTag:
      dBuf = (cType*dLen)(*d)
      lib.ElDiagonalDist_c(A.obj,dLen,dBuf)
    elif A.tag == zTag:
      dBuf = (zType*zLen)(*d)
      lib.ElDiagonalDist_z(A.obj,dLen,dBuf)
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
  if type(A) is Matrix:
    if   A.tag == cTag:
      cPhase = CFUNCTYPE(sType,iType,iType)(phase)
      lib.ElEgorov_c(A.obj,cPhase,n)
    elif A.tag == zTag:
      cPhase = CFUNCTYPE(dType,iType,iType)(phase)
      lib.ElEgorov_z(A.obj,cPhase,n)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == cTag:
      cPhase = CFUNCTYPE(sType,iType,iType)(phase)
      lib.ElEgorovDist_c(A.obj,cPhase,n)
    elif A.tag == zTag:
      cPhase = CFUNCTYPE(dType,iType,iType)(phase)
      lib.ElEgorovDist_z(A.obj,cPhase,n)
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
  if type(A) is Matrix:
    if   A.tag == sTag: 
      cBuf = (sType*cLen)(*c)
      lib.ElFiedler_s(A.obj,cLen,cBuf)
    elif A.tag == dTag:
      cBuf = (dType*cLen)(*c)
      lib.ElFiedler_d(A.obj,cLen,cBuf)
    elif A.tag == cTag:
      cBuf = (cType*cLen)(*c)
      lib.ElFiedler_c(A.obj,cLen,cBuf)
    elif A.tag == zTag:
      cBuf = (zType*cLen)(*c)
      lib.ElFiedler_z(A.obj,cLen,cBuf)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      cBuf = (sType*cLen)(*c)
      lib.ElFiedlerDist_s(A.obj,cLen,cBuf)
    elif A.tag == dTag:
      cBuf = (dType*cLen)(*c)
      lib.ElFiedlerDist_d(A.obj,cLen,cBuf)
    elif A.tag == cTag:
      cBuf = (cType*cLen)(*c)
      lib.ElFiedlerDist_c(A.obj,cLen,cBuf)
    elif A.tag == zTag:
      cBuf = (zType*cLen)(*c)
      lib.ElFiedlerDist_z(A.obj,cLen,cBuf)
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

# GCD matrix
# ----------
# TODO

# Gear matrix
# -----------
# TODO

# GEPP Growth
# -----------
# TODO

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
def Gaussian(A,m,n,mean,stddev):
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
