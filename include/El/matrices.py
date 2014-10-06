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

# Random
# ======
