#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from ..core import *

# BLAS 2
# ======

# Gemv
# ----
lib.ElGemv_s.argtypes = [c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElGemv_s.restype = c_uint
lib.ElGemv_d.argtypes = [c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElGemv_d.restype = c_uint
lib.ElGemv_c.argtypes = [c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElGemv_c.restype = c_uint
lib.ElGemv_z.argtypes = [c_uint,zType,c_void_p,c_void_p,zType,c_void_p]
lib.ElGemv_z.restype = c_uint
lib.ElGemvDist_s.argtypes = [c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElGemvDist_s.restype = c_uint
lib.ElGemvDist_d.argtypes = [c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElGemvDist_d.restype = c_uint
lib.ElGemvDist_c.argtypes = [c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElGemvDist_c.restype = c_uint
lib.ElGemvDist_z.argtypes = [c_uint,zType,c_void_p,c_void_p,zType,c_void_p]
lib.ElGemvDist_z.restype = c_uint
def Gemv(orient,alphaPre,A,x,betaPre,y):
  if type(A) is not type(x) or type(x) is not type(y):
    raise Exception('Types of {A,x,y} must match')
  if A.tag != x.tag or x.tag != y.tag:
    raise Exception('Datatypes of {A,x,y} must match')
  alpha = TagToType(A.tag)(alphaPre)
  beta = TagToType(A.tag)(betaPre)
  args = [orient,alpha,A.obj,x.obj,beta,y.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElGemv_s(*args)
    elif A.tag == dTag: lib.ElGemv_d(*args)
    elif A.tag == cTag: lib.ElGemv_c(*args)
    elif A.tag == zTag: lib.ElGemv_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElGemvDist_s(*args)
    elif A.tag == dTag: lib.ElGemvDist_d(*args)
    elif A.tag == cTag: lib.ElGemvDist_c(*args)
    elif A.tag == zTag: lib.ElGemvDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Ger
# ---
lib.ElGer_s.argtypes = [sType,c_void_p,c_void_p,c_void_p]
lib.ElGer_s.restype = c_uint
lib.ElGer_d.argtypes = [dType,c_void_p,c_void_p,c_void_p]
lib.ElGer_d.restype = c_uint
lib.ElGer_c.argtypes = [cType,c_void_p,c_void_p,c_void_p]
lib.ElGer_c.restype = c_uint
lib.ElGer_z.argtypes = [zType,c_void_p,c_void_p,c_void_p]
lib.ElGer_z.restype = c_uint
lib.ElGerDist_s.argtypes = [sType,c_void_p,c_void_p,c_void_p]
lib.ElGerDist_s.restype = c_uint
lib.ElGerDist_d.argtypes = [dType,c_void_p,c_void_p,c_void_p]
lib.ElGerDist_d.restype = c_uint
lib.ElGerDist_c.argtypes = [cType,c_void_p,c_void_p,c_void_p]
lib.ElGerDist_c.restype = c_uint
lib.ElGerDist_z.argtypes = [zType,c_void_p,c_void_p,c_void_p]
lib.ElGerDist_z.restype = c_uint
def Ger(alphaPre,A,x,y):
  if type(A) is not type(x) or type(x) is not type(y):
    raise Exception('Types of {A,x,y} must match')
  if A.tag != x.tag or x.tag != y.tag:
    raise Exception('Datatypes of {A,x,y} must match')
  alpha = TagToType(A.tag)(alphaPre)
  args = [alpha,A.obj,x.obj,y.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElGer_s(*args)
    elif A.tag == dTag: lib.ElGer_d(*args)
    elif A.tag == cTag: lib.ElGer_c(*args)
    elif A.tag == zTag: lib.ElGer_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElGerDist_s(*args)
    elif A.tag == dTag: lib.ElGerDist_d(*args)
    elif A.tag == cTag: lib.ElGerDist_c(*args)
    elif A.tag == zTag: lib.ElGerDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Geru
# ----
lib.ElGeru_c.argtypes = [cType,c_void_p,c_void_p,c_void_p]
lib.ElGeru_c.restype = c_uint
lib.ElGeru_z.argtypes = [zType,c_void_p,c_void_p,c_void_p]
lib.ElGeru_z.restype = c_uint
lib.ElGeruDist_c.argtypes = [cType,c_void_p,c_void_p,c_void_p]
lib.ElGeruDist_c.restype = c_uint
lib.ElGeruDist_z.argtypes = [zType,c_void_p,c_void_p,c_void_p]
lib.ElGeruDist_z.restype = c_uint
def Geru(alphaPre,A,x,y):
  if type(A) is not type(x) or type(x) is not type(y):
    raise Exception('Types of {A,x,y} must match')
  if A.tag != x.tag or x.tag != y.tag:
    raise Exception('Datatypes of {A,x,y} must match')
  alpha = TagToType(A.tag)(alphaPre)
  args = [alpha,A.obj,x.obj,y.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElGer_s(*args)
    elif A.tag == dTag: lib.ElGer_d(*args)
    elif A.tag == cTag: lib.ElGeru_c(*args)
    elif A.tag == zTag: lib.ElGeru_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElGerDist_s(*args)
    elif A.tag == dTag: lib.ElGerDist_d(*args)
    elif A.tag == cTag: lib.ElGeruDist_c(*args)
    elif A.tag == zTag: lib.ElGeruDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Symv/Hemv
# ---------
lib.ElSymv_s.argtypes = [c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElSymv_s.restype = c_uint
lib.ElSymv_d.argtypes = [c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElSymv_d.restype = c_uint
lib.ElSymv_c.argtypes = [c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElSymv_c.restype = c_uint
lib.ElSymv_z.argtypes = [c_uint,zType,c_void_p,c_void_p,zType,c_void_p]
lib.ElSymv_z.restype = c_uint
lib.ElSymvDist_s.argtypes = [c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElSymvDist_s.restype = c_uint
lib.ElSymvDist_d.argtypes = [c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElSymvDist_d.restype = c_uint
lib.ElSymvDist_c.argtypes = [c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElSymvDist_c.restype = c_uint
lib.ElSymvDist_z.argtypes = [c_uint,zType,c_void_p,c_void_p,zType,c_void_p]
lib.ElSymvDist_z.restype = c_uint

lib.ElHemv_c.argtypes = [c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElHemv_c.restype = c_uint
lib.ElHemv_z.argtypes = [c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElHemv_z.restype = c_uint
lib.ElHemvDist_c.argtypes = [c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElHemvDist_c.restype = c_uint
lib.ElHemvDist_z.argtypes = [c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElHemvDist_z.restype = c_uint

def Symv(uplo,alphaPre,A,x,betaPre,y,conj=False):
  if type(A) is not type(x) or type(x) is not type(y):
    raise Exception('Types of {A,x,y} must match')
  if A.tag != x.tag or x.tag != y.tag:
    raise Exception('Datatypes of {A,x,y} must match')
  alpha = TagToType(A.tag)(alphaPre)
  beta = TagToType(A.tag)(betaPre)
  args = [uplo,alpha,A.obj,x.obj,beta,y.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymv_s(*args)
    elif A.tag == dTag: lib.ElSymv_d(*args)
    elif A.tag == cTag:
      if conj: lib.ElHemv_c(uplo,alpha.real,A.obj,x.obj,beta.real,y.obj)
      else:    lib.ElSymv_c(*args)
    elif A.tag == zTag: 
      if conj: lib.ElHemv_z(uplo,alpha.real,A.obj,x.obj,beta.real,y.obj)
      else:    lib.ElSymv_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymvDist_s(*args)
    elif A.tag == dTag: lib.ElSymvDist_d(*args)
    elif A.tag == cTag:
      if conj: lib.ElHemvDist_c(uplo,alpha.real,A.obj,x.obj,beta.real,y.obj)
      else:    lib.ElSymvDist_c(*args)
    elif A.tag == zTag: 
      if conj: lib.ElHemvDist_z(uplo,alpha.real,A.obj,x.obj,beta.real,y.obj)
      else:    lib.ElSymvDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

def Hemv(uplo,alpha,A,x,beta,y):
  Symv(uplo,alpha,A,x,beta,y,True)

# Syr/Her
# -------
lib.ElSyr_s.argtypes = [c_uint,sType,c_void_p,c_void_p]
lib.ElSyr_s.restype = c_uint
lib.ElSyr_d.argtypes = [c_uint,dType,c_void_p,c_void_p]
lib.ElSyr_d.restype = c_uint
lib.ElSyr_c.argtypes = [c_uint,cType,c_void_p,c_void_p]
lib.ElSyr_c.restype = c_uint
lib.ElSyr_z.argtypes = [c_uint,zType,c_void_p,c_void_p]
lib.ElSyr_z.restype = c_uint
lib.ElSyrDist_s.argtypes = [c_uint,sType,c_void_p,c_void_p]
lib.ElSyrDist_s.restype = c_uint
lib.ElSyrDist_d.argtypes = [c_uint,dType,c_void_p,c_void_p]
lib.ElSyrDist_d.restype = c_uint
lib.ElSyrDist_c.argtypes = [c_uint,cType,c_void_p,c_void_p]
lib.ElSyrDist_c.restype = c_uint
lib.ElSyrDist_z.argtypes = [c_uint,zType,c_void_p,c_void_p]
lib.ElSyrDist_z.restype = c_uint

lib.ElHer_c.argtypes = [c_uint,sType,c_void_p,c_void_p]
lib.ElHer_c.restype = c_uint
lib.ElHer_z.argtypes = [c_uint,dType,c_void_p,c_void_p]
lib.ElHer_z.restype = c_uint
lib.ElHerDist_c.argtypes = [c_uint,sType,c_void_p,c_void_p]
lib.ElHerDist_c.restype = c_uint
lib.ElHerDist_z.argtypes = [c_uint,dType,c_void_p,c_void_p]
lib.ElHerDist_z.restype = c_uint

def Syr(uplo,alphaPre,x,A,conj=False):
  if type(A) is not type(x): raise Exception('Types of A and x must match')
  if A.tag != x.tag: raise Exception('Datatypes of A and x must match')
  alpha = TagToType(A.tag)(alphaPre)
  args = [uplo,alpha,x.obj,A.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSyr_s(*args)
    elif A.tag == dTag: lib.ElSyr_d(*args)
    elif A.tag == cTag: 
      if conj: lib.ElHer_c(uplo,alpha.real,x.obj,A.obj)
      else:    lib.ElSyr_c(*args)
    elif A.tag == zTag: 
      if conj: lib.ElHer_z(uplo,alpha.real,x.obj,A.obj)
      else:    lib.ElSyr_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSyrDist_s(*args)
    elif A.tag == dTag: lib.ElSyrDist_d(*args)
    elif A.tag == cTag: 
      if conj: lib.ElHerDist_c(uplo,alpha.real,x.obj,A.obj)
      else:    lib.ElSyrDist_c(*args)
    elif A.tag == zTag: 
      if conj: lib.ElHerDist_z(uplo,alpha.real,x.obj,A.obj)
      else:    lib.ElSyrDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

def Her(uplo,alpha,x,A):
  Syr(uplo,alpha,x,A,True)

# Syr2/Her2
# ---------
lib.ElSyr2_s.argtypes = [c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElSyr2_s.restype = c_uint
lib.ElSyr2_d.argtypes = [c_uint,dType,c_void_p,c_void_p,c_void_p]
lib.ElSyr2_d.restype = c_uint
lib.ElSyr2_c.argtypes = [c_uint,cType,c_void_p,c_void_p,c_void_p]
lib.ElSyr2_c.restype = c_uint
lib.ElSyr2_z.argtypes = [c_uint,zType,c_void_p,c_void_p,c_void_p]
lib.ElSyr2_z.restype = c_uint
lib.ElSyr2Dist_s.argtypes = [c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElSyr2Dist_s.restype = c_uint
lib.ElSyr2Dist_d.argtypes = [c_uint,dType,c_void_p,c_void_p,c_void_p]
lib.ElSyr2Dist_d.restype = c_uint
lib.ElSyr2Dist_c.argtypes = [c_uint,cType,c_void_p,c_void_p,c_void_p]
lib.ElSyr2Dist_c.restype = c_uint
lib.ElSyr2Dist_z.argtypes = [c_uint,zType,c_void_p,c_void_p,c_void_p]
lib.ElSyr2Dist_z.restype = c_uint

lib.ElHer2_c.argtypes = [c_uint,cType,c_void_p,c_void_p,c_void_p]
lib.ElHer2_c.restype = c_uint
lib.ElHer2_z.argtypes = [c_uint,zType,c_void_p,c_void_p,c_void_p]
lib.ElHer2_z.restype = c_uint
lib.ElHer2Dist_c.argtypes = [c_uint,cType,c_void_p,c_void_p,c_void_p]
lib.ElHer2Dist_c.restype = c_uint
lib.ElHer2Dist_z.argtypes = [c_uint,zType,c_void_p,c_void_p,c_void_p]
lib.ElHer2Dist_z.restype = c_uint

def Syr2(uplo,alphaPre,x,y,A,conj=False):
  if type(A) is not type(x) or type(x) is not type(y):
    raise Exception('Types of {A,x,y} must match')
  if A.tag != x.tag or x.tag != y.tag:
    raise Exception('Datatypes of {A,x,y} must match')
  alpha = TagToType(A.tag)(alphaPre)
  args = [uplo,alpha,x.obj,y.obj,A.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSyr2_s(*args)
    elif A.tag == dTag: lib.ElSyr2_d(*args)
    elif A.tag == cTag:
      if conj: lib.ElHer2_c(*args)
      else:    lib.ElSyr2_c(*args)
    elif A.tag == zTag:
      if conj: lib.ElHer2_z(*args)
      else:    lib.ElSyr2_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSyr2Dist_s(*args)
    elif A.tag == dTag: lib.ElSyr2Dist_d(*args)
    elif A.tag == cTag:
      if conj: lib.ElHer2Dist_c(*args)
      else:    lib.ElSyr2Dist_c(*args)
    elif A.tag == zTag:
      if conj: lib.ElHer2Dist_z(*args)
      else:    lib.ElSyr2Dist_z(*args)
    else: DataExcept()
  else: TypeExcept()

def Her2(uplo,alpha,x,y,A):
  Syr2(uplo,alpha,x,y,A,True)

# Trmv
# ----
lib.ElTrmv_s.argtypes = [c_uint,c_uint,c_uint,c_void_p,c_void_p]
lib.ElTrmv_s.restype = c_uint
lib.ElTrmv_d.argtypes = [c_uint,c_uint,c_uint,c_void_p,c_void_p]
lib.ElTrmv_d.restype = c_uint
lib.ElTrmv_c.argtypes = [c_uint,c_uint,c_uint,c_void_p,c_void_p]
lib.ElTrmv_c.restype = c_uint
lib.ElTrmv_z.argtypes = [c_uint,c_uint,c_uint,c_void_p,c_void_p]
lib.ElTrmv_z.restype = c_uint
#lib.ElTrmvDist_s.argtypes = [c_uint,c_uint,c_uint,c_void_p,c_void_p]
#lib.ElTrmvDist_s.restype = c_uint
#lib.ElTrmvDist_d.argtypes = [c_uint,c_uint,c_uint,c_void_p,c_void_p]
#lib.ElTrmvDist_d.restype = c_uint
#lib.ElTrmvDist_c.argtypes = [c_uint,c_uint,c_uint,c_void_p,c_void_p]
#lib.ElTrmvDist_c.restype = c_uint
#lib.ElTrmvDist_z.argtypes = [c_uint,c_uint,c_uint,c_void_p,c_void_p]
#lib.ElTrmvDist_z.restype = c_uint
def Trmv(uplo,orient,diag,A,x):
  args = [uplo,orient,diag,A.obj,x.obj]
  if type(A) is Matrix: 
    if   A.tag == sTag: lib.ElTrmv_s(*args)  
    elif A.tag == dTag: lib.ElTrmv_d(*args)
    elif A.tag == cTag: lib.ElTrmv_c(*args)
    elif A.tag == zTag: lib.ElTrmv_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    # There is not currently a distributed Trmv implementation
    Trmm(LEFT,uplo,orient,diag,1,A,x)
  else: TypeExcept()

# Trr
# ---
lib.ElTrr_i.argtypes = [c_uint,iType,c_void_p,c_void_p,c_void_p]
lib.ElTrr_i.restype = c_uint
lib.ElTrr_s.argtypes = [c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElTrr_s.restype = c_uint
lib.ElTrr_d.argtypes = [c_uint,dType,c_void_p,c_void_p,c_void_p]
lib.ElTrr_d.restype = c_uint
lib.ElTrr_c.argtypes = [c_uint,cType,c_void_p,c_void_p,c_void_p]
lib.ElTrr_c.restype = c_uint
lib.ElTrr_z.argtypes = [c_uint,zType,c_void_p,c_void_p,c_void_p]
lib.ElTrr_z.restype = c_uint
lib.ElTrrDist_i.argtypes = [c_uint,iType,c_void_p,c_void_p,c_void_p]
lib.ElTrrDist_i.restype = c_uint
lib.ElTrrDist_s.argtypes = [c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElTrrDist_s.restype = c_uint
lib.ElTrrDist_d.argtypes = [c_uint,dType,c_void_p,c_void_p,c_void_p]
lib.ElTrrDist_d.restype = c_uint
lib.ElTrrDist_c.argtypes = [c_uint,cType,c_void_p,c_void_p,c_void_p]
lib.ElTrrDist_c.restype = c_uint
lib.ElTrrDist_z.argtypes = [c_uint,zType,c_void_p,c_void_p,c_void_p]
lib.ElTrrDist_z.restype = c_uint
def Trr(uplo,alphaPre,x,y,A):
  if type(x) is not type(y) or type(y) is not type(A):
    raise Exception('Types of {x,y,A} must be the same')
  if x.tag != y.tag or y.tag != A.tag:
    raise Exception('Datatypes of {x,y,A} must be the same')
  alpha = TagToType(A.tag)(alphaPre)
  args = [uplo,alpha,x.obj,y.obj,A.obj]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElTrr_i(*args)
    elif A.tag == sTag: lib.ElTrr_s(*args)
    elif A.tag == dTag: lib.ElTrr_d(*args)
    elif A.tag == cTag: lib.ElTrr_c(*args)
    elif A.tag == zTag: lib.ElTrr_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElTrrDist_i(*args)
    elif A.tag == sTag: lib.ElTrrDist_s(*args)
    elif A.tag == dTag: lib.ElTrrDist_d(*args)
    elif A.tag == cTag: lib.ElTrrDist_c(*args)
    elif A.tag == zTag: lib.ElTrrDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Trr2
# ----
lib.ElTrr2_i.argtypes = [c_uint,iType,c_void_p,c_void_p,c_void_p]
lib.ElTrr2_i.restype = c_uint
lib.ElTrr2_s.argtypes = [c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElTrr2_s.restype = c_uint
lib.ElTrr2_d.argtypes = [c_uint,dType,c_void_p,c_void_p,c_void_p]
lib.ElTrr2_d.restype = c_uint
lib.ElTrr2_c.argtypes = [c_uint,cType,c_void_p,c_void_p,c_void_p]
lib.ElTrr2_c.restype = c_uint
lib.ElTrr2_z.argtypes = [c_uint,zType,c_void_p,c_void_p,c_void_p]
lib.ElTrr2_z.restype = c_uint
lib.ElTrr2Dist_i.argtypes = [c_uint,iType,c_void_p,c_void_p,c_void_p]
lib.ElTrr2Dist_i.restype = c_uint
lib.ElTrr2Dist_s.argtypes = [c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElTrr2Dist_s.restype = c_uint
lib.ElTrr2Dist_d.argtypes = [c_uint,dType,c_void_p,c_void_p,c_void_p]
lib.ElTrr2Dist_d.restype = c_uint
lib.ElTrr2Dist_c.argtypes = [c_uint,cType,c_void_p,c_void_p,c_void_p]
lib.ElTrr2Dist_c.restype = c_uint
lib.ElTrr2Dist_z.argtypes = [c_uint,zType,c_void_p,c_void_p,c_void_p]
lib.ElTrr2Dist_z.restype = c_uint
def Trr2(uplo,alphaPre,X,Y,A):
  if type(X) is not type(Y) or type(Y) is not type(A):
    raise Exception('Types of {X,Y,A} must be the same')
  if X.tag != Y.tag or Y.tag != A.tag:
    raise Exception('Datatypes of {X,Y,A} must be the same')
  alpha = TagToType(A.tag)(alphaPre)
  args = [uplo,alpha,X.obj,Y.obj,A.obj]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElTrr2_i(*args)
    elif A.tag == sTag: lib.ElTrr2_s(*args)
    elif A.tag == dTag: lib.ElTrr2_d(*args)
    elif A.tag == cTag: lib.ElTrr2_c(*args)
    elif A.tag == zTag: lib.ElTrr2_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElTrr2Dist_i(*args)
    elif A.tag == sTag: lib.ElTrr2Dist_s(*args)
    elif A.tag == dTag: lib.ElTrr2Dist_d(*args)
    elif A.tag == cTag: lib.ElTrr2Dist_c(*args)
    elif A.tag == zTag: lib.ElTrr2Dist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Trsv
# ----
lib.ElTrsv_s.argtypes = [c_uint,c_uint,c_uint,c_void_p,c_void_p]
lib.ElTrsv_s.restype = c_uint
lib.ElTrsv_d.argtypes = [c_uint,c_uint,c_uint,c_void_p,c_void_p]
lib.ElTrsv_d.restype = c_uint
lib.ElTrsv_c.argtypes = [c_uint,c_uint,c_uint,c_void_p,c_void_p]
lib.ElTrsv_c.restype = c_uint
lib.ElTrsv_z.argtypes = [c_uint,c_uint,c_uint,c_void_p,c_void_p]
lib.ElTrsv_z.restype = c_uint
lib.ElTrsvDist_s.argtypes = [c_uint,c_uint,c_uint,c_void_p,c_void_p]
lib.ElTrsvDist_s.restype = c_uint
lib.ElTrsvDist_d.argtypes = [c_uint,c_uint,c_uint,c_void_p,c_void_p]
lib.ElTrsvDist_d.restype = c_uint
lib.ElTrsvDist_c.argtypes = [c_uint,c_uint,c_uint,c_void_p,c_void_p]
lib.ElTrsvDist_c.restype = c_uint
lib.ElTrsvDist_z.argtypes = [c_uint,c_uint,c_uint,c_void_p,c_void_p]
lib.ElTrsvDist_z.restype = c_uint
def Trsv(uplo,orient,diag,A,x):
  args = [uplo,orient,diag,A.obj,x.obj]
  if type(A) is Matrix: 
    if   A.tag == sTag: lib.ElTrsv_s(*args)  
    elif A.tag == dTag: lib.ElTrsv_d(*args)
    elif A.tag == cTag: lib.ElTrsv_c(*args)
    elif A.tag == zTag: lib.ElTrsv_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTrsvDist_s(*args)  
    elif A.tag == dTag: lib.ElTrsvDist_d(*args)
    elif A.tag == cTag: lib.ElTrsvDist_c(*args)
    elif A.tag == zTag: lib.ElTrsvDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
