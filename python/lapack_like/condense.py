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

# Reduction of a general matrix to bidiagonal form
# ================================================
lib.ElBidiag_s.argtypes = \
lib.ElBidiag_d.argtypes = \
lib.ElBidiag_c.argtypes = \
lib.ElBidiag_z.argtypes = \
lib.ElBidiagDist_s.argtypes = \
lib.ElBidiagDist_d.argtypes = \
lib.ElBidiagDist_c.argtypes = \
lib.ElBidiagDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p]

lib.ElBidiagOnly_s.argtypes = \
lib.ElBidiagOnly_d.argtypes = \
lib.ElBidiagOnly_c.argtypes = \
lib.ElBidiagOnly_z.argtypes = \
lib.ElBidiagOnlyDist_s.argtypes = \
lib.ElBidiagOnlyDist_d.argtypes = \
lib.ElBidiagOnlyDist_c.argtypes = \
lib.ElBidiagOnlyDist_z.argtypes = \
  [c_void_p]

def Bidiag(A,bidiagOnly=False):
  if type(A) is Matrix:
    if bidiagOnly:
      args = [A.obj]
      if   A.tag == sTag: lib.ElBidiagOnly_s(*args)
      elif A.tag == dTag: lib.ElBidiagOnly_d(*args)
      elif A.tag == cTag: lib.ElBidiagOnly_c(*args)
      elif A.tag == zTag: lib.ElBidiagOnly_z(*args)
      else: DataExcept()
    else:
      tP = Matrix(A.tag)
      tQ = Matrix(A.tag)
      args = [A.obj,tP.obj,tQ.obj]
      if   A.tag == sTag: lib.ElBidiag_s(*args)
      elif A.tag == dTag: lib.ElBidiag_d(*args)
      elif A.tag == cTag: lib.ElBidiag_c(*args)
      elif A.tag == zTag: lib.ElBidiag_z(*args)
      else: DataExcept()
      return tP, tQ
  elif type(A) is DistMatrix:
    if bidiagOnly:
      args = [A.obj]
      if   A.tag == sTag: lib.ElBidiagOnlyDist_s(*args)
      elif A.tag == dTag: lib.ElBidiagOnlyDist_d(*args)
      elif A.tag == cTag: lib.ElBidiagOnlyDist_c(*args)
      elif A.tag == zTag: lib.ElBidiagOnlyDist_z(*args)
      else: DataExcept()
    else:
      tP = DistMatrix(A.tag,STAR,STAR,A.Grid())
      tQ = DistMatrix(A.tag,STAR,STAR,A.Grid())
      args = [A.obj,tP.obj,tQ.obj]
      if   A.tag == sTag: lib.ElBidiagDist_s(*args)
      elif A.tag == dTag: lib.ElBidiagDist_d(*args)
      elif A.tag == cTag: lib.ElBidiagDist_c(*args)
      elif A.tag == zTag: lib.ElBidiagDist_z(*args)
      else: DataExcept()
      return tP, tQ
  else: TypeExcept()

# Apply Q from B := Q^H A P
lib.ElApplyQAfterBidiag_s.argtypes = \
lib.ElApplyQAfterBidiag_d.argtypes = \
lib.ElApplyQAfterBidiag_c.argtypes = \
lib.ElApplyQAfterBidiag_z.argtypes = \
lib.ElApplyQAfterBidiagDist_s.argtypes = \
lib.ElApplyQAfterBidiagDist_d.argtypes = \
lib.ElApplyQAfterBidiagDist_c.argtypes = \
lib.ElApplyQAfterBidiagDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]

def ApplyQAfterBidiag(side,orient,A,t,B):
  if type(A) is not type(t) or type(t) is not type(B):
    raise Exception('Matrix types of {A,t,B} must match')
  if A.tag != t.tag or t.tag != B.tag:
    raise Exception('Datatypes of {A,t,B} must match')
  args = [side,orient,A.obj,t.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElApplyQAfterBidiag_s(*args)
    elif A.tag == dTag: lib.ElApplyQAfterBidiag_d(*args)
    elif A.tag == cTag: lib.ElApplyQAfterBidiag_c(*args)
    elif A.tag == zTag: lib.ElApplyQAfterBidiag_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElApplyQAfterBidiagDist_s(*args)
    elif A.tag == dTag: lib.ElApplyQAfterBidiagDist_d(*args)
    elif A.tag == cTag: lib.ElApplyQAfterBidiagDist_c(*args)
    elif A.tag == zTag: lib.ElApplyQAfterBidiagDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Apply P from B := Q^H A P
lib.ElApplyPAfterBidiag_s.argtypes = \
lib.ElApplyPAfterBidiag_d.argtypes = \
lib.ElApplyPAfterBidiag_c.argtypes = \
lib.ElApplyPAfterBidiag_z.argtypes = \
lib.ElApplyPAfterBidiagDist_s.argtypes = \
lib.ElApplyPAfterBidiagDist_d.argtypes = \
lib.ElApplyPAfterBidiagDist_c.argtypes = \
lib.ElApplyPAfterBidiagDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]

def ApplyPAfterBidiag(side,orient,A,t,B):
  if type(A) is not type(t) or type(t) is not type(B):
    raise Exception('Matrix types of {A,t,B} must match')
  if A.tag != t.tag or t.tag != B.tag:
    raise Exception('Datatypes of {A,t,B} must match')
  args = [side,orient,A.obj,t.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElApplyPAfterBidiag_s(*args)
    elif A.tag == dTag: lib.ElApplyPAfterBidiag_d(*args)
    elif A.tag == cTag: lib.ElApplyPAfterBidiag_c(*args)
    elif A.tag == zTag: lib.ElApplyPAfterBidiag_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElApplyPAfterBidiagDist_s(*args)
    elif A.tag == dTag: lib.ElApplyPAfterBidiagDist_d(*args)
    elif A.tag == cTag: lib.ElApplyPAfterBidiagDist_c(*args)
    elif A.tag == zTag: lib.ElApplyPAfterBidiagDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Reduction of a Hermitian matrix to real symmetric tridiagonal form
# ==================================================================
(HERMITIAN_TRIDIAG_NORMAL,HERMITIAN_TRIDIAG_SQUARE,HERMITIAN_TRIDIAG_DEFAULT)= \
(0,1,2)

# TODO: Reenable TridiagCtrl

lib.ElHermitianTridiag_s.argtypes = \
lib.ElHermitianTridiag_d.argtypes = \
lib.ElHermitianTridiag_c.argtypes = \
lib.ElHermitianTridiag_z.argtypes = \
lib.ElHermitianTridiagDist_s.argtypes = \
lib.ElHermitianTridiagDist_d.argtypes = \
lib.ElHermitianTridiagDist_c.argtypes = \
lib.ElHermitianTridiagDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p]
lib.ElHermitianTridiagOnly_s.argtypes = \
lib.ElHermitianTridiagOnly_d.argtypes = \
lib.ElHermitianTridiagOnly_c.argtypes = \
lib.ElHermitianTridiagOnly_z.argtypes = \
lib.ElHermitianTridiagOnlyDist_s.argtypes = \
lib.ElHermitianTridiagOnlyDist_d.argtypes = \
lib.ElHermitianTridiagOnlyDist_c.argtypes = \
lib.ElHermitianTridiagOnlyDist_z.argtypes = \
  [c_uint,c_void_p]

def HermitianTridiag(uplo,A,onlyTridiag=False,ctrl=None):
  if type(A) is Matrix:
    if onlyTridiag: 
      args = [uplo,A.obj]
      if   A.tag == sTag: lib.ElHermitianTridiagOnly_s(*args)
      elif A.tag == dTag: lib.ElHermitianTridiagOnly_d(*args)
      elif A.tag == cTag: lib.ElHermitianTridiagOnly_c(*args)
      elif A.tag == zTag: lib.ElHermitianTridiagOnly_z(*args)
      else: DataExcept()
    else:
      t = Matrix(A.tag)
      args = [uplo,A.obj,t.obj]
      if   A.tag == sTag: lib.ElHermitianTridiag_s(*args)
      elif A.tag == dTag: lib.ElHermitianTridiag_d(*args)
      elif A.tag == cTag: lib.ElHermitianTridiag_c(*args)
      elif A.tag == zTag: lib.ElHermitianTridiag_z(*args)
      else: DataExcept()
      return t
  elif type(A) is DistMatrix:
    if onlyTridiag: 
      if ctrl == None:
        args = [uplo,A.obj]
        if   A.tag == sTag: lib.ElHermitianTridiagOnlyDist_s(*args)
        elif A.tag == dTag: lib.ElHermitianTridiagOnlyDist_d(*args)
        elif A.tag == cTag: lib.ElHermitianTridiagOnlyDist_c(*args)
        elif A.tag == zTag: lib.ElHermitianTridiagOnlyDist_z(*args)
        else: DataExcept()
      else: 
        raise Exception('Support for HermitianTridiagCtrl temp. disabled')
    else:
      t = Matrix(A.tag)
      if ctrl == None:
        args = [uplo,A.obj,t.obj]
        if   A.tag == sTag: lib.ElHermitianTridiagDist_s(*args)
        elif A.tag == dTag: lib.ElHermitianTridiagDist_d(*args)
        elif A.tag == cTag: lib.ElHermitianTridiagDist_c(*args)
        elif A.tag == zTag: lib.ElHermitianTridiagDist_z(*args)
        else: DataExcept()
      else:
        raise Exception('Support for HermitianTridiagCtrl temp. disabled')
      return t
  else: TypeExcept()

lib.ElApplyQAfterHermitianTridiag_s.argtypes = \
lib.ElApplyQAfterHermitianTridiag_d.argtypes = \
lib.ElApplyQAfterHermitianTridiag_c.argtypes = \
lib.ElApplyQAfterHermitianTridiag_z.argtypes = \
lib.ElApplyQAfterHermitianTridiagDist_s.argtypes = \
lib.ElApplyQAfterHermitianTridiagDist_d.argtypes = \
lib.ElApplyQAfterHermitianTridiagDist_c.argtypes = \
lib.ElApplyQAfterHermitianTridiagDist_z.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,c_void_p]

def ApplyQAfterHermitianTridiag(side,uplo,orient,A,t,B):
  if type(A) is not type(t) or type(t) is not type(B):
    raise Exception('Matrix types of {A,t,B} must match')
  if A.tag != t.tag or t.tag != B.tag:
    raise Exception('Datatypes of {A,t,B} must match')
  args = [side,uplo,orient,A.obj,t.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElApplyQAfterHermitianTridiag_s(*args)
    elif A.tag == dTag: lib.ElApplyQAfterHermitianTridiag_d(*args)
    elif A.tag == cTag: lib.ElApplyQAfterHermitianTridiag_c(*args)
    elif A.tag == zTag: lib.ElApplyQAfterHermitianTridiag_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElApplyQAfterHermitianTridiagDist_s(*args)
    elif A.tag == dTag: lib.ElApplyQAfterHermitianTridiagDist_d(*args)
    elif A.tag == cTag: lib.ElApplyQAfterHermitianTridiagDist_c(*args)
    elif A.tag == zTag: lib.ElApplyQAfterHermitianTridiagDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Reduction of a square matrix to Hessenberg form by unitary similarity
# =====================================================================

lib.ElHessenbergOnly_s.argtypes = \
lib.ElHessenbergOnly_d.argtypes = \
lib.ElHessenbergOnly_c.argtypes = \
lib.ElHessenbergOnly_z.argtypes = \
lib.ElHessenbergOnlyDist_s.argtypes = \
lib.ElHessenbergOnlyDist_d.argtypes = \
lib.ElHessenbergOnlyDist_c.argtypes = \
lib.ElHessenbergOnlyDist_z.argtypes = \
  [c_uint,c_void_p]
lib.ElHessenberg_s.argtypes = \
lib.ElHessenberg_d.argtypes = \
lib.ElHessenberg_c.argtypes = \
lib.ElHessenberg_z.argtypes = \
lib.ElHessenbergDist_s.argtypes = \
lib.ElHessenbergDist_d.argtypes = \
lib.ElHessenbergDist_c.argtypes = \
lib.ElHessenbergDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p]

def Hessenberg(uplo,A,hessOnly=False):
  if type(A) is Matrix:
    if hessOnly:
      args = [uplo,A.obj]
      if   A.tag == sTag: lib.ElHessenbergOnly_s(*args)
      elif A.tag == dTag: lib.ElHessenbergOnly_d(*args) 
      elif A.tag == cTag: lib.ElHessenbergOnly_c(*args)
      elif A.tag == zTag: lib.ElHessenbergOnly_z(*args)
      else: DataExcept()
    else:
      t = Matrix(A.tag)
      args = [uplo,A.obj,t.obj]
      if   A.tag == sTag: lib.ElHessenberg_s(*args)
      elif A.tag == dTag: lib.ElHessenberg_d(*args)
      elif A.tag == cTag: lib.ElHessenberg_c(*args)
      elif A.tag == zTag: lib.ElHessenberg_z(*args)
      else: DataExcept()
      return t
  elif type(A) is DistMatrix:
    if hessOnly:
      args = [uplo,A.obj]
      if   A.tag == sTag: lib.ElHessenbergOnlyDist_s(*args)
      elif A.tag == dTag: lib.ElHessenbergOnlyDist_d(*args) 
      elif A.tag == cTag: lib.ElHessenbergOnlyDist_c(*args)
      elif A.tag == zTag: lib.ElHessenbergOnlyDist_z(*args)
      else: DataExcept()
    else:
      t = DistMatrix(A.tag,STAR,STAR,A.Grid())
      args = [uplo,A.obj,t.obj]
      if   A.tag == sTag: lib.ElHessenbergDist_s(*args)
      elif A.tag == dTag: lib.ElHessenbergDist_d(*args)
      elif A.tag == cTag: lib.ElHessenbergDist_c(*args)
      elif A.tag == zTag: lib.ElHessenbergDist_z(*args)
      else: DataExcept()
      return t
  else: TypeExcept()

lib.ElApplyQAfterHessenberg_s.argtypes = \
lib.ElApplyQAfterHessenberg_d.argtypes = \
lib.ElApplyQAfterHessenberg_c.argtypes = \
lib.ElApplyQAfterHessenberg_z.argtypes = \
lib.ElApplyQAfterHessenbergDist_s.argtypes = \
lib.ElApplyQAfterHessenbergDist_d.argtypes = \
lib.ElApplyQAfterHessenbergDist_c.argtypes = \
lib.ElApplyQAfterHessenbergDist_z.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,c_void_p]

def ApplyQAfterHessenberg(side,uplo,orient,A,t,B):
  if type(A) is not type(t) or type(t) is not type(B):
    raise Exception('Matrix types of {A,t,B} must match')
  if A.tag != t.tag or t.tag != B.tag:
    raise Exception('Datatypes of {A,t,B} must match')
  args = [side,uplo,orient,A.obj,t.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElApplyQAfterHessenberg_s(*args)
    elif A.tag == dTag: lib.ElApplyQAfterHessenberg_d(*args)
    elif A.tag == cTag: lib.ElApplyQAfterHessenberg_c(*args)
    elif A.tag == zTag: lib.ElApplyQAfterHessenberg_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElApplyQAfterHessenbergDist_s(*args)
    elif A.tag == dTag: lib.ElApplyQAfterHessenbergDist_d(*args)
    elif A.tag == cTag: lib.ElApplyQAfterHessenbergDist_c(*args)
    elif A.tag == zTag: lib.ElApplyQAfterHessenbergDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
