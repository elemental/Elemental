#
#  Copyright (c) 2009-2015, Jack Poulson
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
lib.ElBidiag_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElBidiag_s.restype = c_uint
lib.ElBidiag_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElBidiag_d.restype = c_uint
lib.ElBidiag_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElBidiag_c.restype = c_uint
lib.ElBidiag_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElBidiag_z.restype = c_uint
lib.ElBidiagDist_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElBidiagDist_s.restype = c_uint
lib.ElBidiagDist_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElBidiagDist_d.restype = c_uint
lib.ElBidiagDist_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElBidiagDist_c.restype = c_uint
lib.ElBidiagDist_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElBidiagDist_z.restype = c_uint
lib.ElBidiagOnly_s.argtypes = [c_void_p]
lib.ElBidiagOnly_s.restype = c_uint
lib.ElBidiagOnly_d.argtypes = [c_void_p]
lib.ElBidiagOnly_d.restype = c_uint
lib.ElBidiagOnly_c.argtypes = [c_void_p]
lib.ElBidiagOnly_c.restype = c_uint
lib.ElBidiagOnly_z.argtypes = [c_void_p]
lib.ElBidiagOnly_z.restype = c_uint
lib.ElBidiagOnlyDist_s.argtypes = [c_void_p]
lib.ElBidiagOnlyDist_s.restype = c_uint
lib.ElBidiagOnlyDist_d.argtypes = [c_void_p]
lib.ElBidiagOnlyDist_d.restype = c_uint
lib.ElBidiagOnlyDist_c.argtypes = [c_void_p]
lib.ElBidiagOnlyDist_c.restype = c_uint
lib.ElBidiagOnlyDist_z.argtypes = [c_void_p]
lib.ElBidiagOnlyDist_z.restype = c_uint
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
lib.ElApplyQAfterBidiag_s.argtypes = [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterBidiag_s.restype = c_uint
lib.ElApplyQAfterBidiag_d.argtypes = [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterBidiag_d.restype = c_uint
lib.ElApplyQAfterBidiag_c.argtypes = [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterBidiag_c.restype = c_uint
lib.ElApplyQAfterBidiag_z.argtypes = [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterBidiag_z.restype = c_uint
lib.ElApplyQAfterBidiagDist_s.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterBidiagDist_s.restype = c_uint
lib.ElApplyQAfterBidiagDist_d.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterBidiagDist_d.restype = c_uint
lib.ElApplyQAfterBidiagDist_c.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterBidiagDist_c.restype = c_uint
lib.ElApplyQAfterBidiagDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterBidiagDist_z.restype = c_uint
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
lib.ElApplyPAfterBidiag_s.argtypes = [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyPAfterBidiag_s.restype = c_uint
lib.ElApplyPAfterBidiag_d.argtypes = [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyPAfterBidiag_d.restype = c_uint
lib.ElApplyPAfterBidiag_c.argtypes = [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyPAfterBidiag_c.restype = c_uint
lib.ElApplyPAfterBidiag_z.argtypes = [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyPAfterBidiag_z.restype = c_uint
lib.ElApplyPAfterBidiagDist_s.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyPAfterBidiagDist_s.restype = c_uint
lib.ElApplyPAfterBidiagDist_d.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyPAfterBidiagDist_d.restype = c_uint
lib.ElApplyPAfterBidiagDist_c.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyPAfterBidiagDist_c.restype = c_uint
lib.ElApplyPAfterBidiagDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyPAfterBidiagDist_z.restype = c_uint
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

lib.ElHermitianTridiag_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHermitianTridiag_s.restype = c_uint
lib.ElHermitianTridiag_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHermitianTridiag_d.restype = c_uint
lib.ElHermitianTridiag_c.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHermitianTridiag_c.restype = c_uint
lib.ElHermitianTridiag_z.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHermitianTridiag_z.restype = c_uint
lib.ElHermitianTridiagDist_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHermitianTridiagDist_s.restype = c_uint
lib.ElHermitianTridiagDist_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHermitianTridiagDist_d.restype = c_uint
lib.ElHermitianTridiagDist_c.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHermitianTridiagDist_c.restype = c_uint
lib.ElHermitianTridiagDist_z.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHermitianTridiagDist_z.restype = c_uint
lib.ElHermitianTridiagOnly_s.argtypes = [c_uint,c_void_p]
lib.ElHermitianTridiagOnly_s.restype = c_uint
lib.ElHermitianTridiagOnly_d.argtypes = [c_uint,c_void_p]
lib.ElHermitianTridiagOnly_d.restype = c_uint
lib.ElHermitianTridiagOnly_c.argtypes = [c_uint,c_void_p]
lib.ElHermitianTridiagOnly_c.restype = c_uint
lib.ElHermitianTridiagOnly_z.argtypes = [c_uint,c_void_p]
lib.ElHermitianTridiagOnly_z.restype = c_uint
lib.ElHermitianTridiagOnlyDist_s.argtypes = [c_uint,c_void_p]
lib.ElHermitianTridiagOnlyDist_s.restype = c_uint
lib.ElHermitianTridiagOnlyDist_d.argtypes = [c_uint,c_void_p]
lib.ElHermitianTridiagOnlyDist_d.restype = c_uint
lib.ElHermitianTridiagOnlyDist_c.argtypes = [c_uint,c_void_p]
lib.ElHermitianTridiagOnlyDist_c.restype = c_uint
lib.ElHermitianTridiagOnlyDist_z.argtypes = [c_uint,c_void_p]
lib.ElHermitianTridiagOnlyDist_z.restype = c_uint
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
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterHermitianTridiag_s.restype = c_uint
lib.ElApplyQAfterHermitianTridiag_d.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterHermitianTridiag_d.restype = c_uint
lib.ElApplyQAfterHermitianTridiag_c.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterHermitianTridiag_c.restype = c_uint
lib.ElApplyQAfterHermitianTridiag_z.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterHermitianTridiag_z.restype = c_uint
lib.ElApplyQAfterHermitianTridiagDist_s.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterHermitianTridiagDist_s.restype = c_uint
lib.ElApplyQAfterHermitianTridiagDist_d.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterHermitianTridiagDist_d.restype = c_uint
lib.ElApplyQAfterHermitianTridiagDist_c.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterHermitianTridiagDist_c.restype = c_uint
lib.ElApplyQAfterHermitianTridiagDist_z.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterHermitianTridiagDist_z.restype = c_uint
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
lib.ElHessenbergOnly_s.argtypes = [c_uint,c_void_p]
lib.ElHessenbergOnly_s.restype = c_uint
lib.ElHessenbergOnly_d.argtypes = [c_uint,c_void_p]
lib.ElHessenbergOnly_d.restype = c_uint
lib.ElHessenbergOnly_c.argtypes = [c_uint,c_void_p]
lib.ElHessenbergOnly_c.restype = c_uint
lib.ElHessenbergOnly_z.argtypes = [c_uint,c_void_p]
lib.ElHessenbergOnly_z.restype = c_uint
lib.ElHessenbergOnlyDist_s.argtypes = [c_uint,c_void_p]
lib.ElHessenbergOnlyDist_s.restype = c_uint
lib.ElHessenbergOnlyDist_d.argtypes = [c_uint,c_void_p]
lib.ElHessenbergOnlyDist_d.restype = c_uint
lib.ElHessenbergOnlyDist_c.argtypes = [c_uint,c_void_p]
lib.ElHessenbergOnlyDist_c.restype = c_uint
lib.ElHessenbergOnlyDist_z.argtypes = [c_uint,c_void_p]
lib.ElHessenbergOnlyDist_z.restype = c_uint
lib.ElHessenberg_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHessenberg_s.restype = c_uint
lib.ElHessenberg_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHessenberg_d.restype = c_uint
lib.ElHessenberg_c.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHessenberg_c.restype = c_uint
lib.ElHessenberg_z.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHessenberg_z.restype = c_uint
lib.ElHessenbergDist_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHessenbergDist_s.restype = c_uint
lib.ElHessenbergDist_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHessenbergDist_d.restype = c_uint
lib.ElHessenbergDist_c.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHessenbergDist_c.restype = c_uint
lib.ElHessenbergDist_z.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHessenbergDist_z.restype = c_uint
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
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterHessenberg_s.restype = c_uint
lib.ElApplyQAfterHessenberg_d.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterHessenberg_d.restype = c_uint
lib.ElApplyQAfterHessenberg_c.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterHessenberg_c.restype = c_uint
lib.ElApplyQAfterHessenberg_z.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterHessenberg_z.restype = c_uint
lib.ElApplyQAfterHessenbergDist_s.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterHessenbergDist_s.restype = c_uint
lib.ElApplyQAfterHessenbergDist_d.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterHessenbergDist_d.restype = c_uint
lib.ElApplyQAfterHessenbergDist_c.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterHessenbergDist_c.restype = c_uint
lib.ElApplyQAfterHessenbergDist_z.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterHessenbergDist_z.restype = c_uint
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
