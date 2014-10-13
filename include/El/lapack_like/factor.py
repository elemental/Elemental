#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from ..core import *
import ctypes

# Emulate an enum for LDL pivot types
(BUNCH_KAUFMAN_A,BUNCH_KAUFMAN_C,BUNCH_KAUFMAN_D,BUNCH_KAUFMAN_BOUNDED,
 BUNCH_PARLETT)=(0,1,2,3,4)

class LDLPivot(ctypes.Structure):
  _fields_ = [("nb",iType),("from",(iType*2))]

lib.ElQRCtrlFillDefault_s.argtypes = [c_void_p]
lib.ElQRCtrlFillDefault_s.restype = c_uint
lib.ElQRCtrlFillDefault_d.argtypes = [c_void_p]
lib.ElQRCtrlFillDefault_d.restype = c_uint
class QRCtrl_s(ctypes.Structure):
  _fields_ = [("colPiv",bType),("boundRank",bType),("maxRank",iType),
              ("adaptive",bType),("tol",sType),("alwaysRecomputeNorms",bType)]
  def __init__(self):
    lib.ElQRCtrlFillDefault_s(pointer(self))
class QRCtrl_d(ctypes.Structure):
  _fields_ = [("colPiv",bType),("boundRank",bType),("maxRank",iType),
              ("adaptive",bType),("tol",dType),("alwaysRecomputeNorms",bType)]
  def __init__(self):
    lib.ElQRCtrlFillDefault_d(pointer(self))

# Cholesky
# ========
lib.ElCholesky_s.argtypes = [c_uint,c_void_p]
lib.ElCholesky_s.restype = c_uint
lib.ElCholesky_d.argtypes = [c_uint,c_void_p]
lib.ElCholesky_d.restype = c_uint
lib.ElCholesky_c.argtypes = [c_uint,c_void_p]
lib.ElCholesky_c.restype = c_uint
lib.ElCholesky_z.argtypes = [c_uint,c_void_p]
lib.ElCholesky_z.restype = c_uint
lib.ElCholeskyPiv_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElCholeskyPiv_s.restype = c_uint
lib.ElCholeskyPiv_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElCholeskyPiv_d.restype = c_uint
lib.ElCholeskyPiv_c.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElCholeskyPiv_c.restype = c_uint
lib.ElCholeskyPiv_z.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElCholeskyPiv_z.restype = c_uint
lib.ElCholeskyDist_s.argtypes = [c_uint,c_void_p]
lib.ElCholeskyDist_s.restype = c_uint
lib.ElCholeskyDist_d.argtypes = [c_uint,c_void_p]
lib.ElCholeskyDist_d.restype = c_uint
lib.ElCholeskyDist_c.argtypes = [c_uint,c_void_p]
lib.ElCholeskyDist_c.restype = c_uint
lib.ElCholeskyDist_z.argtypes = [c_uint,c_void_p]
lib.ElCholeskyDist_z.restype = c_uint
lib.ElCholeskyPivDist_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElCholeskyPivDist_s.restype = c_uint
lib.ElCholeskyPivDist_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElCholeskyPivDist_d.restype = c_uint
lib.ElCholeskyPivDist_c.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElCholeskyPivDist_c.restype = c_uint
lib.ElCholeskyPivDist_z.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElCholeskyPivDist_z.restype = c_uint
def Cholesky(uplo,A,p=None):
  if type(A) is Matrix:
    if type(p) is Matrix:
      if p.tag != iTag: raise Exception('p should be integral')
      if   A.tag == sTag: lib.ElCholeskyPiv_s(uplo,A.obj,p.obj)
      elif A.tag == dTag: lib.ElCholeskyPiv_d(uplo,A.obj,p.obj)
      elif A.tag == cTag: lib.ElCholeskyPiv_c(uplo,A.obj,p.obj)
      elif A.tag == zTag: lib.ElCholeskyPiv_z(uplo,A.obj,p.obj)
      else: raise Exception('Unsupported datatype')
    else:
      if   A.tag == sTag: lib.ElCholesky_s(uplo,A.obj)
      elif A.tag == dTag: lib.ElCholesky_d(uplo,A.obj)
      elif A.tag == cTag: lib.ElCholesky_c(uplo,A.obj)
      elif A.tag == zTag: lib.ElCholesky_z(uplo,A.obj)
      else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if type(p) is DistMatrix:
      if p.tag != iTag: raise Exception('p should be integral')
      if   A.tag == sTag: lib.ElCholeskyPivDist_s(uplo,A.obj,p.obj)
      elif A.tag == dTag: lib.ElCholeskyPivDist_d(uplo,A.obj,p.obj)
      elif A.tag == cTag: lib.ElCholeskyPivDist_c(uplo,A.obj,p.obj)
      elif A.tag == zTag: lib.ElCholeskyPivDist_z(uplo,A.obj,p.obj)
      else: raise Exception('Unsupported datatype')
    else:
      if   A.tag == sTag: lib.ElCholeskyDist_s(uplo,A.obj)
      elif A.tag == dTag: lib.ElCholeskyDist_d(uplo,A.obj)
      elif A.tag == cTag: lib.ElCholeskyDist_c(uplo,A.obj)
      elif A.tag == zTag: lib.ElCholeskyDist_z(uplo,A.obj)
      else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElSolveAfterCholesky_s.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSolveAfterCholesky_s.restype = c_uint
lib.ElSolveAfterCholesky_d.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSolveAfterCholesky_d.restype = c_uint
lib.ElSolveAfterCholesky_c.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSolveAfterCholesky_c.restype = c_uint
lib.ElSolveAfterCholesky_z.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSolveAfterCholesky_z.restype = c_uint
lib.ElSolveAfterCholeskyDist_s.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSolveAfterCholeskyDist_s.restype = c_uint
lib.ElSolveAfterCholeskyDist_d.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSolveAfterCholeskyDist_d.restype = c_uint
lib.ElSolveAfterCholeskyDist_c.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSolveAfterCholeskyDist_c.restype = c_uint
lib.ElSolveAfterCholeskyDist_z.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSolveAfterCholeskyDist_z.restype = c_uint
def SolveAfterCholesky(uplo,orient,A,B):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSolveAfterCholesky_s(uplo,orient,A.obj,B.obj)
    elif A.tag == dTag: lib.ElSolveAfterCholesky_d(uplo,orient,A.obj,B.obj)
    elif A.tag == cTag: lib.ElSolveAfterCholesky_c(uplo,orient,A.obj,B.obj)
    elif A.tag == zTag: lib.ElSolveAfterCholesky_z(uplo,orient,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSolveAfterCholeskyDist_s(uplo,orient,A.obj,B.obj)
    elif A.tag == dTag: lib.ElSolveAfterCholeskyDist_d(uplo,orient,A.obj,B.obj)
    elif A.tag == cTag: lib.ElSolveAfterCholeskyDist_c(uplo,orient,A.obj,B.obj)
    elif A.tag == zTag: lib.ElSolveAfterCholeskyDist_z(uplo,orient,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

