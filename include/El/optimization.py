#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *

from ctypes import CFUNCTYPE

# Basis pursuit
# =============
lib.ElBasisPursuit_s.argtypes = \
lib.ElBasisPursuit_d.argtypes = \
lib.ElBasisPursuit_c.argtypes = \
lib.ElBasisPursuit_z.argtypes = \
lib.ElBasisPursuitDist_s.argtypes = \
lib.ElBasisPursuitDist_d.argtypes = \
lib.ElBasisPursuitDist_c.argtypes = \
lib.ElBasisPursuitDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,POINTER(iType)]
lib.ElBasisPursuit_s.restype = \
lib.ElBasisPursuit_d.restype = \
lib.ElBasisPursuit_c.restype = \
lib.ElBasisPursuit_z.restype = \
lib.ElBasisPursuitDist_s.restype = \
lib.ElBasisPursuitDist_d.restype = \
lib.ElBasisPursuitDist_c.restype = \
lib.ElBasisPursuitDist_z.restype = \
  c_uint
def BasisPursuit(A,b):
  if type(A) is not type(b):
    raise Exception('Types of A and b must match')
  if A.tag != b.tag:
    raise Exception('Datatypes of A and b must match')
  numIts = iType()
  if type(A) is Matrix:
    z = Matrix(A.tag)
    args = [A.obj,b.obj,z.obj,pointer(numIts)]
    if   A.tag == sTag: lib.ElBasisPursuit_s(*args)
    elif A.tag == dTag: lib.ElBasisPursuit_d(*args)
    elif A.tag == cTag: lib.ElBasisPursuit_c(*args)
    elif A.tag == zTag: lib.ElBasisPursuit_z(*args)
    else: DataExcept()
    return z, numIts
  elif type(A) is DistMatrix:
    z = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,b.obj,z.obj,pointer(numIts)]
    if   A.tag == sTag: lib.ElBasisPursuitDist_s(*args)
    elif A.tag == dTag: lib.ElBasisPursuitDist_d(*args)
    elif A.tag == cTag: lib.ElBasisPursuitDist_c(*args)
    elif A.tag == zTag: lib.ElBasisPursuitDist_z(*args)
    else: DataExcept()
    return z, numIts
  else: TypeExcept()

# Least Absolute Shrinkage and Selection Operator
# ===============================================
lib.ElLasso_s.argtypes = \
lib.ElLasso_c.argtypes = \
lib.ElLassoDist_s.argtypes = \
lib.ElLassoDist_c.argtypes = \
  [c_void_p,c_void_p,sType,c_void_p,POINTER(iType)]
lib.ElLasso_d.argtypes = \
lib.ElLasso_z.argtypes = \
lib.ElLassoDist_d.argtypes = \
lib.ElLassoDist_z.argtypes = \
  [c_void_p,c_void_p,dType,c_void_p,POINTER(iType)]

lib.ElLasso_s.restype = \
lib.ElLasso_d.restype = \
lib.ElLasso_c.restype = \
lib.ElLasso_z.restype = \
lib.ElLassoDist_s.restype = \
lib.ElLassoDist_d.restype = \
lib.ElLassoDist_c.restype = \
lib.ElLassoDist_z.restype = \
  c_uint

def Lasso(A,b,lamb):
  if type(A) is not type(b): raise Exception('Types of A and b must match')
  if A.tag != b.tag: raise Exception('Datatypes of A and b must match')
  numIts = iType()
  if type(A) is Matrix:
    z = Matrix(A.tag)
    args = [A.obj,b.obj,lamb,z.obj,pointer(numIts)]
    if   A.tag == sTag: lib.ElLasso_s(*args)
    elif A.tag == dTag: lib.ElLasso_d(*args)
    elif A.tag == cTag: lib.ElLasso_c(*args)
    elif A.tag == zTag: lib.ElLasso_z(*args)
    else: DataExcept()
    return z, numIts
  elif type(A) is DistMatrix:
    z = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,b.obj,lamb,z.obj,pointer(numIts)]
    if   A.tag == sTag: lib.ElLassoDist_s(*args)
    elif A.tag == dTag: lib.ElLassoDist_d(*args)
    elif A.tag == cTag: lib.ElLassoDist_c(*args)
    elif A.tag == zTag: lib.ElLassoDist_z(*args)
    else: DataExcept()
    return z, numIts
  else: TypeExcept()

# Linear program
# ==============

(LP_ADMM,LP_IPF,LP_IPF_SELFDUAL,LP_MEHROTRA,LP_MEHROTRA_SELFDUAL)=(0,1,2,3,4)

lib.ElLPIPFLineSearchCtrlDefault_s.argtypes = \
lib.ElLPIPFLineSearchCtrlDefault_d.argtypes = [c_void_p]
lib.ElLPIPFLineSearchCtrlDefault_s.restype = \
lib.ElLPIPFLineSearchCtrlDefault_d.restype = c_uint
class LPIPFLineSearchCtrl_s(ctypes.Structure):
  _fields_ = [("gamma",sType),("beta",sType),("psi",sType),
              ("stepRatio",sType),("progress",bType)]
  def __init__(self):
    lib.ElLPIPFLineSearchCtrlDefault_s(pointer(self))
class LPIPFLineSearchCtrl_d(ctypes.Structure):
  _fields_ = [("gamma",dType),("beta",dType),("psi",dType),
              ("stepRatio",dType),("progress",bType)]
  def __init__(self):
    lib.ElLPIPFLineSearchCtrlDefault_d(pointer(self))

# Direct conic form
# -----------------
lib.ElLPDirectADMMCtrlDefault_s.argtypes = \
lib.ElLPDirectADMMCtrlDefault_d.argtypes = [c_void_p]
lib.ElLPDirectADMMCtrlDefault_s.restype = \
lib.ElLPDirectADMMCtrlDefault_d.restype = c_uint
class LPDirectADMMCtrl_s(ctypes.Structure):
  _fields_ = [("rho",sType),("alpha",sType),
              ("maxIter",iType),
              ("absTol",sType),("relTol",sType),
              ("inv",bType),("progress",bType)]
  def __init__(self):
    lib.ElLPDirectADMMCtrlDefault_s(pointer(self))
class LPDirectADMMCtrl_d(ctypes.Structure):
  _fields_ = [("rho",dType),("alpha",dType),
              ("maxIter",iType),
              ("absTol",dType),("relTol",dType),
              ("inv",bType),("progress",bType)]
  def __init__(self):
    lib.ElLPDirectADMMCtrlDefault_d(pointer(self))

(LP_DIRECT_FULL_KKT,LP_DIRECT_AUGMENTED_KKT,LP_DIRECT_NORMAL_KKT) = (0,1,2)

lib.ElLPDirectIPFCtrlDefault_s.argtypes = \
lib.ElLPDirectIPFCtrlDefault_d.argtypes = [c_void_p,bType]
lib.ElLPDirectIPFCtrlDefault_s.restype = \
lib.ElLPDirectIPFCtrlDefault_d.restype = c_uint
class LPDirectIPFCtrl_s(ctypes.Structure):
  _fields_ = [("primalInitialized",bType),("dualInitialized",bType),
              ("tol",sType),("maxIts",iType),("centering",sType),
              ("system",c_uint),("lineSearchCtrl",LPIPFLineSearchCtrl_s),
              ("progress",bType)]
  def __init__(self,isSparse=True):
    lib.ElLPDirectIPFCtrlDefault_s(pointer(self),isSparse)
class LPDirectIPFCtrl_d(ctypes.Structure):
  _fields_ = [("primalInitialized",bType),("dualInitialized",bType),
              ("tol",dType),("maxIts",iType),("centering",dType),
              ("system",c_uint),("lineSearchCtrl",LPIPFLineSearchCtrl_d),
              ("progress",bType)]
  def __init__(self,isSparse=True):
    lib.ElLPDirectIPFCtrlDefault_d(pointer(self),isSparse)

lib.ElLPDirectMehrotraCtrlDefault_s.argtypes = \
lib.ElLPDirectMehrotraCtrlDefault_d.argtypes = [c_void_p,bType]
lib.ElLPDirectMehrotraCtrlDefault_s.restype = \
lib.ElLPDirectMehrotraCtrlDefault_d.restype = c_uint
class LPDirectMehrotraCtrl_s(ctypes.Structure):
  _fields_ = [("primalInitialized",bType),("dualInitialized",bType),
              ("tol",sType),("maxIts",iType),("maxStepRatio",sType),
              ("system",c_uint),("progress",bType)]
  def __init__(self,isSparse=True):
    lib.ElLPDirectMehrotraCtrlDefault_s(pointer(self),isSparse)
class LPDirectMehrotraCtrl_d(ctypes.Structure):
  _fields_ = [("primalInitialized",bType),("dualInitialized",bType),
              ("tol",dType),("maxIts",iType),("maxStepRatio",dType),
              ("system",c_uint),("progress",bType)]
  def __init__(self,isSparse=True):
    lib.ElLPDirectMehrotraCtrlDefault_d(pointer(self),isSparse)

lib.ElLPDirectCtrlDefault_s.argtypes = \
lib.ElLPDirectCtrlDefault_d.argtypes = [c_void_p,bType]
lib.ElLPDirectCtrlDefault_s.restype = \
lib.ElLPDirectCtrlDefault_d.restype = c_uint
class LPDirectCtrl_s(ctypes.Structure):
  _fields_ = [("approach",c_uint),("admmCtrl",LPDirectADMMCtrl_s),
              ("ipfCtrl",LPDirectIPFCtrl_s),
              ("mehrotraCtrl",LPDirectMehrotraCtrl_s)]
  def __init__(self,isSparse=True):
    lib.ElLPDirectCtrlDefault_s(pointer(self),isSparse)
class LPDirectCtrl_d(ctypes.Structure):
  _fields_ = [("approach",c_uint),("admmCtrl",LPDirectADMMCtrl_d),
              ("ipfCtrl",LPDirectIPFCtrl_d),
              ("mehrotraCtrl",LPDirectMehrotraCtrl_d)]
  def __init__(self,isSparse=True):
    lib.ElLPDirectCtrlDefault_d(pointer(self),isSparse)

lib.ElLPDirect_s.argtypes = \
lib.ElLPDirect_d.argtypes = \
lib.ElLPDirectDist_s.argtypes = \
lib.ElLPDirectDist_d.argtypes = \
lib.ElLPDirectSparse_s.argtypes = \
lib.ElLPDirectSparse_d.argtypes = \
lib.ElLPDirectDistSparse_s.argtypes = \
lib.ElLPDirectDistSparse_d.argtypes = \
  [c_void_p,
   c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p]
lib.ElLPDirect_s.restype = \
lib.ElLPDirect_d.restype = \
lib.ElLPDirectDist_s.restype = \
lib.ElLPDirectDist_d.restype = \
lib.ElLPDirectSparse_s.restype = \
lib.ElLPDirectSparse_d.restype = \
lib.ElLPDirectDistSparse_s.restype = \
lib.ElLPDirectDistSparse_d.restype = \
  c_uint

lib.ElLPDirectX_s.argtypes = \
lib.ElLPDirectXSparse_s.argtypes = \
lib.ElLPDirectXDist_s.argtypes = \
lib.ElLPDirectXDistSparse_s.argtypes = \
  [c_void_p,
   c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   LPDirectCtrl_s]
lib.ElLPDirectX_d.argtypes = \
lib.ElLPDirectXSparse_d.argtypes = \
lib.ElLPDirectXDist_d.argtypes = \
lib.ElLPDirectXDistSparse_d.argtypes = \
  [c_void_p,
   c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   LPDirectCtrl_d]
lib.ElLPDirectX_s.restype = \
lib.ElLPDirectX_d.restype = \
lib.ElLPDirectXSparse_s.restype = \
lib.ElLPDirectXSparse_d.restype = \
lib.ElLPDirectXDist_s.restype = \
lib.ElLPDirectXDist_d.restype = \
lib.ElLPDirectXDistSparse_s.restype = \
lib.ElLPDirectXDistSparse_d.restype = \
  c_uint

def LPDirect(A,b,c,x,y,z,ctrl=None):
  if type(b) is not type(c) or type(b) is not type(x) or \
     type(b) is not type(y) or type(b) is not type(z):
    raise Exception('{b,c,x,y,z} must be of the same type')
  args = [A.obj,b.obj,c.obj,x.obj,y.obj,z.obj]
  argsCtrl = [A.obj,b.obj,c.obj,x.obj,y.obj,z.obj,ctrl]
  if type(A) is Matrix:
    if type(b) is not Matrix: raise Exception('b must be a Matrix')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElLPDirect_s(*args)
      else:            lib.ElLPDirectX_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElLPDirect_d(*args)
      else:            lib.ElLPDirectX_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if type(b) is not DistMatrix: raise Exception('b must be a DistMatrix')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElLPDirectDist_s(*args)
      else:            lib.ElLPDirectXDist_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElLPDirectDist_d(*args)
      else:            lib.ElLPDirectXDist_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if type(b) is not Matrix: raise Exception('b must be a Matrix')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElLPDirectSparse_s(*args)
      else:            lib.ElLPDirectXSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElLPDirectSparse_d(*args)
      else:            lib.ElLPDirectXSparse_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if type(b) is not DistMultiVec: raise Exception('b must be a DistMultiVec')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElLPDirectDistSparse_s(*args)
      else:            lib.ElLPDirectXDistSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElLPDirectDistSparse_d(*args)
      else:            lib.ElLPDirectXDistSparse_d(*argsCtrl)
    else: DataExcept()
  else: TypeExcept()

# Affine conic form
# -----------------
lib.ElLPAffineIPFCtrlDefault_s.argtypes = \
lib.ElLPAffineIPFCtrlDefault_d.argtypes = [c_void_p]
lib.ElLPAffineIPFCtrlDefault_s.restype = \
lib.ElLPAffineIPFCtrlDefault_d.restype = c_uint
class LPAffineIPFCtrl_s(ctypes.Structure):
  _fields_ = [("primalInitialized",bType),("dualInitialized",bType),
              ("tol",sType),("maxIts",iType),("centering",sType),
              ("lineSearchCtrl",LPIPFLineSearchCtrl_s),
              ("progress",bType)]
  def __init__(self):
    lib.ElLPAffineIPFCtrlDefault_s(pointer(self))
class LPAffineIPFCtrl_d(ctypes.Structure):
  _fields_ = [("primalInitialized",bType),("dualInitialized",bType),
              ("tol",dType),("maxIts",iType),("centering",dType),
              ("lineSearchCtrl",LPIPFLineSearchCtrl_d),
              ("progress",bType)]
  def __init__(self):
    lib.ElLPAffineIPFCtrlDefault_d(pointer(self))

lib.ElLPAffineMehrotraCtrlDefault_s.argtypes = \
lib.ElLPAffineMehrotraCtrlDefault_d.argtypes = [c_void_p]
lib.ElLPAffineMehrotraCtrlDefault_s.restype = \
lib.ElLPAffineMehrotraCtrlDefault_d.restype = c_uint
class LPAffineMehrotraCtrl_s(ctypes.Structure):
  _fields_ = [("primalInitialized",bType),("dualInitialized",bType),
              ("tol",sType),("maxIts",iType),("maxStepRatio",sType),
              ("progress",bType)]
  def __init__(self):
    lib.ElLPAffineMehrotraCtrlDefault_s(pointer(self))
class LPAffineMehrotraCtrl_d(ctypes.Structure):
  _fields_ = [("primalInitialized",bType),("dualInitialized",bType),
              ("tol",dType),("maxIts",iType),("maxStepRatio",dType),
              ("progress",bType)]
  def __init__(self):
    lib.ElLPAffineMehrotraCtrlDefault_d(pointer(self))

lib.ElLPAffineCtrlDefault_s.argtypes = \
lib.ElLPAffineCtrlDefault_d.argtypes = [c_void_p]
lib.ElLPAffineCtrlDefault_s.restype = \
lib.ElLPAffineCtrlDefault_d.restype = c_uint
class LPAffineCtrl_s(ctypes.Structure):
  _fields_ = [("approach",c_uint),
              ("ipfCtrl",LPAffineIPFCtrl_s),
              ("mehrotraCtrl",LPAffineMehrotraCtrl_s)]
  def __init__(self):
    lib.ElLPAffineCtrlDefault_s(pointer(self))
class LPAffineCtrl_d(ctypes.Structure):
  _fields_ = [("approach",c_uint),
              ("ipfCtrl",LPAffineIPFCtrl_d),
              ("mehrotraCtrl",LPAffineMehrotraCtrl_d)]
  def __init__(self):
    lib.ElLPAffineCtrlDefault_d(pointer(self))

lib.ElLPAffine_s.argtypes = \
lib.ElLPAffine_d.argtypes = \
lib.ElLPAffineDist_s.argtypes = \
lib.ElLPAffineDist_d.argtypes = \
lib.ElLPAffineSparse_s.argtypes = \
lib.ElLPAffineSparse_d.argtypes = \
lib.ElLPAffineDistSparse_s.argtypes = \
lib.ElLPAffineDistSparse_d.argtypes = \
  [c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLPAffine_s.restype = \
lib.ElLPAffine_d.restype = \
lib.ElLPAffineDist_s.restype = \
lib.ElLPAffineDist_d.restype = \
lib.ElLPAffineSparse_s.restype = \
lib.ElLPAffineSparse_d.restype = \
lib.ElLPAffineDistSparse_s.restype = \
lib.ElLPAffineDistSparse_d.restype = \
  c_uint

lib.ElLPAffineX_s.argtypes = \
lib.ElLPAffineXDist_s.argtypes = \
lib.ElLPAffineXSparse_s.argtypes = \
lib.ElLPAffineXDistSparse_s.argtypes = \
  [c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,c_void_p,
   LPAffineCtrl_s]
lib.ElLPAffineX_d.argtypes = \
lib.ElLPAffineXDist_d.argtypes = \
lib.ElLPAffineXSparse_d.argtypes = \
lib.ElLPAffineXDistSparse_d.argtypes = \
  [c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,c_void_p,
   LPAffineCtrl_d]
lib.ElLPAffineX_s.restype = \
lib.ElLPAffineX_d.restype = \
lib.ElLPAffineXDist_s.restype = \
lib.ElLPAffineXDist_d.restype = \
lib.ElLPAffineXSparse_s.restype = \
lib.ElLPAffineXSparse_d.restype = \
lib.ElLPAffineXDistSparse_s.restype = \
lib.ElLPAffineXDistSparse_d.restype = \
  c_uint

def LPAffine(A,G,b,c,h,x,y,z,s,ctrl=None):
  if type(A) is not type(G):
    raise Exception('A and G must be of the same type')
  if type(b) is not type(c) or type(b) is not type(c) or \
     type(b) is not type(h) or type(b) is not type(x) or \
     type(b) is not type(y) or type(b) is not type(z) or \
     type(b) is not type(s):
    raise Exception('{b,c,h,x,y,z,s} must be of the same type')
  args = [A.obj,G.obj,b.obj,c.obj,h.obj,x.obj,y.obj,z.obj,s.obj]
  argsCtrl = [A.obj,G.obj,b.obj,c.obj,h.obj,x.obj,y.obj,z.obj,s.obj,ctrl]
  if type(A) is Matrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElLPAffine_s(*args)
      else:            lib.ElLPAffineX_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElLPAffine_d(*args)
      else:            lib.ElLPAffineX_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if type(b) is not DistMatrix:
      raise Exception('b must be a DistMatrix')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElLPAffineDist_s(*args)
      else:            lib.ElLPAffineXDist_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElLPAffineDist_d(*args)
      else:            lib.ElLPAffineXDist_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElLPAffineSparse_s(*args)
      else:            lib.ElLPAffineXSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElLPAffineSparse_d(*args)
      else:            lib.ElLPAffineXSparse_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if type(b) is not DistMultiVec:
      raise Exception('b must be a DistMultiVec')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElLPAffineDistSparse_s(*args)
      else:            lib.ElLPAffineXDistSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElLPAffineDistSparse_d(*args)
      else:            lib.ElLPAffineXDistSparse_d(*argsCtrl)
    else: DataExcept()
  else: TypeExcept()

# Logistic regression
# ===================
(NO_PENALTY,L1_PENALTY,L2_PENALTY)=(0,1,2)

lib.ElLogisticRegression_s.argtypes = \
lib.ElLogisticRegressionDist_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,sType,c_uint,POINTER(iType)]
lib.ElLogisticRegression_d.argtypes = \
lib.ElLogisticRegressionDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,dType,c_uint,POINTER(iType)]

lib.ElLogisticRegression_s.restype = \
lib.ElLogisticRegression_d.restype = \
lib.ElLogisticRegressionDist_s.restype = \
lib.ElLogisticRegressionDist_d.restype = \
  c_uint

def LogisticRegression(G,q,gamma,penalty=L1_PENALTY):
  if type(G) is not type(q):
    raise Exception('Types of G and q must match')
  if G.tag != q.tag:
    raise Exception('Datatypes of G and q must match')
  numIts = iType()
  if type(G) is Matrix:
    z = Matrix(G.tag)
    args = [G.obj,q.obj,z.obj,gamma,penalty,pointer(numIts)]
    if   G.tag == sTag: lib.ElLogisticRegression_s(*args)
    elif G.tag == dTag: lib.ElLogisticRegression_d(*args)
    else: DataExcept()
    return z, numIts
  elif type(G) is DistMatrix:
    z = DistMatrix(G.tag,MC,MR,G.Grid())
    args = [G.obj,q.obj,z.obj,gamma,penalty,pointer(numIts)]
    if   G.tag == sTag: lib.ElLogisticRegressionDist_s(*args)
    elif G.tag == dTag: lib.ElLogisticRegressionDist_d(*args)
    else: DataExcept()
    return z, numIts
  else: TypeExcept()

# Model fit
# =========
lib.ElModelFit_s.argtypes = \
lib.ElModelFitDist_s.argtypes = \
  [CFUNCTYPE(None,c_void_p,sType),CFUNCTYPE(None,c_void_p,sType),
   c_void_p,c_void_p,c_void_p,sType,POINTER(iType)]
lib.ElModelFit_d.argtypes = \
lib.ElModelFitDist_d.argtypes = \
  [CFUNCTYPE(None,c_void_p,dType),CFUNCTYPE(None,c_void_p,dType),
   c_void_p,c_void_p,c_void_p,dType,POINTER(iType)]
lib.ElModelFit_s.restype = \
lib.ElModelFit_d.restype = \
lib.ElModelFitDist_s.restype = \
lib.ElModelFitDist_d.restype = \
  c_uint

def ModelFit(lossProx,regProx,A,b,rho=1.2):
  if type(A) is not type(b):
    raise Exception('Types of A and b must match')
  if A.tag != b.tag:
    raise Exception('Datatypes of A and b must match')
  numIts = iType()
  cLoss = CFUNCTYPE(None,c_void_p,TagToType(A.tag))(lossProx)
  cReg = CFUNCTYPE(None,c_void_p,TagToType(A.tag))(regProx)
  if type(A) is Matrix:
    w = Matrix(A.tag)
    args = [cLoss,cReg,A.obj,b.obj,w.obj,rho,pointer(numIts)]
    if   A.tag == sTag: lib.ElModelFit_s(*args)
    elif A.tag == dTag: lib.ElModelFit_d(*args)
    else: DataExcept()
    return w, numIts
  elif type(A) is DistMatrix:
    w = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [cLoss,cReg,A.obj,b.obj,w.obj,rho,pointer(numIts)]
    if   A.tag == sTag: lib.ElModelFitDist_s(*args)
    elif A.tag == dTag: lib.ElModelFitDist_d(*args)
    else: DataExcept()
    return w, numIts
  else: TypeExcept()

# Non-negative matrix factorization
# =================================
lib.ElNMF_s.argtypes = \
lib.ElNMF_d.argtypes = \
lib.ElNMFDist_s.argtypes = \
lib.ElNMFDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p]
lib.ElNMF_s.restype = \
lib.ElNMF_d.restype = \
lib.ElNMFDist_s.restype = \
lib.ElNMFDist_d.restype = \
  c_uint

def NMF(A):
  args = [A.obj,X.obj,Y.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElNMF_s(*args)
    elif A.tag == dTag: lib.ElNMF_d(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElNMFDist_s(*args)
    elif A.tag == dTag: lib.ElNMFDist_d(*args)
    else: DataExcept()
  else: TypeExcept()

# Non-negative least squares
# ==========================
lib.ElNonNegativeLeastSquares_s.argtypes = \
lib.ElNonNegativeLeastSquares_d.argtypes = \
lib.ElNonNegativeLeastSquaresDist_s.argtypes = \
lib.ElNonNegativeLeastSquaresDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,POINTER(iType)]
lib.ElNonNegativeLeastSquares_s.restype = \
lib.ElNonNegativeLeastSquares_d.restype = \
lib.ElNonNegativeLeastSquaresDist_s.restype = \
lib.ElNonNegativeLeastSquaresDist_d.restype = \
  c_uint

def NonNegativeLeastSquares(A,Y):
  if type(A) is not type(Y):
    raise Exception('Types of A and Y must match')
  if A.tag != Y.tag:
    raise Exception('Datatypes of A and Y must match')
  numIts = iType()
  if type(A) is Matrix:
    Z = Matrix(A.tag)
    args = [A.obj,Y.obj,Z.obj,pointer(numIts)]
    if   A.tag == sTag: lib.ElNonNegativeLeastSquares_s(*args)
    elif A.tag == dTag: lib.ElNonNegativeLeastSquares_d(*args)
    else: DataExcept()
    return Z, numIts
  elif type(A) is DistMatrix:
    Z = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,Y.obj,Z.obj,pointer(numIts)]
    if   A.tag == sTag: lib.ElNonNegativeLeastSquaresDist_s(*args)
    elif A.tag == dTag: lib.ElNonNegativeLeastSquaresDist_d(*args)
    else: DataExcept()
    return Z, numIts
  else: TypeExcept()

# Quadratic program
# =================
(QP_ADMM,QP_IPF,QP_IPF_SELFDUAL,QP_MEHROTRA,QP_MEHROTRA_SELFDUAL)=(0,1,2,3,4)

lib.ElQPIPFLineSearchCtrlDefault_s.argtypes = \
lib.ElQPIPFLineSearchCtrlDefault_d.argtypes = [c_void_p]
lib.ElQPIPFLineSearchCtrlDefault_s.restype = \
lib.ElQPIPFLineSearchCtrlDefault_d.restype = c_uint
class QPIPFLineSearchCtrl_s(ctypes.Structure):
  _fields_ = [("gamma",sType),("beta",sType),("psi",sType),
              ("stepRatio",sType),("progress",bType)]
  def __init__(self):
    lib.ElQPIPFLineSearchCtrlDefault_s(pointer(self))
class QPIPFLineSearchCtrl_d(ctypes.Structure):
  _fields_ = [("gamma",dType),("beta",dType),("psi",dType),
              ("stepRatio",dType),("progress",bType)]
  def __init__(self):
    lib.ElQPIPFLineSearchCtrlDefault_d(pointer(self))

# Direct conic form
# -----------------
(QP_DIRECT_FULL_KKT,QP_DIRECT_AUGMENTED_KKT) = (0,1)

lib.ElQPDirectIPFCtrlDefault_s.argtypes = \
lib.ElQPDirectIPFCtrlDefault_d.argtypes = [c_void_p]
lib.ElQPDirectIPFCtrlDefault_s.restype = \
lib.ElQPDirectIPFCtrlDefault_d.restype = c_uint
class QPDirectIPFCtrl_s(ctypes.Structure):
  _fields_ = [("primalInitialized",bType),("dualInitialized",bType),
              ("tol",sType),("maxIts",iType),("centering",sType),
              ("system",c_uint),("lineSearchCtrl",QPIPFLineSearchCtrl_s),
              ("progress",bType)]
  def __init__(self):
    lib.ElQPDirectIPFCtrlDefault_s(pointer(self))
class QPDirectIPFCtrl_d(ctypes.Structure):
  _fields_ = [("primalInitialized",bType),("dualInitialized",bType),
              ("tol",dType),("maxIts",iType),("centering",dType),
              ("system",c_uint),("lineSearchCtrl",QPIPFLineSearchCtrl_d),
              ("progress",bType)]
  def __init__(self):
    lib.ElQPDirectIPFCtrlDefault_d(pointer(self))

lib.ElQPDirectMehrotraCtrlDefault_s.argtypes = \
lib.ElQPDirectMehrotraCtrlDefault_d.argtypes = [c_void_p]
lib.ElQPDirectMehrotraCtrlDefault_s.restype = \
lib.ElQPDirectMehrotraCtrlDefault_d.restype = c_uint
class QPDirectMehrotraCtrl_s(ctypes.Structure):
  _fields_ = [("primalInitialized",bType),("dualInitialized",bType),
              ("tol",sType),("maxIts",iType),("maxStepRatio",sType),
              ("system",c_uint),("progress",bType)]
  def __init__(self):
    lib.ElQPDirectMehrotraCtrlDefault_s(pointer(self))
class QPDirectMehrotraCtrl_d(ctypes.Structure):
  _fields_ = [("primalInitialized",bType),("dualInitialized",bType),
              ("tol",dType),("maxIts",iType),("maxStepRatio",dType),
              ("system",c_uint),("progress",bType)]
  def __init__(self):
    lib.ElQPDirectMehrotraCtrlDefault_d(pointer(self))

lib.ElQPDirectCtrlDefault_s.argtypes = \
lib.ElQPDirectCtrlDefault_d.argtypes = [c_void_p]
lib.ElQPDirectCtrlDefault_s.restype = \
lib.ElQPDirectCtrlDefault_d.restype = c_uint
class QPDirectCtrl_s(ctypes.Structure):
  _fields_ = [("approach",c_uint),
              ("ipfCtrl",QPDirectIPFCtrl_s),
              ("mehrotraCtrl",QPDirectMehrotraCtrl_s)]
  def __init__(self):
    lib.ElQPDirectCtrlDefault_s(pointer(self))
class QPDirectCtrl_d(ctypes.Structure):
  _fields_ = [("approach",c_uint),
              ("ipfCtrl",QPDirectIPFCtrl_d),
              ("mehrotraCtrl",QPDirectMehrotraCtrl_d)]
  def __init__(self):
    lib.ElQPDirectCtrlDefault_d(pointer(self))

lib.ElQPDirect_s.argtypes = \
lib.ElQPDirect_d.argtypes = \
lib.ElQPDirectDist_s.argtypes = \
lib.ElQPDirectDist_d.argtypes = \
lib.ElQPDirectSparse_s.argtypes = \
lib.ElQPDirectSparse_d.argtypes = \
lib.ElQPDirectDistSparse_s.argtypes = \
lib.ElQPDirectDistSparse_d.argtypes = \
  [c_void_p,c_void_p,
   c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p]
lib.ElQPDirect_s.restype = \
lib.ElQPDirect_d.restype = \
lib.ElQPDirectDist_s.restype = \
lib.ElQPDirectDist_d.restype = \
lib.ElQPDirectSparse_s.restype = \
lib.ElQPDirectSparse_d.restype = \
lib.ElQPDirectDistSparse_s.restype = \
lib.ElQPDirectDistSparse_d.restype = \
  c_uint

lib.ElQPDirectX_s.argtypes = \
lib.ElQPDirectXSparse_s.argtypes = \
lib.ElQPDirectXDist_s.argtypes = \
lib.ElQPDirectXDistSparse_s.argtypes = \
  [c_void_p,c_void_p,
   c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   QPDirectCtrl_s]
lib.ElQPDirectX_d.argtypes = \
lib.ElQPDirectXSparse_d.argtypes = \
lib.ElQPDirectXDist_d.argtypes = \
lib.ElQPDirectXDistSparse_d.argtypes = \
  [c_void_p,c_void_p,
   c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   QPDirectCtrl_d]
lib.ElQPDirectX_s.restype = \
lib.ElQPDirectX_d.restype = \
lib.ElQPDirectXSparse_s.restype = \
lib.ElQPDirectXSparse_d.restype = \
lib.ElQPDirectXDist_s.restype = \
lib.ElQPDirectXDist_d.restype = \
lib.ElQPDirectXDistSparse_s.restype = \
lib.ElQPDirectXDistSparse_d.restype = \
  c_uint

def QPDirect(Q,A,b,c,x,y,z,ctrl=None):
  if type(Q) is not type(A):
    raise Exception('A and Q must be of the same type')
  if type(b) is not type(c) or type(b) is not type(x) or \
     type(b) is not type(y) or type(b) is not type(z):
    raise Exception('{b,c,x,y,z} must be of the same type')
  args = [Q.obj,A.obj,b.obj,c.obj,x.obj,y.obj,z.obj]
  argsCtrl = [Q.obj,A.obj,b.obj,c.obj,x.obj,y.obj,z.obj,ctrl]
  if type(A) is Matrix:
    if type(b) is not Matrix: raise Exception('b must be a Matrix')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElQPDirect_s(*args)
      else:            lib.ElQPDirectX_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElQPDirect_d(*args)
      else:            lib.ElQPDirectX_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if type(b) is not DistMatrix: raise Exception('b must be a DistMatrix')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElQPDirectDist_s(*args)
      else:            lib.ElQPDirectXDist_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElQPDirectDist_d(*args)
      else:            lib.ElQPDirectXDist_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if type(b) is not Matrix: raise Exception('b must be a Matrix')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElQPDirectSparse_s(*args)
      else:            lib.ElQPDirectXSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElQPDirectSparse_d(*args)
      else:            lib.ElQPDirectXSparse_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if type(b) is not DistMultiVec: raise Exception('b must be a DistMultiVec')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElQPDirectDistSparse_s(*args)
      else:            lib.ElQPDirectXDistSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElQPDirectDistSparse_d(*args)
      else:            lib.ElQPDirectXDistSparse_d(*argsCtrl)
    else: DataExcept()
  else: TypeExcept()

# Affine conic form
# -----------------
lib.ElQPAffineIPFCtrlDefault_s.argtypes = \
lib.ElQPAffineIPFCtrlDefault_d.argtypes = [c_void_p]
lib.ElQPAffineIPFCtrlDefault_s.restype = \
lib.ElQPAffineIPFCtrlDefault_d.restype = c_uint
class QPAffineIPFCtrl_s(ctypes.Structure):
  _fields_ = [("primalInitialized",bType),("dualInitialized",bType),
              ("tol",sType),("maxIts",iType),("centering",sType),
              ("lineSearchCtrl",QPIPFLineSearchCtrl_s),
              ("progress",bType)]
  def __init__(self):
    lib.ElQPAffineIPFCtrlDefault_s(pointer(self))
class QPAffineIPFCtrl_d(ctypes.Structure):
  _fields_ = [("primalInitialized",bType),("dualInitialized",bType),
              ("tol",dType),("maxIts",iType),("centering",dType),
              ("lineSearchCtrl",QPIPFLineSearchCtrl_d),
              ("progress",bType)]
  def __init__(self):
    lib.ElQPAffineIPFCtrlDefault_d(pointer(self))

lib.ElQPAffineMehrotraCtrlDefault_s.argtypes = \
lib.ElQPAffineMehrotraCtrlDefault_d.argtypes = [c_void_p]
lib.ElQPAffineMehrotraCtrlDefault_s.restype = \
lib.ElQPAffineMehrotraCtrlDefault_d.restype = c_uint
class QPAffineMehrotraCtrl_s(ctypes.Structure):
  _fields_ = [("primalInitialized",bType),("dualInitialized",bType),
              ("tol",sType),("maxIts",iType),("maxStepRatio",sType),
              ("progress",bType)]
  def __init__(self):
    lib.ElQPAffineMehrotraCtrlDefault_s(pointer(self))
class QPAffineMehrotraCtrl_d(ctypes.Structure):
  _fields_ = [("primalInitialized",bType),("dualInitialized",bType),
              ("tol",dType),("maxIts",iType),("maxStepRatio",dType),
              ("progress",bType)]
  def __init__(self):
    lib.ElQPAffineMehrotraCtrlDefault_d(pointer(self))

lib.ElQPAffineCtrlDefault_s.argtypes = \
lib.ElQPAffineCtrlDefault_d.argtypes = [c_void_p]
lib.ElQPAffineCtrlDefault_s.restype = \
lib.ElQPAffineCtrlDefault_d.restype = c_uint
class QPAffineCtrl_s(ctypes.Structure):
  _fields_ = [("approach",c_uint),
              ("ipfCtrl",QPAffineIPFCtrl_s),
              ("mehrotraCtrl",QPAffineMehrotraCtrl_s)]
  def __init__(self):
    lib.ElQPAffineCtrlDefault_s(pointer(self))
class QPAffineCtrl_d(ctypes.Structure):
  _fields_ = [("approach",c_uint),
              ("ipfCtrl",QPAffineIPFCtrl_d),
              ("mehrotraCtrl",QPAffineMehrotraCtrl_d)]
  def __init__(self):
    lib.ElQPAffineCtrlDefault_d(pointer(self))

lib.ElQPAffine_s.argtypes = \
lib.ElQPAffine_d.argtypes = \
lib.ElQPAffineDist_s.argtypes = \
lib.ElQPAffineDist_d.argtypes = \
lib.ElQPAffineSparse_s.argtypes = \
lib.ElQPAffineSparse_d.argtypes = \
lib.ElQPAffineDistSparse_s.argtypes = \
lib.ElQPAffineDistSparse_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElQPAffine_s.restype = \
lib.ElQPAffine_d.restype = \
lib.ElQPAffineDist_s.restype = \
lib.ElQPAffineDist_d.restype = \
lib.ElQPAffineSparse_s.restype = \
lib.ElQPAffineSparse_d.restype = \
lib.ElQPAffineDistSparse_s.restype = \
lib.ElQPAffineDistSparse_d.restype = \
  c_uint

lib.ElQPAffineX_s.argtypes = \
lib.ElQPAffineXDist_s.argtypes = \
lib.ElQPAffineXSparse_s.argtypes = \
lib.ElQPAffineXDistSparse_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,c_void_p,
   QPAffineCtrl_s]
lib.ElQPAffineX_d.argtypes = \
lib.ElQPAffineXDist_d.argtypes = \
lib.ElQPAffineXSparse_d.argtypes = \
lib.ElQPAffineXDistSparse_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,c_void_p,
   QPAffineCtrl_d]
lib.ElQPAffineX_s.restype = \
lib.ElQPAffineX_d.restype = \
lib.ElQPAffineXDist_s.restype = \
lib.ElQPAffineXDist_d.restype = \
lib.ElQPAffineXSparse_s.restype = \
lib.ElQPAffineXSparse_d.restype = \
lib.ElQPAffineXDistSparse_s.restype = \
lib.ElQPAffineXDistSparse_d.restype = \
  c_uint

def QPAffine(Q,A,G,b,c,h,x,y,z,s,ctrl=None):
  if type(Q) is not type(A) or type(A) is not type(G):
    raise Exception('{Q,A,G} must be of the same type')
  if type(b) is not type(c) or type(b) is not type(c) or \
     type(b) is not type(h) or type(b) is not type(x) or \
     type(b) is not type(y) or type(b) is not type(z) or \
     type(b) is not type(s):
    raise Exception('{b,c,h,x,y,z,s} must be of the same type')
  args = [Q.obj,A.obj,G.obj,b.obj,c.obj,h.obj,x.obj,y.obj,z.obj,s.obj]
  argsCtrl = [Q.obj,A.obj,G.obj,b.obj,c.obj,h.obj,x.obj,y.obj,z.obj,s.obj,ctrl]
  if type(A) is Matrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElQPAffine_s(*args)
      else:            lib.ElQPAffineX_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElQPAffine_d(*args)
      else:            lib.ElQPAffineX_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if type(b) is not DistMatrix:
      raise Exception('b must be a DistMatrix')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElQPAffineDist_s(*args)
      else:            lib.ElQPAffineXDist_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElQPAffineDist_d(*args)
      else:            lib.ElQPAffineXDist_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElQPAffineSparse_s(*args)
      else:            lib.ElQPAffineXSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElQPAffineSparse_d(*args)
      else:            lib.ElQPAffineXSparse_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if type(b) is not DistMultiVec:
      raise Exception('b must be a DistMultiVec')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElQPAffineDistSparse_s(*args)
      else:            lib.ElQPAffineXDistSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElQPAffineDistSparse_d(*args)
      else:            lib.ElQPAffineXDistSparse_d(*argsCtrl)
    else: DataExcept()
  else: TypeExcept()

# Box form
# --------
lib.ElQPBoxADMM_s.argtypes = \
lib.ElQPBoxADMMDist_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,c_void_p,POINTER(iType)]
lib.ElQPBoxADMM_d.argtypes = \
lib.ElQPBoxADMMDist_d.argtypes = \
  [c_void_p,c_void_p,dType,dType,c_void_p,POINTER(iType)]
lib.ElQPBoxADMM_s.restype = \
lib.ElQPBoxADMM_d.restype = \
lib.ElQPBoxADMMDist_s.restype = \
lib.ElQPBoxADMMDist_d.restype = \
  c_uint

def QPBoxADMM(Q,C,lb,ub):
  if type(Q) is not type(C):
    raise Exception('Types of Q and C must match')
  if Q.tag != C.tag:
    raise Exception('Datatypes of Q and C must match')
  numIts = iType()
  if type(Q) is Matrix:
    Z = Matrix(Q.tag)
    args = [Q.obj,C.obj,lb,ub,Z.obj,pointer(numIts)]
    if   Q.tag == sTag: lib.ElQPBoxADMM_s(*args)
    elif Q.tag == dTag: lib.ElQPBoxADMM_d(*args)
    else: DataExcept()
    return Z, numIts
  elif type(Q) is DistMatrix:
    Z = DistMatrix(Q.tag,MC,MR,Q.Grid())
    args = [Q.obj,C.obj,lb,ub,Z.obj,pointer(numIts)]
    if   Q.tag == sTag: lib.ElQPBoxADMMDist_s(*args)
    elif Q.tag == dTag: lib.ElQPBoxADMMDist_d(*args)
    else: DataExcept()
    return Z, numIts
  else: TypeExcept()

# Robust Principal Component Analysis
# ===================================
lib.ElRPCA_s.argtypes = \
lib.ElRPCA_d.argtypes = \
lib.ElRPCA_c.argtypes = \
lib.ElRPCA_z.argtypes = \
lib.ElRPCADist_s.argtypes = \
lib.ElRPCADist_d.argtypes = \
lib.ElRPCADist_c.argtypes = \
lib.ElRPCADist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p]
lib.ElRPCA_s.restype = \
lib.ElRPCA_d.restype = \
lib.ElRPCA_c.restype = \
lib.ElRPCA_z.restype = \
lib.ElRPCADist_s.restype = \
lib.ElRPCADist_d.restype = \
lib.ElRPCADist_c.restype = \
lib.ElRPCADist_z.restype = \
  c_uint

def RPCA(M):
  if type(M) is Matrix:
    L = Matrix(M.tag)
    S = Matrix(M.tag)
    args = [M.obj,L.obj,S.obj]
    if   M.tag == sTag: lib.ElRPCA_s(*args)
    elif M.tag == dTag: lib.ElRPCA_d(*args)
    elif M.tag == cTag: lib.ElRPCA_c(*args)
    elif M.tag == zTag: lib.ElRPCA_z(*args)
    return L, S
  elif type(M) is DistMatrix:
    L = DistMatrix(M.tag,MC,MR,M.Grid())
    S = DistMatrix(M.tag,MC,MR,M.Grid())
    args = [M.obj,L.obj,S.obj]
    if   M.tag == sTag: lib.ElRPCADist_s(*args)
    elif M.tag == dTag: lib.ElRPCADist_d(*args)
    elif M.tag == cTag: lib.ElRPCADist_c(*args)
    elif M.tag == zTag: lib.ElRPCADist_z(*args)
    return L, S
  else: TypeExcept()

# Sparse inverse covariance selection
# ===================================
lib.ElSparseInvCov_s.argtypes = \
lib.ElSparseInvCov_c.argtypes = \
lib.ElSparseInvCovDist_s.argtypes = \
lib.ElSparseInvCovDist_c.argtypes = \
  [c_void_p,sType,c_void_p,POINTER(iType)]
lib.ElSparseInvCov_d.argtypes = \
lib.ElSparseInvCov_z.argtypes = \
lib.ElSparseInvCovDist_d.argtypes = \
lib.ElSparseInvCovDist_z.argtypes = \
  [c_void_p,dType,c_void_p,POINTER(iType)]

lib.ElSparseInvCov_s.restype = \
lib.ElSparseInvCov_d.restype = \
lib.ElSparseInvCov_c.restype = \
lib.ElSparseInvCov_z.restype = \
lib.ElSparseInvCovDist_s.restype = \
lib.ElSparseInvCovDist_d.restype = \
lib.ElSparseInvCovDist_c.restype = \
lib.ElSparseInvCovDist_z.restype = \
  c_uint

def SparseInvCov(D,lamb):
  numIts = iType()
  if type(D) is Matrix:
    Z = Matrix(D.tag)
    args = [D.obj,lamb,Z.obj,pointer(numIts)]
    if   D.tag == sTag: lib.ElSparseInvCov_s(*args)
    elif D.tag == dTag: lib.ElSparseInvCov_d(*args)
    elif D.tag == cTag: lib.ElSparseInvCov_c(*args)
    elif D.tag == zTag: lib.ElSparseInvCov_z(*args)
    else: DataExcept()
    return Z, numIts
  elif type(D) is DistMatrix:
    Z = DistMatrix(D.tag,MC,MR,D.Grid())
    args = [D.obj,lamb,Z.obj,pointer(numIts)]
    if   D.tag == sTag: lib.ElSparseInvCovDist_s(*args)
    elif D.tag == dTag: lib.ElSparseInvCovDist_d(*args)
    elif D.tag == cTag: lib.ElSparseInvCovDist_c(*args)
    elif D.tag == zTag: lib.ElSparseInvCovDist_z(*args)
    else: DataExcept()
    return Z, numIts
  else: TypeExcept()

# Support Vector Machine
# ======================
lib.ElSVM_s.argtypes = \
lib.ElSVMDist_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,sType,POINTER(iType)]
lib.ElSVM_d.argtypes = \
lib.ElSVMDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,dType,POINTER(iType)]
lib.ElSVM_s.restype = \
lib.ElSVM_d.restype = \
lib.ElSVMDist_s.restype = \
lib.ElSVMDist_d.restype = \
  c_uint

def SVM(G,q,gamma):
  if type(G) is not type(q):
    raise Exception('Types of G and q must match')
  if G.tag != q.tag:
    raise Exception('Datatypes of G and q must match')
  numIts = iType()
  if type(G) is Matrix:
    z = Matrix(G.tag)
    args = [G.obj,q.obj,z.obj,gamma,pointer(numIts)]
    if   G.tag == sTag: lib.ElSVM_s(*args)
    elif G.tag == dTag: lib.ElSVM_d(*args)
    else: DataExcept()
    return z, numIts
  elif type(G) is DistMatrix:
    z = DistMatrix(G.tag,MC,MR,G.Grid())
    args = [G.obj,q.obj,z.obj,gamma,pointer(numIts)]
    if   G.tag == sTag: lib.ElSVMDist_s(*args)
    elif G.tag == dTag: lib.ElSVMDist_d(*args)
    else: DataExcept()
    return z, numIts
  else: TypeExcept()

# Utilities
# =========

# Clipping
# --------
lib.ElLowerClip_s.argtypes = \
lib.ElLowerClipDist_s.argtypes = \
  [c_void_p,sType]
lib.ElLowerClip_d.argtypes = \
lib.ElLowerClipDist_d.argtypes = \
  [c_void_p,dType]
lib.ElLowerClip_s.restype = \
lib.ElLowerClip_d.restype = \
lib.ElLowerClipDist_s.restype = \
lib.ElLowerClipDist_d.restype = \
  c_uint

def LowerClip(X,lowerBound=0):
  args = [X.obj,lowerBound]
  if type(X) is Matrix:
    if   X.tag == sTag: lib.ElLowerClip_s(*args)
    elif X.tag == dTag: lib.ElLowerClip_d(*args)
    else: DataExcept()
  elif type(X) is DistMatrix:
    if   X.tag == sTag: lib.ElLowerClipDist_s(*args)
    elif X.tag == dTag: lib.ElLowerClipDist_d(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElUpperClip_s.argtypes = \
lib.ElUpperClipDist_s.argtypes = \
  [c_void_p,sType]
lib.ElUpperClip_d.argtypes = \
lib.ElUpperClipDist_d.argtypes = \
  [c_void_p,dType]
lib.ElUpperClip_s.restype = \
lib.ElUpperClip_d.restype = \
lib.ElUpperClipDist_s.restype = \
lib.ElUpperClipDist_d.restype = \
  c_uint

def UpperClip(X,upperBound=1):
  args = [X.obj,upperBound]
  if type(X) is Matrix:
    if   X.tag == sTag: lib.ElUpperClip_s(*args)
    elif X.tag == dTag: lib.ElUpperClip_d(*args)
    else: DataExcept()
  elif type(X) is DistMatrix:
    if   X.tag == sTag: lib.ElUpperClipDist_s(*args)
    elif X.tag == dTag: lib.ElUpperClipDist_d(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElClip_s.argtypes = \
lib.ElClipDist_s.argtypes = \
  [c_void_p,sType,sType]
lib.ElClip_d.argtypes = \
lib.ElClipDist_d.argtypes = \
  [c_void_p,dType,dType]
lib.ElClip_s.restype = \
lib.ElClip_d.restype = \
lib.ElClipDist_s.restype = \
lib.ElClipDist_d.restype = \
  c_uint

def Clip(X,lowerBound=0,upperBound=1):
  args = [X.obj,lowerBound,upperBound]
  if type(X) is Matrix:
    if   X.tag == sTag: lib.ElClip_s(*args)
    elif X.tag == dTag: lib.ElClip_d(*args)
    else: DataExcept()
  elif type(X) is DistMatrix:
    if   X.tag == sTag: lib.ElClipDist_s(*args)
    elif X.tag == dTag: lib.ElClipDist_d(*args)
    else: DataExcept()
  else: TypeExcept()

# Coherence
# ---------
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
# ----------
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

# Frobenius-norm proximal map
# ---------------------------
lib.ElFrobeniusProx_s.argtypes = \
lib.ElFrobeniusProx_c.argtypes = \
lib.ElFrobeniusProxDist_s.argtypes = \
lib.ElFrobeniusProxDist_c.argtypes = \
  [c_void_p,sType]
lib.ElFrobeniusProx_d.argtypes = \
lib.ElFrobeniusProx_z.argtypes = \
lib.ElFrobeniusProxDist_d.argtypes = \
lib.ElFrobeniusProxDist_z.argtypes = \
  [c_void_p,dType]

lib.ElFrobeniusProx_s.restype = \
lib.ElFrobeniusProx_d.restype = \
lib.ElFrobeniusProx_c.restype = \
lib.ElFrobeniusProx_z.restype = \
lib.ElFrobeniusProxDist_s.restype = \
lib.ElFrobeniusProxDist_d.restype = \
lib.ElFrobeniusProxDist_c.restype = \
lib.ElFrobeniusProxDist_z.restype = \
  c_uint

def FrobeniusProx(A,rho):
  args = [A.obj,rho]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElFrobeniusProx_s(*args)
    elif A.tag == dTag: lib.ElFrobeniusProx_d(*args)
    elif A.tag == cTag: lib.ElFrobeniusProx_c(*args)
    elif A.tag == zTag: lib.ElFrobeniusProx_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElFrobeniusProxDist_s(*args)
    elif A.tag == dTag: lib.ElFrobeniusProxDist_d(*args)
    elif A.tag == cTag: lib.ElFrobeniusProxDist_c(*args)
    elif A.tag == zTag: lib.ElFrobeniusProxDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Hinge-loss proximal map
# -----------------------
lib.ElHingeLossProx_s.argtypes = \
lib.ElHingeLossProxDist_s.argtypes = \
  [c_void_p,sType]
lib.ElHingeLossProx_d.argtypes = \
lib.ElHingeLossProxDist_d.argtypes = \
  [c_void_p,dType]

lib.ElHingeLossProx_s.restype = \
lib.ElHingeLossProx_d.restype = \
lib.ElHingeLossProxDist_s.restype = \
lib.ElHingeLossProxDist_d.restype = \
  c_uint

def HingeLossProx(A,rho):
  args = [A.obj,rho]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHingeLossProx_s(*args)
    elif A.tag == dTag: lib.ElHingeLossProx_d(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHingeLossProxDist_s(*args)
    elif A.tag == dTag: lib.ElHingeLossProxDist_d(*args)
    else: DataExcept()
  else: TypeExcept()

# Log barrier
# -----------
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
# ------------------
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

# Logistic proximal map
# ---------------------
lib.ElLogisticProx_s.argtypes = \
lib.ElLogisticProxDist_s.argtypes = \
  [c_void_p,sType]
lib.ElLogisticProx_d.argtypes = \
lib.ElLogisticProxDist_d.argtypes = \
  [c_void_p,dType]

lib.ElLogisticProx_s.restype = \
lib.ElLogisticProx_d.restype = \
lib.ElLogisticProxDist_s.restype = \
lib.ElLogisticProxDist_d.restype = \
  c_uint

def LogisticProx(A,rho):
  args = [A.obj,rho]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElLogisticProx_s(*args)
    elif A.tag == dTag: lib.ElLogisticProx_d(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElLogisticProxDist_s(*args)
    elif A.tag == dTag: lib.ElLogisticProxDist_d(*args)
    else: DataExcept()
  else: TypeExcept()

# Singular-value soft-thresholding
# --------------------------------
lib.ElSVT_s.argtypes = \
lib.ElSVT_c.argtypes = \
lib.ElSVTDist_s.argtypes = \
lib.ElSVTDist_c.argtypes = \
  [c_void_p,sType,bType]
lib.ElSVT_d.argtypes = \
lib.ElSVT_z.argtypes = \
lib.ElSVTDist_d.argtypes = \
lib.ElSVTDist_z.argtypes = \
  [c_void_p,dType,bType]

lib.ElSVT_s.restype = \
lib.ElSVT_d.restype = \
lib.ElSVT_c.restype = \
lib.ElSVT_z.restype = \
lib.ElSVTDist_s.restype = \
lib.ElSVTDist_d.restype = \
lib.ElSVTDist_c.restype = \
lib.ElSVTDist_z.restype = \
  c_uint

def SVT(A,rho,relative=False):
  args = [A.obj,rho,relative]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSVT_s(*args)
    elif A.tag == dTag: lib.ElSVT_d(*args)
    elif A.tag == cTag: lib.ElSVT_c(*args)
    elif A.tag == zTag: lib.ElSVT_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSVTDist_s(*args)
    elif A.tag == dTag: lib.ElSVTDist_d(*args)
    elif A.tag == cTag: lib.ElSVTDist_c(*args)
    elif A.tag == zTag: lib.ElSVTDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Soft-thresholding
# -----------------
lib.ElSoftThreshold_s.argtypes = \
lib.ElSoftThreshold_c.argtypes = \
lib.ElSoftThresholdDist_s.argtypes = \
lib.ElSoftThresholdDist_c.argtypes = \
  [c_void_p,sType,bType]
lib.ElSoftThreshold_d.argtypes = \
lib.ElSoftThreshold_z.argtypes = \
lib.ElSoftThresholdDist_d.argtypes = \
lib.ElSoftThresholdDist_z.argtypes = \
  [c_void_p,dType,bType]

lib.ElSoftThreshold_s.restype = \
lib.ElSoftThreshold_d.restype = \
lib.ElSoftThreshold_c.restype = \
lib.ElSoftThreshold_z.restype = \
lib.ElSoftThresholdDist_s.restype = \
lib.ElSoftThresholdDist_d.restype = \
lib.ElSoftThresholdDist_c.restype = \
lib.ElSoftThresholdDist_z.restype = \
  c_uint

def SoftThreshold(A,rho,relative=False):
  args = [A.obj,rho,relative]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSoftThreshold_s(*args)
    elif A.tag == dTag: lib.ElSoftThreshold_d(*args)
    elif A.tag == cTag: lib.ElSoftThreshold_c(*args)
    elif A.tag == zTag: lib.ElSoftThreshold_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSoftThresholdDist_s(*args)
    elif A.tag == dTag: lib.ElSoftThresholdDist_d(*args)
    elif A.tag == cTag: lib.ElSoftThresholdDist_c(*args)
    elif A.tag == zTag: lib.ElSoftThresholdDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
