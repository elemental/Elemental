#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *

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
  if A.tag != b.tag or b.tag != c.tag or c.tag != x.tag or \
     x.tag != y.tag or y.tag != z.tag:
    raise Exception('Datatypes of {A,b,c,x,y,z} must match')
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

