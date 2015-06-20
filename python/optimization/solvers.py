#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *
from El.lapack_like.factor import *

(FULL_KKT,AUGMENTED_KKT,NORMAL_KKT) = (0,1,2)

lib.ElIPFLineSearchCtrlDefault_s.argtypes = \
lib.ElIPFLineSearchCtrlDefault_d.argtypes = \
  [c_void_p]
class IPFLineSearchCtrl_s(ctypes.Structure):
  _fields_ = [("gamma",sType),("beta",sType),("psi",sType),
              ("stepRatio",sType),("progress",bType)]
  def __init__(self):
    lib.ElIPFLineSearchCtrlDefault_s(pointer(self))
class IPFLineSearchCtrl_d(ctypes.Structure):
  _fields_ = [("gamma",dType),("beta",dType),("psi",dType),
              ("stepRatio",dType),("progress",bType)]
  def __init__(self):
    lib.ElIPFLineSearchCtrlDefault_d(pointer(self))

# Linear program
# ==============

(LP_ADMM,LP_IPF,LP_IPF_SELFDUAL,LP_MEHROTRA,LP_MEHROTRA_SELFDUAL)=(0,1,2,3,4)

# Direct conic form
# -----------------
lib.ElLPDirectADMMCtrlDefault_s.argtypes = \
lib.ElLPDirectADMMCtrlDefault_d.argtypes = \
  [c_void_p]
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

lib.ElLPDirectIPFCtrlDefault_s.argtypes = \
lib.ElLPDirectIPFCtrlDefault_d.argtypes = \
  [c_void_p,bType]
class LPDirectIPFCtrl_s(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",sType),("targetTol",sType),
              ("maxIts",iType),("centering",sType),
              ("system",c_uint),("qsdCtrl",RegQSDCtrl_s),
              ("lineSearchCtrl",IPFLineSearchCtrl_s),
              ("equilibrate",bType),("progress",bType),("time",bType)]
  def __init__(self,isSparse=True):
    lib.ElLPDirectIPFCtrlDefault_s(pointer(self),isSparse)
class LPDirectIPFCtrl_d(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",dType),("targetTol",dType),
              ("maxIts",iType),("centering",dType),
              ("system",c_uint),("qsdCtrl",RegQSDCtrl_d),
              ("lineSearchCtrl",IPFLineSearchCtrl_d),
              ("equilibrate",bType),("progress",bType),("time",bType)]
  def __init__(self,isSparse=True):
    lib.ElLPDirectIPFCtrlDefault_d(pointer(self),isSparse)

lib.ElLPDirectMehrotraCtrlDefault_s.argtypes = \
lib.ElLPDirectMehrotraCtrlDefault_d.argtypes = \
  [c_void_p,bType]
class LPDirectMehrotraCtrl_s(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",sType),("targetTol",sType),
              ("maxIts",iType),("maxStepRatio",sType),
              ("system",c_uint),("qsdCtrl",RegQSDCtrl_s),
              ("outerEquil",bType),("innerEquil",bType),
              ("scaleTwoNorm",bType),("basisSize",iType),
              ("progress",bType),("time",bType)]
  def __init__(self,isSparse=True):
    lib.ElLPDirectMehrotraCtrlDefault_s(pointer(self),isSparse)
class LPDirectMehrotraCtrl_d(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",dType),("targetTol",dType),
              ("maxIts",iType),("maxStepRatio",dType),
              ("system",c_uint),("qsdCtrl",RegQSDCtrl_d),
              ("outerEquil",bType),("innerEquil",bType),
              ("scaleTwoNorm",bType),("basisSize",iType),
              ("progress",bType),("time",bType)]
  def __init__(self,isSparse=True):
    lib.ElLPDirectMehrotraCtrlDefault_d(pointer(self),isSparse)

lib.ElLPDirectCtrlDefault_s.argtypes = \
lib.ElLPDirectCtrlDefault_d.argtypes = \
  [c_void_p,bType]
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
lib.ElLPAffineIPFCtrlDefault_d.argtypes = \
  [c_void_p]
class LPAffineIPFCtrl_s(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",sType),("targetTol",sType),
              ("maxIts",iType),("centering",sType),
              ("qsdCtrl",RegQSDCtrl_s),
              ("lineSearchCtrl",IPFLineSearchCtrl_s),
              ("equilibrate",bType),("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElLPAffineIPFCtrlDefault_s(pointer(self))
class LPAffineIPFCtrl_d(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",dType),("targetTol",dType),
              ("maxIts",iType),("centering",dType),
              ("qsdCtrl",RegQSDCtrl_d),
              ("lineSearchCtrl",IPFLineSearchCtrl_d),
              ("equilibrate",bType),("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElLPAffineIPFCtrlDefault_d(pointer(self))

lib.ElLPAffineMehrotraCtrlDefault_s.argtypes = \
lib.ElLPAffineMehrotraCtrlDefault_d.argtypes = \
  [c_void_p]
class LPAffineMehrotraCtrl_s(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",sType),("targetTol",sType),
              ("maxIts",iType),("maxStepRatio",sType),
              ("qsdCtrl",RegQSDCtrl_s),
              ("outerEquil",bType),("innerEquil",bType),
              ("scaleTwoNorm",bType),("basisSize",iType),
              ("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElLPAffineMehrotraCtrlDefault_s(pointer(self))
class LPAffineMehrotraCtrl_d(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",dType),("targetTol",dType),
              ("maxIts",iType),("maxStepRatio",dType),
              ("qsdCtrl",RegQSDCtrl_d),
              ("outerEquil",bType),("innerEquil",bType),
              ("scaleTwoNorm",bType),("basisSize",iType),
              ("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElLPAffineMehrotraCtrlDefault_d(pointer(self))

lib.ElLPAffineCtrlDefault_s.argtypes = \
lib.ElLPAffineCtrlDefault_d.argtypes = \
  [c_void_p]
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

# Direct conic form
# -----------------
lib.ElQPDirectIPFCtrlDefault_s.argtypes = \
lib.ElQPDirectIPFCtrlDefault_d.argtypes = \
  [c_void_p]
class QPDirectIPFCtrl_s(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",sType),("targetTol",sType),
              ("maxIts",iType),("centering",sType),
              ("system",c_uint),("qsdCtrl",RegQSDCtrl_s),
              ("lineSearchCtrl",IPFLineSearchCtrl_s),("equilibrate",bType),
              ("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElQPDirectIPFCtrlDefault_s(pointer(self))
class QPDirectIPFCtrl_d(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",dType),("targetTol",dType),
              ("maxIts",iType),("centering",dType),
              ("system",c_uint),("qsdCtrl",RegQSDCtrl_d),
              ("lineSearchCtrl",IPFLineSearchCtrl_d),("equilibrate",bType),
              ("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElQPDirectIPFCtrlDefault_d(pointer(self))

lib.ElQPDirectMehrotraCtrlDefault_s.argtypes = \
lib.ElQPDirectMehrotraCtrlDefault_d.argtypes = \
  [c_void_p]
class QPDirectMehrotraCtrl_s(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",sType),("targetTol",sType),
              ("maxIts",iType),("maxStepRatio",sType),
              ("system",c_uint),("qsdCtrl",RegQSDCtrl_s),
              ("outerEquil",bType),("innerEquil",bType),
              ("scaleTwoNorm",bType),("basisSize",iType),
              ("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElQPDirectMehrotraCtrlDefault_s(pointer(self))
class QPDirectMehrotraCtrl_d(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",dType),("targetTol",dType),
              ("maxIts",iType),("maxStepRatio",dType),
              ("system",c_uint),("qsdCtrl",RegQSDCtrl_d),
              ("outerEquil",bType),("innerEquil",bType),
              ("scaleTwoNorm",bType),("basisSize",iType),
              ("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElQPDirectMehrotraCtrlDefault_d(pointer(self))

lib.ElQPDirectCtrlDefault_s.argtypes = \
lib.ElQPDirectCtrlDefault_d.argtypes = \
  [c_void_p]
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
lib.ElQPAffineIPFCtrlDefault_d.argtypes = \
  [c_void_p]
class QPAffineIPFCtrl_s(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",sType),("targetTol",sType),
              ("maxIts",iType),("centering",sType),
              ("qsdCtrl",RegQSDCtrl_s),
              ("lineSearchCtrl",IPFLineSearchCtrl_s),
              ("equilibrate",bType),("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElQPAffineIPFCtrlDefault_s(pointer(self))
class QPAffineIPFCtrl_d(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",dType),("targetTol",dType),
              ("maxIts",iType),("centering",dType),
              ("qsdCtrl",RegQSDCtrl_d),
              ("lineSearchCtrl",IPFLineSearchCtrl_d),
              ("equilibrate",bType),("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElQPAffineIPFCtrlDefault_d(pointer(self))

lib.ElQPAffineMehrotraCtrlDefault_s.argtypes = \
lib.ElQPAffineMehrotraCtrlDefault_d.argtypes = \
  [c_void_p]
class QPAffineMehrotraCtrl_s(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",sType),("targetTol",sType),
              ("maxIts",iType),("maxStepRatio",sType),
              ("qsdCtrl",RegQSDCtrl_s),
              ("outerEquil",bType),("innerEquil",bType),
              ("scaleTwoNorm",bType),("basisSize",iType),
              ("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElQPAffineMehrotraCtrlDefault_s(pointer(self))
class QPAffineMehrotraCtrl_d(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",dType),("targetTol",dType),
              ("maxIts",iType),("maxStepRatio",dType),
              ("qsdCtrl",RegQSDCtrl_d),
              ("outerEquil",bType),("innerEquil",bType),
              ("scaleTwoNorm",bType),("basisSize",iType),
              ("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElQPAffineMehrotraCtrlDefault_d(pointer(self))

lib.ElQPAffineCtrlDefault_s.argtypes = \
lib.ElQPAffineCtrlDefault_d.argtypes = \
  [c_void_p]
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
lib.ElQPBoxADMMCtrlDefault_s.argtypes = \
lib.ElQPBoxADMMCtrlDefault_d.argtypes = \
  [c_void_p]
class QPBoxADMMCtrl_s(ctypes.Structure):
  _fields_ = [("rho",sType),("alpha",sType),
              ("maxIter",iType),
              ("absTol",sType),("relTol",sType),
              ("inv",bType),("progress",bType)]
  def __init__(self):
    lib.ElQPBoxADMMCtrlDefault_s(pointer(self))
class QPBoxADMMCtrl_d(ctypes.Structure):
  _fields_ = [("rho",dType),("alpha",dType),
              ("maxIter",iType),
              ("absTol",dType),("relTol",dType),
              ("inv",bType),("progress",bType)]
  def __init__(self):
    lib.ElQPBoxADMMCtrlDefault_d(pointer(self))

lib.ElQPBoxADMM_s.argtypes = \
lib.ElQPBoxADMMDist_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,c_void_p,POINTER(iType)]
lib.ElQPBoxADMM_d.argtypes = \
lib.ElQPBoxADMMDist_d.argtypes = \
  [c_void_p,c_void_p,dType,dType,c_void_p,POINTER(iType)]

lib.ElQPBoxADMMX_s.argtypes = \
lib.ElQPBoxADMMXDist_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,c_void_p,QPBoxADMMCtrl_s,POINTER(iType)]
lib.ElQPBoxADMMX_d.argtypes = \
lib.ElQPBoxADMMXDist_d.argtypes = \
  [c_void_p,c_void_p,dType,dType,c_void_p,QPBoxADMMCtrl_d,POINTER(iType)]

def QPBoxADMM(Q,C,lb,ub,ctrl=None):
  if type(Q) is not type(C):
    raise Exception('Types of Q and C must match')
  if Q.tag != C.tag:
    raise Exception('Datatypes of Q and C must match')
  numIts = iType()
  if type(Q) is Matrix:
    Z = Matrix(Q.tag)
    args = [Q.obj,C.obj,lb,ub,Z.obj,pointer(numIts)]
    argsCtrl = [Q.obj,C.obj,lb,ub,Z.obj,ctrl,pointer(numIts)]
    if   Q.tag == sTag: 
      if ctrl==None: lib.ElQPBoxADMM_s(*args)
      else:          lib.ElQPBoxADMMX_s(*argsCtrl)
    elif Q.tag == dTag: 
      if ctrl==None: lib.ElQPBoxADMM_d(*args)
      else:          lib.ElQPBoxADMMX_d(*argsCtrl)
    else: DataExcept()
    return Z, numIts
  elif type(Q) is DistMatrix:
    Z = DistMatrix(Q.tag,MC,MR,Q.Grid())
    args = [Q.obj,C.obj,lb,ub,Z.obj,pointer(numIts)]
    argsCtrl = [Q.obj,C.obj,lb,ub,Z.obj,ctrl,pointer(numIts)]
    if   Q.tag == sTag: 
      if ctrl==None: lib.ElQPBoxADMMDist_s(*args)
      else:          lib.ElQPBoxADMMXDist_s(*argsCtrl)
    elif Q.tag == dTag: 
      if ctrl==None: lib.ElQPBoxADMMDist_d(*args)
      else:          lib.ElQPBoxADMMXDist_d(*argsCtrl)
    else: DataExcept()
    return Z, numIts
  else: TypeExcept()

# Second-order cone programs
# ==========================
(SOCP_ADMM,SOCP_IPF,SOCP_IPF_SELFDUAL,SOCP_MEHROTRA,SOCP_MEHROTRA_SELFDUAL)=\
  (0,1,2,3,4)

# Direct conic form
# -----------------
lib.ElSOCPDirectMehrotraCtrlDefault_s.argtypes = \
lib.ElSOCPDirectMehrotraCtrlDefault_d.argtypes = \
  [c_void_p,bType]
class SOCPDirectMehrotraCtrl_s(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",sType),("targetTol",sType),
              ("maxIts",iType),("maxStepRatio",sType),
              ("qsdCtrl",RegQSDCtrl_s),
              ("outerEquil",bType),("innerEquil",bType),
              ("scaleTwoNorm",bType),("basisSize",iType),
              ("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElSOCPDirectMehrotraCtrlDefault_s(pointer(self))
class SOCPDirectMehrotraCtrl_d(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",dType),("targetTol",dType),
              ("maxIts",iType),("maxStepRatio",dType),
              ("qsdCtrl",RegQSDCtrl_d),
              ("outerEquil",bType),("innerEquil",bType),
              ("scaleTwoNorm",bType),("basisSize",iType),
              ("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElSOCPDirectMehrotraCtrlDefault_d(pointer(self))

lib.ElSOCPDirectCtrlDefault_s.argtypes = \
lib.ElSOCPDirectCtrlDefault_d.argtypes = \
  [c_void_p,bType]
class SOCPDirectCtrl_s(ctypes.Structure):
  _fields_ = [("approach",c_uint),
              ("mehrotraCtrl",SOCPDirectMehrotraCtrl_s)]
  def __init__(self):
    lib.ElSOCPDirectCtrlDefault_s(pointer(self))
class SOCPDirectCtrl_d(ctypes.Structure):
  _fields_ = [("approach",c_uint),
              ("mehrotraCtrl",SOCPDirectMehrotraCtrl_d)]
  def __init__(self):
    lib.ElSOCPDirectCtrlDefault_d(pointer(self))

lib.ElSOCPDirect_s.argtypes = \
lib.ElSOCPDirect_d.argtypes = \
lib.ElSOCPDirectDist_s.argtypes = \
lib.ElSOCPDirectDist_d.argtypes = \
lib.ElSOCPDirectSparse_s.argtypes = \
lib.ElSOCPDirectSparse_d.argtypes = \
lib.ElSOCPDirectDistSparse_s.argtypes = \
lib.ElSOCPDirectDistSparse_d.argtypes = \
  [c_void_p,
   c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p]

lib.ElSOCPDirectX_s.argtypes = \
lib.ElSOCPDirectXSparse_s.argtypes = \
lib.ElSOCPDirectXDist_s.argtypes = \
lib.ElSOCPDirectXDistSparse_s.argtypes = \
  [c_void_p,
   c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   SOCPDirectCtrl_s]
lib.ElSOCPDirectX_d.argtypes = \
lib.ElSOCPDirectXSparse_d.argtypes = \
lib.ElSOCPDirectXDist_d.argtypes = \
lib.ElSOCPDirectXDistSparse_d.argtypes = \
  [c_void_p,
   c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   SOCPDirectCtrl_d]

def SOCPDirect(A,b,c,orders,firstInds,labels,x,y,z,ctrl=None):
  if A.tag != b.tag or b.tag != c.tag or c.tag != x.tag or \
     x.tag != y.tag or y.tag != z.tag:
    raise Exception('Datatypes of {A,b,c,x,y,z} must match')
  if orders.tag != iTag or firstInds.tag != iTag or labels.tag != iTag:
    raise Exception('Datatypes of conic descriptions should be integers')
  if type(b) is not type(c) or type(b) is not type(x) or \
     type(b) is not type(y) or type(b) is not type(z) or \
     type(b) is not type(orders) or type(b) is not type(firstInds) or \
     type(b) is not type(labels):
    raise Exception('{b,c,x,y,z,orders,firstInds,labels} must have same type')
  args = [A.obj,b.obj,c.obj,orders.obj,firstInds.obj,labels.obj, 
          x.obj,y.obj,z.obj]
  argsCtrl = [A.obj,b.obj,c.obj,orders.obj,firstInds.obj,labels.obj,
              x.obj,y.obj,z.obj,ctrl]
  if type(A) is Matrix:
    if type(b) is not Matrix: raise Exception('b must be a Matrix')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElSOCPDirect_s(*args)
      else:            lib.ElSOCPDirectX_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElSOCPDirect_d(*args)
      else:            lib.ElSOCPDirectX_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if type(b) is not DistMatrix: raise Exception('b must be a DistMatrix')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElSOCPDirectDist_s(*args)
      else:            lib.ElSOCPDirectXDist_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElSOCPDirectDist_d(*args)
      else:            lib.ElSOCPDirectXDist_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if type(b) is not Matrix: raise Exception('b must be a Matrix')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElSOCPDirectSparse_s(*args)
      else:            lib.ElSOCPDirectXSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElSOCPDirectSparse_d(*args)
      else:            lib.ElSOCPDirectXSparse_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if type(b) is not DistMultiVec: raise Exception('b must be a DistMultiVec')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElSOCPDirectDistSparse_s(*args)
      else:            lib.ElSOCPDirectXDistSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElSOCPDirectDistSparse_d(*args)
      else:            lib.ElSOCPDirectXDistSparse_d(*argsCtrl)
    else: DataExcept()
  else: TypeExcept()

# Affine conic form
# -----------------
lib.ElSOCPAffineMehrotraCtrlDefault_s.argtypes = \
lib.ElSOCPAffineMehrotraCtrlDefault_d.argtypes = \
  [c_void_p]
class SOCPAffineMehrotraCtrl_s(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",sType),("targetTol",sType),
              ("maxIts",iType),("maxStepRatio",sType),
              ("qsdCtrl",RegQSDCtrl_s),
              ("outerEquil",bType),("innerEquil",bType),
              ("scaleTwoNorm",bType),("basisSize",iType),
              ("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElSOCPAffineMehrotraCtrlDefault_s(pointer(self))
class SOCPAffineMehrotraCtrl_d(ctypes.Structure):
  _fields_ = [("primalInit",bType),("dualInit",bType),
              ("minTol",dType),("targetTol",dType),
              ("maxIts",iType),("maxStepRatio",dType),
              ("qsdCtrl",RegQSDCtrl_d),
              ("outerEquil",bType),("innerEquil",bType),
              ("scaleTwoNorm",bType),("basisSize",iType),
              ("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElSOCPAffineMehrotraCtrlDefault_d(pointer(self))

lib.ElSOCPAffineCtrlDefault_s.argtypes = \
lib.ElSOCPAffineCtrlDefault_d.argtypes = \
  [c_void_p]
class SOCPAffineCtrl_s(ctypes.Structure):
  _fields_ = [("approach",c_uint),
              ("mehrotraCtrl",SOCPAffineMehrotraCtrl_s)]
  def __init__(self):
    lib.ElSOCPAffineCtrlDefault_s(pointer(self))
class SOCPAffineCtrl_d(ctypes.Structure):
  _fields_ = [("approach",c_uint),
              ("mehrotraCtrl",SOCPAffineMehrotraCtrl_d)]
  def __init__(self):
    lib.ElSOCPAffineCtrlDefault_d(pointer(self))

lib.ElSOCPAffine_s.argtypes = \
lib.ElSOCPAffine_d.argtypes = \
lib.ElSOCPAffineDist_s.argtypes = \
lib.ElSOCPAffineDist_d.argtypes = \
lib.ElSOCPAffineSparse_s.argtypes = \
lib.ElSOCPAffineSparse_d.argtypes = \
lib.ElSOCPAffineDistSparse_s.argtypes = \
lib.ElSOCPAffineDistSparse_d.argtypes = \
  [c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,c_void_p]

lib.ElSOCPAffineX_s.argtypes = \
lib.ElSOCPAffineXDist_s.argtypes = \
lib.ElSOCPAffineXSparse_s.argtypes = \
lib.ElSOCPAffineXDistSparse_s.argtypes = \
  [c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,c_void_p,
   SOCPAffineCtrl_s]
lib.ElSOCPAffineX_d.argtypes = \
lib.ElSOCPAffineXDist_d.argtypes = \
lib.ElSOCPAffineXSparse_d.argtypes = \
lib.ElSOCPAffineXDistSparse_d.argtypes = \
  [c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,
   c_void_p,c_void_p,c_void_p,c_void_p,
   SOCPAffineCtrl_d]

# TODO
def SOCPAffine(A,G,b,c,h,orders,firstInds,labels,x,y,z,s,ctrl=None):
  if type(A) is not type(G):
    raise Exception('A and G must be of the same type')
  if orders.tag != iTag or firstInds.tag != iTag or labels.tag != iTag:
    raise Exception('cone descriptions must have integer datatypes')
  if type(b) is not type(c) or type(b) is not type(c) or \
     type(b) is not type(h) or type(b) is not type(x) or \
     type(b) is not type(y) or type(b) is not type(z) or \
     type(b) is not type(s) or type(b) is not type(orders) or \
     type(b) is not type(firstInds) or type(b) is not type(labels):
    raise Exception('vectors must be of the same type')
  args = [A.obj,G.obj,b.obj,c.obj,h.obj,orders.obj,firstInds.obj,labels.obj,
          x.obj,y.obj,z.obj,s.obj]
  argsCtrl = [A.obj,G.obj,b.obj,c.obj,h.obj,orders.obj,firstInds.obj,labels.obj,
              x.obj,y.obj,z.obj,s.obj,ctrl]
  if type(A) is Matrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElSOCPAffine_s(*args)
      else:            lib.ElSOCPAffineX_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElSOCPAffine_d(*args)
      else:            lib.ElSOCPAffineX_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if type(b) is not DistMatrix:
      raise Exception('b must be a DistMatrix')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElSOCPAffineDist_s(*args)
      else:            lib.ElSOCPAffineXDist_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElSOCPAffineDist_d(*args)
      else:            lib.ElSOCPAffineXDist_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElSOCPAffineSparse_s(*args)
      else:            lib.ElSOCPAffineXSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElSOCPAffineSparse_d(*args)
      else:            lib.ElSOCPAffineXSparse_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if type(b) is not DistMultiVec:
      raise Exception('b must be a DistMultiVec')
    if   A.tag == sTag: 
      if ctrl == None: lib.ElSOCPAffineDistSparse_s(*args)
      else:            lib.ElSOCPAffineXDistSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElSOCPAffineDistSparse_d(*args)
      else:            lib.ElSOCPAffineXDistSparse_d(*argsCtrl)
    else: DataExcept()
  else: TypeExcept()
