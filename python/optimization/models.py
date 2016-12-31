#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *
from solvers import *

from ctypes import CFUNCTYPE

# Basis pursuit
# =============
lib.ElBPADMMCtrlDefault_s.argtypes = \
lib.ElBPADMMCtrlDefault_d.argtypes = \
  [c_void_p]
class BPADMMCtrl_s(ctypes.Structure):
  _fields_ = [("rho",sType),("alpha",sType),("maxIter",iType),
              ("absTol",sType),("relTol",sType),
              ("usePinv",bType),("pinvTol",sType),("progress",bType)]
  def __init__(self):
    lib.ElBPADMMCtrlDefault_s(pointer(self))
class BPADMMCtrl_d(ctypes.Structure):
  _fields_ = [("rho",dType),("alpha",dType),("maxIter",iType),
              ("absTol",dType),("relTol",dType),
              ("usePinv",bType),("pinvTol",dType),("progress",bType)]
  def __init__(self):
    lib.ElBPADMMCtrlDefault_d(pointer(self))

lib.ElBPCtrlDefault_s.argtypes = \
lib.ElBPCtrlDefault_d.argtypes = \
  [c_void_p,bType]
class BPCtrl_s(ctypes.Structure):
  _fields_ = [("useIPM",bType),("useSOCP",bType),
              ("admmCtrl",BPADMMCtrl_s),
              ("lpIPMCtrl",LPDirectCtrl_s),
              ("socpIPMCtrl",SOCPDirectCtrl_s)]
  def __init__(self,isSparse):
    lib.ElBPCtrlDefault_s(pointer(self),isSparse)
class BPCtrl_d(ctypes.Structure):
  _fields_ = [("useIPM",bType),("useSOCP",bType),
              ("admmCtrl",BPADMMCtrl_d),
              ("lpIPMCtrl",LPDirectCtrl_d),
              ("socpIPMCtrl",SOCPDirectCtrl_d)]
  def __init__(self,isSparse):
    lib.ElBPCtrlDefault_d(pointer(self),isSparse)

lib.ElBPCtrlDefault_c.argtypes = \
lib.ElBPCtrlDefault_z.argtypes = \
  [c_void_p]
class BPCtrl_c(ctypes.Structure):
  _fields_ = [("ipmCtrl",SOCPDirectCtrl_s)]
  def __init__(self):
    lib.ElBPCtrlDefault_c(pointer(self))
class BPCtrl_z(ctypes.Structure):
  _fields_ = [("ipmCtrl",SOCPDirectCtrl_d)]
  def __init__(self):
    lib.ElBPCtrlDefault_z(pointer(self))

lib.ElBP_s.argtypes = \
lib.ElBP_d.argtypes = \
lib.ElBP_c.argtypes = \
lib.ElBP_z.argtypes = \
lib.ElBPDist_s.argtypes = \
lib.ElBPDist_d.argtypes = \
lib.ElBPDist_c.argtypes = \
lib.ElBPDist_z.argtypes = \
lib.ElBPSparse_s.argtypes = \
lib.ElBPSparse_d.argtypes = \
lib.ElBPSparse_c.argtypes = \
lib.ElBPSparse_z.argtypes = \
lib.ElBPDistSparse_s.argtypes = \
lib.ElBPDistSparse_d.argtypes = \
lib.ElBPDistSparse_c.argtypes = \
lib.ElBPDistSparse_z.argtypes = \
  [c_void_p,c_void_p,c_void_p]

lib.ElBPX_s.argtypes = \
lib.ElBPXDist_s.argtypes = \
lib.ElBPXSparse_s.argtypes = \
lib.ElBPXDistSparse_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,
   BPCtrl_s]
lib.ElBPX_d.argtypes = \
lib.ElBPXDist_d.argtypes = \
lib.ElBPXSparse_d.argtypes = \
lib.ElBPXDistSparse_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,
   BPCtrl_d]
lib.ElBPX_c.argtypes = \
lib.ElBPXDist_c.argtypes = \
lib.ElBPXSparse_c.argtypes = \
lib.ElBPXDistSparse_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,
   BPCtrl_c]
lib.ElBPX_z.argtypes = \
lib.ElBPXDist_z.argtypes = \
lib.ElBPXSparse_z.argtypes = \
lib.ElBPXDistSparse_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,
   BPCtrl_z]

def BP(A,b,ctrl=None):
  if A.tag != b.tag:
    raise Exception('Datatypes of A and b must match')
  if type(A) is Matrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    x = Matrix(A.tag)
    args = [A.obj,b.obj,x.obj]
    argsCtrl = [A.obj,b.obj,x.obj,ctrl] 
    if   A.tag == sTag: 
      if ctrl == None: lib.ElBP_s(*args)
      else:            lib.ElBPX_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElBP_d(*args)
      else:            lib.ElBPX_d(*argsCtrl)
    elif A.tag == cTag:
      if ctrl == None: lib.ElBP_c(*args)
      else:            lib.ElBPX_c(*argsCtrl)
    elif A.tag == zTag:
      if ctrl == None: lib.ElBP_z(*args)
      else:            lib.ElBPX_z(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistMatrix:
    if type(b) is not DistMatrix:
      raise Exception('b must be a DistMatrix')
    x = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,b.obj,x.obj]
    argsCtrl = [A.obj,b.obj,x.obj,ctrl] 
    if   A.tag == sTag: 
      if ctrl == None: lib.ElBPDist_s(*args)
      else:            lib.ElBPXDist_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElBPDist_d(*args)
      else:            lib.ElBPXDist_d(*argsCtrl)
    elif A.tag == cTag: 
      if ctrl == None: lib.ElBPDist_c(*args)
      else:            lib.ElBPXDist_c(*argsCtrl)
    elif A.tag == zTag: 
      if ctrl == None: lib.ElBPDist_z(*args)
      else:            lib.ElBPXDist_z(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is SparseMatrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    x = Matrix(A.tag)
    args = [A.obj,b.obj,x.obj]
    argsCtrl = [A.obj,b.obj,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElBPSparse_s(*args)
      else:            lib.ElBPXSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElBPSparse_d(*args)
      else:            lib.ElBPXSparse_d(*argsCtrl)
    elif A.tag == cTag: 
      if ctrl == None: lib.ElBPSparse_c(*args)
      else:            lib.ElBPXSparse_c(*argsCtrl)
    elif A.tag == zTag: 
      if ctrl == None: lib.ElBPSparse_z(*args)
      else:            lib.ElBPXSparse_z(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistSparseMatrix:
    if type(b) is not DistMultiVec:
      raise Exception('b must be a DistMultiVec')
    x = DistMultiVec(A.tag,A.Grid())
    args = [A.obj,b.obj,x.obj]
    argsCtrl = [A.obj,b.obj,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElBPDistSparse_s(*args)
      else:            lib.ElBPXDistSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElBPDistSparse_d(*args)
      else:            lib.ElBPXDistSparse_d(*argsCtrl)
    elif A.tag == cTag: 
      if ctrl == None: lib.ElBPDistSparse_c(*args)
      else:            lib.ElBPXDistSparse_c(*argsCtrl)
    elif A.tag == zTag: 
      if ctrl == None: lib.ElBPDistSparse_z(*args)
      else:            lib.ElBPXDistSparse_z(*argsCtrl)
    else: DataExcept()
    return x
  else: TypeExcept()

# Chebyshev point
# ===============
lib.ElCP_s.argtypes = \
lib.ElCP_d.argtypes = \
lib.ElCPDist_s.argtypes = \
lib.ElCPDist_d.argtypes = \
lib.ElCPSparse_s.argtypes = \
lib.ElCPSparse_d.argtypes = \
lib.ElCPDistSparse_s.argtypes = \
lib.ElCPDistSparse_d.argtypes = \
  [c_void_p,c_void_p,c_void_p]

lib.ElCPX_s.argtypes = \
lib.ElCPXDist_s.argtypes = \
lib.ElCPXSparse_s.argtypes = \
lib.ElCPXDistSparse_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,
   LPAffineCtrl_s]
lib.ElCPX_d.argtypes = \
lib.ElCPXDist_d.argtypes = \
lib.ElCPXSparse_d.argtypes = \
lib.ElCPXDistSparse_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,
   LPAffineCtrl_d]

def CP(A,b,ctrl=None):
  if A.tag != b.tag:
    raise Exception('Datatypes of A and b must match')
  if type(A) is Matrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    x = Matrix(A.tag)
    args = [A.obj,b.obj,x.obj]
    argsCtrl = [A.obj,b.obj,x.obj,ctrl] 
    if   A.tag == sTag: 
      if ctrl == None: lib.ElCP_s(*args)
      else:            lib.ElCPX_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElCP_d(*args)
      else:            lib.ElCPX_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistMatrix:
    if type(b) is not DistMatrix:
      raise Exception('b must be a DistMatrix')
    x = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,b.obj,x.obj]
    argsCtrl = [A.obj,b.obj,x.obj,ctrl] 
    if   A.tag == sTag: 
      if ctrl == None: lib.ElCPDist_s(*args)
      else:            lib.ElCPXDist_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElCPDist_d(*args)
      else:            lib.ElCPXDist_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is SparseMatrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    x = Matrix(A.tag)
    args = [A.obj,b.obj,x.obj]
    argsCtrl = [A.obj,b.obj,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElCPSparse_s(*args)
      else:            lib.ElCPXSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElCPSparse_d(*args)
      else:            lib.ElCPXSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistSparseMatrix:
    if type(b) is not DistMultiVec:
      raise Exception('b must be a DistMultiVec')
    x = DistMultiVec(A.tag,A.Grid())
    args = [A.obj,b.obj,x.obj]
    argsCtrl = [A.obj,b.obj,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElCPDistSparse_s(*args)
      else:            lib.ElCPXDistSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElCPDistSparse_d(*args)
      else:            lib.ElCPXDistSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  else: TypeExcept()

# Dantzig selector
# ================
lib.ElDS_s.argtypes = \
lib.ElDSDist_s.argtypes = \
lib.ElDSSparse_s.argtypes = \
lib.ElDSDistSparse_s.argtypes = \
  [c_void_p,c_void_p,sType,c_void_p]
lib.ElDS_d.argtypes = \
lib.ElDSDist_d.argtypes = \
lib.ElDSSparse_d.argtypes = \
lib.ElDSDistSparse_d.argtypes = \
  [c_void_p,c_void_p,dType,c_void_p]

lib.ElDSX_s.argtypes = \
lib.ElDSXDist_s.argtypes = \
lib.ElDSXSparse_s.argtypes = \
lib.ElDSXDistSparse_s.argtypes = \
  [c_void_p,c_void_p,sType,c_void_p,
   LPAffineCtrl_s]
lib.ElDSX_d.argtypes = \
lib.ElDSXDist_d.argtypes = \
lib.ElDSXSparse_d.argtypes = \
lib.ElDSXDistSparse_d.argtypes = \
  [c_void_p,c_void_p,dType,c_void_p,
   LPAffineCtrl_d]

def DS(A,b,lambdaPre,ctrl=None):
  if A.tag != b.tag:
    raise Exception('Datatypes of A and b must match')
  lambd = TagToType(A.tag)(lambdaPre)
  if type(A) is Matrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    x = Matrix(A.tag)
    args = [A.obj,b.obj,lambd,x.obj]
    argsCtrl = [A.obj,b.obj,lambd,x.obj,ctrl] 
    if   A.tag == sTag: 
      if ctrl == None: lib.ElDS_s(*args)
      else:            lib.ElDSX_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElDS_d(*args)
      else:            lib.ElDSX_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistMatrix:
    if type(b) is not DistMatrix:
      raise Exception('b must be a DistMatrix')
    x = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,b.obj,lambd,x.obj]
    argsCtrl = [A.obj,b.obj,lambd,x.obj,ctrl] 
    if   A.tag == sTag: 
      if ctrl == None: lib.ElDSDist_s(*args)
      else:            lib.ElDSXDist_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElDSDist_d(*args)
      else:            lib.ElDSXDist_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is SparseMatrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    x = Matrix(A.tag)
    args = [A.obj,b.obj,lambd,x.obj]
    argsCtrl = [A.obj,b.obj,lambd,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElDSSparse_s(*args)
      else:            lib.ElDSXSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElDSSparse_d(*args)
      else:            lib.ElDSXSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistSparseMatrix:
    if type(b) is not DistMultiVec:
      raise Exception('b must be a DistMultiVec')
    x = DistMultiVec(A.tag,A.Grid())
    args = [A.obj,b.obj,lambd,x.obj]
    argsCtrl = [A.obj,b.obj,lambd,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElDSDistSparse_s(*args)
      else:            lib.ElDSXDistSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElDSDistSparse_d(*args)
      else:            lib.ElDSXDistSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  else: TypeExcept()

# Least Absolute Value regression
# ===============================
lib.ElLAV_s.argtypes = \
lib.ElLAV_d.argtypes = \
lib.ElLAVDist_s.argtypes = \
lib.ElLAVDist_d.argtypes = \
lib.ElLAVSparse_s.argtypes = \
lib.ElLAVSparse_d.argtypes = \
lib.ElLAVDistSparse_s.argtypes = \
lib.ElLAVDistSparse_d.argtypes = \
  [c_void_p,c_void_p,c_void_p]

lib.ElLAVX_s.argtypes = \
lib.ElLAVXDist_s.argtypes = \
lib.ElLAVXSparse_s.argtypes = \
lib.ElLAVXDistSparse_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,
   LPAffineCtrl_s]
lib.ElLAVX_d.argtypes = \
lib.ElLAVXDist_d.argtypes = \
lib.ElLAVXSparse_d.argtypes = \
lib.ElLAVXDistSparse_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,
   LPAffineCtrl_d]

def LAV(A,b,ctrl=None):
  if A.tag != b.tag:
    raise Exception('Datatypes of A and b must match')
  if type(A) is Matrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    x = Matrix(A.tag)
    args = [A.obj,b.obj,x.obj]
    argsCtrl = [A.obj,b.obj,x.obj,ctrl] 
    if   A.tag == sTag: 
      if ctrl == None: lib.ElLAV_s(*args)
      else:            lib.ElLAVX_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElLAV_d(*args)
      else:            lib.ElLAVX_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistMatrix:
    if type(b) is not DistMatrix:
      raise Exception('b must be a DistMatrix')
    x = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,b.obj,x.obj]
    argsCtrl = [A.obj,b.obj,x.obj,ctrl] 
    if   A.tag == sTag: 
      if ctrl == None: lib.ElLAVDist_s(*args)
      else:            lib.ElLAVXDist_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElLAVDist_d(*args)
      else:            lib.ElLAVXDist_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is SparseMatrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    x = Matrix(A.tag)
    args = [A.obj,b.obj,x.obj]
    argsCtrl = [A.obj,b.obj,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElLAVSparse_s(*args)
      else:            lib.ElLAVXSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElLAVSparse_d(*args)
      else:            lib.ElLAVXSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistSparseMatrix:
    if type(b) is not DistMultiVec:
      raise Exception('b must be a DistMultiVec')
    x = DistMultiVec(A.tag,A.Grid())
    args = [A.obj,b.obj,x.obj]
    argsCtrl = [A.obj,b.obj,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElLAVDistSparse_s(*args)
      else:            lib.ElLAVXDistSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElLAVDistSparse_d(*args)
      else:            lib.ElLAVXDistSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  else: TypeExcept()

# Robust least squares
# ====================
lib.ElRLS_s.argtypes = \
lib.ElRLSDist_s.argtypes = \
lib.ElRLSSparse_s.argtypes = \
lib.ElRLSDistSparse_s.argtypes = \
  [c_void_p,c_void_p,sType,c_void_p]
lib.ElRLS_d.argtypes = \
lib.ElRLSDist_d.argtypes = \
lib.ElRLSSparse_d.argtypes = \
lib.ElRLSDistSparse_d.argtypes = \
  [c_void_p,c_void_p,dType,c_void_p]

lib.ElRLSX_s.argtypes = \
lib.ElRLSXDist_s.argtypes = \
lib.ElRLSXSparse_s.argtypes = \
lib.ElRLSXDistSparse_s.argtypes = \
  [c_void_p,c_void_p,sType,c_void_p,SOCPAffineCtrl_s]
lib.ElRLSX_d.argtypes = \
lib.ElRLSXDist_d.argtypes = \
lib.ElRLSXSparse_d.argtypes = \
lib.ElRLSXDistSparse_d.argtypes = \
  [c_void_p,c_void_p,dType,c_void_p,SOCPAffineCtrl_d]

def RLS(A,b,rho,ctrl=None):
  if A.tag != b.tag:
    raise Exception('Datatypes of A and b must match')
  if type(A) is Matrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    x = Matrix(A.tag)
    args = [A.obj,b.obj,rho,x.obj]
    argsCtrl = [A.obj,b.obj,rho,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl==None: lib.ElRLS_s(*args)
      else:          lib.ElRLSX_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl==None: lib.ElRLS_d(*args)
      else:          lib.ElRLSX_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistMatrix:
    if type(b) is not DistMatrix:
      raise Exception('b must be a DistMatrix')
    x = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,b.obj,rho,x.obj]
    argsCtrl = [A.obj,b.obj,rho,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl==None: lib.ElRLSDist_s(*args)
      else:          lib.ElRLSXDist_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl==None: lib.ElRLSDist_d(*args)
      else:          lib.ElRLSXDist_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is SparseMatrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    x = SparseMatrix(A.tag)
    args = [A.obj,b.obj,rho,x.obj]
    argsCtrl = [A.obj,b.obj,rho,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl==None: lib.ElRLSSparse_s(*args)
      else:          lib.ElRLSXSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl==None: lib.ElRLSSparse_d(*args)
      else:          lib.ElRLSXSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistSparseMatrix:
    if type(b) is not DistMultiVec:
      raise Exception('b must be a DistMultiVec')
    x = DistMultiVec(A.tag,A.Grid())
    args = [A.obj,b.obj,rho,x.obj]
    argsCtrl = [A.obj,b.obj,rho,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl==None: lib.ElRLSDistSparse_s(*args)
      else:          lib.ElRLSXDistSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl==None: lib.ElRLSDistSparse_d(*args)
      else:          lib.ElRLSXDistSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  else: TypeExcept()

# Robust non-negative least squares
# =================================
lib.ElRNNLS_s.argtypes = \
lib.ElRNNLSDist_s.argtypes = \
lib.ElRNNLSSparse_s.argtypes = \
lib.ElRNNLSDistSparse_s.argtypes = \
  [c_void_p,c_void_p,sType,c_void_p]
lib.ElRNNLS_d.argtypes = \
lib.ElRNNLSDist_d.argtypes = \
lib.ElRNNLSSparse_d.argtypes = \
lib.ElRNNLSDistSparse_d.argtypes = \
  [c_void_p,c_void_p,dType,c_void_p]

lib.ElRNNLSX_s.argtypes = \
lib.ElRNNLSXDist_s.argtypes = \
lib.ElRNNLSXSparse_s.argtypes = \
lib.ElRNNLSXDistSparse_s.argtypes = \
  [c_void_p,c_void_p,sType,c_void_p,SOCPAffineCtrl_s]
lib.ElRNNLSX_d.argtypes = \
lib.ElRNNLSXDist_d.argtypes = \
lib.ElRNNLSXSparse_d.argtypes = \
lib.ElRNNLSXDistSparse_d.argtypes = \
  [c_void_p,c_void_p,dType,c_void_p,SOCPAffineCtrl_d]

def RNNLS(A,b,rho,ctrl=None):
  if A.tag != b.tag:
    raise Exception('Datatypes of A and b must match')
  if type(A) is Matrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    x = Matrix(A.tag)
    args = [A.obj,b.obj,rho,x.obj]
    argsCtrl = [A.obj,b.obj,rho,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl==None: lib.ElRNNLS_s(*args)
      else:          lib.ElRNNLSX_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl==None: lib.ElRNNLS_d(*args)
      else:          lib.ElRNNLSX_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistMatrix:
    if type(b) is not DistMatrix:
      raise Exception('b must be a DistMatrix')
    x = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,b.obj,rho,x.obj]
    argsCtrl = [A.obj,b.obj,rho,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl==None: lib.ElRNNLSDist_s(*args)
      else:          lib.ElRNNLSXDist_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl==None: lib.ElRNNLSDist_d(*args)
      else:          lib.ElRNNLSXDist_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is SparseMatrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    x = SparseMatrix(A.tag)
    args = [A.obj,b.obj,rho,x.obj]
    argsCtrl = [A.obj,b.obj,rho,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl==None: lib.ElRNNLSSparse_s(*args)
      else:          lib.ElRNNLSXSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl==None: lib.ElRNNLSSparse_d(*args)
      else:          lib.ElRNNLSXSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistSparseMatrix:
    if type(b) is not DistMultiVec:
      raise Exception('b must be a DistMultiVec')
    x = DistMultiVec(A.tag,A.Grid())
    args = [A.obj,b.obj,rho,x.obj]
    argsCtrl = [A.obj,b.obj,rho,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl==None: lib.ElRNNLSDistSparse_s(*args)
      else:          lib.ElRNNLSXDistSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl==None: lib.ElRNNLSDistSparse_d(*args)
      else:          lib.ElRNNLSXDistSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  else: TypeExcept()

# Non-negative least squares
# ==========================
lib.ElNNLSCtrlDefault_s.argtypes = \
lib.ElNNLSCtrlDefault_d.argtypes = \
  [c_void_p]
(NNLS_ADMM,NNLS_QP,NNLS_SOCP)=(0,1,2)
class NNLSCtrl_s(ctypes.Structure):
  _fields_ = [("approach",c_uint),
              ("admmCtrl",ADMMCtrl_s),
              ("qpCtrl",QPDirectCtrl_s),
              ("socpCtrl",SOCPAffineCtrl_s)]
  def __init__(self):
    lib.ElNNLSCtrlDefault_s(pointer(self))
class NNLSCtrl_d(ctypes.Structure):
  _fields_ = [("approach",c_uint),
              ("admmCtrl",ADMMCtrl_d),
              ("qpCtrl",QPDirectCtrl_d),
              ("socpCtrl",SOCPAffineCtrl_d)]
  def __init__(self):
    lib.ElNNLSCtrlDefault_d(pointer(self))

lib.ElNNLS_s.argtypes = \
lib.ElNNLS_d.argtypes = \
lib.ElNNLSDist_s.argtypes = \
lib.ElNNLSDist_d.argtypes = \
lib.ElNNLSSparse_s.argtypes = \
lib.ElNNLSSparse_d.argtypes = \
lib.ElNNLSDistSparse_s.argtypes = \
lib.ElNNLSDistSparse_d.argtypes = \
  [c_void_p,c_void_p,c_void_p]
lib.ElNNLSX_s.argtypes = \
lib.ElNNLSXDist_s.argtypes = \
lib.ElNNLSXSparse_s.argtypes = \
lib.ElNNLSXDistSparse_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,NNLSCtrl_s]
lib.ElNNLSX_d.argtypes = \
lib.ElNNLSXDist_d.argtypes = \
lib.ElNNLSXSparse_d.argtypes = \
lib.ElNNLSXDistSparse_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,NNLSCtrl_d]

def NNLS(A,b,ctrl=None):
  if A.tag != b.tag:
    raise Exception('Datatypes of A and b must match')
  if type(A) is Matrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    x = Matrix(A.tag)
    args = [A.obj,b.obj,x.obj]
    argsCtrl = [A.obj,b.obj,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl==None: lib.ElNNLS_s(*args)
      else:          lib.ElNNLSX_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl==None: lib.ElNNLS_d(*args)
      else:          lib.ElNNLSX_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistMatrix:
    if type(b) is not DistMatrix:
      raise Exception('b must be a DistMatrix')
    x = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,b.obj,x.obj]
    argsCtrl = [A.obj,b.obj,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl==None: lib.ElNNLSDist_s(*args)
      else:          lib.ElNNLSXDist_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl==None: lib.ElNNLSDist_d(*args)
      else:          lib.ElNNLSXDist_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is SparseMatrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    x = SparseMatrix(A.tag)
    args = [A.obj,b.obj,x.obj]
    argsCtrl = [A.obj,b.obj,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl==None: lib.ElNNLSSparse_s(*args)
      else:          lib.ElNNLSXSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl==None: lib.ElNNLSSparse_d(*args)
      else:          lib.ElNNLSXSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistSparseMatrix:
    if type(b) is not DistMultiVec:
      raise Exception('b must be a DistMultiVec')
    x = DistMultiVec(A.tag,A.Grid())
    args = [A.obj,b.obj,x.obj]
    argsCtrl = [A.obj,b.obj,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl==None: lib.ElNNLSDistSparse_s(*args)
      else:          lib.ElNNLSXDistSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl==None: lib.ElNNLSDistSparse_d(*args)
      else:          lib.ElNNLSXDistSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  else: TypeExcept()

# Non-negative matrix factorization
# =================================
lib.ElNMFCtrlDefault_s.argtypes = \
lib.ElNMFCtrlDefault_d.argtypes = \
  [c_void_p]
class NMFCtrl_s(ctypes.Structure):
  _fields_ = [("nnlsCtrl",NNLSCtrl_s),("maxIter",iType)]
  def __init__(self):
    lib.ElNMFCtrlDefault_s(pointer(self))
class NMFCtrl_d(ctypes.Structure):
  _fields_ = [("nnlsCtrl",NNLSCtrl_d),("maxIter",iType)]
  def __init__(self):
    lib.ElNMFCtrlDefault_d(pointer(self))

lib.ElNMF_s.argtypes = \
lib.ElNMF_d.argtypes = \
lib.ElNMFDist_s.argtypes = \
lib.ElNMFDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p]
lib.ElNMFX_s.argtypes = \
lib.ElNMFXDist_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,NMFCtrl_s]
lib.ElNMFX_d.argtypes = \
lib.ElNMFXDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,NMFCtrl_d]

def NMF(A,ctrl=None):
  args = [A.obj,X.obj,Y.obj]
  argsCtrl = [A.obj,X.obj,Y.obj,ctrl]
  if type(A) is Matrix:
    if   A.tag == sTag: 
      if ctrl==None: lib.ElNMF_s(*args)
      else:          lib.ElNMFX_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl==None: lib.ElNMF_d(*args)
      else:          lib.ElNMFX_d(*argsCtrl)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      if ctrl==None: lib.ElNMFDist_s(*args)
      else:          lib.ElNMFXDist_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl==None: lib.ElNMFDist_d(*args)
      else:          lib.ElNMFXDist_d(*argsCtrl)
    else: DataExcept()
  else: TypeExcept()

# Basis pursuit denoising
# =======================
lib.ElBPDNADMMCtrlDefault_s.argtypes = \
lib.ElBPDNADMMCtrlDefault_d.argtypes = \
  [c_void_p]
class BPDNADMMCtrl_s(ctypes.Structure):
  _fields_ = [("rho",sType),("alpha",sType),("maxIter",iType),
              ("absTol",sType),("relTol",sType),
              ("inv",bType),("progress",bType)]
  def __init__(self):
    lib.ElBPDNADMMCtrlDefault_s(pointer(self))
class BPDNADMMCtrl_d(ctypes.Structure):
  _fields_ = [("rho",dType),("alpha",dType),("maxIter",iType),
              ("absTol",dType),("relTol",dType),
              ("inv",bType),("progress",bType)]
  def __init__(self):
    lib.ElBPDNADMMCtrlDefault_d(pointer(self))

lib.ElBPDNCtrlDefault_s.argtypes = \
lib.ElBPDNCtrlDefault_d.argtypes = \
  [c_void_p]
class BPDNCtrl_s(ctypes.Structure):
  _fields_ = [("useIPM",bType),
              ("admmCtrl",BPDNADMMCtrl_s),("ipmCtrl",QPAffineCtrl_s)]
  def __init__(self):
    lib.ElBPDNCtrlDefault_s(pointer(self))
class BPDNCtrl_d(ctypes.Structure):
  _fields_ = [("useIPM",bType),
              ("admmCtrl",BPDNADMMCtrl_d),("ipmCtrl",QPAffineCtrl_d)]
  def __init__(self):
    lib.ElBPDNCtrlDefault_d(pointer(self))

lib.ElBPDN_s.argtypes = \
lib.ElBPDNDist_s.argtypes = \
lib.ElBPDNSparse_s.argtypes = \
lib.ElBPDNDistSparse_s.argtypes = \
  [c_void_p,c_void_p,sType,c_void_p]
lib.ElBPDN_d.argtypes = \
lib.ElBPDNDist_d.argtypes = \
lib.ElBPDNSparse_d.argtypes = \
lib.ElBPDNDistSparse_d.argtypes = \
  [c_void_p,c_void_p,dType,c_void_p]

lib.ElBPDNX_s.argtypes = \
lib.ElBPDNXDist_s.argtypes = \
lib.ElBPDNXSparse_s.argtypes = \
lib.ElBPDNXDistSparse_s.argtypes = \
  [c_void_p,c_void_p,sType,c_void_p,BPDNCtrl_s]
lib.ElBPDNX_d.argtypes = \
lib.ElBPDNXDist_d.argtypes = \
lib.ElBPDNXSparse_d.argtypes = \
lib.ElBPDNXDistSparse_d.argtypes = \
  [c_void_p,c_void_p,dType,c_void_p,BPDNCtrl_d]

def BPDN(A,b,lambdPre,ctrl=None):
  if A.tag != b.tag:
    raise Exception('Datatypes of A and b must match')
  lambd = TagToType(A.tag)(lambdPre)
  if type(A) is Matrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    x = Matrix(A.tag)
    args = [A.obj,b.obj,lambd,x.obj]
    argsCtrl = [A.obj,b.obj,lambd,x.obj,ctrl] 
    if   A.tag == sTag: 
      if ctrl == None: lib.ElBPDN_s(*args)
      else:            lib.ElBPDNX_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElBPDN_d(*args)
      else:            lib.ElBPDNX_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistMatrix:
    if type(b) is not DistMatrix:
      raise Exception('b must be a DistMatrix')
    x = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,b.obj,lambd,x.obj]
    argsCtrl = [A.obj,b.obj,lambd,x.obj,ctrl] 
    if   A.tag == sTag: 
      if ctrl == None: lib.ElBPDNDist_s(*args)
      else:            lib.ElBPDNXDist_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElBPDNDist_d(*args)
      else:            lib.ElBPDNXDist_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is SparseMatrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    x = Matrix(A.tag)
    args = [A.obj,b.obj,lambd,x.obj]
    argsCtrl = [A.obj,b.obj,lambd,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElBPDNSparse_s(*args)
      else:            lib.ElBPDNXSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElBPDNSparse_d(*args)
      else:            lib.ElBPDNXSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistSparseMatrix:
    if type(b) is not DistMultiVec:
      raise Exception('b must be a DistMultiVec')
    x = DistMultiVec(A.tag,A.Grid())
    args = [A.obj,b.obj,lambd,x.obj]
    argsCtrl = [A.obj,b.obj,lambd,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElBPDNDistSparse_s(*args)
      else:            lib.ElBPDNXDistSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElBPDNDistSparse_d(*args)
      else:            lib.ElBPDNXDistSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  else: TypeExcept()

# Elastic net
# ===========
lib.ElEN_s.argtypes = \
lib.ElENDist_s.argtypes = \
lib.ElENSparse_s.argtypes = \
lib.ElENDistSparse_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,c_void_p]
lib.ElEN_d.argtypes = \
lib.ElENDist_d.argtypes = \
lib.ElENSparse_d.argtypes = \
lib.ElENDistSparse_d.argtypes = \
  [c_void_p,c_void_p,dType,dType,c_void_p]

lib.ElENX_s.argtypes = \
lib.ElENXDist_s.argtypes = \
lib.ElENXSparse_s.argtypes = \
lib.ElENXDistSparse_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,c_void_p,
   QPAffineCtrl_s]
lib.ElENX_d.argtypes = \
lib.ElENXDist_d.argtypes = \
lib.ElENXSparse_d.argtypes = \
lib.ElENXDistSparse_d.argtypes = \
  [c_void_p,c_void_p,dType,dType,c_void_p,
   QPAffineCtrl_d]

def EN(A,b,lambda1Pre,lambda2Pre,ctrl=None):
  if A.tag != b.tag:
    raise Exception('Datatypes of A and b must match')
  lambda1 = TagToType(A.tag)(lambda1Pre)
  lambda2 = TagToType(A.tag)(lambda2Pre)
  if type(A) is Matrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    x = Matrix(A.tag)
    args = [A.obj,b.obj,lambda1,lambda2,x.obj]
    argsCtrl = [A.obj,b.obj,lambda1,lambda2,x.obj,ctrl] 
    if   A.tag == sTag: 
      if ctrl == None: lib.ElEN_s(*args)
      else:            lib.ElENX_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElEN_d(*args)
      else:            lib.ElENX_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistMatrix:
    if type(b) is not DistMatrix:
      raise Exception('b must be a DistMatrix')
    x = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,b.obj,lambda1,lambda2,x.obj]
    argsCtrl = [A.obj,b.obj,lambda1,lambda2,x.obj,ctrl] 
    if   A.tag == sTag: 
      if ctrl == None: lib.ElENDist_s(*args)
      else:            lib.ElENXDist_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElENDist_d(*args)
      else:            lib.ElENXDist_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is SparseMatrix:
    if type(b) is not Matrix:
      raise Exception('b must be a Matrix')
    x = Matrix(A.tag)
    args = [A.obj,b.obj,lambda1,lambda2,x.obj]
    argsCtrl = [A.obj,b.obj,lambda1,lambda2,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElENSparse_s(*args)
      else:            lib.ElENXSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElENSparse_d(*args)
      else:            lib.ElENXSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistSparseMatrix:
    if type(b) is not DistMultiVec:
      raise Exception('b must be a DistMultiVec')
    x = DistMultiVec(A.tag,A.Grid())
    args = [A.obj,b.obj,lambda1,lambda2,x.obj]
    argsCtrl = [A.obj,b.obj,lambda1,lambda2,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElENDistSparse_s(*args)
      else:            lib.ElENXDistSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElENDistSparse_d(*args)
      else:            lib.ElENXDistSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  else: TypeExcept()

# Robust Principal Component Analysis
# ===================================

lib.ElRPCACtrlDefault_s.argtypes = \
lib.ElRPCACtrlDefault_d.argtypes = \
  [c_void_p]
class RPCACtrl_s(ctypes.Structure):
  _fields_ = [("useALM",bType),("usePivQR",bType),("progress",bType),
              ("numPivSteps",iType),("maxIts",iType),
              ("tau",sType),("beta",sType),("rho",sType),("tol",sType)]
  def __init__(self):
    lib.ElRPCACtrlDefault_s(pointer(self))
class RPCACtrl_d(ctypes.Structure):
  _fields_ = [("useALM",bType),("usePivQR",bType),("progress",bType),
              ("numPivSteps",iType),("maxIts",iType),
              ("tau",dType),("beta",dType),("rho",dType),("tol",dType)]
  def __init__(self):
    lib.ElRPCACtrlDefault_d(pointer(self))

lib.ElRPCA_s.argtypes = \
lib.ElRPCA_d.argtypes = \
lib.ElRPCA_c.argtypes = \
lib.ElRPCA_z.argtypes = \
lib.ElRPCADist_s.argtypes = \
lib.ElRPCADist_d.argtypes = \
lib.ElRPCADist_c.argtypes = \
lib.ElRPCADist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p]

lib.ElRPCAX_s.argtypes = \
lib.ElRPCAX_c.argtypes = \
lib.ElRPCAXDist_s.argtypes = \
lib.ElRPCAXDist_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,RPCACtrl_s]

lib.ElRPCAX_d.argtypes = \
lib.ElRPCAX_z.argtypes = \
lib.ElRPCAXDist_d.argtypes = \
lib.ElRPCAXDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,RPCACtrl_d]

def RPCA(M,ctrl=None):
  if type(M) is Matrix:
    L = Matrix(M.tag)
    S = Matrix(M.tag)
    args = [M.obj,L.obj,S.obj]
    argsCtrl = [M.obj,L.obj,S.obj,ctrl]
    if   M.tag == sTag: 
      if ctrl==None: lib.ElRPCA_s(*args)
      else:          lib.ElRPCAX_s(*argsCtrl)
    elif M.tag == dTag: 
      if ctrl==None: lib.ElRPCA_d(*args)
      else:          lib.ElRPCAX_d(*argsCtrl)
    elif M.tag == cTag: 
      if ctrl==None: lib.ElRPCA_c(*args)
      else:          lib.ElRPCAX_c(*argsCtrl)
    elif M.tag == zTag: 
      if ctrl==None: lib.ElRPCA_z(*args)
      else:          lib.ElRPCAX_z(*argsCtrl)
    return L, S
  elif type(M) is DistMatrix:
    L = DistMatrix(M.tag,MC,MR,M.Grid())
    S = DistMatrix(M.tag,MC,MR,M.Grid())
    args = [M.obj,L.obj,S.obj]
    argsCtrl = [M.obj,L.obj,S.obj,ctrl]
    if   M.tag == sTag: 
      if ctrl==None: lib.ElRPCADist_s(*args)
      else:          lib.ElRPCAXDist_s(*argsCtrl)
    elif M.tag == dTag: 
      if ctrl==None: lib.ElRPCADist_d(*args)
      else:          lib.ElRPCAXDist_d(*argsCtrl)
    elif M.tag == cTag: 
      if ctrl==None: lib.ElRPCADist_c(*args)
      else:          lib.ElRPCAXDist_c(*argsCtrl)
    elif M.tag == zTag: 
      if ctrl==None: lib.ElRPCADist_z(*args)
      else:          lib.ElRPCAXDist_z(*argsCtrl)
    return L, S
  else: TypeExcept()

# Sparse inverse covariance selection
# ===================================

lib.ElSparseInvCovCtrlDefault_s.argtypes = \
lib.ElSparseInvCovCtrlDefault_d.argtypes = \
  [c_void_p]
class SparseInvCovCtrl_s(ctypes.Structure):
  _fields_ = [("rho",sType),("alpha",sType),
              ("maxIter",iType),
              ("absTol",sType),("relTol",sType),
              ("progress",bType)]
  def __init__(self):
    lib.ElSparseInvCovCtrlDefault_s(pointer(self))
class SparseInvCovCtrl_d(ctypes.Structure):
  _fields_ = [("rho",dType),("alpha",dType),
              ("maxIter",iType),
              ("absTol",dType),("relTol",dType),
              ("progress",bType)]
  def __init__(self):
    lib.ElSparseInvCovCtrlDefault_d(pointer(self))

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

lib.ElSparseInvCovX_s.argtypes = \
lib.ElSparseInvCovX_c.argtypes = \
lib.ElSparseInvCovXDist_s.argtypes = \
lib.ElSparseInvCovXDist_c.argtypes = \
  [c_void_p,sType,c_void_p,SparseInvCovCtrl_s,POINTER(iType)]
lib.ElSparseInvCovX_d.argtypes = \
lib.ElSparseInvCovX_z.argtypes = \
lib.ElSparseInvCovXDist_d.argtypes = \
lib.ElSparseInvCovXDist_z.argtypes = \
  [c_void_p,dType,c_void_p,SparseInvCovCtrl_d,POINTER(iType)]

def SparseInvCov(D,lambdaPre,ctrl=None):
  numIts = iType()
  lambd = TagToType(Base(D.tag))(lambdaPre) 
  if type(D) is Matrix:
    Z = Matrix(D.tag)
    args = [D.obj,lambd,Z.obj,pointer(numIts)]
    argsCtrl = [D.obj,lambd,Z.obj,ctrl,pointer(numIts)]
    if   D.tag == sTag: 
      if ctrl==None: lib.ElSparseInvCov_s(*args)
      else:          lib.ElSparseInvCovX_s(*argsCtrl)
    elif D.tag == dTag: 
      if ctrl==None: lib.ElSparseInvCov_d(*args)
      else:          lib.ElSparseInvCovX_d(*argsCtrl)
    elif D.tag == cTag:
      if ctrl==None: lib.ElSparseInvCov_c(*args)
      else:          lib.ElSparseInvCovX_c(*argsCtrl)
    elif D.tag == zTag:
      if ctrl==None: lib.ElSparseInvCov_z(*args)
      else:          lib.ElSparseInvCovX_z(*argsCtrl)
    else: DataExcept()
    return Z, numIts
  elif type(D) is DistMatrix:
    Z = DistMatrix(D.tag,MC,MR,D.Grid())
    args = [D.obj,lambd,Z.obj,pointer(numIts)]
    argsCtrl = [D.obj,lambd,Z.obj,ctrl,pointer(numIts)]
    if   D.tag == sTag: 
      if ctrl==None: lib.ElSparseInvCovDist_s(*args)
      else:          lib.ElSparseInvCovXDist_s(*argsCtrl)
    elif D.tag == dTag: 
      if ctrl==None: lib.ElSparseInvCovDist_d(*args)
      else:          lib.ElSparseInvCovXDist_d(*argsCtrl)
    elif D.tag == cTag:
      if ctrl==None: lib.ElSparseInvCovDist_c(*args)
      else:          lib.ElSparseInvCovXDist_c(*argsCtrl)
    elif D.tag == zTag:
      if ctrl==None: lib.ElSparseInvCovDist_z(*args)
      else:          lib.ElSparseInvCovXDist_z(*argsCtrl)
    else: DataExcept()
    return Z, numIts
  else: TypeExcept()

# Support Vector Machine
# ======================
lib.ElSVMCtrlDefault_s.argtypes = \
lib.ElSVMCtrlDefault_d.argtypes = \
  [c_void_p]
class SVMCtrl_s(ctypes.Structure):
  _fields_ = [("ipmCtrl",QPAffineCtrl_s)]
  def __init__(self):
    lib.ElSVMCtrlDefault_s(pointer(self))
class SVMCtrl_d(ctypes.Structure):
  _fields_ = [("ipmCtrl",QPAffineCtrl_d)]
  def __init__(self):
    lib.ElSVMCtrlDefault_d(pointer(self))

lib.ElSVM_s.argtypes = \
lib.ElSVMDist_s.argtypes = \
lib.ElSVMSparse_s.argtypes = \
lib.ElSVMDistSparse_s.argtypes = \
  [c_void_p,c_void_p,sType,c_void_p]
lib.ElSVM_d.argtypes = \
lib.ElSVMDist_d.argtypes = \
lib.ElSVMSparse_d.argtypes = \
lib.ElSVMDistSparse_d.argtypes = \
  [c_void_p,c_void_p,dType,c_void_p]

lib.ElSVMX_s.argtypes = \
lib.ElSVMXDist_s.argtypes = \
lib.ElSVMXSparse_s.argtypes = \
lib.ElSVMXDistSparse_s.argtypes = \
  [c_void_p,c_void_p,sType,c_void_p,
   SVMCtrl_s]
lib.ElSVMX_d.argtypes = \
lib.ElSVMXDist_d.argtypes = \
lib.ElSVMXSparse_d.argtypes = \
lib.ElSVMXDistSparse_d.argtypes = \
  [c_void_p,c_void_p,dType,c_void_p,
   SVMCtrl_d]

def SVM(A,d,lambdPre,ctrl=None):
  if A.tag != d.tag:
    raise Exception('Datatypes of A and d must match')
  lambd = TagToType(A.tag)(lambdPre)
  if type(A) is Matrix:
    if type(d) is not Matrix:
      raise Exception('d must be a Matrix')
    x = Matrix(A.tag)
    args = [A.obj,d.obj,lambd,x.obj]
    argsCtrl = [A.obj,d.obj,lambd,x.obj,ctrl] 
    if   A.tag == sTag: 
      if ctrl == None: lib.ElSVM_s(*args)
      else:            lib.ElSVMX_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElSVM_d(*args)
      else:            lib.ElSVMX_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistMatrix:
    if type(d) is not DistMatrix:
      raise Exception('d must be a DistMatrix')
    x = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,d.obj,lambd,x.obj]
    argsCtrl = [A.obj,d.obj,lambd,x.obj,ctrl] 
    if   A.tag == sTag: 
      if ctrl == None: lib.ElSVMDist_s(*args)
      else:            lib.ElSVMXDist_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElSVMDist_d(*args)
      else:            lib.ElSVMXDist_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is SparseMatrix:
    if type(d) is not Matrix:
      raise Exception('d must be a Matrix')
    x = Matrix(A.tag)
    args = [A.obj,d.obj,lambd,x.obj]
    argsCtrl = [A.obj,d.obj,lambd,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElSVMSparse_s(*args)
      else:            lib.ElSVMXSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElSVMSparse_d(*args)
      else:            lib.ElSVMXSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(A) is DistSparseMatrix:
    if type(d) is not DistMultiVec:
      raise Exception('d must be a DistMultiVec')
    x = DistMultiVec(A.tag,A.Grid())
    args = [A.obj,d.obj,lambd,x.obj]
    argsCtrl = [A.obj,d.obj,lambd,x.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElSVMDistSparse_s(*args)
      else:            lib.ElSVMXDistSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl == None: lib.ElSVMDistSparse_d(*args)
      else:            lib.ElSVMXDistSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  else: TypeExcept()

# Total variation denoising
# =========================
lib.ElTV_s.argtypes = \
lib.ElTVDist_s.argtypes = \
lib.ElTVDistSparse_s.argtypes = \
  [c_void_p,sType,c_void_p]
lib.ElTV_d.argtypes = \
lib.ElTVDist_d.argtypes = \
lib.ElTVDistSparse_d.argtypes = \
  [c_void_p,dType,c_void_p]

lib.ElTVX_s.argtypes = \
lib.ElTVXDist_s.argtypes = \
lib.ElTVXDistSparse_s.argtypes = \
  [c_void_p,sType,c_void_p,QPAffineCtrl_s]
lib.ElTVX_d.argtypes = \
lib.ElTVXDist_d.argtypes = \
lib.ElTVXDistSparse_d.argtypes = \
  [c_void_p,dType,c_void_p,QPAffineCtrl_d]

def TV(b,lambdPre,ctrl=None):
  lambd = TagToType(b.tag)(lambdPre)
  if type(b) is Matrix:
    x = Matrix(b.tag)
    args = [b.obj,lambd,x.obj]
    argsCtrl = [b.obj,lambd,x.obj,ctrl] 
    if   b.tag == sTag: 
      if ctrl == None: lib.ElTV_s(*args)
      else:            lib.ElTVX_s(*argsCtrl)
    elif b.tag == dTag: 
      if ctrl == None: lib.ElTV_d(*args)
      else:            lib.ElTVX_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(b) is DistMatrix:
    x = DistMatrix(b.tag,MC,MR,b.Grid())
    args = [b.obj,lambd,x.obj]
    argsCtrl = [b.obj,lambd,x.obj,ctrl] 
    if   b.tag == sTag: 
      if ctrl == None: lib.ElTVDist_s(*args)
      else:            lib.ElTVXDist_s(*argsCtrl)
    elif b.tag == dTag: 
      if ctrl == None: lib.ElTVDist_d(*args)
      else:            lib.ElTVXDist_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(b) is DistMultiVec:
    x = DistMultiVec(b.tag,b.Grid())
    args = [b.obj,lambd,x.obj]
    argsCtrl = [b.obj,lambd,x.obj,ctrl]
    if   b.tag == sTag: 
      if ctrl == None: lib.ElTVDistSparse_s(*args)
      else:            lib.ElTVXDistSparse_s(*argsCtrl)
    elif b.tag == dTag: 
      if ctrl == None: lib.ElTVDistSparse_d(*args)
      else:            lib.ElTVXDistSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  else: TypeExcept()

# Long-only portfolio
# ===================
lib.ElLongOnlyPortfolioSparse_s.argtypes = \
lib.ElLongOnlyPortfolioDistSparse_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,sType,c_void_p]

lib.ElLongOnlyPortfolioSparse_d.argtypes = \
lib.ElLongOnlyPortfolioDistSparse_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,dType,c_void_p]

lib.ElLongOnlyPortfolioXSparse_s.argtypes = \
lib.ElLongOnlyPortfolioXDistSparse_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,sType,c_void_p,SOCPAffineCtrl_s]

lib.ElLongOnlyPortfolioXSparse_d.argtypes = \
lib.ElLongOnlyPortfolioXDistSparse_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,dType,c_void_p,SOCPAffineCtrl_d]

lib.ElLongOnlyPortfolioSparse_s.restype = \
lib.ElLongOnlyPortfolioSparse_d.restype = \
lib.ElLongOnlyPortfolioDistSparse_s.restype = \
lib.ElLongOnlyPortfolioDistSparse_d.restype = \
lib.ElLongOnlyPortfolioXSparse_s.restype = \
lib.ElLongOnlyPortfolioXSparse_d.restype = \
lib.ElLongOnlyPortfolioXDistSparse_s.restype = \
lib.ElLongOnlyPortfolioXDistSparse_d.restype = \
  c_uint

def LongOnlyPortfolio(d,F,c,gammaPre,ctrl=None):
  gamma = TagToType(d.tag)(gammaPre)
  if type(F) is SparseMatrix:
    x = Matrix(d.tag)
    args = [d.obj,F.obj,c.obj,gamma,x.obj]
    argsCtrl = [d.obj,F.obj,c.obj,gamma,x.obj,ctrl]
    if   d.tag == sTag:
      if ctrl==None: lib.ElLongOnlyPortfolioSparse_s(*args)
      else:          lib.ElLongOnlyPortfolioXSparse_s(*argsCtrl)
    elif d.tag == dTag:
      if ctrl==None: lib.ElLongOnlyPortfolioSparse_d(*args)
      else:          lib.ElLongOnlyPortfolioXSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  elif type(F) is DistSparseMatrix:
    x = DistMultiVec(d.tag,d.Grid())
    args = [d.obj,F.obj,c.obj,gamma,x.obj]
    argsCtrl = [d.obj,F.obj,c.obj,gamma,x.obj,ctrl]
    if   d.tag == sTag:
      if ctrl==None: lib.ElLongOnlyPortfolioDistSparse_s(*args)
      else:          lib.ElLongOnlyPortfolioXDistSparse_s(*argsCtrl)
    elif d.tag == dTag:
      if ctrl==None: lib.ElLongOnlyPortfolioDistSparse_d(*args)
      else:          lib.ElLongOnlyPortfolioXDistSparse_d(*argsCtrl)
    else: DataExcept()
    return x
  else: TypeExcept()
