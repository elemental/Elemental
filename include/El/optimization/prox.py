#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *

# Clipping
# ========
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

# Frobenius-norm proximal map
# ===========================
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
# =======================
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

# Logistic proximal map
# =====================
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
# ================================
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
# =================
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
