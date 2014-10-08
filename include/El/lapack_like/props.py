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

# Emulate an enum for the sign scaling
(SIGN_SCALE_NONE,SIGN_SCALE_DET,SIGN_SCALE_FROB)=(0,1,2)

lib.ElSignCtrlDefault_s.argtypes = [c_void_p]
lib.ElSignCtrlDefault_s.restype = c_uint
class SignCtrl_s(ctypes.Structure):
  _fields_ = [("maxIts",iType),
              ("tol",sType),
              ("power",sType),
              ("scaling",c_uint)]
  def __init__(self):
    lib.ElSignCtrlDefault_s(pointer(self))

lib.ElSignCtrlDefault_d.argtypes = [c_void_p]
lib.ElSignCtrlDefault_d.restype = c_uint
class SignCtrl_d(ctypes.Structure):
  _fields_ = [("maxIts",iType),
              ("tol",dType),
              ("power",dType),
              ("scaling",c_uint)]
  def __init__(self):
    lib.ElSignCtrlDefault_d(pointer(self))

lib.ElHessQRCtrlDefault.argtypes = [c_void_p]
lib.ElHessQRCtrlDefault.restype = c_uint
class HessQRCtrl(ctypes.Structure):
  _fields_ = [("distAED",bType),
              ("blockHeight",iType),("blockWidth",iType)]
  def __init__(self):
    lib.ElHessQRCtrlDefault(pointer(self))

lib.ElSDCCtrlDefault_s.argtypes = [c_void_p]
lib.ElSDCCtrlDefault_s.restype = c_uint
class SDCCtrl_s(ctypes.Structure):
  _fields_ = [("cutoff",iType),
              ("maxInnerIts",iType),("maxOuterIts",iType),
              ("tol",sType),
              ("spreadFactor",sType),
              ("random",bType),
              ("progress",bType),
              ("signCtrl",SignCtrl_s)]
  def __init__(self):
    lib.ElSDCCtrlDefault_s(pointer(self))

lib.ElSDCCtrlDefault_d.argtypes = [c_void_p]
lib.ElSDCCtrlDefault_d.restype = c_uint
class SDCCtrl_d(ctypes.Structure):
  _fields_ = [("cutoff",iType),
              ("maxInnerIts",iType),("maxOuterIts",iType),
              ("tol",dType),
              ("spreadFactor",dType),
              ("random",bType),
              ("progress",bType),
              ("signCtrl",SignCtrl_d)]
  def __init__(self):
    lib.ElSDCCtrlDefault_d(pointer(self))

lib.ElSchurCtrlDefault_s.argtypes = [c_void_p]
lib.ElSchurCtrlDefault_s.restype = c_uint
class SchurCtrl_s(ctypes.Structure):
  _fields_ = [("useSDC",bType),
              ("qrCtrl",HessQRCtrl),
              ("sdcCtrl",SDCCtrl_s)]
  def __init__(self):
    lib.ElSchurCtrlDefault_s(pointer(self))

lib.ElSchurCtrlDefault_d.argtypes = [c_void_p]
lib.ElSchurCtrlDefault_d.restype = c_uint
class SchurCtrl_d(ctypes.Structure):
  _fields_ = [("useSDC",bType),
              ("qrCtrl",HessQRCtrl),
              ("sdcCtrl",SDCCtrl_d)]
  def __init__(self):
    lib.ElSchurCtrlDefault_d(pointer(self))

lib.ElSnapshotCtrlDefault.argtypes = [c_void_p]
lib.ElSnapshotCtrlDefault.restype = c_uint
lib.ElSnapshotCtrlDestroy.argtypes = [c_void_p]
lib.ElSnapshotCtrlDestroy.restype = c_uint
class SnapshotCtrl(ctypes.Structure):
  _fields_ = [("realSize",iType),("imagSize",iType),
              ("imgSaveFreq",iType),("numSaveFreq",iType),
              ("imgDispFreq",iType),
              ("imgSaveCount",iType),("numSaveCount",iType),
              ("imgDispCount",iType),
              ("imgBase",c_char_p),("numBase",c_char_p),
              ("imgFormat",c_uint),("numFormat",c_uint),
              ("itCounts",bType)]
  def __init__(self):
    lib.ElSnaphsotCtrlDefault(pointer(self)) 
  def Destroy(self):
    lib.ElSnapshotCtrlDestroy(pointer(self))

# Emulate an enum for the pseudospectral norm
(PS_TWO_NORM,PS_ONE_NORM)=(0,1)

lib.ElPseudospecCtrlDefault_s.argtypes = [c_void_p]
lib.ElPseudospecCtrlDefault_s.restype = c_uint
lib.ElPseudospecCtrlDestroy_s.argtypes = [c_void_p]
lib.ElPseudospecCtrlDestroy_s.restype = c_uint
class PseudospecCtrl_s(ctypes.Structure):
  _fields_ = [("norm",c_uint),
              ("blockWidth",iType),
              ("schur",bType),
              ("forceComplexSchur",bType),
              ("forceComplexPs",bType),
              ("schurCtrl",SchurCtrl_s),
              ("maxIts",iType),
              ("tol",sType),
              ("deflate",bType),
              ("arnoldi",bType),
              ("basisSize",iType),
              ("reorthog",bType),
              ("progress",bType),
              ("snapCtrl",SnapshotCtrl)]
  def __init__(self):
    lib.ElPseudospecCtrlDefault_s(pointer(self))
  def Destroy(self):
    lib.ElPseudospecCtrlDestroy_s(pointer(self))

lib.ElPseudospecCtrlDefault_d.argtypes = [c_void_p]
lib.ElPseudospecCtrlDefault_d.restype = c_uint
lib.ElPseudospecCtrlDestroy_d.argtypes = [c_void_p]
lib.ElPseudospecCtrlDestroy_d.restype = c_uint
class PseudospecCtrl_d(ctypes.Structure):
  _fields_ = [("norm",c_uint),
              ("blockWidth",iType),
              ("schur",bType),
              ("forceComplexSchur",bType),
              ("forceComplexPs",bType),
              ("schurCtrl",SchurCtrl_d),
              ("maxIts",iType),
              ("tol",dType),
              ("deflate",bType),
              ("arnoldi",bType),
              ("basisSize",iType),
              ("reorthog",bType),
              ("progress",bType),
              ("snapCtrl",SnapshotCtrl)]
  def __init__(self):
    lib.ElPseudospecCtrlDefault_d(pointer(self))
  def Destroy(self):
    lib.ElPseudospecCtrlDestroy_d(pointer(self))

lib.ElPseudospectralAutoWindowX_s.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_s]
lib.ElPseudospectralAutoWindowX_s.restype = c_uint
lib.ElPseudospectralAutoWindowX_d.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_d]
lib.ElPseudospectralAutoWindowX_d.restype = c_uint
lib.ElPseudospectralAutoWindowX_c.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_s]
lib.ElPseudospectralAutoWindowX_c.restype = c_uint
lib.ElPseudospectralAutoWindowX_z.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_d]
lib.ElPseudospectralAutoWindowX_z.restype = c_uint
lib.ElPseudospectralAutoWindowXDist_s.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_s]
lib.ElPseudospectralAutoWindowXDist_s.restype = c_uint
lib.ElPseudospectralAutoWindowXDist_d.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_d]
lib.ElPseudospectralAutoWindowXDist_d.restype = c_uint
lib.ElPseudospectralAutoWindowXDist_c.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_s]
lib.ElPseudospectralAutoWindowXDist_c.restype = c_uint
lib.ElPseudospectralAutoWindowXDist_z.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_d]
lib.ElPseudospectralAutoWindowXDist_z.restype = c_uint
def PseudospectralAutoWindow(A,invNormMap,realSize,imagSize,ctrl):
  if type(A) is not type(invNormMap): 
    raise Exception('Types of matrices must match')
  if invNormMap.tag != Base(A.tag): 
    raise Exception('Datatype of invNormMap must be base of datatype of A')
  if type(A) is Matrix:
    if   A.tag == sTag: 
      lib.ElPseudospectralAutoWindowX_s \
      (A.obj,invNormMap.obj,realSize,imagSize,ctrl)
    elif A.tag == dTag:
      lib.ElPseudospectralAutoWindowX_d \
      (A.obj,invNormMap.obj,realSize,imagSize,ctrl)
    elif A.tag == cTag:
      lib.ElPseudospectralAutoWindowX_c \
      (A.obj,invNormMap.obj,realSize,imagSize,ctrl)
    elif A.tag == zTag:
      lib.ElPseudospectralAutoWindowX_z \
      (A.obj,invNormMap.obj,realSize,imagSize,ctrl)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElPseudospectralAutoWindowXDist_s \
      (A.obj,invNormMap.obj,realSize,imagSize,ctrl)
    elif A.tag == dTag:
      lib.ElPseudospectralAutoWindowXDist_d \
      (A.obj,invNormMap.obj,realSize,imagSize,ctrl)
    elif A.tag == cTag:
      lib.ElPseudospectralAutoWindowXDist_c \
      (A.obj,invNormMap.obj,realSize,imagSize,ctrl)
    elif A.tag == zTag:
      lib.ElPseudospectralAutoWindowXDist_z \
      (A.obj,invNormMap.obj,realSize,imagSize,ctrl)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElPseudospectralWindowX_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,sType,iType,iType,PseudospecCtrl_s]
lib.ElPseudospectralWindowX_s.restype = c_uint
lib.ElPseudospectralWindowX_d.argtypes = \
  [c_void_p,c_void_p,dType,dType,dType,iType,iType,PseudospecCtrl_d]
lib.ElPseudospectralWindowX_d.restype = c_uint
lib.ElPseudospectralWindowX_c.argtypes = \
  [c_void_p,c_void_p,cType,sType,sType,iType,iType,PseudospecCtrl_s]
lib.ElPseudospectralWindowX_c.restype = c_uint
lib.ElPseudospectralWindowX_z.argtypes = \
  [c_void_p,c_void_p,zType,dType,dType,iType,iType,PseudospecCtrl_d]
lib.ElPseudospectralWindowX_z.restype = c_uint
lib.ElPseudospectralWindowXDist_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,sType,iType,iType,PseudospecCtrl_s]
lib.ElPseudospectralWindowXDist_s.restype = c_uint
lib.ElPseudospectralWindowXDist_d.argtypes = \
  [c_void_p,c_void_p,dType,dType,dType,iType,iType,PseudospecCtrl_d]
lib.ElPseudospectralWindowXDist_d.restype = c_uint
lib.ElPseudospectralWindowXDist_c.argtypes = \
  [c_void_p,c_void_p,cType,sType,sType,iType,iType,PseudospecCtrl_s]
lib.ElPseudospectralWindowXDist_c.restype = c_uint
lib.ElPseudospectralWindowXDist_z.argtypes = \
  [c_void_p,c_void_p,zType,dType,dType,iType,iType,PseudospecCtrl_d]
lib.ElPseudospectralWindowXDist_z.restype = c_uint
def PseudospectralWindow \
    (A,invNormMap,centerPre,realWidth,imagWidth,realSize,imagSize,ctrl):
  center = TagToType(A.tag)(centerPre)
  if type(A) is not type(invNormMap): 
    raise Exception('Types of matrices must match')
  if invNormMap.tag != Base(A.tag): 
    raise Exception('Datatype of invNormMap must be base of datatype of A')
  if type(A) is Matrix:
    if   A.tag == sTag: 
      lib.ElPseudospectralWindowX_s \
      (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize,ctrl)
    elif A.tag == dTag:
      lib.ElPseudospectralWindowX_d \
      (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize,ctrl)
    elif A.tag == cTag:
      lib.ElPseudospectralWindowX_c \
      (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize,ctrl)
    elif A.tag == zTag:
      lib.ElPseudospectralWindowX_z \
      (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize,ctrl)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElPseudospectralWindowXDist_s \
      (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize,ctrl)
    elif A.tag == dTag:
      lib.ElPseudospectralWindowXDist_d \
      (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize,ctrl)
    elif A.tag == cTag:
      lib.ElPseudospectralWindowXDist_c \
      (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize,ctrl)
    elif A.tag == zTag:
      lib.ElPseudospectralWindowXDist_z \
      (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize,ctrl)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElPseudospectralCloudX_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_s]
lib.ElPseudospectralCloudX_s.restype = c_uint
lib.ElPseudospectralCloudX_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_d]
lib.ElPseudospectralCloudX_d.restype = c_uint
lib.ElPseudospectralCloudX_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_s]
lib.ElPseudospectralCloudX_c.restype = c_uint
lib.ElPseudospectralCloudX_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_d]
lib.ElPseudospectralCloudX_z.restype = c_uint
lib.ElPseudospectralCloudXDist_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_s]
lib.ElPseudospectralCloudXDist_s.restype = c_uint
lib.ElPseudospectralCloudXDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_d]
lib.ElPseudospectralCloudXDist_d.restype = c_uint
lib.ElPseudospectralCloudXDist_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_s]
lib.ElPseudospectralCloudXDist_c.restype = c_uint
lib.ElPseudospectralCloudXDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_d]
lib.ElPseudospectralCloudXDist_z.restype = c_uint
def PseudospectralCloud(A,shifts,invNorms,ctrl):
  if type(A) is not type(invNorms): 
    raise Exception('Types of matrices must match')
  if invNorms.tag != Base(A.tag): 
    raise Exception('Datatype of invNorms must be base of datatype of A')
  if type(A) is Matrix:
    if   A.tag == sTag: 
      lib.ElPseudospectralCloudX_s(A.obj,shifts.obj,invNorms.obj,ctrl)
    elif A.tag == dTag:
      lib.ElPseudospectralCloudX_d(A.obj,shifts.obj,invNorms.obj,ctrl)
    elif A.tag == cTag:
      lib.ElPseudospectralCloudX_c(A.obj,shifts.obj,invNorms.obj,ctrl)
    elif A.tag == zTag:
      lib.ElPseudospectralCloudX_z(A.obj,shifts.obj,invNorms.obj,ctrl)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElPseudospectralCloudXDist_s(A.obj,shifts.obj,invNorms.obj,ctrl)
    elif A.tag == dTag:
      lib.ElPseudospectralCloudXDist_d(A.obj,shifts.obj,invNorms.obj,ctrl)
    elif A.tag == cTag:
      lib.ElPseudospectralCloudXDist_c(A.obj,shifts.obj,invNorms.obj,ctrl)
    elif A.tag == zTag:
      lib.ElPseudospectralCloudXDist_z(A.obj,shifts.obj,invNorms.obj,ctrl)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')
