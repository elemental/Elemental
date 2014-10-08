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

class SignCtrl_s(ctypes.Structure):
  _fields_ = [("maxIts",iType),
              ("tol",sType),
              ("power",sType),
              ("scaling",c_uint)]
lib.ElSignCtrlDefault_s.argtypes = [POINTER(SignCtrl_s)]
lib.ElSignCtrlDefault_s.restype = c_uint
def SignCtrlDefault_s(ctrl):
  lib.ElSignCtrlDefault_s(pointer(ctrl))

class SignCtrl_d(ctypes.Structure):
  _fields_ = [("maxIts",iType),
              ("tol",dType),
              ("power",dType),
              ("scaling",c_uint)]
lib.ElSignCtrlDefault_d.argtypes = [POINTER(SignCtrl_d)]
lib.ElSignCtrlDefault_d.restype = c_uint
def SignCtrlDefault_d(ctrl):
  lib.ElSignCtrlDefault_d(pointer(ctrl))

class HessQRCtrl(ctypes.Structure):
  _fields_ = [("distAED",bType),
              ("blockHeight",iType),("blockWidth",iType)]
lib.ElHessQRCtrlDefault.argtypes = [POINTER(HessQRCtrl)]
lib.ElHessQRCtrlDefault.restype = c_uint
def HessQRCtrlDefault(ctrl):
  lib.ElHessQRCtrlDefault(pointer(ctrl))

class SDCCtrl_s(ctypes.Structure):
  _fields_ = [("cutoff",iType),
              ("maxInnerIts",iType),("maxOuterIts",iType),
              ("tol",sType),
              ("spreadFactor",sType),
              ("random",bType),
              ("progress",bType),
              ("signCtrl",SignCtrl_s)]
lib.ElSDCCtrlDefault_s.argtypes = [POINTER(SDCCtrl_s)]
lib.ElSDCCtrlDefault_s.restype = c_uint
def SDCCtrlDefault_s(ctrl):
  lib.ElSDCCtrlDefault_s(pointer(ctrl))

class SDCCtrl_d(ctypes.Structure):
  _fields_ = [("cutoff",iType),
              ("maxInnerIts",iType),("maxOuterIts",iType),
              ("tol",dType),
              ("spreadFactor",dType),
              ("random",bType),
              ("progress",bType),
              ("signCtrl",SignCtrl_d)]
lib.ElSDCCtrlDefault_d.argtypes = [POINTER(SDCCtrl_d)]
lib.ElSDCCtrlDefault_d.restype = c_uint
def SDCCtrlDefault_d(ctrl):
  lib.ElSDCCtrlDefault_d(pointer(ctrl))

class SchurCtrl_s(ctypes.Structure):
  _fields_ = [("useSDC",bType),
              ("qrCtrl",HessQRCtrl),
              ("sdcCtrl",SDCCtrl_s)]
lib.ElSchurCtrlDefault_s.argtypes = [POINTER(SchurCtrl_s)]
lib.ElSchurCtrlDefault_s.restype = c_uint
def SchurCtrlDefault_s(ctrl):
  lib.ElSchurCtrlDefault_s(pointer(ctrl))

class SchurCtrl_d(ctypes.Structure):
  _fields_ = [("useSDC",bType),
              ("qrCtrl",HessQRCtrl),
              ("sdcCtrl",SDCCtrl_d)]
lib.ElSchurCtrlDefault_d.argtypes = [POINTER(SchurCtrl_d)]
lib.ElSchurCtrlDefault_d.restype = c_uint
def SchurCtrlDefault_d(ctrl):
  lib.ElSchurCtrlDefault_d(pointer(ctrl))

class SnapshotCtrl(ctypes.Structure):
  _fields_ = [("realSize",iType),("imagSize",iType),
              ("imgSaveFreq",iType),("numSaveFreq",iType),
              ("imgDispFreq",iType),
              ("imgSaveCount",iType),("numSaveCount",iType),
              ("imgDispCount",iType),
              ("imgBase",c_char_p),("numBase",c_char_p),
              ("imgFormat",c_uint),("numFormat",c_uint),
              ("itCounts",bType)]
lib.ElSnapshotCtrlDefault.argtypes = [POINTER(SnapshotCtrl)]
lib.ElSnapshotCtrlDefault.restype = c_uint
lib.ElSnapshotCtrlDestroy.argtypes = [POINTER(SnapshotCtrl)]
lib.ElSnapshotCtrlDestroy.restype = c_uint
def SnapshotCtrlDefault(ctrl):
  lib.ElSnaphsotCtrlDefault(pointer(ctrl)) 
def SnapshotCtrlDestroy(ctrl):
  lib.ElSnapshotCtrlDestroy(pointer(ctrl))

# Emulate an enum for the pseudospectral norm
(PS_TWO_NORM,PS_ONE_NORM)=(0,1)

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
lib.ElPseudospecCtrlDefault_s.argtypes = [POINTER(PseudospecCtrl_s)]
lib.ElPseudospecCtrlDefault_s.restype = c_uint
lib.ElPseudospecCtrlDestroy_s.argtypes = [POINTER(PseudospecCtrl_s)]
lib.ElPseudospecCtrlDestroy_s.restype = c_uint
def PseudospecCtrlDefault_s(ctrl):
  lib.ElPseudospecCtrlDefault_s(pointer(ctrl))
def PseudospecCtrlDestroy_s(ctrl):
  lib.ElPseudospecCtrlDestroy_s(pointer(ctrl))

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
lib.ElPseudospecCtrlDefault_d.argtypes = [POINTER(PseudospecCtrl_d)]
lib.ElPseudospecCtrlDefault_d.restype = c_uint
lib.ElPseudospecCtrlDestroy_d.argtypes = [POINTER(PseudospecCtrl_d)]
lib.ElPseudospecCtrlDestroy_d.restype = c_uint
def PseudospecCtrlDefault_d(ctrl):
  lib.ElPseudospecCtrlDefault_d(pointer(ctrl))
def PseudospecCtrlDestroy_d(ctrl):
  lib.ElPseudospecCtrlDestroy_d(pointer(ctrl))

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
