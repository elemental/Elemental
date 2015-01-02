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

# Apply packed reflectors
# =======================
lib.ElApplyPackedReflectors_s.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,iType,c_void_p,c_void_p,c_void_p]
lib.ElApplyPackedReflectors_s.restype = c_uint
lib.ElApplyPackedReflectors_d.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,iType,c_void_p,c_void_p,c_void_p]
lib.ElApplyPackedReflectors_d.restype = c_uint
lib.ElApplyPackedReflectors_c.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,iType,c_void_p,c_void_p,c_void_p]
lib.ElApplyPackedReflectors_c.restype = c_uint
lib.ElApplyPackedReflectors_z.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,iType,c_void_p,c_void_p,c_void_p]
lib.ElApplyPackedReflectors_z.restype = c_uint
lib.ElApplyPackedReflectorsDist_s.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,iType,c_void_p,c_void_p,c_void_p]
lib.ElApplyPackedReflectorsDist_s.restype = c_uint
lib.ElApplyPackedReflectorsDist_d.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,iType,c_void_p,c_void_p,c_void_p]
lib.ElApplyPackedReflectorsDist_d.restype = c_uint
lib.ElApplyPackedReflectorsDist_c.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,iType,c_void_p,c_void_p,c_void_p]
lib.ElApplyPackedReflectorsDist_c.restype = c_uint
lib.ElApplyPackedReflectorsDist_z.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,iType,c_void_p,c_void_p,c_void_p]
lib.ElApplyPackedReflectorsDist_z.restype = c_uint
def ApplyPackedReflectors(side,uplo,dir,order,offset,H,t,A):
  if type(H) is not type(t) or type(t) is not type(A):
    raise Exception('Matrix types of {H,t,A} must match')
  if H.tag != t.tag or t.tag != A.tag:
    raise Exception('Datatypes of {H,t,A} must match')
  args = [side,uplo,dir,order,offset,H.obj,t.obj,A.obj]
  if type(H) is Matrix:
    if   H.tag == sTag: lib.ElApplyPackedReflectors_s(*args)
    elif H.tag == dTag: lib.ElApplyPackedReflectors_d(*args)
    elif H.tag == cTag: lib.ElApplyPackedReflectors_c(*args)
    elif H.tag == zTag: lib.ElApplyPackedReflectors_z(*args)
    else: DataExcept()
  elif type(H) is DistMatrix:
    if   H.tag == sTag: lib.ElApplyPackedReflectorsDist_s(*args)
    elif H.tag == dTag: lib.ElApplyPackedReflectorsDist_d(*args)
    elif H.tag == cTag: lib.ElApplyPackedReflectorsDist_c(*args)
    elif H.tag == zTag: lib.ElApplyPackedReflectorsDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Expand packed reflectors
# ========================
lib.ElExpandPackedReflectors_s.argtypes = \
  [c_uint,c_uint,iType,c_void_p,c_void_p]
lib.ElExpandPackedReflectors_s.restype = c_uint
lib.ElExpandPackedReflectors_d.argtypes = \
  [c_uint,c_uint,iType,c_void_p,c_void_p]
lib.ElExpandPackedReflectors_d.restype = c_uint
lib.ElExpandPackedReflectors_c.argtypes = \
  [c_uint,c_uint,iType,c_void_p,c_void_p]
lib.ElExpandPackedReflectors_c.restype = c_uint
lib.ElExpandPackedReflectors_z.argtypes = \
  [c_uint,c_uint,iType,c_void_p,c_void_p]
lib.ElExpandPackedReflectors_z.restype = c_uint
lib.ElExpandPackedReflectorsDist_s.argtypes = \
  [c_uint,c_uint,iType,c_void_p,c_void_p]
lib.ElExpandPackedReflectorsDist_s.restype = c_uint
lib.ElExpandPackedReflectorsDist_d.argtypes = \
  [c_uint,c_uint,iType,c_void_p,c_void_p]
lib.ElExpandPackedReflectorsDist_d.restype = c_uint
lib.ElExpandPackedReflectorsDist_c.argtypes = \
  [c_uint,c_uint,iType,c_void_p,c_void_p]
lib.ElExpandPackedReflectorsDist_c.restype = c_uint
lib.ElExpandPackedReflectorsDist_z.argtypes = \
  [c_uint,c_uint,iType,c_void_p,c_void_p]
lib.ElExpandPackedReflectorsDist_z.restype = c_uint
def ExpandPackedReflectors(uplo,dir,offset,H,t):
  if type(H) is not type(t):
    raise Exception('Types of H and t must match')
  if H.tag != t.tag:
    raise Exception('Datatypes of H and t must match')
  args = [uplo,dir,offset,H.obj,t.obj]
  if type(H) is Matrix:
    if   H.tag == sTag: lib.ElExpandPackedReflectors_s(*args)
    elif H.tag == dTag: lib.ElExpandPackedReflectors_d(*args)
    elif H.tag == cTag: lib.ElExpandPackedReflectors_c(*args)
    elif H.tag == zTag: lib.ElExpandPackedReflectors_z(*args)
    else: DataExcept()
  elif type(H) is DistMatrix:
    if   H.tag == sTag: lib.ElExpandPackedReflectorsDist_s(*args)
    elif H.tag == dTag: lib.ElExpandPackedReflectorsDist_d(*args)
    elif H.tag == cTag: lib.ElExpandPackedReflectorsDist_c(*args)
    elif H.tag == zTag: lib.ElExpandPackedReflectorsDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Hyperbolic reflector
# ====================

# Left application
# ----------------
lib.ElLeftHyperbolicReflector_s.argtypes = \
  [POINTER(sType),c_void_p,POINTER(sType)]
lib.ElLeftHyperbolicReflector_s.restype = c_uint
lib.ElLeftHyperbolicReflector_d.argtypes = \
  [POINTER(dType),c_void_p,POINTER(dType)]
lib.ElLeftHyperbolicReflector_d.restype = c_uint
lib.ElLeftHyperbolicReflector_c.argtypes = \
  [POINTER(cType),c_void_p,POINTER(cType)]
lib.ElLeftHyperbolicReflector_c.restype = c_uint
lib.ElLeftHyperbolicReflector_z.argtypes = \
  [POINTER(zType),c_void_p,POINTER(zType)]
lib.ElLeftHyperbolicReflector_z.restype = c_uint
lib.ElLeftHyperbolicReflectorDist_s.argtypes = \
  [POINTER(sType),c_void_p,POINTER(sType)]
lib.ElLeftHyperbolicReflectorDist_s.restype = c_uint
lib.ElLeftHyperbolicReflectorDist_d.argtypes = \
  [POINTER(dType),c_void_p,POINTER(dType)]
lib.ElLeftHyperbolicReflectorDist_d.restype = c_uint
lib.ElLeftHyperbolicReflectorDist_c.argtypes = \
  [POINTER(cType),c_void_p,POINTER(cType)]
lib.ElLeftHyperbolicReflectorDist_c.restype = c_uint
lib.ElLeftHyperbolicReflectorDist_z.argtypes = \
  [POINTER(zType),c_void_p,POINTER(zType)]
lib.ElLeftHyperbolicReflectorDist_z.restype = c_uint
def LeftHyperbolicReflector(chi,x):
  alpha = TagToType(x.tag)(chi)
  tau = TagToType(x.tag)()
  args = [pointer(alpha),x.obj,pointer(tau)]
  if type(x) is Matrix:
    if   x.tag == sTag: lib.ElLeftHyperbolicReflector_s(*args)
    elif x.tag == dTag: lib.ElLeftHyperbolicReflector_d(*args)
    elif x.tag == cTag: lib.ElLeftHyperbolicReflector_c(*args)
    elif x.tag == zTag: lib.ElLeftHyperbolicReflector_z(*args)
    else: DataExcept()
  elif type(x) is DistMatrix:
    if   x.tag == sTag: lib.ElLeftHyperbolicReflectorDist_s(*args)
    elif x.tag == dTag: lib.ElLeftHyperbolicReflectorDist_d(*args)
    elif x.tag == cTag: lib.ElLeftHyperbolicReflectorDist_c(*args)
    elif x.tag == zTag: lib.ElLeftHyperbolicReflectorDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return alpha, tau

# Right application
# -----------------
lib.ElRightHyperbolicReflector_s.argtypes = \
  [POINTER(sType),c_void_p,POINTER(sType)]
lib.ElRightHyperbolicReflector_s.restype = c_uint
lib.ElRightHyperbolicReflector_d.argtypes = \
  [POINTER(dType),c_void_p,POINTER(dType)]
lib.ElRightHyperbolicReflector_d.restype = c_uint
lib.ElRightHyperbolicReflector_c.argtypes = \
  [POINTER(cType),c_void_p,POINTER(cType)]
lib.ElRightHyperbolicReflector_c.restype = c_uint
lib.ElRightHyperbolicReflector_z.argtypes = \
  [POINTER(zType),c_void_p,POINTER(zType)]
lib.ElRightHyperbolicReflector_z.restype = c_uint
lib.ElRightHyperbolicReflectorDist_s.argtypes = \
  [POINTER(sType),c_void_p,POINTER(sType)]
lib.ElRightHyperbolicReflectorDist_s.restype = c_uint
lib.ElRightHyperbolicReflectorDist_d.argtypes = \
  [POINTER(dType),c_void_p,POINTER(dType)]
lib.ElRightHyperbolicReflectorDist_d.restype = c_uint
lib.ElRightHyperbolicReflectorDist_c.argtypes = \
  [POINTER(cType),c_void_p,POINTER(cType)]
lib.ElRightHyperbolicReflectorDist_c.restype = c_uint
lib.ElRightHyperbolicReflectorDist_z.argtypes = \
  [POINTER(zType),c_void_p,POINTER(zType)]
lib.ElRightHyperbolicReflectorDist_z.restype = c_uint
def RightHyperbolicReflector(chi,x):
  alpha = TagToType(x.tag)(chi)
  tau = TagToType(x.tag)()
  args = [pointer(alpha),x.obj,pointer(tau)]
  if type(x) is Matrix:
    if   x.tag == sTag: lib.ElRightHyperbolicReflector_s(*args)
    elif x.tag == dTag: lib.ElRightHyperbolicReflector_d(*args)
    elif x.tag == cTag: lib.ElRightHyperbolicReflector_c(*args)
    elif x.tag == zTag: lib.ElRightHyperbolicReflector_z(*args)
    else: DataExcept()
  elif type(x) is DistMatrix:
    if   x.tag == sTag: lib.ElRightHyperbolicReflectorDist_s(*args)
    elif x.tag == dTag: lib.ElRightHyperbolicReflectorDist_d(*args)
    elif x.tag == cTag: lib.ElRightHyperbolicReflectorDist_c(*args)
    elif x.tag == zTag: lib.ElRightHyperbolicReflectorDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return alpha, tau

# Householder reflector
# =====================

# Left application
# ----------------
lib.ElLeftReflector_s.argtypes = [POINTER(sType),c_void_p,POINTER(sType)]
lib.ElLeftReflector_s.restype = c_uint
lib.ElLeftReflector_d.argtypes = [POINTER(dType),c_void_p,POINTER(dType)]
lib.ElLeftReflector_d.restype = c_uint
lib.ElLeftReflector_c.argtypes = [POINTER(cType),c_void_p,POINTER(cType)]
lib.ElLeftReflector_c.restype = c_uint
lib.ElLeftReflector_z.argtypes = [POINTER(zType),c_void_p,POINTER(zType)]
lib.ElLeftReflector_z.restype = c_uint
lib.ElLeftReflectorDist_s.argtypes = [POINTER(sType),c_void_p,POINTER(sType)]
lib.ElLeftReflectorDist_s.restype = c_uint
lib.ElLeftReflectorDist_d.argtypes = [POINTER(dType),c_void_p,POINTER(dType)]
lib.ElLeftReflectorDist_d.restype = c_uint
lib.ElLeftReflectorDist_c.argtypes = [POINTER(cType),c_void_p,POINTER(cType)]
lib.ElLeftReflectorDist_c.restype = c_uint
lib.ElLeftReflectorDist_z.argtypes = [POINTER(zType),c_void_p,POINTER(zType)]
lib.ElLeftReflectorDist_z.restype = c_uint
def LeftReflector(chi,x):
  alpha = TagToType(x.tag)(chi)
  tau = TagToType(x.tag)()
  args = [pointer(alpha),x.obj,pointer(tau)]
  if type(x) is Matrix:
    if   x.tag == sTag: lib.ElLeftReflector_s(*args)
    elif x.tag == dTag: lib.ElLeftReflector_d(*args)
    elif x.tag == cTag: lib.ElLeftReflector_c(*args)
    elif x.tag == zTag: lib.ElLeftReflector_z(*args)
    else: DataExcept()
  elif type(x) is DistMatrix:
    if   x.tag == sTag: lib.ElLeftReflectorDist_s(*args)
    elif x.tag == dTag: lib.ElLeftReflectorDist_d(*args)
    elif x.tag == cTag: lib.ElLeftReflectorDist_c(*args)
    elif x.tag == zTag: lib.ElLeftReflectorDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return alpha, tau

# Right application
# -----------------
lib.ElRightReflector_s.argtypes = [POINTER(sType),c_void_p,POINTER(sType)]
lib.ElRightReflector_s.restype = c_uint
lib.ElRightReflector_d.argtypes = [POINTER(dType),c_void_p,POINTER(dType)]
lib.ElRightReflector_d.restype = c_uint
lib.ElRightReflector_c.argtypes = [POINTER(cType),c_void_p,POINTER(cType)]
lib.ElRightReflector_c.restype = c_uint
lib.ElRightReflector_z.argtypes = [POINTER(zType),c_void_p,POINTER(zType)]
lib.ElRightReflector_z.restype = c_uint
lib.ElRightReflectorDist_s.argtypes = [POINTER(sType),c_void_p,POINTER(sType)]
lib.ElRightReflectorDist_s.restype = c_uint
lib.ElRightReflectorDist_d.argtypes = [POINTER(dType),c_void_p,POINTER(dType)]
lib.ElRightReflectorDist_d.restype = c_uint
lib.ElRightReflectorDist_c.argtypes = [POINTER(cType),c_void_p,POINTER(cType)]
lib.ElRightReflectorDist_c.restype = c_uint
lib.ElRightReflectorDist_z.argtypes = [POINTER(zType),c_void_p,POINTER(zType)]
lib.ElRightReflectorDist_z.restype = c_uint
def RightReflector(chi,x):
  alpha = TagToType(x.tag)(chi)
  tau = TagToType(x.tag)()
  args = [pointer(alpha),x.obj,pointer(tau)]
  if type(x) is Matrix:
    if   x.tag == sTag: lib.ElRightReflector_s(*args)
    elif x.tag == dTag: lib.ElRightReflector_d(*args)
    elif x.tag == cTag: lib.ElRightReflector_c(*args)
    elif x.tag == zTag: lib.ElRightReflector_z(*args)
    else: DataExcept()
  elif type(x) is DistMatrix:
    if   x.tag == sTag: lib.ElRightReflectorDist_s(*args)
    elif x.tag == dTag: lib.ElRightReflectorDist_d(*args)
    elif x.tag == cTag: lib.ElRightReflectorDist_c(*args)
    elif x.tag == zTag: lib.ElRightReflectorDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return alpha, tau
