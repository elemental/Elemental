#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *

# Coherence
# =========
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
# ==========
lib.ElCovariance_s.argtypes = \
lib.ElCovariance_d.argtypes = \
lib.ElCovariance_c.argtypes = \
lib.ElCovariance_z.argtypes = \
lib.ElCovarianceDist_s.argtypes = \
lib.ElCovarianceDist_d.argtypes = \
lib.ElCovarianceDist_c.argtypes = \
lib.ElCovarianceDist_z.argtypes = \
  [c_void_p,c_void_p]

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

# Log barrier
# ===========
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
# ==================
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

# SOC Identity
# ============
# TODO

# SOC dots
# ========
lib.ElSOCDots_s.argtypes = \
lib.ElSOCDots_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSOCDotsDist_s.argtypes = \
lib.ElSOCDotsDist_d.argtypes = \
lib.ElSOCDotsDistMultiVec_s.argtypes = \
lib.ElSOCDotsDistMultiVec_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_int]

def SOCDots(x,y,orders,firstInds,cutoff=1000):
  # TODO: Sanity checking
  if type(x) is Matrix:
    z = Matrix(x.tag)
    args = [x.obj,y.obj,z.obj,orders.obj,firstInds.obj]
    if   x.tag == sTag: lib.ElSOCDots_s(*args)
    elif x.tag == dTag: lib.ElSOCDots_d(*args)
    else: DataExcept()
    return z
  elif type(x) is DistMatrix:
    z = DistMatrix(x.tag,VC,STAR,x.Grid())
    args = [x.obj,y.obj,z.obj,orders.obj,firstInds.obj,cutoff]
    if   x.tag == sTag: lib.ElSOCDotsDist_s(*args)
    elif x.tag == dTag: lib.ElSOCDotsDist_d(*args)
    else: DataExcept()
    return z
  elif type(x) is DistMultiVec:
    z = DistMultiVec(x.tag,x.Comm())
    args = [x.obj,y.obj,z.obj,orders.obj,firstInds.obj,cutoff]
    if   x.tag == sTag: lib.ElSOCDotsDistMultiVec_s(*args)
    elif x.tag == dTag: lib.ElSOCDotsDistMultiVec_d(*args)
    else: DataExcept()
    return z
  else: TypeExcept()

# Cone Broadcast
# ==============
lib.ElConeBroadcast_s.argtypes = \
lib.ElConeBroadcast_d.argtypes = \
  [c_void_p,c_void_p,c_void_p]
lib.ElConeBroadcastDist_s.argtypes = \
lib.ElConeBroadcastDist_d.argtypes = \
lib.ElConeBroadcastDistMultiVec_s.argtypes = \
lib.ElConeBroadcastDistMultiVec_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_int]

def ConeBroadcast(x,orders,firstInds,cutoff=1000):
  # TODO: Sanity checking
  if type(x) is Matrix:
    args = [x.obj,orders.obj,firstInds.obj]
    if   x.tag == sTag: lib.ElConeBroadcast_s(*args)
    elif x.tag == dTag: lib.ElConeBroadcast_d(*args)
    else: DataExcept()
  elif type(x) is DistMatrix: 
    args = [x.obj,orders.obj,firstInds.obj,cutoff]
    if   x.tag == sTag: lib.ElConeBroadcastDist_s(*args)
    elif x.tag == dTag: lib.ElConeBroadcastDist_d(*args)
    else: DataExcept()
  elif type(x) is DistMultiVec:
    args = [x.obj,orders.obj,firstInds.obj,cutoff]
    if   x.tag == sTag: lib.ElConeBroadcastDistMultiVec_s(*args)
    elif x.tag == dTag: lib.ElConeBroadcastDistMultiVec_d(*args)
    else: DataExcept()
  else: TypeExcept()

# SOC Identity
# ============
lib.ElSOCIdentity_s.argtypes = \
lib.ElSOCIdentity_d.argtypes = \
lib.ElSOCIdentityDist_s.argtypes = \
lib.ElSOCIdentityDist_d.argtypes = \
lib.ElSOCIdentityDistMultiVec_s.argtypes = \
lib.ElSOCIdentityDistMultiVec_d.argtypes = \
  [c_void_p,c_void_p,c_void_p]

def SOCIdentity(x,orders,firstInds):
  # TODO: Sanity checking
  args = [x.obj,orders.obj,firstInds.obj]
  if type(x) is Matrix:
    if   x.tag == sTag: lib.ElSOCIdentity_s(*args)
    elif x.tag == dTag: lib.ElSOCIdentity_d(*args)
    else: DataExcept()
  elif type(x) is DistMatrix:
    if   x.tag == sTag: lib.ElSOCIdentityDist_s(*args)
    elif x.tag == dTag: lib.ElSOCIdentityDist_d(*args)
    else: DataExcept()
  elif type(x) is DistMultiVec:
    if   x.tag == sTag: lib.ElSOCIdentityDistMultiVec_s(*args)
    elif x.tag == dTag: lib.ElSOCIdentityDistMultiVec_d(*args)
    else: DataExcept()
  else: TypeExcept()

# SOC Reflect
# ===========
lib.ElSOCReflect_s.argtypes = \
lib.ElSOCReflect_d.argtypes = \
lib.ElSOCReflectDist_s.argtypes = \
lib.ElSOCReflectDist_d.argtypes = \
lib.ElSOCReflectDistMultiVec_s.argtypes = \
lib.ElSOCReflectDistMultiVec_d.argtypes = \
  [c_void_p,c_void_p,c_void_p]

def SOCReflect(x,orders,firstInds):
  # TODO: Sanity checking
  args = [x.obj,orders.obj,firstInds.obj]
  if type(x) is Matrix:
    if   x.tag == sTag: lib.ElSOCReflect_s(*args)
    elif x.tag == dTag: lib.ElSOCReflect_d(*args)
    else: DataExcept()
  elif type(x) is DistMatrix:
    if   x.tag == sTag: lib.ElSOCReflectDist_s(*args)
    elif x.tag == dTag: lib.ElSOCReflectDist_d(*args)
    else: DataExcept()
  elif type(x) is DistMultiVec:
    if   x.tag == sTag: lib.ElSOCReflectDistMultiVec_s(*args)
    elif x.tag == dTag: lib.ElSOCReflectDistMultiVec_d(*args)
    else: DataExcept()
  else: TypeExcept()

# SOC Determinants
# ================
lib.ElSOCDets_s.argtypes = \
lib.ElSOCDets_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSOCDetsDist_s.argtypes = \
lib.ElSOCDetsDist_d.argtypes = \
lib.ElSOCDetsDistMultiVec_s.argtypes = \
lib.ElSOCDetsDistMultiVec_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_int]

def SOCDets(x,orders,firstInds,cutoff=1000):
  # TODO: Sanity checking
  if type(x) is Matrix:
    d = Matrix(x.tag)
    args = [x.obj,d.obj,orders.obj,firstInds.obj]
    if   x.tag == sTag: lib.ElSOCDets_s(*args)
    elif x.tag == dTag: lib.ElSOCDets_d(*args)
    else: DataExcept()
    return d
  elif type(x) is DistMatrix:
    d = DistMatrix(x.tag,VC,STAR,x.Grid())
    args = [x.obj,d.obj,orders.obj,firstInds.obj,cutoff]
    if   x.tag == sTag: lib.ElSOCDetsDist_s(*args)
    elif x.tag == dTag: lib.ElSOCDetsDist_d(*args)
    else: DataExcept()
    return d
  elif type(x) is DistMultiVec:
    d = DistMultiVec(x.tag,x.Comm())
    args = [x.obj,d.obj,orders.obj,firstInds.obj,cutoff]
    if   x.tag == sTag: lib.ElSOCDetsDistMultiVec_s(*args)
    elif x.tag == dTag: lib.ElSOCDetsDistMultiVec_d(*args)
    else: DataExcept()
    return d
  else: TypeExcept()

# Num non-SOC
# ===========
lib.ElNumNonSOC_s.argtypes = \
lib.ElNumNonSOC_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,POINTER(c_int)]
lib.ElNumNonSOCDist_s.argtypes = \
lib.ElNumNonSOCDist_d.argtypes = \
lib.ElNumNonSOCDistMultiVec_s.argtypes = \
lib.ElNumNonSOCDistMultiVec_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_int,POINTER(c_int)]

def NumNonSOC(x,orders,firstInds,cutoff=1000):
  # TODO: Sanity checking
  numNonSOC = iType()
  if type(x) is Matrix:
    args = [x.obj,orders.obj,firstInds.obj,pointer(numNonSOC)]
    if   x.tag == sTag: lib.ElNumNonSOC_s(*args)
    elif x.tag == dTag: lib.ElNumNonSOC_d(*args)
    else: DataExcept()
  elif type(x) is DistMatrix: 
    args = [x.obj,orders.obj,firstInds.obj,cutoff,pointer(numNonSOC)]
    if   x.tag == sTag: lib.ElNumNonSOCDist_s(*args)
    elif x.tag == dTag: lib.ElNumNonSOCDist_d(*args)
    else: DataExcept()
  elif type(x) is DistMultiVec:
    args = [x.obj,orders.obj,firstInds.obj,cutoff,pointer(numNonSOC)]
    if   x.tag == sTag: lib.ElNumNonSOCDistMultiVec_s(*args)
    elif x.tag == dTag: lib.ElNumNonSOCDistMultiVec_d(*args)
    else: DataExcept()
  else: TypeExcept()
  return numNonSOC.value

# SOC Apply
# =========
lib.ElSOCApply_s.argtypes = \
lib.ElSOCApply_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSOCApplyDist_s.argtypes = \
lib.ElSOCApplyDist_d.argtypes = \
lib.ElSOCApplyDistMultiVec_s.argtypes = \
lib.ElSOCApplyDistMultiVec_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_int]

def SOCApply(x,y,orders,firstInds,cutoff=1000):
  # TODO: Sanity checking
  if type(x) is Matrix:
    z = Matrix(x.tag)
    args = [x.obj,y.obj,z.obj,orders.obj,firstInds.obj]
    if   x.tag == sTag: lib.ElSOCApply_s(*args)
    elif x.tag == dTag: lib.ElSOCApply_d(*args)
    else: DataExcept()
    return z
  elif type(x) is DistMatrix:
    z = DistMatrix(x.tag,VC,STAR,x.Grid())
    args = [x.obj,y.obj,z.obj,orders.obj,firstInds.obj,cutoff]
    if   x.tag == sTag: lib.ElSOCApplyDist_s(*args)
    elif x.tag == dTag: lib.ElSOCApplyDist_d(*args)
    else: DataExcept()
    return z
  elif type(x) is DistMultiVec:
    z = DistMultiVec(x.tag,x.Comm())
    args = [x.obj,y.obj,z.obj,orders.obj,firstInds.obj,cutoff]
    if   x.tag == sTag: lib.ElSOCApplyDistMultiVec_s(*args)
    elif x.tag == dTag: lib.ElSOCApplyDistMultiVec_d(*args)
    else: DataExcept()
    return z
  else: TypeExcept()

# SOC Apply quadratic
# ===================
lib.ElSOCApplyQuadratic_s.argtypes = \
lib.ElSOCApplyQuadratic_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSOCApplyQuadraticDist_s.argtypes = \
lib.ElSOCApplyQuadraticDist_d.argtypes = \
lib.ElSOCApplyQuadraticDistMultiVec_s.argtypes = \
lib.ElSOCApplyQuadraticDistMultiVec_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_int]

def SOCApplyQuadratic(x,y,orders,firstInds,cutoff=1000):
  # TODO: Sanity checking
  if type(x) is Matrix:
    z = Matrix(x.tag)
    args = [x.obj,y.obj,z.obj,orders.obj,firstInds.obj]
    if   x.tag == sTag: lib.ElSOCApplyQuadratic_s(*args)
    elif x.tag == dTag: lib.ElSOCApplyQuadratic_d(*args)
    else: DataExcept()
    return z
  elif type(x) is DistMatrix:
    z = DistMatrix(x.tag,VC,STAR,x.Grid())
    args = [x.obj,y.obj,z.obj,orders.obj,firstInds.obj,cutoff]
    if   x.tag == sTag: lib.ElSOCApplyQuadraticDist_s(*args)
    elif x.tag == dTag: lib.ElSOCApplyQuadraticDist_d(*args)
    else: DataExcept()
    return z
  elif type(x) is DistMultiVec:
    z = DistMultiVec(x.tag,x.Comm())
    args = [x.obj,y.obj,z.obj,orders.obj,firstInds.obj,cutoff]
    if   x.tag == sTag: lib.ElSOCApplyQuadraticDistMultiVec_s(*args)
    elif x.tag == dTag: lib.ElSOCApplyQuadraticDistMultiVec_d(*args)
    else: DataExcept()
    return z
  else: TypeExcept()

# SOC Inverse
# ===========
lib.ElSOCInverse_s.argtypes = \
lib.ElSOCInverse_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSOCInverseDist_s.argtypes = \
lib.ElSOCInverseDist_d.argtypes = \
lib.ElSOCInverseDistMultiVec_s.argtypes = \
lib.ElSOCInverseDistMultiVec_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_int]

def SOCInverse(x,orders,firstInds,cutoff=1000):
  # TODO: Sanity checking
  if type(x) is Matrix:
    xInv = Matrix(x.tag)
    args = [x.obj,xInv.obj,orders.obj,firstInds.obj]
    if   x.tag == sTag: lib.ElSOCInverse_s(*args)
    elif x.tag == dTag: lib.ElSOCInverse_d(*args)
    else: DataExcept()
    return xInv
  elif type(x) is DistMatrix:
    xInv = DistMatrix(x.tag,VC,STAR,x.Grid())
    args = [x.obj,xInv.obj,orders.obj,firstInds.obj,cutoff]
    if   x.tag == sTag: lib.ElSOCInverseDist_s(*args)
    elif x.tag == dTag: lib.ElSOCInverseDist_d(*args)
    else: DataExcept()
    return xInv
  elif type(x) is DistMultiVec:
    xInv = DistMultiVec(x.tag,x.Comm())
    args = [x.obj,xInv.obj,orders.obj,firstInds.obj,cutoff]
    if   x.tag == sTag: lib.ElSOCInverseDistMultiVec_s(*args)
    elif x.tag == dTag: lib.ElSOCInverseDistMultiVec_d(*args)
    else: DataExcept()
    return xInv
  else: TypeExcept()

# SOC Square-root
# ===============
lib.ElSOCSquareRoot_s.argtypes = \
lib.ElSOCSquareRoot_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSOCSquareRootDist_s.argtypes = \
lib.ElSOCSquareRootDist_d.argtypes = \
lib.ElSOCSquareRootDistMultiVec_s.argtypes = \
lib.ElSOCSquareRootDistMultiVec_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_int]

def SOCSquareRoot(x,orders,firstInds,cutoff=1000):
  # TODO: Sanity checking
  if type(x) is Matrix:
    xRoot = Matrix(x.tag)
    args = [x.obj,xRoot.obj,orders.obj,firstInds.obj]
    if   x.tag == sTag: lib.ElSOCSquareRoot_s(*args)
    elif x.tag == dTag: lib.ElSOCSquareRoot_d(*args)
    else: DataExcept()
    return xRoot
  elif type(x) is DistMatrix:
    xRoot = DistMatrix(x.tag,VC,STAR,x.Grid())
    args = [x.obj,xRoot.obj,orders.obj,firstInds.obj,cutoff]
    if   x.tag == sTag: lib.ElSOCSquareRootDist_s(*args)
    elif x.tag == dTag: lib.ElSOCSquareRootDist_d(*args)
    else: DataExcept()
    return xRoot
  elif type(x) is DistMultiVec:
    xRoot = DistMultiVec(x.tag,x.Comm())
    args = [x.obj,xRoot.obj,orders.obj,firstInds.obj,cutoff]
    if   x.tag == sTag: lib.ElSOCSquareRootDistMultiVec_s(*args)
    elif x.tag == dTag: lib.ElSOCSquareRootDistMultiVec_d(*args)
    else: DataExcept()
    return xRoot
  else: TypeExcept()

# SOC Nesterov-Todd
# =================
lib.ElSOCNesterovTodd_s.argtypes = \
lib.ElSOCNesterovTodd_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSOCNesterovToddDist_s.argtypes = \
lib.ElSOCNesterovToddDist_d.argtypes = \
lib.ElSOCNesterovToddDistMultiVec_s.argtypes = \
lib.ElSOCNesterovToddDistMultiVec_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_int]

def SOCNesterovTodd(s,z,orders,firstInds,cutoff=1000):
  # TODO: Sanity checking
  if type(s) is Matrix:
    w = Matrix(s.tag)
    args = [s.obj,z.obj,w.obj,orders.obj,firstInds.obj]
    if   s.tag == sTag: lib.ElSOCNesterovTodd_s(*args)
    elif s.tag == dTag: lib.ElSOCNesterovTodd_d(*args)
    else: DataExcept()
    return w
  elif type(s) is DistMatrix:
    w = DistMatrix(s.tag,VC,STAR,s.Grid())
    args = [s.obj,z.obj,w.obj,orders.obj,firstInds.obj,cutoff]
    if   s.tag == sTag: lib.ElSOCNesterovToddDist_s(*args)
    elif s.tag == dTag: lib.ElSOCNesterovToddDist_d(*args)
    else: DataExcept()
    return w
  elif type(s) is DistMultiVec:
    w = DistMultiVec(s.tag,s.Comm())
    args = [s.obj,z.obj,w.obj,orders.obj,firstInds.obj,cutoff]
    if   s.tag == sTag: lib.ElSOCNesterovToddDistMultiVec_s(*args)
    elif s.tag == dTag: lib.ElSOCNesterovToddDistMultiVec_d(*args)
    else: DataExcept()
    return w
  else: TypeExcept()

# Max step in SOC
# ===============
lib.ElMaxStepInSOC_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,sType,POINTER(sType)]
lib.ElMaxStepInSOC_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,dType,POINTER(dType)]
lib.ElMaxStepInSOCDist_s.argtypes = \
lib.ElMaxStepInSOCDistMultiVec_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,sType,iType,POINTER(sType)]
lib.ElMaxStepInSOCDist_d.argtypes = \
lib.ElMaxStepInSOCDistMultiVec_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,dType,iType,POINTER(dType)]

def MaxStepInSOC(x,y,orders,firstInds,upperBound,cutoff=1000):
  # TODO: Sanity checking
  alpha = TagToType(x.tag)()
  if type(x) is Matrix:
    args = [x.obj,y.obj,orders.obj,firstInds.obj,upperBound,pointer(alpha)]
    if   x.tag == sTag: lib.ElMaxStepInSOC_s(*args)
    elif x.tag == dTag: lib.ElMaxStepInSOC_d(*args)
    else: DataExcept()
  elif type(x) is DistMatrix:
    args = [x.obj,y.obj,orders.obj,firstInds.obj,
            upperBound,cutoff,pointer(alpha)]
    if   x.tag == sTag: lib.ElMaxStepInSOCDist_s(*args)
    elif x.tag == dTag: lib.ElMaxStepInSOCDist_d(*args)
    else: DataExcept()
  elif type(x) is DistMultiVec:
    args = [x.obj,y.obj,orders.obj,firstInds.obj,
            upperBound,cutoff,pointer(alpha)]
    if   x.tag == sTag: lib.ElMaxStepInSOCDistMultiVec_s(*args)
    elif x.tag == dTag: lib.ElMaxStepInSOCDistMultiVec_d(*args)
    else: DataExcept()
  else: TypeExcept()
  return alpha.value
