#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from environment import *

# Basic element manipulation
# ==========================

# Return the complex argument of a scalar
# ---------------------------------------
lib.ElArg_s.argtypes = [sType,POINTER(sType)]
lib.ElArg_d.argtypes = [dType,POINTER(dType)]
lib.ElArg_c.argtypes = [cType,POINTER(sType)]
lib.ElArg_z.argtypes = [zType,POINTER(dType)]
def Arg(alpha):
  if   type(alpha) is sType:
    result = sType()
    lib.ElArg_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElArg_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = sType()
    lib.ElArg_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = dType()
    lib.ElArg_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

# Construct a complex number from its polar coordinates
# -----------------------------------------------------
lib.ElComplexFromPolar_c.argtypes = [sType,sType,POINTER(cType)]
lib.ElComplexFromPolar_z.argtyped = [dType,dType,POINTER(zType)]
def ComplexFromPolar(r,theta):
  if   type(r) is sType:
    result = cType()
    lib.ElComplexFromPolar_c(r,theta,pointer(result))
    return result
  elif type(r) is dType:
    result = zType()
    lib.ElComplexFromPolar_z(r,theta,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

# Magnitude and sign
# ==================

lib.ElAbs_i.argtypes = [iType,POINTER(iType)]
lib.ElAbs_s.argtypes = [sType,POINTER(sType)]
lib.ElAbs_d.argtypes = [dType,POINTER(dType)]
lib.ElAbs_c.argtypes = [cType,POINTER(sType)]
lib.ElAbs_z.argtypes = [zType,POINTER(dType)]
def Abs(alpha):
  if   type(alpha) is iType:
    result = iType()
    lib.ElAbs_i(alpha,pointer(result))
    return result
  elif type(alpha) is sType:
    result = sType()
    lib.ElAbs_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElAbs_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = sType()
    lib.ElAbs_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = dType()
    lib.ElAbs_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

lib.ElSafeAbs_c.argtypes = [cType,POINTER(sType)]
lib.ElSafeAbs_z.argtypes = [zType,POINTER(dType)]
def SafeAbs(alpha):
  if   type(alpha) is iType:
    result = iType()
    lib.ElAbs_i(alpha,pointer(result))
    return result
  elif type(alpha) is sType:
    result = sType()
    lib.ElAbs_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElAbs_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = sType()
    lib.ElSafeAbs_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = dType()
    lib.ElSafeAbs_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

lib.ElFastAbs_c.argtypes = [cType,POINTER(sType)]
lib.ElFastAbs_z.argtypes = [zType,POINTER(dType)]
def FastAbs(alpha):
  if   type(alpha) is iType:
    result = iType()
    lib.ElAbs_i(alpha,pointer(result))
    return result
  elif type(alpha) is sType:
    result = sType()
    lib.ElAbs_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElAbs_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = sType()
    lib.ElFastAbs_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = dType()
    lib.ElFastAbs_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

lib.ElSgn_i.argtypes = [iType,bType,POINTER(iType)]
lib.ElSgn_s.argtypes = [sType,bType,POINTER(sType)]
lib.ElSgn_d.argtypes = [dType,bType,POINTER(dType)]
def Sgn(alpha,symm=True):
  if   type(alpha) is iType:
    result = iType()
    lib.ElSgn_i(alpha,symm,pointer(result))
    return result
  elif type(alpha) is sType:
    result = sType()
    lib.ElSgn_s(alpha,symm,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElSgn_d(alpha,symm,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

# Exponentiation
# ==============
lib.ElExp_s.argtypes = [sType,POINTER(sType)]
lib.ElExp_d.argtypes = [dType,POINTER(dType)]
lib.ElExp_c.argtypes = [cType,POINTER(cType)]
lib.ElExp_z.argtypes = [zType,POINTER(zType)]
def Exp(alpha):
  if   type(alpha) is sType:
    result = sType()
    lib.ElExp_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElExp_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = cType()
    lib.ElExp_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = zType()
    lib.ElExp_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

lib.ElPow_s.argtypes = [sType,sType,POINTER(sType)]
lib.ElPow_d.argtypes = [dType,dType,POINTER(dType)]
lib.ElPow_c.argtypes = [cType,cType,POINTER(cType)]
lib.ElPow_z.argtypes = [zType,zType,POINTER(zType)]
def Pow(alpha,beta):
  if   type(alpha) is sType:
    result = sType()
    lib.ElPow_s(alpha,beta,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElPow_d(alpha,beta,pointer(result))
    return result
  elif type(alpha) is cType:
    result = cType()
    lib.ElPow_c(alpha,beta,pointer(result))
    return result
  elif type(alpha) is zType:
    result = zType()
    lib.ElPow_z(alpha,beta,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

lib.ElLog_s.argtypes = [sType,POINTER(sType)]
lib.ElLog_d.argtypes = [dType,POINTER(dType)]
lib.ElLog_c.argtypes = [cType,POINTER(cType)]
lib.ElLog_z.argtypes = [zType,POINTER(zType)]
def Log(alpha):
  if   type(alpha) is sType:
    result = sType()
    lib.ElLog_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElLog_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = cType()
    lib.ElLog_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = zType()
    lib.ElLog_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

lib.ElSqrt_s.argtypes = [sType,POINTER(sType)]
lib.ElSqrt_d.argtypes = [dType,POINTER(dType)]
lib.ElSqrt_c.argtypes = [cType,POINTER(cType)]
lib.ElSqrt_z.argtypes = [zType,POINTER(zType)]
def Sqrt(alpha):
  if   type(alpha) is sType:
    result = sType()
    lib.ElSqrt_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElSqrt_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = cType()
    lib.ElSqrt_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = zType()
    lib.ElSqrt_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

# Trigonometric functions
# =======================
lib.ElCos_s.argtypes = [sType,POINTER(sType)]
lib.ElCos_d.argtypes = [dType,POINTER(dType)]
lib.ElCos_c.argtypes = [cType,POINTER(cType)]
lib.ElCos_z.argtypes = [zType,POINTER(zType)]
def Cos(alpha):
  if   type(alpha) is sType:
    result = sType()
    lib.ElCos_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElCos_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = cType()
    lib.ElCos_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = zType()
    lib.ElCos_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

lib.ElSin_s.argtypes = [sType,POINTER(sType)]
lib.ElSin_d.argtypes = [dType,POINTER(dType)]
lib.ElSin_c.argtypes = [cType,POINTER(cType)]
lib.ElSin_z.argtypes = [zType,POINTER(zType)]
def Sin(alpha):
  if   type(alpha) is sType:
    result = sType()
    lib.ElSin_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElSin_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = cType()
    lib.ElSin_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = zType()
    lib.ElSin_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

lib.ElTan_s.argtypes = [sType,POINTER(sType)]
lib.ElTan_d.argtypes = [dType,POINTER(dType)]
lib.ElTan_c.argtypes = [cType,POINTER(cType)]
lib.ElTan_z.argtypes = [zType,POINTER(zType)]
def Tan(alpha):
  if   type(alpha) is sType:
    result = sType()
    lib.ElTan_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElTan_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = cType()
    lib.ElTan_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = zType()
    lib.ElTan_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

lib.ElAcos_s.argtypes = [sType,POINTER(sType)]
lib.ElAcos_d.argtypes = [dType,POINTER(dType)]
lib.ElAcos_c.argtypes = [cType,POINTER(cType)]
lib.ElAcos_z.argtypes = [zType,POINTER(zType)]
def Acos(alpha):
  if   type(alpha) is sType:
    result = sType()
    lib.ElAcos_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElAcos_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = cType()
    lib.ElAcos_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = zType()
    lib.ElAcos_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

lib.ElAsin_s.argtypes = [sType,POINTER(sType)]
lib.ElAsin_d.argtypes = [dType,POINTER(dType)]
lib.ElAsin_c.argtypes = [cType,POINTER(cType)]
lib.ElAsin_z.argtypes = [zType,POINTER(zType)]
def Asin(alpha):
  if   type(alpha) is sType:
    result = sType()
    lib.ElAsin_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElAsin_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = cType()
    lib.ElAsin_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = zType()
    lib.ElAsin_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

lib.ElAtan_s.argtypes = [sType,POINTER(sType)]
lib.ElAtan_d.argtypes = [dType,POINTER(dType)]
lib.ElAtan_c.argtypes = [cType,POINTER(cType)]
lib.ElAtan_z.argtypes = [zType,POINTER(zType)]
def Atan(alpha):
  if   type(alpha) is sType:
    result = sType()
    lib.ElAtan_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElAtan_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = cType()
    lib.ElAtan_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = zType()
    lib.ElAtan_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

lib.ElAtan2_s.argtypes = [sType,sType,POINTER(sType)]
lib.ElAtan2_d.argtypes = [dType,dType,POINTER(dType)]
def Atan2(y,x):
  if   type(y) is sType:
    result = sType()
    lib.ElAtan2_s(y,x,pointer(result))
    return result
  elif type(y) is dType:
    result = dType()
    lib.ElAtan2_d(y,x,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

# Hyperbolic functions
# ====================
lib.ElCosh_s.argtypes = [sType,POINTER(sType)]
lib.ElCosh_d.argtypes = [dType,POINTER(dType)]
lib.ElCosh_c.argtypes = [cType,POINTER(cType)]
lib.ElCosh_z.argtypes = [zType,POINTER(zType)]
def Cosh(alpha):
  if   type(alpha) is sType:
    result = sType()
    lib.ElCosh_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElCosh_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = cType()
    lib.ElCosh_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = zType()
    lib.ElCosh_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

lib.ElSinh_s.argtypes = [sType,POINTER(sType)]
lib.ElSinh_d.argtypes = [dType,POINTER(dType)]
lib.ElSinh_c.argtypes = [cType,POINTER(cType)]
lib.ElSinh_z.argtypes = [zType,POINTER(zType)]
def Sinh(alpha):
  if   type(alpha) is sType:
    result = sType()
    lib.ElSinh_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElSinh_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = cType()
    lib.ElSinh_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = zType()
    lib.ElSinh_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

lib.ElTanh_s.argtypes = [sType,POINTER(sType)]
lib.ElTanh_d.argtypes = [dType,POINTER(dType)]
lib.ElTanh_c.argtypes = [cType,POINTER(cType)]
lib.ElTanh_z.argtypes = [zType,POINTER(zType)]
def Tanh(alpha):
  if   type(alpha) is sType:
    result = sType()
    lib.ElTanh_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElTanh_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = cType()
    lib.ElTanh_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = zType()
    lib.ElTanh_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

lib.ElAcosh_s.argtypes = [sType,POINTER(sType)]
lib.ElAcosh_d.argtypes = [dType,POINTER(dType)]
lib.ElAcosh_c.argtypes = [cType,POINTER(cType)]
lib.ElAcosh_z.argtypes = [zType,POINTER(zType)]
def Acosh(alpha):
  if   type(alpha) is sType:
    result = sType()
    lib.ElAcosh_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElAcosh_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = cType()
    lib.ElAcosh_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = zType()
    lib.ElAcosh_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

lib.ElAsinh_s.argtypes = [sType,POINTER(sType)]
lib.ElAsinh_d.argtypes = [dType,POINTER(dType)]
lib.ElAsinh_c.argtypes = [cType,POINTER(cType)]
lib.ElAsinh_z.argtypes = [zType,POINTER(zType)]
def Asinh(alpha):
  if   type(alpha) is sType:
    result = sType()
    lib.ElAsinh_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElAsinh_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = cType()
    lib.ElAsinh_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = zType()
    lib.ElAsinh_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')

lib.ElAtanh_s.argtypes = [sType,POINTER(sType)]
lib.ElAtanh_d.argtypes = [dType,POINTER(dType)]
lib.ElAtanh_c.argtypes = [cType,POINTER(cType)]
lib.ElAtanh_z.argtypes = [zType,POINTER(zType)]
def Atanh(alpha):
  if   type(alpha) is sType:
    result = sType()
    lib.ElAtanh_s(alpha,pointer(result))
    return result
  elif type(alpha) is dType:
    result = dType()
    lib.ElAtanh_d(alpha,pointer(result))
    return result
  elif type(alpha) is cType:
    result = cType()
    lib.ElAtanh_c(alpha,pointer(result))
    return result
  elif type(alpha) is zType:
    result = zType()
    lib.ElAtanh_z(alpha,pointer(result))
    return result
  else: raise Exception('Unsupported datatype')
