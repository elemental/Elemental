#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from ..core import *

# BLAS 1
# ======

lib.ElCopy_s.argtypes = [c_void_p,c_void_p]
lib.ElCopy_s.restype = c_uint
lib.ElCopy_d.argtypes = [c_void_p,c_void_p]
lib.ElCopy_d.restype = c_uint
lib.ElCopy_c.argtypes = [c_void_p,c_void_p]
lib.ElCopy_c.restype = c_uint
lib.ElCopy_z.argtypes = [c_void_p,c_void_p]
lib.ElCopy_z.restype = c_uint
lib.ElCopyDist_s.argtypes = [c_void_p,c_void_p]
lib.ElCopyDist_s.restype = c_uint
lib.ElCopyDist_d.argtypes = [c_void_p,c_void_p]
lib.ElCopyDist_d.restype = c_uint
lib.ElCopyDist_c.argtypes = [c_void_p,c_void_p]
lib.ElCopyDist_c.restype = c_uint
lib.ElCopyDist_z.argtypes = [c_void_p,c_void_p]
lib.ElCopyDist_z.restype = c_uint
def Copy(A,B):
  if A.tag != B.tag:
    raise Exception('Copying between datatypes is not yet supported in Python')
  if type(A) is not type(B): raise Exception('Matrix types must match')
  if type(A) is Matrix:
    if   B.tag == iTag: lib.ElCopy_i(A.obj,B.obj)
    elif B.tag == sTag: lib.ElCopy_s(A.obj,B.obj)
    elif B.tag == dTag: lib.ElCopy_d(A.obj,B.obj)
    elif B.tag == cTag: lib.ElCopy_c(A.obj,B.obj)
    elif B.tag == zTag: lib.ElCopy_z(A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   B.tag == iTag: lib.ElCopyDist_i(A.obj,B.obj)
    elif B.tag == sTag: lib.ElCopyDist_s(A.obj,B.obj)
    elif B.tag == dTag: lib.ElCopyDist_d(A.obj,B.obj)
    elif B.tag == cTag: lib.ElCopyDist_c(A.obj,B.obj)
    elif B.tag == zTag: lib.ElCopyDist_z(A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix types')

lib.ElTranspose_s.argtypes = [c_void_p,c_void_p]
lib.ElTranspose_s.restype = c_uint
lib.ElTranspose_d.argtypes = [c_void_p,c_void_p]
lib.ElTranspose_d.restype = c_uint
lib.ElTranspose_c.argtypes = [c_void_p,c_void_p]
lib.ElTranspose_c.restype = c_uint
lib.ElTranspose_z.argtypes = [c_void_p,c_void_p]
lib.ElTranspose_z.restype = c_uint
lib.ElTransposeDist_s.argtypes = [c_void_p,c_void_p]
lib.ElTransposeDist_s.restype = c_uint
lib.ElTransposeDist_d.argtypes = [c_void_p,c_void_p]
lib.ElTransposeDist_d.restype = c_uint
lib.ElTransposeDist_c.argtypes = [c_void_p,c_void_p]
lib.ElTransposeDist_c.restype = c_uint
lib.ElTransposeDist_z.argtypes = [c_void_p,c_void_p]
lib.ElTransposeDist_z.restype = c_uint
def Transpose(A,B):
  if A.tag != B.tag:
    raise Exception('Transposing between datatypes not yet supported in Python')
  if type(A) is not type(B): raise Exception('Matrix types must match')
  if type(A) is Matrix:
    if   B.tag == iTag: lib.ElTranspose_i(A.obj,B.obj)
    elif B.tag == sTag: lib.ElTranspose_s(A.obj,B.obj)
    elif B.tag == dTag: lib.ElTranspose_d(A.obj,B.obj)
    elif B.tag == cTag: lib.ElTranspose_c(A.obj,B.obj)
    elif B.tag == zTag: lib.ElTranspose_z(A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   B.tag == iTag: lib.ElTransposeDist_i(A.obj,B.obj)
    elif B.tag == sTag: lib.ElTransposeDist_s(A.obj,B.obj)
    elif B.tag == dTag: lib.ElTransposeDist_d(A.obj,B.obj)
    elif B.tag == cTag: lib.ElTransposeDist_c(A.obj,B.obj)
    elif B.tag == zTag: lib.ElTransposeDist_z(A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix types')

lib.ElAdjoint_c.argtypes = [c_void_p,c_void_p]
lib.ElAdjoint_c.restype = c_uint
lib.ElAdjoint_z.argtypes = [c_void_p,c_void_p]
lib.ElAdjoint_z.restype = c_uint
lib.ElAdjointDist_c.argtypes = [c_void_p,c_void_p]
lib.ElAdjointDist_c.restype = c_uint
lib.ElAdjointDist_z.argtypes = [c_void_p,c_void_p]
lib.ElAdjointDist_z.restype = c_uint
def Adjoint(A,B):
  if A.tag != B.tag:
    raise Exception('Transposing between datatypes not yet supported in Python')
  if type(A) is not type(B): raise Exception('Matrix types must match')
  if type(A) is Matrix:
    if   B.tag == cTag: lib.ElAdjoint_c(A.obj,B.obj)
    elif B.tag == zTag: lib.ElAdjoint_z(A.obj,B.obj)
    else: Transpose(A,B)
  elif type(A) is DistMatrix:
    if   B.tag == cTag: lib.ElAdjointDist_c(A.obj,B.obj)
    elif B.tag == zTag: lib.ElAdjointDist_z(A.obj,B.obj)
    else: Transpose(A,B)
  else: raise Exception('Unsupported matrix types')

# TODO: Version which returns the result instead?
lib.ElRealPart_i.argtypes = [c_void_p,c_void_p]
lib.ElRealPart_i.restype = c_uint
lib.ElRealPart_s.argtypes = [c_void_p,c_void_p]
lib.ElRealPart_s.restype = c_uint
lib.ElRealPart_d.argtypes = [c_void_p,c_void_p]
lib.ElRealPart_d.restype = c_uint
lib.ElRealPart_c.argtypes = [c_void_p,c_void_p]
lib.ElRealPart_c.restype = c_uint
lib.ElRealPart_z.argtypes = [c_void_p,c_void_p]
lib.ElRealPart_z.restype = c_uint
lib.ElRealPartDist_i.argtypes = [c_void_p,c_void_p]
lib.ElRealPartDist_i.restype = c_uint
lib.ElRealPartDist_s.argtypes = [c_void_p,c_void_p]
lib.ElRealPartDist_s.restype = c_uint
lib.ElRealPartDist_d.argtypes = [c_void_p,c_void_p]
lib.ElRealPartDist_d.restype = c_uint
lib.ElRealPartDist_c.argtypes = [c_void_p,c_void_p]
lib.ElRealPartDist_c.restype = c_uint
lib.ElRealPartDist_z.argtypes = [c_void_p,c_void_p]
lib.ElRealPartDist_z.restype = c_uint
def RealPart(A,AReal):
  if AReal.tag != Base(A.tag):
    raise Exception('AReal must have the base datatype of A')
  if type(A) is not type(AReal): raise Exception('Matrix types must match')
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElRealPart_i(A.obj,AReal.obj)
    elif A.tag == sTag: lib.ElRealPart_s(A.obj,AReal.obj)
    elif A.tag == dTag: lib.ElRealPart_d(A.obj,AReal.obj)
    elif A.tag == cTag: lib.ElRealPart_c(A.obj,AReal.obj)
    elif A.tag == zTag: lib.ElRealPart_z(A.obj,AReal.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElRealPartDist_i(A.obj,AReal.obj)
    elif A.tag == sTag: lib.ElRealPartDist_s(A.obj,AReal.obj)
    elif A.tag == dTag: lib.ElRealPartDist_d(A.obj,AReal.obj)
    elif A.tag == cTag: lib.ElRealPartDist_c(A.obj,AReal.obj)
    elif A.tag == zTag: lib.ElRealPartDist_z(A.obj,AReal.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix types')

# TODO: Version which returns the result instead?
lib.ElImagPart_i.argtypes = [c_void_p,c_void_p]
lib.ElImagPart_i.restype = c_uint
lib.ElImagPart_s.argtypes = [c_void_p,c_void_p]
lib.ElImagPart_s.restype = c_uint
lib.ElImagPart_d.argtypes = [c_void_p,c_void_p]
lib.ElImagPart_d.restype = c_uint
lib.ElImagPart_c.argtypes = [c_void_p,c_void_p]
lib.ElImagPart_c.restype = c_uint
lib.ElImagPart_z.argtypes = [c_void_p,c_void_p]
lib.ElImagPart_z.restype = c_uint
lib.ElImagPartDist_i.argtypes = [c_void_p,c_void_p]
lib.ElImagPartDist_i.restype = c_uint
lib.ElImagPartDist_s.argtypes = [c_void_p,c_void_p]
lib.ElImagPartDist_s.restype = c_uint
lib.ElImagPartDist_d.argtypes = [c_void_p,c_void_p]
lib.ElImagPartDist_d.restype = c_uint
lib.ElImagPartDist_c.argtypes = [c_void_p,c_void_p]
lib.ElImagPartDist_c.restype = c_uint
lib.ElImagPartDist_z.argtypes = [c_void_p,c_void_p]
lib.ElImagPartDist_z.restype = c_uint
def ImagPart(A,AImag):
  if AImag.tag != Base(A.tag):
    raise Exception('AImag must have the base datatype of A')
  if type(A) is not type(AImag): raise Exception('Matrix types must match')
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElImagPart_i(A.obj,AImag.obj)
    elif A.tag == sTag: lib.ElImagPart_s(A.obj,AImag.obj)
    elif A.tag == dTag: lib.ElImagPart_d(A.obj,AImag.obj)
    elif A.tag == cTag: lib.ElImagPart_c(A.obj,AImag.obj)
    elif A.tag == zTag: lib.ElImagPart_z(A.obj,AImag.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElImagPartDist_i(A.obj,AImag.obj)
    elif A.tag == sTag: lib.ElImagPartDist_s(A.obj,AImag.obj)
    elif A.tag == dTag: lib.ElImagPartDist_d(A.obj,AImag.obj)
    elif A.tag == cTag: lib.ElImagPartDist_c(A.obj,AImag.obj)
    elif A.tag == zTag: lib.ElImagPartDist_z(A.obj,AImag.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix types')
