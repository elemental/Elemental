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

def Copy(A,B):
  CheckTag(A.tag)
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

def Transpose(A,B):
  CheckTag(A.tag)
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

def Adjoint(A,B):
  CheckTag(A.tag)
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
