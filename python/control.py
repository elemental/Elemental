#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *

lib.ElLyapunov_s.argtypes = \
lib.ElLyapunov_d.argtypes = \
lib.ElLyapunov_c.argtypes = \
lib.ElLyapunov_z.argtypes = \
lib.ElLyapunovDist_s.argtypes = \
lib.ElLyapunovDist_d.argtypes = \
lib.ElLyapunovDist_c.argtypes = \
lib.ElLyapunovDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p]

def Lyapunov(A,C):
  if type(A) is not type(C):
    raise Exception('Matrix types must match')
  if A.tag != C.tag:
    raise Exception('Datatypes must match')
  if type(A) is Matrix:
    X = Matrix(A.tag)
    args = [A.obj,C.obj,X.obj]
    if   A.tag == sTag: lib.ElLyapunov_s(*args)
    elif A.tag == dTag: lib.ElLyapunov_d(*args)
    elif A.tag == cTag: lib.ElLyapunov_c(*args)
    elif A.tag == zTag: lib.ElLyapunov_z(*args)
    else: DataExcept()
    return X
  elif type(A) is DistMatrix:
    X = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,C.obj,X.obj]
    if   A.tag == sTag: lib.ElLyapunovDist_s(*args)
    elif A.tag == dTag: lib.ElLyapunovDist_d(*args)
    elif A.tag == cTag: lib.ElLyapunovDist_c(*args)
    elif A.tag == zTag: lib.ElLyapunovDist_z(*args)
    else: DataExcept()
    return X
  else: TypeExcept()

lib.ElRiccati_s.argtypes = \
lib.ElRiccati_d.argtypes = \
lib.ElRiccati_c.argtypes = \
lib.ElRiccati_z.argtypes = \
lib.ElRiccatiDist_s.argtypes = \
lib.ElRiccatiDist_d.argtypes = \
lib.ElRiccatiDist_c.argtypes = \
lib.ElRiccatiDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]

def Riccati(uplo,A,K,L):
  if type(A) is Matrix:
    X = Matrix(A.tag)
    args = [uplo,A.obj,K.obj,L.obj,X.obj]
    if   A.tag == sTag: lib.ElRiccati_s(*args)
    elif A.tag == dTag: lib.ElRiccati_d(*args)
    elif A.tag == cTag: lib.ElRiccati_c(*args)
    elif A.tag == zTag: lib.ElRiccati_z(*args)
    else: DataExcept()
    return X
  elif type(A) is DistMatrix:
    X = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [uplo,A.obj,K.obj,L.obj,X.obj]
    if   A.tag == sTag: lib.ElRiccatiDist_s(*args)
    elif A.tag == dTag: lib.ElRiccatiDist_d(*args)
    elif A.tag == cTag: lib.ElRiccatiDist_c(*args)
    elif A.tag == zTag: lib.ElRiccatiDist_z(*args)
    else: DataExcept()
    return X
  else: TypeExcept()

lib.ElRiccatiPreformed_s.argtypes = \
lib.ElRiccatiPreformed_d.argtypes = \
lib.ElRiccatiPreformed_c.argtypes = \
lib.ElRiccatiPreformed_z.argtypes = \
lib.ElRiccatiPreformedDist_s.argtypes = \
lib.ElRiccatiPreformedDist_d.argtypes = \
lib.ElRiccatiPreformedDist_c.argtypes = \
lib.ElRiccatiPreformedDist_z.argtypes = \
  [c_void_p,c_void_p]

def RiccatiPreformed(W):
  if type(W) is Matrix:
    X = Matrix(W.tag)
    args = [W.obj,X.obj]
    if   W.tag == sTag: lib.ElRiccatiPreformed_s(*args)
    elif W.tag == dTag: lib.ElRiccatiPreformed_d(*args)
    elif W.tag == cTag: lib.ElRiccatiPreformed_c(*args)
    elif W.tag == zTag: lib.ElRiccatiPreformed_z(*args)
    else: DataExcept()
    return X
  elif type(W) is DistMatrix:
    X = DistMatrix(W.tag,MC,MR,W.Grid())
    args = [W.obj,X.obj]
    if   W.tag == sTag: lib.ElRiccatiPreformedDist_s(*args)
    elif W.tag == dTag: lib.ElRiccatiPreformedDist_d(*args)
    elif W.tag == cTag: lib.ElRiccatiPreformedDist_c(*args)
    elif W.tag == zTag: lib.ElRiccatiPreformedDist_z(*args)
    else: DataExcept()
    return X
  else: TypeExcept()

lib.ElSylvester_s.argtypes = \
lib.ElSylvester_d.argtypes = \
lib.ElSylvester_c.argtypes = \
lib.ElSylvester_z.argtypes = \
lib.ElSylvesterDist_s.argtypes = \
lib.ElSylvesterDist_d.argtypes = \
lib.ElSylvesterDist_c.argtypes = \
lib.ElSylvesterDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p]

def Sylvester(A,B,C):
  if type(A) is not type(B) or type(B) is not type(C):
    raise Exception('Matrix types of {A,B,C} must match')
  if A.tag != B.tag or B.tag != C.tag:
    raise Exception('Datatypes of {A,B,C} must match')
  if type(A) is Matrix:
    X = Matrix(A.tag)
    args = [A.obj,B.obj,C.obj,X.obj]
    if   A.tag == sTag: lib.ElSylvester_s(*args)
    elif A.tag == dTag: lib.ElSylvester_d(*args)
    elif A.tag == cTag: lib.ElSylvester_c(*args)
    elif A.tag == zTag: lib.ElSylvester_z(*args)
    else: DataExcept()
    return X
  elif type(A) is DistMatrix:
    X = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,B.obj,C.obj,X.obj]
    if   A.tag == sTag: lib.ElSylvesterDist_s(*args)
    elif A.tag == dTag: lib.ElSylvesterDist_d(*args)
    elif A.tag == cTag: lib.ElSylvesterDist_c(*args)
    elif A.tag == zTag: lib.ElSylvesterDist_z(*args)
    else: DataExcept()
    return X
  else: TypeExcept()

lib.ElSylvesterPreformed_s.argtypes = \
lib.ElSylvesterPreformed_d.argtypes = \
lib.ElSylvesterPreformed_c.argtypes = \
lib.ElSylvesterPreformed_z.argtypes = \
lib.ElSylvesterPreformedDist_s.argtypes = \
lib.ElSylvesterPreformedDist_d.argtypes = \
lib.ElSylvesterPreformedDist_c.argtypes = \
lib.ElSylvesterPreformedDist_z.argtypes = \
  [iType,c_void_p,c_void_p]

def SylvesterPreformed(m,W):
  if type(W) is Matrix:
    X = Matrix(W.tag)
    args = [m,W.obj,X.obj]
    if   W.tag == sTag: lib.ElSylvesterPreformed_s(*args)
    elif W.tag == dTag: lib.ElSylvesterPreformed_d(*args)
    elif W.tag == cTag: lib.ElSylvesterPreformed_c(*args)
    elif W.tag == zTag: lib.ElSylvesterPreformed_z(*args)
    else: DataExcept()
    return X
  elif type(W) is DistMatrix:
    X = DistMatrix(W.tag,MC,MR,W.Grid())
    args = [m,W.obj,X.obj]
    if   W.tag == sTag: lib.ElSylvesterPreformedDist_s(*args)
    elif W.tag == dTag: lib.ElSylvesterPreformedDist_d(*args)
    elif W.tag == cTag: lib.ElSylvesterPreformedDist_c(*args)
    elif W.tag == zTag: lib.ElSylvesterPreformedDist_z(*args)
    else: DataExcept()
    return X
  else: TypeExcept()
