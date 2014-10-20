#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *

lib.ElLyapunov_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLyapunov_s.restype = c_uint
lib.ElLyapunov_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLyapunov_d.restype = c_uint
lib.ElLyapunov_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLyapunov_c.restype = c_uint
lib.ElLyapunov_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLyapunov_z.restype = c_uint
lib.ElLyapunovDist_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLyapunovDist_s.restype = c_uint
lib.ElLyapunovDist_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLyapunovDist_d.restype = c_uint
lib.ElLyapunovDist_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLyapunovDist_c.restype = c_uint
lib.ElLyapunovDist_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLyapunovDist_z.restype = c_uint
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

lib.ElRicatti_s.argtypes = [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElRicatti_s.restype = c_uint
lib.ElRicatti_d.argtypes = [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElRicatti_d.restype = c_uint
lib.ElRicatti_c.argtypes = [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElRicatti_c.restype = c_uint
lib.ElRicatti_z.argtypes = [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElRicatti_z.restype = c_uint
lib.ElRicattiDist_s.argtypes = [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElRicattiDist_s.restype = c_uint
lib.ElRicattiDist_d.argtypes = [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElRicattiDist_d.restype = c_uint
lib.ElRicattiDist_c.argtypes = [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElRicattiDist_c.restype = c_uint
lib.ElRicattiDist_z.argtypes = [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElRicattiDist_z.restype = c_uint
def Ricatti(uplo,A,K,L):
  if type(A) is Matrix:
    X = Matrix(A.tag)
    args = [uplo,A.obj,K.obj,L.obj,X.obj]
    if   A.tag == sTag: lib.ElRicatti_s(*args)
    elif A.tag == dTag: lib.ElRicatti_d(*args)
    elif A.tag == cTag: lib.ElRicatti_c(*args)
    elif A.tag == zTag: lib.ElRicatti_z(*args)
    else: DataExcept()
    return X
  elif type(A) is DistMatrix:
    X = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [uplo,A.obj,K.obj,L.obj,X.obj]
    if   A.tag == sTag: lib.ElRicattiDist_s(*args)
    elif A.tag == dTag: lib.ElRicattiDist_d(*args)
    elif A.tag == cTag: lib.ElRicattiDist_c(*args)
    elif A.tag == zTag: lib.ElRicattiDist_z(*args)
    else: DataExcept()
    return X
  else: TypeExcept()

lib.ElRicattiPreformed_s.argtypes = [c_void_p,c_void_p]
lib.ElRicattiPreformed_s.restype = c_uint
lib.ElRicattiPreformed_d.argtypes = [c_void_p,c_void_p]
lib.ElRicattiPreformed_d.restype = c_uint
lib.ElRicattiPreformed_c.argtypes = [c_void_p,c_void_p]
lib.ElRicattiPreformed_c.restype = c_uint
lib.ElRicattiPreformed_z.argtypes = [c_void_p,c_void_p]
lib.ElRicattiPreformed_z.restype = c_uint
lib.ElRicattiPreformedDist_s.argtypes = [c_void_p,c_void_p]
lib.ElRicattiPreformedDist_s.restype = c_uint
lib.ElRicattiPreformedDist_d.argtypes = [c_void_p,c_void_p]
lib.ElRicattiPreformedDist_d.restype = c_uint
lib.ElRicattiPreformedDist_c.argtypes = [c_void_p,c_void_p]
lib.ElRicattiPreformedDist_c.restype = c_uint
lib.ElRicattiPreformedDist_z.argtypes = [c_void_p,c_void_p]
lib.ElRicattiPreformedDist_z.restype = c_uint
def RicattiPreformed(W):
  if type(W) is Matrix:
    X = Matrix(W.tag)
    args = [W.obj,X.obj]
    if   W.tag == sTag: lib.ElRicattiPreformed_s(*args)
    elif W.tag == dTag: lib.ElRicattiPreformed_d(*args)
    elif W.tag == cTag: lib.ElRicattiPreformed_c(*args)
    elif W.tag == zTag: lib.ElRicattiPreformed_z(*args)
    else: DataExcept()
    return X
  elif type(W) is DistMatrix:
    X = DistMatrix(W.tag,MC,MR,W.Grid())
    args = [W.obj,X.obj]
    if   W.tag == sTag: lib.ElRicattiPreformedDist_s(*args)
    elif W.tag == dTag: lib.ElRicattiPreformedDist_d(*args)
    elif W.tag == cTag: lib.ElRicattiPreformedDist_c(*args)
    elif W.tag == zTag: lib.ElRicattiPreformedDist_z(*args)
    else: DataExcept()
    return X
  else: TypeExcept()

lib.ElSylvester_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSylvester_s.restype = c_uint
lib.ElSylvester_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSylvester_d.restype = c_uint
lib.ElSylvester_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSylvester_c.restype = c_uint
lib.ElSylvester_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSylvester_z.restype = c_uint
lib.ElSylvesterDist_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSylvesterDist_s.restype = c_uint
lib.ElSylvesterDist_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSylvesterDist_d.restype = c_uint
lib.ElSylvesterDist_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSylvesterDist_c.restype = c_uint
lib.ElSylvesterDist_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSylvesterDist_z.restype = c_uint
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

lib.ElSylvesterPreformed_s.argtypes = [iType,c_void_p,c_void_p]
lib.ElSylvesterPreformed_s.restype = c_uint
lib.ElSylvesterPreformed_d.argtypes = [iType,c_void_p,c_void_p]
lib.ElSylvesterPreformed_d.restype = c_uint
lib.ElSylvesterPreformed_c.argtypes = [iType,c_void_p,c_void_p]
lib.ElSylvesterPreformed_c.restype = c_uint
lib.ElSylvesterPreformed_z.argtypes = [iType,c_void_p,c_void_p]
lib.ElSylvesterPreformed_z.restype = c_uint
lib.ElSylvesterPreformedDist_s.argtypes = [iType,c_void_p,c_void_p]
lib.ElSylvesterPreformedDist_s.restype = c_uint
lib.ElSylvesterPreformedDist_d.argtypes = [iType,c_void_p,c_void_p]
lib.ElSylvesterPreformedDist_d.restype = c_uint
lib.ElSylvesterPreformedDist_c.argtypes = [iType,c_void_p,c_void_p]
lib.ElSylvesterPreformedDist_c.restype = c_uint
lib.ElSylvesterPreformedDist_z.argtypes = [iType,c_void_p,c_void_p]
lib.ElSylvesterPreformedDist_z.restype = c_uint
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
