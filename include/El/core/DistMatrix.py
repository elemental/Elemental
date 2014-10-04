#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from environment import *
import ctypes, numpy

from Grid import DefaultGrid

# (Abstract)DistMatrix
# ====================

class DistData(ctypes.Structure):
  _fields_ = [('colDist',ctypes.c_uint),
              ('rowDist',ctypes.c_uint), 
              ('colAlign',iType),
              ('rowAlign',iType),
              ('root',iType),
              ('grid',ctypes.c_void_p)]

class DistMatrix(object):
  def __init__(self,tag=dTag,colDist=MC,rowDist=MR,grid=DefaultGrid()):
    self.obj = ctypes.c_void_p()
    CheckTag(tag)
    colDistVal = ctypes.c_uint(colDist)
    rowDistVal = ctypes.c_uint(rowDist)
    if   tag == iTag: 
      lib.ElDistMatrixCreateSpecific_i \
      (colDistVal,rowDistVal,grid,pointer(self.obj))
    elif tag == sTag: 
      lib.ElDistMatrixCreateSpecific_s \
      (colDistVal,rowDistVal,grid,pointer(self.obj))
    elif tag == dTag: 
      lib.ElDistMatrixCreateSpecific_d \
      (colDistVal,rowDistVal,grid,pointer(self.obj))
    elif tag == cTag: 
      lib.ElDistMatrixCreateSpecific_c \
      (colDistVal,rowDistVal,grid,pointer(self.obj))
    elif tag == zTag: 
      lib.ElDistMatrixCreateSpecific_z \
      (colDistVal,rowDistVal,grid,pointer(self.obj))
    self.tag = tag
