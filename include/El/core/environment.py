#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#

import ctypes

# TODO: Greatly improve this search functionality. At the moment, it is likely
#       necessary to manually specify the path
from ctypes.util import find_library
libPath = find_library('El')
if libPath == None:
  raise Exception("Could not find Elemental library")
lib = ctypes.cdll.LoadLibrary(libPath)
#lib = ctypes.cdll.LoadLibrary('/home/poulson/Source/Internal/Elemental/build/libEl.so')

# Environment
# ===========

# Basic types
# -----------

from ctypes import c_void_p, c_char_p, c_uint, c_int, c_float, c_double
from ctypes import pointer, POINTER

# TODO: Switch to a different boolean type if appropriate
bType = c_int
# TODO: Switch from c_int if Elemental was configured for 64-bit integers
iType = c_int
sType = c_float
dType = c_double

class ComplexFloat(ctypes.Structure):
  _fields_ = [("real",sType),("imag",sType)]
cType = ComplexFloat
class ComplexDouble(ctypes.Structure):
  _fields_ = [("real",dType),("imag",dType)]
zType = ComplexDouble

# Query Elemental to determine whether MPI_Comm is an 'int' or a void pointer
commIsVoidP = bType()
lib.ElMPICommIsVoidPointer(pointer(commIsVoidP))
if commIsVoidP:
  MPI_Comm = c_void_p
else:
  MPI_Comm = c_int

# Query Elemental to determine whether MPI_Group is an 'int' or a void pointer
groupIsVoidP = bType()
lib.ElMPIGroupIsVoidPointer(pointer(groupIsVoidP))
if groupIsVoidP:
  MPI_Group = c_void_p
else:
  MPI_Group = c_int

# Create a simple enum for the supported datatypes
(iTag,sTag,dTag,cTag,zTag)=(0,1,2,3,4)
def CheckTag(tag):
  if tag != iTag and \
     tag != sTag and tag != dTag and \
     tag != cTag and tag != zTag:
    print 'Unsupported datatype'

def Base(tag):
  if   tag == iTag: return iTag
  elif tag == sTag: return sTag
  elif tag == dTag: return dTag
  elif tag == cTag: return sTag
  elif tag == zTag: return dTag
  else: raise Exception('Invalid tag')

def TagToType(tag):
  if   tag == iTag: return iType
  elif tag == sTag: return sType
  elif tag == dTag: return dType
  elif tag == cTag: return cType
  elif tag == zTag: return zType
  else: raise Exception('Invalid tag')

# Emulate an enum for matrix distributions
(MC,MD,MR,VC,VR,STAR,CIRC)=(0,1,2,3,4,5,6)

# Emulate an enum for grid ordering
(ROW_MAJOR,COL_MAJOR)=(0,1)

# Emulate an enum for upper or lower
(LOWER,UPPER)=(0,1)

# Emulate an enum for left or right
(LEFT,RIGHT)=(0,1)

# Emulate the file format enum
(AUTO,ASCII,ASCII_MATLAB,BINARY,BINARY_FLAT,BMP,JPG,JPEG,MATRIX_MARKET,
 PNG,PPM,XBM,XPM)=(0,1,2,3,4,5,6,7,8,9,10,11,12)

# Emulate a colormap enum
(GRAYSCALE,GRAYSCALE_DISCRETE,RED_BLACK_GREEN,BLUE_RED)=(0,1,2,3)

# TODO: Many more enums

# Initialization
# --------------

lib.ElInitialize.argtypes = [POINTER(c_int),POINTER(POINTER(c_char_p))]
lib.ElInitialize.restype = c_uint
def Initialize():
  argc = c_int()
  argv = POINTER(c_char_p)()
  ctypes.pythonapi.Py_GetArgcArgv(ctypes.byref(argc),ctypes.byref(argv))
  lib.ElInitialize(pointer(argc),pointer(argv))

lib.ElFinalize.argtypes = []
lib.ElFinalize.restype = c_uint
def Finalize():
  lib.ElFinalize()

lib.ElInitialized.argtypes = [bType]
lib.ElInitialized.restype = c_uint
def Initialized():
  # NOTE: This is not expected to be portable and should be fixed
  active = bType()
  activeP = pointer(active)
  lib.ElInitialized( activeP )
  return active

# Initialize MPI
Initialize()
