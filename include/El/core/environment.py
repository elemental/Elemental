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

# TODO: Switch to a different boolean type if appropriate
from ctypes import c_int    as bType
# TODO: Switch from c_int if Elemental was configured for 64-bit integers
from ctypes import c_int    as iType
from ctypes import c_float  as sType
from ctypes import c_double as dType
from ctypes import pointer
from ctypes import POINTER

# Query Elemental to determine whether MPI_Comm is an 'int' or a void pointer
commIsVoidP = bType()
lib.ElMPICommIsVoidPointer(pointer(commIsVoidP))
if commIsVoidP:
  MPI_Comm = ctypes.c_void_p
else:
  MPI_Comm = ctypes.c_int

# Query Elemental to determine whether MPI_Group is an 'int' or a void pointer
groupIsVoidP = bType()
lib.ElMPIGroupIsVoidPointer(pointer(groupIsVoidP))
if groupIsVoidP:
  MPI_Group = ctypes.c_void_p
else:
  MPI_Group = ctypes.c_int

# Create a simple enum for the supported datatypes
(iTag,sTag,dTag,cTag,zTag)=(0,1,2,3,4)
def CheckTag(tag):
  if tag != iTag and \
     tag != sTag and tag != dTag and \
     tag != cTag and tag != zTag:
    print 'Unsupported datatype'

# Emulate an enum for matrix distributions
(MC,MD,MR,VC,VR,STAR,CIRC)=(0,1,2,3,4,5,6)

# Emulate an enum for grid ordering
(ROW_MAJOR,COL_MAJOR)=(0,1)

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

import sys
def Initialize():
  argc = ctypes.c_int(len(sys.argv))
  _argv = ""
  for arg in sys.argv:
    _argv += arg + ' '
  argv = pointer(ctypes.c_char_p(_argv))
  lib.ElInitialize(pointer(argc),pointer(argv))

def Finalize():
  lib.ElFinalize()

def Initialized():
  # NOTE: This is not expected to be portable and should be fixed
  active = bType()
  activeP = pointer(active)
  lib.ElInitialized( activeP )
  return active

class ComplexFloat(ctypes.Structure):
  _fields_ = [("real",ctypes.c_float),("imag",ctypes.c_float)]
cType = ComplexFloat
class ComplexDouble(ctypes.Structure):
  _fields_ = [("real",ctypes.c_double),("imag",ctypes.c_double)]
zType = ComplexDouble

# Initialize MPI
Initialize()
