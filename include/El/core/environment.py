#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#

import ctypes
#from ctypes.util import find_library
#libPath = find_library('El')
#if libPath == None:
#  raise Exception("Could not find Elemental library")
#lib = ctypes.cdll.LoadLibrary(libPath)
lib = ctypes.cdll.LoadLibrary('../build/libEl.so')

# TODO: Switch to a different boolean type if appropriate
from ctypes import c_int    as bType
# TODO: Switch from c_int if Elemental was configured for 64-bit integers
from ctypes import c_int    as iType
from ctypes import c_float  as sType
from ctypes import c_double as dType
from ctypes import pointer

# Environment
# ===========
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

# Create a simple enum for the supported datatypes
(iTag,sTag,dTag,cTag,zTag)=(0,1,2,3,4)
def CheckTag(tag):
  if tag != iTag and \
     tag != sTag and tag != dTag and \
     tag != cTag and tag != zTag:
    print 'Unsupported datatype'

class ComplexFloat(ctypes.Structure):
  _fields_ = [("real",ctypes.c_float),("imag",ctypes.c_float)]
cType = ComplexFloat
class ComplexDouble(ctypes.Structure):
  _fields_ = [("real",ctypes.c_double),("imag",ctypes.c_double)]
zType = ComplexDouble
