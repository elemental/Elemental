#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import ctypes, numpy, sys
from ctypes import byref, cdll, c_double, pointer, POINTER
lib = cdll.LoadLibrary('libEl.so')

def Initialize():
  argc = ctypes.c_int(len(sys.argv))
  _argv = ""
  for arg in sys.argv:
    _argv += arg + ' '
  argv = pointer(ctypes.c_char_p(_argv))
  lib.ElInitialize( pointer(argc), pointer(argv) )

def Finalize():
  lib.ElFinalize()

def Initialized():
  # NOTE: This is not expected to be portable and should be fixed
  active = ctypes.c_int()
  activeP = pointer(active) 
  lib.ElInitialized( activeP )
  return active

