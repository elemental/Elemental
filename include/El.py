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

# Core functionality
# ******************

# Environment
# ===========

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

# Matrix
# ======

class Matrix_i(object):
  def __init__(self):
    self.obj = ctypes.c_void_p()
    lib.ElMatrixCreate_i(pointer(self.obj))
  def Destroy(self):
    lib.ElMatrixDestroy_i(self.obj)
  def Resize(self,height,width):
    lib.ElMatrixResize_i(self.obj,height,width)
  def ResizeWithLDim(self,height,width,ldim):
    lib.ElMatrixResize_i(self.obj,height,width,ldim)
  def Empty(self):
    lib.ElMatrixEmpty_i(self.obj)

class Matrix_s(object):
  def __init__(self):
    self.obj = ctypes.c_void_p()
    lib.ElMatrixCreate_s(pointer(self.obj))
  def Destroy(self):
    lib.ElMatrixDestroy_s(self.obj)
  def Resize(self,height,width):
    lib.ElMatrixResize_s(self.obj,height,width)
  def ResizeWithLDim(self,height,width,ldim):
    lib.ElMatrixResize_s(self.obj,height,width,ldim)
  def Empty(self):
    lib.ElMatrixEmpty_s(self.obj)

class Matrix_d(object):
  def __init__(self):
    self.obj = ctypes.c_void_p()
    lib.ElMatrixCreate_d(pointer(self.obj))
  def Destroy(self):
    lib.ElMatrixDestroy_d(self.obj)
  def Resize(self,height,width):
    lib.ElMatrixResize_d(self.obj,height,width)
  def ResizeWithLDim(self,height,width,ldim):
    lib.ElMatrixResize_d(self.obj,height,width,ldim)
  def Empty(self):
    lib.ElMatrixEmpty_d(self.obj)

class Matrix_c(object):
  def __init__(self):
    self.obj = ctypes.c_void_p()
    lib.ElMatrixCreate_c(pointer(self.obj))
  def Destroy(self):
    lib.ElMatrixDestroy_c(self.obj)
  def Resize(self,height,width):
    lib.ElMatrixResize_c(self.obj,height,width)
  def ResizeWithLDim(self,height,width,ldim):
    lib.ElMatrixResize_c(self.obj,height,width,ldim)
  def Empty(self):
    lib.ElMatrixEmpty_c(self.obj)

class Matrix_z(object):
  def __init__(self):
    self.obj = ctypes.c_void_p()
    lib.ElMatrixCreate_z(pointer(self.obj))
  def Destroy(self):
    lib.ElMatrixDestroy_z(self.obj)
  def Resize(self,height,width):
    lib.ElMatrixResize_z(self.obj,height,width)
  def ResizeWithLDim(self,height,width,ldim):
    lib.ElMatrixResize_z(self.obj,height,width,ldim)
  def Empty(self):
    lib.ElMatrixEmpty_z(self.obj)

# Input/Output
# ************
def Print_i(A,s):
  lib.ElPrint_i(A.obj,ctypes.c_char_p(s))
def Print_s(matrix_s,s):
  lib.ElPrint_s(A.obj,ctypes.c_char_p(s))
def Print_d(A,s):
  lib.ElPrint_d(A.obj,ctypes.c_char_p(s))
def Print_c(A,s):
  lib.ElPrint_c(A.obj,ctypes.c_char_p(s))
def Print_z(A,s):
  lib.ElPrint_z(A.obj,ctypes.c_char_p(s))
