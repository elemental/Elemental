#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import ctypes, numpy, sys
from ctypes import pointer
lib = ctypes.cdll.LoadLibrary('libEl.so')

buffer_from_memory = ctypes.pythonapi.PyBuffer_FromReadWriteMemory
buffer_from_memory.restype = ctypes.py_object

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
  lib.ElInitialize(pointer(argc),pointer(argv))

def Finalize():
  lib.ElFinalize()

def Initialized():
  # NOTE: This is not expected to be portable and should be fixed
  active = ctypes.c_int()
  activeP = pointer(active) 
  lib.ElInitialized( activeP )
  return active

# Create a simple enum for the supported datatypes
(iType,sType,dType,cType,zType)=(0,1,2,3,4)
def CheckType(datatype):
  if datatype != iType and \
     datatype != sType and datatype != dType and \
     datatype != cType and datatype != zType:
    print 'Unsupported datatype'

class ComplexFloat(ctypes.Structure):
  _fields_ = [("real",ctypes.c_float),("imag",ctypes.c_float)]
class ComplexDouble(ctypes.Structure):
  _fields_ = [("real",ctypes.c_double),("imag",ctypes.c_double)]

# Matrix
# ======

class Matrix(object):
  def __init__(self,initType):
    self.obj = ctypes.c_void_p()
    CheckType(initType)
    if   initType == iType: lib.ElMatrixCreate_i(pointer(self.obj))
    elif initType == sType: lib.ElMatrixCreate_s(pointer(self.obj))
    elif initType == dType: lib.ElMatrixCreate_d(pointer(self.obj))
    elif initType == cType: lib.ElMatrixCreate_c(pointer(self.obj))
    elif initType == zType: lib.ElMatrixCreate_z(pointer(self.obj))
    self.datatype = initType
  def Destroy(self):
    if   self.datatype == iType: lib.ElMatrixDestroy_i(self.obj)
    elif self.datatype == sType: lib.ElMatrixDestroy_s(self.obj)
    elif self.datatype == dType: lib.ElMatrixDestroy_d(self.obj)
    elif self.datatype == cType: lib.ElMatrixDestroy_c(self.obj)
    elif self.datatype == zType: lib.ElMatrixDestroy_z(self.obj)
  def SetType(self,newType):
    self.Destroy()
    CheckType(newType)
    if   newType == iType: lib.ElMatrixCreate_i(pointer(self.obj))
    elif newType == sType: lib.ElMatrixCreate_s(pointer(self.obj))
    elif newType == dType: lib.ElMatrixCreate_d(pointer(self.obj))
    elif newType == cType: lib.ElMatrixCreate_c(pointer(self.obj))
    elif newType == zType: lib.ElMatrixCreate_z(pointer(self.obj))
    self.datatype = newType
  def Resize(self,height,width):
    if   self.datatype == iType: lib.ElMatrixResize_i(self.obj,height,width)
    elif self.datatype == sType: lib.ElMatrixResize_s(self.obj,height,width)
    elif self.datatype == dType: lib.ElMatrixResize_d(self.obj,height,width)
    elif self.datatype == cType: lib.ElMatrixResize_c(self.obj,height,width)
    elif self.datatype == zType: lib.ElMatrixResize_z(self.obj,height,width)
  def ResizeWithLDim(self,height,width,ldim):
    if   self.datatype == iType: 
      lib.ElMatrixResizeWithLDim_i(self.obj,height,width,ldim)
    elif self.datatype == sType: 
      lib.ElMatrixResizeWithLDim_s(self.obj,height,width,ldim)
    elif self.datatype == dType: 
      lib.ElMatrixResizeWithLDim_d(self.obj,height,width,ldim)
    elif self.datatype == cType: 
      lib.ElMatrixResizeWithLDim_c(self.obj,height,width,ldim)
    elif self.datatype == zType: 
      lib.ElMatrixResizeWithLDim_z(self.obj,height,width,ldim)
  def Empty(self):
    if   self.datatype == iType: lib.ElMatrixEmpty_i(self.obj)
    elif self.datatype == sType: lib.ElMatrixEmpty_s(self.obj)
    elif self.datatype == dType: lib.ElMatrixEmpty_d(self.obj)
    elif self.datatype == cType: lib.ElMatrixEmpty_c(self.obj)
    elif self.datatype == zType: lib.ElMatrixEmpty_z(self.obj)
  def Attach(self,height,width,buffer,ldim):
    print 'Attaching buffers is not yet supported in Python'
  def LockedAttach(self,height,width,buffer,ldim):
    print 'Attaching buffers is not yet supported in Python'
  def Control(self,height,width,buffer,ldim):
    print 'Controlling buffers is not yet supported in Python'
  def Height(self):
    height = ctypes.c_int()
    if   self.datatype == iType: lib.ElMatrixHeight_i(self.obj,pointer(height))
    elif self.datatype == sType: lib.ElMatrixHeight_s(self.obj,pointer(height))
    elif self.datatype == dType: lib.ElMatrixHeight_d(self.obj,pointer(height))
    elif self.datatype == cType: lib.ElMatrixHeight_c(self.obj,pointer(height))
    elif self.datatype == zType: lib.ElMatrixHeight_z(self.obj,pointer(height))
    return height
  def Width(self):
    # TODO: Switch to 64-bit based upon Elemental's configuration
    width = ctypes.c_int()
    if   self.datatype == iType: lib.ElMatrixWidth_i(self.obj,pointer(width))
    elif self.datatype == sType: lib.ElMatrixWidth_s(self.obj,pointer(width))
    elif self.datatype == dType: lib.ElMatrixWidth_d(self.obj,pointer(width))
    elif self.datatype == cType: lib.ElMatrixWidth_c(self.obj,pointer(width))
    elif self.datatype == zType: lib.ElMatrixWidth_z(self.obj,pointer(width))
    return width
  def LDim(self):
    # TODO: Switch to 64-bit based upon Elemental's configuration
    ldim = ctypes.c_int()
    if   self.datatype == iType: lib.ElMatrixLDim_i(self.obj,pointer(ldim))
    elif self.datatype == sType: lib.ElMatrixLDim_s(self.obj,pointer(ldim))
    elif self.datatype == dType: lib.ElMatrixLDim_d(self.obj,pointer(ldim))
    elif self.datatype == cType: lib.ElMatrixLDim_c(self.obj,pointer(ldim))
    elif self.datatype == zType: lib.ElMatrixLDim_z(self.obj,pointer(ldim))
    return ldim 
  def MemorySize(self):
    # TODO: Switch to 64-bit based upon Elemental's configuration
    memSize = ctypes.c_int()
    if   self.datatype == iType: 
      lib.ElMatrixMemorySize_i(self.obj,pointer(memSize))
    elif self.datatype == sType:
      lib.ElMatrixMemorySize_s(self.obj,pointer(memSize))
    elif self.datatype == dType:
      lib.ElMatrixMemorySize_d(self.obj,pointer(memSize))
    elif self.datatype == cType:
      lib.ElMatrixMemorySize_c(self.obj,pointer(memSize))
    elif self.datatype == zType:
      lib.ElMatrixMemorySize_z(self.obj,pointer(memSize))
    return memSize
  def DiagonalLength(self):
    # TODO: Switch to 64-bit based upon Elemental's configuration
    length = ctypes.c_int()
    if   self.datatype == iType: 
      lib.ElMatrixDiagonalLength_i(self.obj,pointer(length))
    if   self.datatype == sType: 
      lib.ElMatrixDiagonalLength_s(self.obj,pointer(length))
    if   self.datatype == dType: 
      lib.ElMatrixDiagonalLength_d(self.obj,pointer(length))
    if   self.datatype == cType: 
      lib.ElMatrixDiagonalLength_c(self.obj,pointer(length))
    if   self.datatype == zType: 
      lib.ElMatrixDiagonalLength_z(self.obj,pointer(length))
    return length
  def Get(self,i,j):
    if   self.datatype == iType:
      # TODO: Switch to 64-bit based upon Elemental's configuration
      value = ctypes.c_int()
      lib.ElMatrixGet_i(self.obj,i,j,pointer(value))
      return value
    elif self.datatype == sType:
      value = ctypes.c_float()
      lib.ElMatrixGet_s(self.obj,i,j,pointer(value))
      return value
    elif self.datatype == dType:
      value = ctypes.c_double()
      lib.ElMatrixGet_d(self.obj,i,j,pointer(value))
    elif self.datatype == cType:
      value = ComplexFloat()
      lib.ElMatrixGet_c(self.obj,i,j,pointer(value))
      return value
    elif self.datatype == zType:
      value = ComplexDouble()
      lib.ElMatrixGet_z(self.obj,i,j,pointer(value))
      return value
  def Set(self,i,j,value):
    # TODO: Ensure that the input type is compatible; convert otherwise
    if   self.datatype == iType: lib.ElMatrixSet_i(self.obj,i,j,value)
    elif self.datatype == sType: lib.ElMatrixSet_s(self.obj,i,j,value)
    elif self.datatype == dType: lib.ElMatrixSet_d(self.obj,i,j,value)
    elif self.datatype == cType: lib.ElMatrixSet_c(self.obj,i,j,value)
    elif self.datatype == zType: lib.ElMatrixSet_z(self.obj,i,j,value)
  def ToNumPy(self):
    if   self.datatype == iType:
      # TODO: Switch to 64-bit based upon Elemental's configuration
      buf = buffer_from_memory(Buffer(),4*LDim()*Width())
      return numpy.frombuffer(buf,numpy.int32)
    elif self.datatype == sType:
      buf = buffer_from_memory(Buffer(),4*LDim()*Width())
      return buffer_from_memory(buf,numpy.float32)
    elif self.datatype == dType:
      buf = buffer_from_memory(Buffer(),8*LDim()*Width())
      return numpy.frombuffer(buf,numpy.float64)
    elif self.datatype == cType: 
      buf = buffer_from_memory(Buffer(),8*LDim()*Width())
      return numpy.frombuffer(buf,numpy.complex64)
    elif self.datatype == zType:
      buf = buffer_from_memory(Buffer(),16*LDim()*Width())
      return numpy.frombuffer(buf,numpy.complex128)

def Copy(A,B):
  CheckType(A.datatype)
  if A.datatype != B.datatype:
    print 'Copying between datatypes is not yet supported in Python'
    return
  if   B.datatype == iType: lib.ElMatrixCopy_i(A.obj,B.obj)
  elif B.datatype == sType: lib.ElMatrixCopy_s(A.obj,B.obj)
  elif B.datatype == dType: lib.ElMatrixCopy_d(A.obj,B.obj)
  elif B.datatype == cType: lib.ElMatrixCopy_c(A.obj,B.obj)
  elif B.datatype == zType: lib.ElMatrixCopy_z(A.obj,B.obj)

# Input/Output
# ************
def Print(A,s):
  if   A.datatype == iType: lib.ElPrint_i(A.obj,ctypes.c_char_p(s))
  elif A.datatype == sType: lib.ElPrint_s(A.obj,ctypes.c_char_p(s))
  elif A.datatype == dType: lib.ElPrint_d(A.obj,ctypes.c_char_p(s))
  elif A.datatype == cType: lib.ElPrint_c(A.obj,ctypes.c_char_p(s))
  elif A.datatype == zType: lib.ElPrint_z(A.obj,ctypes.c_char_p(s))
