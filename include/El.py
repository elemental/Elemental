#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import ctypes, numpy, sys
lib = ctypes.cdll.LoadLibrary('libEl.so')

# TODO: Switch to a different boolean type if appropriate
from ctypes import c_int    as bType
# TODO: Switch from c_int if Elemental was configured for 64-bit integers
from ctypes import c_int    as iType
from ctypes import c_float  as sType
from ctypes import c_double as dType
from ctypes import pointer
from ctypes import POINTER

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

# Matrix
# ======

class Matrix(object):
  def __init__(self,initTag=dTag):
    self.obj = ctypes.c_void_p()
    CheckTag(initTag)
    if   initTag == iTag: lib.ElMatrixCreate_i(pointer(self.obj))
    elif initTag == sTag: lib.ElMatrixCreate_s(pointer(self.obj))
    elif initTag == dTag: lib.ElMatrixCreate_d(pointer(self.obj))
    elif initTag == cTag: lib.ElMatrixCreate_c(pointer(self.obj))
    elif initTag == zTag: lib.ElMatrixCreate_z(pointer(self.obj))
    self.tag = initTag
  def Destroy(self):
    if   self.tag == iTag: lib.ElMatrixDestroy_i(self.obj)
    elif self.tag == sTag: lib.ElMatrixDestroy_s(self.obj)
    elif self.tag == dTag: lib.ElMatrixDestroy_d(self.obj)
    elif self.tag == cTag: lib.ElMatrixDestroy_c(self.obj)
    elif self.tag == zTag: lib.ElMatrixDestroy_z(self.obj)
  def SetType(self,newTag=dTag):
    self.Destroy()
    CheckTag(newTag)
    if   newTag == iTag: lib.ElMatrixCreate_i(pointer(self.obj))
    elif newTag == sTag: lib.ElMatrixCreate_s(pointer(self.obj))
    elif newTag == dTag: lib.ElMatrixCreate_d(pointer(self.obj))
    elif newTag == cTag: lib.ElMatrixCreate_c(pointer(self.obj))
    elif newTag == zTag: lib.ElMatrixCreate_z(pointer(self.obj))
    self.tag = newTag
  def Resize(self,height,width):
    if   self.tag == iTag: lib.ElMatrixResize_i(self.obj,height,width)
    elif self.tag == sTag: lib.ElMatrixResize_s(self.obj,height,width)
    elif self.tag == dTag: lib.ElMatrixResize_d(self.obj,height,width)
    elif self.tag == cTag: lib.ElMatrixResize_c(self.obj,height,width)
    elif self.tag == zTag: lib.ElMatrixResize_z(self.obj,height,width)
  def ResizeWithLDim(self,height,width,ldim):
    if   self.tag == iTag: 
      lib.ElMatrixResizeWithLDim_i(self.obj,height,width,ldim)
    elif self.tag == sTag: 
      lib.ElMatrixResizeWithLDim_s(self.obj,height,width,ldim)
    elif self.tag == dTag: 
      lib.ElMatrixResizeWithLDim_d(self.obj,height,width,ldim)
    elif self.tag == cTag: 
      lib.ElMatrixResizeWithLDim_c(self.obj,height,width,ldim)
    elif self.tag == zTag: 
      lib.ElMatrixResizeWithLDim_z(self.obj,height,width,ldim)
  def Empty(self):
    if   self.tag == iTag: lib.ElMatrixEmpty_i(self.obj)
    elif self.tag == sTag: lib.ElMatrixEmpty_s(self.obj)
    elif self.tag == dTag: lib.ElMatrixEmpty_d(self.obj)
    elif self.tag == cTag: lib.ElMatrixEmpty_c(self.obj)
    elif self.tag == zTag: lib.ElMatrixEmpty_z(self.obj)
  def Attach(self,height,width,buf,ldim):
    if   self.tag == iTag:
      if type(buf) is not POINTER(iType):
        print 'Invalid buffer type'
        return
      lib.ElMatrixAttach_i(self.obj,height,width,buf,ldim)
    elif self.tag == sTag:
      if type(buf) is not POINTER(sType):
        print 'Invalid buffer type'
        return
      lib.ElMatrixAttach_s(self.obj,height,width,buf,ldim)
    elif self.tag == dTag:
      if type(buf) is not POINTER(dType):
        print 'Invalid buffer type'
        return
      lib.ElMatrixAttach_d(self.obj,height,width,buf,ldim)
    elif self.tag == cTag:
      if type(buf) is not POINTER(cType):
        print 'Invalid buffer type'
        return
      lib.ElMatrixAttach_c(self.obj,height,width,buf,ldim)
    elif self.tag == zTag:
      if type(buf) is not POINTER(zType):
        print 'Invalid buffer type'
        return
      lib.ElMatrixAttach_z(self.obj,height,width,buf,ldim)
  def LockedAttach(self,height,width,buf,ldim):
    if   self.tag == iTag:
      if type(buf) is not POINTER(iType):
        print 'Invalid buffer type'
        return
      lib.ElMatrixLockedAttach_i(self.obj,height,width,buf,ldim)
    elif self.tag == sTag:
      if type(buf) is not POINTER(sType):
        print 'Invalid buffer type'
        return
      lib.ElMatrixLockedAttach_s(self.obj,height,width,buf,ldim)
    elif self.tag == dTag:
      if type(buf) is not POINTER(dType):
        print 'Invalid buffer type'
        return
      lib.ElMatrixLockedAttach_d(self.obj,height,width,buf,ldim)
    elif self.tag == cTag:
      if type(buf) is not POINTER(cType):
        print 'Invalid buffer type'
        return
      lib.ElMatrixLockedAttach_c(self.obj,height,width,buf,ldim)
    elif self.tag == zTag:
      if type(buf) is not POINTER(zType):
        print 'Invalid buffer type'
        return
      lib.ElMatrixLockedAttach_z(self.obj,height,width,buf,ldim)
  def Control(self,height,width,buf,ldim):
    if   self.tag == iTag:
      if type(buf) is not POINTER(iType):
        print 'Invalid buffer type'
        return
      lib.ElMatrixControl_i(self.obj,height,width,buf,ldim)
    elif self.tag == sTag:
      if type(buf) is not POINTER(sType):
        print 'Invalid buffer type'
        return
      lib.ElMatrixControl_s(self.obj,height,width,buf,ldim)
    elif self.tag == dTag:
      if type(buf) is not POINTER(dType):
        print 'Invalid buffer type'
        return
      lib.ElMatrixControl_d(self.obj,height,width,buf,ldim)
    elif self.tag == cTag:
      if type(buf) is not POINTER(cType):
        print 'Invalid buffer type'
        return
      lib.ElMatrixControl_c(self.obj,height,width,buf,ldim)
    elif self.tag == zTag:
      if type(buf) is not POINTER(zType):
        print 'Invalid buffer type'
        return
      lib.ElMatrixControl_z(self.obj,height,width,buf,ldim)
  def Height(self):
    height = iType()
    if   self.tag == iTag: lib.ElMatrixHeight_i(self.obj,pointer(height))
    elif self.tag == sTag: lib.ElMatrixHeight_s(self.obj,pointer(height))
    elif self.tag == dTag: lib.ElMatrixHeight_d(self.obj,pointer(height))
    elif self.tag == cTag: lib.ElMatrixHeight_c(self.obj,pointer(height))
    elif self.tag == zTag: lib.ElMatrixHeight_z(self.obj,pointer(height))
    return height
  def Width(self):
    width = iType()
    if   self.tag == iTag: lib.ElMatrixWidth_i(self.obj,pointer(width))
    elif self.tag == sTag: lib.ElMatrixWidth_s(self.obj,pointer(width))
    elif self.tag == dTag: lib.ElMatrixWidth_d(self.obj,pointer(width))
    elif self.tag == cTag: lib.ElMatrixWidth_c(self.obj,pointer(width))
    elif self.tag == zTag: lib.ElMatrixWidth_z(self.obj,pointer(width))
    return width
  def LDim(self):
    ldim = iType()
    if   self.tag == iTag: lib.ElMatrixLDim_i(self.obj,pointer(ldim))
    elif self.tag == sTag: lib.ElMatrixLDim_s(self.obj,pointer(ldim))
    elif self.tag == dTag: lib.ElMatrixLDim_d(self.obj,pointer(ldim))
    elif self.tag == cTag: lib.ElMatrixLDim_c(self.obj,pointer(ldim))
    elif self.tag == zTag: lib.ElMatrixLDim_z(self.obj,pointer(ldim))
    return ldim 
  def MemorySize(self):
    memSize = iType()
    if   self.tag == iTag: 
      lib.ElMatrixMemorySize_i(self.obj,pointer(memSize))
    elif self.tag == sTag:
      lib.ElMatrixMemorySize_s(self.obj,pointer(memSize))
    elif self.tag == dTag:
      lib.ElMatrixMemorySize_d(self.obj,pointer(memSize))
    elif self.tag == cTag:
      lib.ElMatrixMemorySize_c(self.obj,pointer(memSize))
    elif self.tag == zTag:
      lib.ElMatrixMemorySize_z(self.obj,pointer(memSize))
    return memSize
  def DiagonalLength(self,offset=iType(0)):
    length = iType()
    if   self.tag == iTag: 
      lib.ElMatrixDiagonalLength_i(self.obj,offset,pointer(length))
    if   self.tag == sTag: 
      lib.ElMatrixDiagonalLength_s(self.obj,offset,pointer(length))
    if   self.tag == dTag: 
      lib.ElMatrixDiagonalLength_d(self.obj,offset,pointer(length))
    if   self.tag == cTag: 
      lib.ElMatrixDiagonalLength_c(self.obj,offset,pointer(length))
    if   self.tag == zTag: 
      lib.ElMatrixDiagonalLength_z(self.obj,offset,pointer(length))
    return length
  def Viewing(self):
    viewing = bType()
    if   self.tag == iTag: lib.ElMatrixViewing_i(self.obj,pointer(viewing))
    elif self.tag == sTag: lib.ElMatrixViewing_s(self.obj,pointer(viewing))
    elif self.tag == dTag: lib.ElMatrixViewing_d(self.obj,pointer(viewing))
    elif self.tag == cTag: lib.ElMatrixViewing_c(self.obj,pointer(viewing))
    elif self.tag == zTag: lib.ElMatrixViewing_z(self.obj,pointer(viewing))
    return viewing 
  def FixedSize(self):
    fixed = bType()
    if   self.tag == iTag: lib.ElMatrixFixedSize_i(self.obj,pointer(fixed))
    elif self.tag == sTag: lib.ElMatrixFixedSize_s(self.obj,pointer(fixed))
    elif self.tag == dTag: lib.ElMatrixFixedSize_d(self.obj,pointer(fixed))
    elif self.tag == cTag: lib.ElMatrixFixedSize_c(self.obj,pointer(fixed))
    elif self.tag == zTag: lib.ElMatrixFixedSize_z(self.obj,pointer(fixed))
    return fixed
  def Locked(self):
    locked = bType()
    if   self.tag == iTag: lib.ElMatrixLocked_i(self.obj,pointer(locked))
    elif self.tag == sTag: lib.ElMatrixLocked_s(self.obj,pointer(locked))
    elif self.tag == dTag: lib.ElMatrixLocked_d(self.obj,pointer(locked))
    elif self.tag == cTag: lib.ElMatrixLocked_c(self.obj,pointer(locked))
    elif self.tag == zTag: lib.ElMatrixLocked_z(self.obj,pointer(locked))
    return locked
  def Buffer(self):
    if   self.tag == iTag:
      buf = ctypes.POINTER(iType)()
      lib.ElMatrixBuffer_i(self.obj,pointer(buf)) 
      return buf
    elif self.tag == sTag:
      buf = ctypes.POINTER(sType)()
      lib.ElMatrixBuffer_s(self.obj,pointer(buf))
      return buf
    elif self.tag == dTag: 
      buf = ctypes.POINTER(dType)()
      lib.ElMatrixBuffer_d(self.obj,pointer(buf))
      return buf
    elif self.tag == cTag:
      buf = ctypes.POINTER(cType)()
      lib.ElMatrixBuffer_c(self.obj,pointer(buf))
      return buf
    elif self.tag == zTag:
      buf = ctypes.POINTER(zType)()
      lib.ElMatrixBuffer_z(self.obj,pointer(buf))
      return buf
  def LockedBuffer(self):
    if   self.tag == iTag:
      buf = ctypes.POINTER(iType)()
      lib.ElMatrixLockedBuffer_i(self.obj,pointer(buf)) 
      return buf
    elif self.tag == sTag:
      buf = ctypes.POINTER(sType)()
      lib.ElMatrixLockedBuffer_s(self.obj,pointer(buf))
      return buf
    elif self.tag == dTag: 
      buf = ctypes.POINTER(dType)()
      lib.ElMatrixLockedBuffer_d(self.obj,pointer(buf))
      return buf
    elif self.tag == cTag:
      buf = ctypes.POINTER(cType)()
      lib.ElMatrixLockedBuffer_c(self.obj,pointer(buf))
      return buf
    elif self.tag == zTag:
      buf = ctypes.POINTER(zType)()
      lib.ElMatrixLockedBuffer_z(self.obj,pointer(buf))
      return buf
  def Get(self,i,j):
    if   self.tag == iTag:
      value = iType()
      lib.ElMatrixGet_i(self.obj,i,j,pointer(value))
      return value
    elif self.tag == sTag:
      value = sType()
      lib.ElMatrixGet_s(self.obj,i,j,pointer(value))
      return value
    elif self.tag == dTag:
      value = dType()
      lib.ElMatrixGet_d(self.obj,i,j,pointer(value))
    elif self.tag == cTag:
      value = cType()
      lib.ElMatrixGet_c(self.obj,i,j,pointer(value))
      return value
    elif self.tag == zTag:
      value = zType()
      lib.ElMatrixGet_z(self.obj,i,j,pointer(value))
      return value
  def GetRealPart(self,i,j):
    if self.tag == cTag:
      value = sType()
      lib.ElMatrixGetRealPart_c(self.obj,i,j,pointer(value))
      return value
    elif self.tag == zTag:
      value = dType()
      lib.ElMatrixGetRealPart_z(self.obj,i,j,pointer(value))
      return value
    else: return Get(i,j)
  def GetImagPart(self,i,j):
    if   self.tag == iTag: return iType(0)
    elif self.tag == sTag: return sType(0)
    elif self.tag == dTag: return dType(0)
    elif self.tag == cTag:
      value = ctypes.c_float()
      lib.ElMatrixGetImagPart_c(self.obj,i,j,pointer(value))
      return value
    elif self.tag == zTag:
      value = ctypes.c_double()
      lib.ElMatrixGetImagPart_z(self.obj,i,j,pointer(value))
      return value
  def Set(self,i,j,value):
    if   self.tag == iTag: 
      if type(value) is not iType:
        print 'Invalid scalar type'
        return
      lib.ElMatrixSet_i(self.obj,i,j,iType(value))
    elif self.tag == sTag: 
      if type(value) is not sType:
        print 'Invalid scalar type'
        return
      lib.ElMatrixSet_s(self.obj,i,j,value)
    elif self.tag == dTag: 
      if type(value) is not dType: 
        print 'Invalid scalar type'
        return
      lib.ElMatrixSet_d(self.obj,i,j,value)
    elif self.tag == cTag: 
      if type(value) is not cType: 
        print 'Invalid scalar type'
        return
      lib.ElMatrixSet_c(self.obj,i,j,value)
    elif self.tag == zTag: 
      if type(value) is not zType: 
        print 'Invalid scalar type'
        return
      lib.ElMatrixSet_z(self.obj,i,j,value)
  def SetRealPart(self,i,j,value):
    if self.tag == cTag: 
      lib.ElMatrixSetRealPart_c(self.obj,i,j,sType(value))
    elif self.tag == zTag: 
      lib.ElMatrixSetRealPart_z(self.obj,i,j,dType(value))
    else: Set(i,j,value)
  def SetImagPart(self,i,j,value):
    if self.tag == cTag: 
      lib.ElMatrixSetImagPart_c(self.obj,i,j,sType(value))
    elif self.tag == zTag: 
      lib.ElMatrixSetImagPart_z(self.obj,i,j,dType(value))
    else: print 'Datatype does not have an imaginary component'
  def Update(self,i,j,value):
    if   self.tag == iTag: 
      if type(value) is not iType:
        print 'Invalid scalar type'
        return
      lib.ElMatrixUpdate_i(self.obj,i,j,iType(value))
    elif self.tag == sTag: 
      if type(value) is not sType:
        print 'Invalid scalar type'
        return
      lib.ElMatrixUpdate_s(self.obj,i,j,value)
    elif self.tag == dTag: 
      if type(value) is not dType: 
        print 'Invalid scalar type'
        return
      lib.ElMatrixUpdate_d(self.obj,i,j,value)
    elif self.tag == cTag: 
      if type(value) is not cType: 
        print 'Invalid scalar type'
        return
      lib.ElMatrixUpdate_c(self.obj,i,j,value)
    elif self.tag == zTag: 
      if type(value) is not zType: 
        print 'Invalid scalar type'
        return
      lib.ElMatrixUpdate_z(self.obj,i,j,value)
  def UpdateRealPart(self,i,j,value):
    if self.tag == cTag: 
      lib.ElMatrixUpdateRealPart_c(self.obj,i,j,sType(value))
    elif self.tag == zTag: 
      lib.ElMatrixUpdateRealPart_z(self.obj,i,j,dType(value))
    else: Update(i,j,value)
  def UpdateImagPart(self,i,j,value):
    if self.tag == cTag: 
      lib.ElMatrixUpdateImagPart_c(self.obj,i,j,sType(value))
    elif self.tag == zTag: 
      lib.ElMatrixUpdateImagPart_z(self.obj,i,j,dType(value))
    else: print 'Datatype does not have an imaginary component'
  def MakeReal(self,i,j):
    if   self.tag == cTag: lib.ElMatrixMakeReal_c(self.obj,i,j)
    elif self.tag == zTag: lib.ElMatrixMakeReal_z(self.obj,i,j) 
  def Conjugate(self,i,j):
    if   self.tag == cTag: lib.ElMatrixConjugate_c(self.obj,i,j)
    elif self.tag == zTag: lib.ElMatrixConjugate_z(self.obj,i,j)
  def GetDiagonal(self,offset=iType(0)):
    print 'GetDiagonal not yet supported by Python interface'
  def GetRealPartOfDiagonal(self,offset=iType(0)):
    print 'GetRealPartOfDiagonal not yet supported by Python interface'
  def GetImagPartOfDiagonal(self,offset=iType(0)):
    print 'GetImagPartOfDiagonal not yet supported by Python interface'
  def SetDiagonal(self,diag,offset=iType(0)):
    print 'SetDiagonal not yet supported by Python interface'
  def SetRealPartOfDiagonal(self,diag,offset=iType(0)):
    print 'SetRealPartOfDiagonal not yet supported by Python interface'
  def SetImagPartOfDiagonal(self,diag,offset=iType(0)):
    print 'SetImagPartOfDiagonal not yet supported by Python interface'
  def UpdateDiagonal(self,diag,offset=iType(0)):
    print 'UpdateDiagonal not yet supported by Python interface'
  def UpdateRealPartOfDiagonal(self,diag,offset=iType(0)):
    print 'UpdateRealPartOfDiagonal not yet supported by Python interface'
  def UpdateImagPartOfDiagonal(self,diag,offset=iType(0)):
    print 'UpdateImagPartOfDiagonal not yet supported by Python interface'
  def MakeDiagonalReal(self,offset=iType(0)):
    if   self.tag == iTag: lib.ElMatrixMakeDiagonalReal_i(self.obj,offset)
    elif self.tag == sTag: lib.ElMatrixMakeDiagonalReal_s(self.obj,offset)
    elif self.tag == dTag: lib.ElMatrixMakeDiagonalReal_d(self.obj,offset)
    elif self.tag == cTag: lib.ElMatrixMakeDiagonalReal_c(self.obj,offset)
    elif self.tag == zTag: lib.ElMatrixMakeDiagonalReal_z(self.obj,offset)
  def ConjugateDiagonal(self,offset=iType(0)):
    if   self.tag == iTag: lib.ElMatrixConjugateDiagonal_i(self.obj,offset)
    elif self.tag == sTag: lib.ElMatrixConjugateDiagonal_s(self.obj,offset)
    elif self.tag == dTag: lib.ElMatrixConjugateDiagonal_d(self.obj,offset)
    elif self.tag == cTag: lib.ElMatrixConjugateDiagonal_c(self.obj,offset)
    elif self.tag == zTag: lib.ElMatrixConjugateDiagonal_z(self.obj,offset)
  def GetSubmatrix(self,I,J):
    print 'GetSubmatrix not yet supported by Python interface'
  def GetRealPartOfSubmatrix(self,I,J):
    print 'GetRealPartOfSubmatrix not yet supported by Python interface'
  def GetImagPartOfSubmatrix(self,I,J):
    print 'GetImagPartOfSubmatrix not yet supported by Python interface'
  def SetSubmatrix(self,I,J,ASub):
    print 'SetSubmatrix not yet supported by Python interface'
  def SetRealPartOfSubmatrix(self,I,J,ASub):
    print 'SetRealPartOfSubmatrix not yet supported by Python interface'
  def SetImagPartOfSubmatrix(self,I,J,ASub):
    print 'SetImagPartOfSubmatrix not yet supported by Python interface'
  def UpdateSubmatrix(self,I,J,ASub):
    print 'UpdateSubmatrix not yet supported by Python interface'
  def UpdateRealPartOfSubmatrix(self,I,J,ASub):
    print 'UpdateRealPartOfSubmatrix not yet supported by Python interface'
  def UpdateImagPartOfSubmatrix(self,I,J,ASub):
    print 'UpdateImagPartOfSubmatrix not yet supported by Python interface'
  def MakeSubmatrixReal(self,I,J):
    print 'MakeSubmatrixReal not yet supported by Python interface'
  def ConjugateSubmatrix(self,I,J):
    print 'ConjugateSubmatrix not yet supported by Python interface'
  def ToNumPy(self):
    if   self.tag == iTag:
      # TODO: Switch to 64-bit based upon Elemental's configuration
      buf = buffer_from_memory(Buffer(),4*LDim()*Width())
      return numpy.frombuffer(buf,numpy.int32)
    elif self.tag == sTag:
      buf = buffer_from_memory(Buffer(),4*LDim()*Width())
      return buffer_from_memory(buf,numpy.float32)
    elif self.tag == dTag:
      buf = buffer_from_memory(Buffer(),8*LDim()*Width())
      return numpy.frombuffer(buf,numpy.float64)
    elif self.tag == cTag: 
      buf = buffer_from_memory(Buffer(),8*LDim()*Width())
      return numpy.frombuffer(buf,numpy.complex64)
    elif self.tag == zTag:
      buf = buffer_from_memory(Buffer(),16*LDim()*Width())
      return numpy.frombuffer(buf,numpy.complex128)

def Copy(A,B):
  CheckTag(A.tag)
  if A.tag != B.tag:
    print 'Copying between datatypes is not yet supported in Python'
    return
  if type(A) is Matrix and type(B) is Matrix:
    if   B.tag == iTag: lib.ElMatrixCopy_i(A.obj,B.obj)
    elif B.tag == sTag: lib.ElMatrixCopy_s(A.obj,B.obj)
    elif B.tag == dTag: lib.ElMatrixCopy_d(A.obj,B.obj)
    elif B.tag == cTag: lib.ElMatrixCopy_c(A.obj,B.obj)
    elif B.tag == zTag: lib.ElMatrixCopy_z(A.obj,B.obj)
  else: print 'Unsupported matrix types'

# Input/Output
# ************
def Print(A,s):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElPrint_i(A.obj,ctypes.c_char_p(s))
    elif A.tag == sTag: lib.ElPrint_s(A.obj,ctypes.c_char_p(s))
    elif A.tag == dTag: lib.ElPrint_d(A.obj,ctypes.c_char_p(s))
    elif A.tag == cTag: lib.ElPrint_c(A.obj,ctypes.c_char_p(s))
    elif A.tag == zTag: lib.ElPrint_z(A.obj,ctypes.c_char_p(s))
  else: print 'Unsupported matrix type'
