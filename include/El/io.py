#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *

from ctypes import c_void_p, c_char_p, c_uint, c_int

# Input/Output
# ************

lib.ElPrint_i.argtypes = [c_void_p,c_char_p]
lib.ElPrint_i.restype = c_uint
lib.ElPrint_s.argtypes = [c_void_p,c_char_p]
lib.ElPrint_s.restype = c_uint
lib.ElPrint_d.argtypes = [c_void_p,c_char_p]
lib.ElPrint_d.restype = c_uint
lib.ElPrint_c.argtypes = [c_void_p,c_char_p]
lib.ElPrint_c.restype = c_uint
lib.ElPrint_z.argtypes = [c_void_p,c_char_p]
lib.ElPrint_z.restype = c_uint
lib.ElPrintDist_i.argtypes = [c_void_p,c_char_p]
lib.ElPrintDist_i.restype = c_uint
lib.ElPrintDist_s.argtypes = [c_void_p,c_char_p]
lib.ElPrintDist_s.restype = c_uint
lib.ElPrintDist_d.argtypes = [c_void_p,c_char_p]
lib.ElPrintDist_d.restype = c_uint
lib.ElPrintDist_c.argtypes = [c_void_p,c_char_p]
lib.ElPrintDist_c.restype = c_uint
lib.ElPrintDist_z.argtypes = [c_void_p,c_char_p]
lib.ElPrintDist_z.restype = c_uint
def Print(A,title):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElPrint_i(A.obj,title)
    elif A.tag == sTag: lib.ElPrint_s(A.obj,title)
    elif A.tag == dTag: lib.ElPrint_d(A.obj,title)
    elif A.tag == cTag: lib.ElPrint_c(A.obj,title)
    elif A.tag == zTag: lib.ElPrint_z(A.obj,title)
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElPrintDist_i(A.obj,title)
    elif A.tag == sTag: lib.ElPrintDist_s(A.obj,title)
    elif A.tag == dTag: lib.ElPrintDist_d(A.obj,title)
    elif A.tag == cTag: lib.ElPrintDist_c(A.obj,title)
    elif A.tag == zTag: lib.ElPrintDist_z(A.obj,title)
  else: raise Exception('Unsupported matrix type')

lib.ElSetColorMap.argtypes = [c_uint]
lib.ElSetColorMap.restype = c_uint
def SetColorMap(colorMap):
  lib.ElSetColorMap(colorMap)

lib.ElGetColorMap.argtypes = [POINTER(c_uint)]
lib.ElGetColorMap.restype = c_uint
def ColorMap():
  colorMap = c_uint()
  lib.ElGetColorMap(pointer(colorMap))
  return colorMap

lib.ElSetNumDiscreteColors.argtypes = [iType]
lib.ElSetNumDiscreteColors.restype = c_uint
def SetNumDiscreteColors(numColors):
  lib.ElSetNumDiscreteColors(numColors)

lib.ElNumDiscreteColors.argtypes = [POINTER(iType)]
lib.ElNumDiscreteColors.restype = c_uint
def NumDiscreteColors():
  numDiscrete = iType()
  lib.ElNumDiscreteColors(pointer(numDiscrete))
  return numDiscrete

lib.ElProcessEvents.argtypes = [c_int]
lib.ElProcessEvents.restype = c_uint
def ProcessEvents(numMsecs):
  lib.ElProcessEvents(numMsecs)

lib.ElDisplay_i.argtypes = [c_void_p,c_char_p]
lib.ElDisplay_i.restype = c_uint
lib.ElDisplay_s.argtypes = [c_void_p,c_char_p]
lib.ElDisplay_s.restype = c_uint
lib.ElDisplay_d.argtypes = [c_void_p,c_char_p]
lib.ElDisplay_d.restype = c_uint
lib.ElDisplay_c.argtypes = [c_void_p,c_char_p]
lib.ElDisplay_c.restype = c_uint
lib.ElDisplay_z.argtypes = [c_void_p,c_char_p]
lib.ElDisplay_z.restype = c_uint
lib.ElDisplayDist_i.argtypes = [c_void_p,c_char_p]
lib.ElDisplayDist_i.restype = c_uint
lib.ElDisplayDist_s.argtypes = [c_void_p,c_char_p]
lib.ElDisplayDist_s.restype = c_uint
lib.ElDisplayDist_d.argtypes = [c_void_p,c_char_p]
lib.ElDisplayDist_d.restype = c_uint
lib.ElDisplayDist_c.argtypes = [c_void_p,c_char_p]
lib.ElDisplayDist_c.restype = c_uint
lib.ElDisplayDist_z.argtypes = [c_void_p,c_char_p]
lib.ElDisplayDist_z.restype = c_uint
def Display(A,title):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElDisplay_i(A.obj,title)
    elif A.tag == sTag: lib.ElDisplay_s(A.obj,title)
    elif A.tag == dTag: lib.ElDisplay_d(A.obj,title)
    elif A.tag == cTag: lib.ElDisplay_c(A.obj,title)
    elif A.tag == zTag: lib.ElDisplay_z(A.obj,title)
    # Process an extra 200 milliseconds
    ProcessEvents(200)
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElDisplayDist_i(A.obj,title)
    elif A.tag == sTag: lib.ElDisplayDist_s(A.obj,title)
    elif A.tag == dTag: lib.ElDisplayDist_d(A.obj,title)
    elif A.tag == cTag: lib.ElDisplayDist_c(A.obj,title)
    elif A.tag == zTag: lib.ElDisplayDist_z(A.obj,title)
    # Process an extra 200 milliseconds
    ProcessEvents(200)
  else: raise Exception('Unsupported matrix type')

lib.ElSpy_i.argtypes = [c_void_p,c_char_p,iType]
lib.ElSpy_i.restype = c_uint
lib.ElSpy_s.argtypes = [c_void_p,c_char_p,sType]
lib.ElSpy_s.restype = c_uint
lib.ElSpy_d.argtypes = [c_void_p,c_char_p,dType]
lib.ElSpy_d.restype = c_uint
lib.ElSpy_c.argtypes = [c_void_p,c_char_p,sType]
lib.ElSpy_c.restype = c_uint
lib.ElSpy_z.argtypes = [c_void_p,c_char_p,dType]
lib.ElSpy_z.restype = c_uint
lib.ElSpyDist_i.argtypes = [c_void_p,c_char_p,iType]
lib.ElSpyDist_i.restype = c_uint
lib.ElSpyDist_s.argtypes = [c_void_p,c_char_p,sType]
lib.ElSpyDist_s.restype = c_uint
lib.ElSpyDist_d.argtypes = [c_void_p,c_char_p,dType]
lib.ElSpyDist_d.restype = c_uint
lib.ElSpyDist_c.argtypes = [c_void_p,c_char_p,sType]
lib.ElSpyDist_c.restype = c_uint
lib.ElSpyDist_z.argtypes = [c_void_p,c_char_p,dType]
lib.ElSpyDist_z.restype = c_uint
def Spy(A,title,tol):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElSpy_i(A.obj,title,tol)
    elif A.tag == sTag: lib.ElSpy_s(A.obj,title,tol)
    elif A.tag == dTag: lib.ElSpy_d(A.obj,title,tol)
    elif A.tag == cTag: lib.ElSpy_c(A.obj,title,tol)
    elif A.tag == zTag: lib.ElSpy_z(A.obj,title,tol)
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElSpyDist_i(A.obj,title,tol)
    elif A.tag == sTag: lib.ElSpyDist_s(A.obj,title,tol)
    elif A.tag == dTag: lib.ElSpyDist_d(A.obj,title,tol)
    elif A.tag == cTag: lib.ElSpyDist_c(A.obj,title,tol)
    elif A.tag == zTag: lib.ElSpyDist_z(A.obj,title,tol)
  else: raise Exception('Unsupported matrix type')

lib.ElRead_i.argtypes = [c_void_p,c_char_p,c_uint]
lib.ElRead_i.restype = c_uint
lib.ElRead_s.argtypes = [c_void_p,c_char_p,c_uint]
lib.ElRead_s.restype = c_uint
lib.ElRead_d.argtypes = [c_void_p,c_char_p,c_uint]
lib.ElRead_d.restype = c_uint
lib.ElRead_c.argtypes = [c_void_p,c_char_p,c_uint]
lib.ElRead_c.restype = c_uint
lib.ElRead_z.argtypes = [c_void_p,c_char_p,c_uint]
lib.ElRead_z.restype = c_uint
lib.ElReadDist_i.argtypes = [c_void_p,c_char_p,c_uint]
lib.ElReadDist_i.restype = c_uint
lib.ElReadDist_s.argtypes = [c_void_p,c_char_p,c_uint]
lib.ElReadDist_s.restype = c_uint
lib.ElReadDist_d.argtypes = [c_void_p,c_char_p,c_uint]
lib.ElReadDist_d.restype = c_uint
lib.ElReadDist_c.argtypes = [c_void_p,c_char_p,c_uint]
lib.ElReadDist_c.restype = c_uint
lib.ElReadDist_z.argtypes = [c_void_p,c_char_p,c_uint]
lib.ElReadDist_z.restype = c_uint
def Read(A,filename,fileFormat):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElRead_i(A.obj,filename,fileFormat)
    elif A.tag == sTag: lib.ElRead_s(A.obj,filename,fileFormat)
    elif A.tag == dTag: lib.ElRead_d(A.obj,filename,fileFormat)
    elif A.tag == cTag: lib.ElRead_c(A.obj,filename,fileFormat)
    elif A.tag == zTag: lib.ElRead_z(A.obj,filename,fileFormat)
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElReadDist_i(A.obj,filename,fileFormat)
    elif A.tag == sTag: lib.ElReadDist_s(A.obj,filename,fileFormat)
    elif A.tag == dTag: lib.ElReadDist_d(A.obj,filename,fileFormat)
    elif A.tag == cTag: lib.ElReadDist_c(A.obj,filename,fileFormat)
    elif A.tag == zTag: lib.ElReadDist_z(A.obj,filename,fileFormat)
  else: raise Exception('Unsupported matrix type')

lib.ElWrite_i.argtypes = [c_void_p,c_char_p,c_uint,c_char_p]
lib.ElWrite_i.restype = c_uint
lib.ElWrite_s.argtypes = [c_void_p,c_char_p,c_uint,c_char_p]
lib.ElWrite_s.restype = c_uint
lib.ElWrite_d.argtypes = [c_void_p,c_char_p,c_uint,c_char_p]
lib.ElWrite_d.restype = c_uint
lib.ElWrite_c.argtypes = [c_void_p,c_char_p,c_uint,c_char_p]
lib.ElWrite_c.restype = c_uint
lib.ElWrite_z.argtypes = [c_void_p,c_char_p,c_uint,c_char_p]
lib.ElWrite_z.restype = c_uint
lib.ElWriteDist_i.argtypes = [c_void_p,c_char_p,c_uint,c_char_p]
lib.ElWriteDist_i.restype = c_uint
lib.ElWriteDist_s.argtypes = [c_void_p,c_char_p,c_uint,c_char_p]
lib.ElWriteDist_s.restype = c_uint
lib.ElWriteDist_d.argtypes = [c_void_p,c_char_p,c_uint,c_char_p]
lib.ElWriteDist_d.restype = c_uint
lib.ElWriteDist_c.argtypes = [c_void_p,c_char_p,c_uint,c_char_p]
lib.ElWriteDist_c.restype = c_uint
lib.ElWriteDist_z.argtypes = [c_void_p,c_char_p,c_uint,c_char_p]
lib.ElWriteDist_z.restype = c_uint
def Write(A,basename,fileFormat,title):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElWrite_i(A.obj,basename,fileFormat,title)
    elif A.tag == sTag: lib.ElWrite_s(A.obj,basename,fileFormat,title)
    elif A.tag == dTag: lib.ElWrite_d(A.obj,basename,fileFormat,title)
    elif A.tag == cTag: lib.ElWrite_c(A.obj,basename,fileFormat,title)
    elif A.tag == zTag: lib.ElWrite_z(A.obj,basename,fileFormat,title)
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElWriteDist_i(A.obj,basename,fileFormat,title)
    elif A.tag == sTag: lib.ElWriteDist_s(A.obj,basename,fileFormat,title)
    elif A.tag == dTag: lib.ElWriteDist_d(A.obj,basename,fileFormat,title)
    elif A.tag == cTag: lib.ElWriteDist_c(A.obj,basename,fileFormat,title)
    elif A.tag == zTag: lib.ElWriteDist_z(A.obj,basename,fileFormat,title)
  else: raise Exception('Unsupported matrix type')
