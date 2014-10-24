#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *
from El.blas_like import Copy, RealPart, ImagPart

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
lib.ElPrintGraph.argtypes = [c_void_p,c_char_p]
lib.ElPrintGraph.restype = c_uint
lib.ElPrintDistGraph.argtypes = [c_void_p,c_char_p]
lib.ElPrintDistGraph.restype = c_uint
lib.ElPrintSparse_i.argtypes = [c_void_p,c_char_p]
lib.ElPrintSparse_i.restype = c_uint
lib.ElPrintSparse_s.argtypes = [c_void_p,c_char_p]
lib.ElPrintSparse_s.restype = c_uint
lib.ElPrintSparse_d.argtypes = [c_void_p,c_char_p]
lib.ElPrintSparse_d.restype = c_uint
lib.ElPrintSparse_c.argtypes = [c_void_p,c_char_p]
lib.ElPrintSparse_c.restype = c_uint
lib.ElPrintSparse_z.argtypes = [c_void_p,c_char_p]
lib.ElPrintSparse_z.restype = c_uint
lib.ElPrintDistSparse_i.argtypes = [c_void_p,c_char_p]
lib.ElPrintDistSparse_i.restype = c_uint
lib.ElPrintDistSparse_s.argtypes = [c_void_p,c_char_p]
lib.ElPrintDistSparse_s.restype = c_uint
lib.ElPrintDistSparse_d.argtypes = [c_void_p,c_char_p]
lib.ElPrintDistSparse_d.restype = c_uint
lib.ElPrintDistSparse_c.argtypes = [c_void_p,c_char_p]
lib.ElPrintDistSparse_c.restype = c_uint
lib.ElPrintDistSparse_z.argtypes = [c_void_p,c_char_p]
lib.ElPrintDistSparse_z.restype = c_uint
def Print(A,title=''):
  args = [A.obj,title]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElPrint_i(*args)
    elif A.tag == sTag: lib.ElPrint_s(*args)
    elif A.tag == dTag: lib.ElPrint_d(*args)
    elif A.tag == cTag: lib.ElPrint_c(*args)
    elif A.tag == zTag: lib.ElPrint_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElPrintDist_i(*args)
    elif A.tag == sTag: lib.ElPrintDist_s(*args)
    elif A.tag == dTag: lib.ElPrintDist_d(*args)
    elif A.tag == cTag: lib.ElPrintDist_c(*args)
    elif A.tag == zTag: lib.ElPrintDist_z(*args)
    else: DataExcept()
  elif type(A) is Graph:
    lib.ElPrintGraph(*args)
  elif type(A) is DistGraph:
    lib.ElPrintDistGraph(*args)
  elif type(A) is SparseMatrix:
    if   A.tag == iTag: lib.ElPrintSparse_i(*args)
    elif A.tag == sTag: lib.ElPrintSparse_s(*args)
    elif A.tag == dTag: lib.ElPrintSparse_d(*args)
    elif A.tag == cTag: lib.ElPrintSparse_c(*args)
    elif A.tag == zTag: lib.ElPrintSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == iTag: lib.ElPrintDistSparse_i(*args)
    elif A.tag == sTag: lib.ElPrintDistSparse_s(*args)
    elif A.tag == dTag: lib.ElPrintDistSparse_d(*args)
    elif A.tag == cTag: lib.ElPrintDistSparse_c(*args)
    elif A.tag == zTag: lib.ElPrintDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

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
lib.ElDisplayGraph.argtypes = [c_void_p,c_char_p]
lib.ElDisplayGraph.restype = c_uint
lib.ElDisplayDistGraph.argtypes = [c_void_p,c_char_p]
lib.ElDisplayDistGraph.restype = c_uint
lib.ElDisplaySparse_i.argtypes = [c_void_p,c_char_p]
lib.ElDisplaySparse_i.restype = c_uint
lib.ElDisplaySparse_s.argtypes = [c_void_p,c_char_p]
lib.ElDisplaySparse_s.restype = c_uint
lib.ElDisplaySparse_d.argtypes = [c_void_p,c_char_p]
lib.ElDisplaySparse_d.restype = c_uint
lib.ElDisplaySparse_c.argtypes = [c_void_p,c_char_p]
lib.ElDisplaySparse_c.restype = c_uint
lib.ElDisplaySparse_z.argtypes = [c_void_p,c_char_p]
lib.ElDisplaySparse_z.restype = c_uint
lib.ElDisplayDistSparse_i.argtypes = [c_void_p,c_char_p]
lib.ElDisplayDistSparse_i.restype = c_uint
lib.ElDisplayDistSparse_s.argtypes = [c_void_p,c_char_p]
lib.ElDisplayDistSparse_s.restype = c_uint
lib.ElDisplayDistSparse_d.argtypes = [c_void_p,c_char_p]
lib.ElDisplayDistSparse_d.restype = c_uint
lib.ElDisplayDistSparse_c.argtypes = [c_void_p,c_char_p]
lib.ElDisplayDistSparse_c.restype = c_uint
lib.ElDisplayDistSparse_z.argtypes = [c_void_p,c_char_p]
lib.ElDisplayDistSparse_z.restype = c_uint
def Display(A,title='',tryPython=True):
  if tryPython: 
    if type(A) is Matrix or type(A) is DistMatrix:
      try:  
        import numpy as np
        import matplotlib.pyplot as plt
        isVec = min(A.Height(),A.Width()) == 1
        if type(A) is Matrix:
          if A.tag == cTag or A.tag == zTag:
            AReal = Matrix(Base(A.tag))
            AImag = Matrix(Base(A.tag))
            RealPart(A,AReal)
            ImagPart(A,AImag)
            fig, (ax1,ax2) = plt.subplots(1,2)
            ax1.set_title('Real part')
            ax2.set_title('Imag part')
            if isVec:
              ax1.plot(np.squeeze(AReal.ToNumPy()),'bo-')
              ax2.plot(np.squeeze(AImag.ToNumPy()),'bo-')
            else:
              imReal = ax1.imshow(AReal.ToNumPy())
              cBarReal = fig.colorbar(imReal,ax=ax1)
              imImag = ax2.imshow(AImag.ToNumPy())
              cBarImag = fig.colorbar(imImag,ax=ax2)
            plt.suptitle(title)
            plt.tight_layout()
          else:
            fig = plt.figure()
            axis = fig.add_axes([0.1,0.1,0.8,0.8])
            if isVec:
              axis.plot(np.squeeze(A.ToNumPy()),'bo-')
            else:
              im = axis.imshow(A.ToNumPy())
              fig.colorbar(im,ax=axis)
            plt.title(title)
          plt.draw()
          plt.show(block=False)
        elif type(A) is DistMatrix:
          A_CIRC_CIRC = DistMatrix(A.tag,CIRC,CIRC,A.Grid())
          Copy(A,A_CIRC_CIRC)
          if A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root():
            Display(A_CIRC_CIRC.Matrix(),title,tryPython)
            return
        else: raise Exception('Unsupported matrix type')
        return
      except: 
        print 'Could not import matplotlib.pyplot'
    # TODO: Gather distributed graph on a single process
    elif type(A) is Graph:
      try:  
        import matplotlib.pyplot as plt
        import networkx as nx
        numEdges = A.NumEdges() 
        G = nx.DiGraph()
        for edge in xrange(0,numEdges):
          source = A.Source(edge)
          target = A.Target(edge)
          G.add_edge(source,target)
        fig = plt.figure()
        plt.title(title)
        nx.draw(G)
        plt.draw()
        plt.show(block=False)
        return
      except:
        print 'Could not import networkx and matplotlib.pyplot'
  # Fall back to the built-in Display if we have not succeeded
  args = [A.obj,title]
  numMsExtra = 200
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElDisplay_i(*args)
    elif A.tag == sTag: lib.ElDisplay_s(*args)
    elif A.tag == dTag: lib.ElDisplay_d(*args)
    elif A.tag == cTag: lib.ElDisplay_c(*args)
    elif A.tag == zTag: lib.ElDisplay_z(*args)
    else: DataExcept()
    ProcessEvents(numMsExtra)
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElDisplayDist_i(*args)
    elif A.tag == sTag: lib.ElDisplayDist_s(*args)
    elif A.tag == dTag: lib.ElDisplayDist_d(*args)
    elif A.tag == cTag: lib.ElDisplayDist_c(*args)
    elif A.tag == zTag: lib.ElDisplayDist_z(*args)
    else: DataExcept()
    ProcessEvents(numMsExtra)
  elif type(A) is Graph:
    lib.ElDisplayGraph(*args)
    ProcessEvents(numMsExtra)
  elif type(A) is DistGraph:
    lib.ElDisplayDistGraph(*args)
    ProcessEvents(numMsExtra)
  elif type(A) is SparseMatrix:
    if   A.tag == iTag: lib.ElDisplaySparse_i(*args)
    elif A.tag == sTag: lib.ElDisplaySparse_s(*args)
    elif A.tag == dTag: lib.ElDisplaySparse_d(*args)
    elif A.tag == cTag: lib.ElDisplaySparse_c(*args)
    elif A.tag == zTag: lib.ElDisplaySparse_z(*args)
    else: DataExcept()
    ProcessEvents(numMsExtra)
  elif type(A) is DistSparseMatrix:
    if   A.tag == iTag: lib.ElDisplayDistSparse_i(*args)
    elif A.tag == sTag: lib.ElDisplayDistSparse_s(*args)
    elif A.tag == dTag: lib.ElDisplayDistSparse_d(*args)
    elif A.tag == cTag: lib.ElDisplayDistSparse_c(*args)
    elif A.tag == zTag: lib.ElDisplayDistSparse_z(*args)
    else: DataExcept()
    ProcessEvents(numMsExtra)
  else: TypeExcept()

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
def Spy(A,title='',tol=0):
  args = [A.obj,title,tol]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElSpy_i(*args)
    elif A.tag == sTag: lib.ElSpy_s(*args)
    elif A.tag == dTag: lib.ElSpy_d(*args)
    elif A.tag == cTag: lib.ElSpy_c(*args)
    elif A.tag == zTag: lib.ElSpy_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElSpyDist_i(*args)
    elif A.tag == sTag: lib.ElSpyDist_s(*args)
    elif A.tag == dTag: lib.ElSpyDist_d(*args)
    elif A.tag == cTag: lib.ElSpyDist_c(*args)
    elif A.tag == zTag: lib.ElSpyDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

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
def Read(A,filename,fileFormat=AUTO):
  args = [A.obj,filename,fileFormat]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElRead_i(*args)
    elif A.tag == sTag: lib.ElRead_s(*args)
    elif A.tag == dTag: lib.ElRead_d(*args)
    elif A.tag == cTag: lib.ElRead_c(*args)
    elif A.tag == zTag: lib.ElRead_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElReadDist_i(*args)
    elif A.tag == sTag: lib.ElReadDist_s(*args)
    elif A.tag == dTag: lib.ElReadDist_d(*args)
    elif A.tag == cTag: lib.ElReadDist_c(*args)
    elif A.tag == zTag: lib.ElReadDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

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
def Write(A,basename,fileFormat,title=''):
  args = [A.obj,basename,fileFormat,title]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElWrite_i(*args)
    elif A.tag == sTag: lib.ElWrite_s(*args)
    elif A.tag == dTag: lib.ElWrite_d(*args)
    elif A.tag == cTag: lib.ElWrite_c(*args)
    elif A.tag == zTag: lib.ElWrite_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElWriteDist_i(*args)
    elif A.tag == sTag: lib.ElWriteDist_s(*args)
    elif A.tag == dTag: lib.ElWriteDist_d(*args)
    elif A.tag == cTag: lib.ElWriteDist_c(*args)
    elif A.tag == zTag: lib.ElWriteDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
