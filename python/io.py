#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *
from El.blas_like import Copy, CopyFromRoot, CopyFromNonRoot, RealPart, ImagPart

# Input/Output
# ************

lib.ElPrint_i.argtypes = \
lib.ElPrint_s.argtypes = \
lib.ElPrint_d.argtypes = \
lib.ElPrint_c.argtypes = \
lib.ElPrint_z.argtypes = \
lib.ElPrintDist_i.argtypes = \
lib.ElPrintDist_s.argtypes = \
lib.ElPrintDist_d.argtypes = \
lib.ElPrintDist_c.argtypes = \
lib.ElPrintDist_z.argtypes = \
lib.ElPrintDistMultiVec_i.argtypes = \
lib.ElPrintDistMultiVec_s.argtypes = \
lib.ElPrintDistMultiVec_d.argtypes = \
lib.ElPrintDistMultiVec_c.argtypes = \
lib.ElPrintDistMultiVec_z.argtypes = \
lib.ElPrintGraph.argtypes = \
lib.ElPrintDistGraph.argtypes = \
lib.ElPrintSparse_i.argtypes = \
lib.ElPrintSparse_s.argtypes = \
lib.ElPrintSparse_d.argtypes = \
lib.ElPrintSparse_c.argtypes = \
lib.ElPrintSparse_z.argtypes = \
lib.ElPrintDistSparse_i.argtypes = \
lib.ElPrintDistSparse_s.argtypes = \
lib.ElPrintDistSparse_d.argtypes = \
lib.ElPrintDistSparse_c.argtypes = \
lib.ElPrintDistSparse_z.argtypes = \
  [c_void_p,c_char_p]

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
  elif type(A) is DistMultiVec:
    if   A.tag == iTag: lib.ElPrintDistMultiVec_i(*args)
    elif A.tag == sTag: lib.ElPrintDistMultiVec_s(*args)
    elif A.tag == dTag: lib.ElPrintDistMultiVec_d(*args)
    elif A.tag == cTag: lib.ElPrintDistMultiVec_c(*args)
    elif A.tag == zTag: lib.ElPrintDistMultiVec_z(*args)
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
def SetColorMap(colorMap):
  lib.ElSetColorMap(colorMap)

lib.ElGetColorMap.argtypes = [POINTER(c_uint)]
def ColorMap():
  colorMap = c_uint()
  lib.ElGetColorMap(pointer(colorMap))
  return colorMap

lib.ElSetNumDiscreteColors.argtypes = [iType]
def SetNumDiscreteColors(numColors):
  lib.ElSetNumDiscreteColors(numColors)

lib.ElNumDiscreteColors.argtypes = [POINTER(iType)]
def NumDiscreteColors():
  numDiscrete = iType()
  lib.ElNumDiscreteColors(pointer(numDiscrete))
  return numDiscrete

lib.ElProcessEvents.argtypes = [c_int]
def ProcessEvents(numMsecs):
  lib.ElProcessEvents(numMsecs)

lib.ElDisplay_i.argtypes = \
lib.ElDisplay_s.argtypes = \
lib.ElDisplay_d.argtypes = \
lib.ElDisplay_c.argtypes = \
lib.ElDisplay_z.argtypes = \
lib.ElDisplayDist_i.argtypes = \
lib.ElDisplayDist_s.argtypes = \
lib.ElDisplayDist_d.argtypes = \
lib.ElDisplayDist_c.argtypes = \
lib.ElDisplayDist_z.argtypes = \
lib.ElDisplayDistMultiVec_i.argtypes = \
lib.ElDisplayDistMultiVec_s.argtypes = \
lib.ElDisplayDistMultiVec_d.argtypes = \
lib.ElDisplayDistMultiVec_c.argtypes = \
lib.ElDisplayDistMultiVec_z.argtypes = \
lib.ElDisplayGraph.argtypes = \
lib.ElDisplayDistGraph.argtypes = \
lib.ElDisplaySparse_i.argtypes = \
lib.ElDisplaySparse_s.argtypes = \
lib.ElDisplaySparse_d.argtypes = \
lib.ElDisplaySparse_c.argtypes = \
lib.ElDisplaySparse_z.argtypes = \
lib.ElDisplayDistSparse_i.argtypes = \
lib.ElDisplayDistSparse_s.argtypes = \
lib.ElDisplayDistSparse_d.argtypes = \
lib.ElDisplayDistSparse_c.argtypes = \
lib.ElDisplayDistSparse_z.argtypes = \
  [c_void_p,c_char_p]

def DisplayPyPlot(A,title=''):
  isInline = 'inline' in mpl.get_backend()
  isVec = min(A.Height(),A.Width()) == 1
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
  if not isInline:
    plt.show(block=False)

def DisplayNetworkX(A,title=''):
  numEdges = A.NumEdges() 
  G = nx.DiGraph()
  for edge in xrange(numEdges):
    source = A.Source(edge)
    target = A.Target(edge)
    G.add_edge(source,target)
  fig = plt.figure()
  plt.title(title)
  nx.draw(G)
  plt.draw()
  if not isInline:
    plt.show(block=False)

def DisplayCxx(A,title=''):
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
  elif type(A) is DistMultiVec:
    if   A.tag == iTag: lib.ElDisplayDistMultiVec_i(*args)
    elif A.tag == sTag: lib.ElDisplayDistMultiVec_s(*args)
    elif A.tag == dTag: lib.ElDisplayDistMultiVec_d(*args)
    elif A.tag == cTag: lib.ElDisplayDistMultiVec_c(*args)
    elif A.tag == zTag: lib.ElDisplayDistMultiVec_z(*args)
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

def Display(A,title='',tryPython=True):
  if tryPython: 
    if type(A) is Matrix:
      if havePyPlot:
        DisplayPyPlot(A,title)
        return
    elif type(A) is DistMatrix:
      A_CIRC_CIRC = DistMatrix(A.tag,CIRC,CIRC,A.Grid())
      Copy(A,A_CIRC_CIRC)
      if A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root():
        Display(A_CIRC_CIRC.Matrix(),title,True)
      return
    elif type(A) is DistMultiVec:
      if mpi.Rank(A.Comm()) == 0:
        ASeq = Matrix(A.tag)
        CopyFromRoot(A,ASeq)
        Display(ASeq,title,True)
      else:
        CopyFromNonRoot(A)
      return
    elif type(A) is Graph:
      if haveNetworkX:
        DisplayNetworkX(A,title)
        return
    elif type(A) is DistGraph:
      if mpi.Rank(A.Comm()) == 0:
        ASeq = Graph()
        CopyFromRoot(A,ASeq)
        Display(ASeq,title,True)
      else:
        CopyFromNonRoot(A)
      return
    elif type(A) is SparseMatrix:
      ADense = Matrix(A.tag)
      Copy(A,ADense)
      Display(ADense,title,True)
      return
    elif type(A) is DistSparseMatrix:
      grid = Grid(A.Comm())
      ADense = DistMatrix(A.tag,MC,MR,grid)
      Copy(A,ADense)
      Display(ADense,title,True)
      return

  # Fall back to the internal Display routine
  DisplayCxx(A,title)

lib.ElSpy_i.argtypes = \
lib.ElSpyDist_i.argtypes = \
  [c_void_p,c_char_p,iType]
lib.ElSpy_s.argtypes = \
lib.ElSpyDist_s.argtypes = \
  [c_void_p,c_char_p,sType]
lib.ElSpy_d.argtypes = \
lib.ElSpyDist_d.argtypes = \
  [c_void_p,c_char_p,dType]
lib.ElSpy_c.argtypes = \
lib.ElSpyDist_c.argtypes = \
  [c_void_p,c_char_p,sType]
lib.ElSpy_z.argtypes = \
lib.ElSpyDist_z.argtypes = \
  [c_void_p,c_char_p,dType]

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

lib.ElRead_i.argtypes = \
lib.ElRead_s.argtypes = \
lib.ElRead_d.argtypes = \
lib.ElRead_c.argtypes = \
lib.ElRead_z.argtypes = \
lib.ElReadDist_i.argtypes = \
lib.ElReadDist_s.argtypes = \
lib.ElReadDist_d.argtypes = \
lib.ElReadDist_c.argtypes = \
lib.ElReadDist_z.argtypes = \
  [c_void_p,c_char_p,c_uint]

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

lib.ElWrite_i.argtypes = \
lib.ElWrite_s.argtypes = \
lib.ElWrite_d.argtypes = \
lib.ElWrite_c.argtypes = \
lib.ElWrite_z.argtypes = \
lib.ElWriteDist_i.argtypes = \
lib.ElWriteDist_s.argtypes = \
lib.ElWriteDist_d.argtypes = \
lib.ElWriteDist_c.argtypes = \
lib.ElWriteDist_z.argtypes = \
  [c_void_p,c_char_p,c_uint,c_char_p]

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
