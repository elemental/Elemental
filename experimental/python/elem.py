import ctypes, sys
from ctypes import cdll
lib = cdll.LoadLibrary('./libElem.so')

#
# Core components
#

def Initialize():
  argc = ctypes.c_int(len(sys.argv))
  _argv = ""
  for arg in sys.argv:
    _argv += arg + ' '
  argv = ctypes.pointer(ctypes.c_char_p(_argv))
  lib.Initialize( ctypes.byref(argc), ctypes.byref(argv) )

def Finalize():
  lib.Finalize()

class Grid(object):
  def __init__(self):
    self.obj = lib.DefaultGrid()
  # TODO: Figure out how to wrap MPI_Comm
  @classmethod
  def FromPointer(cls,gridPtr):
    grid = Grid()
    grid.obj = gridPtr
    return grid
  def Free(self):
    lib.FreeGrid( self.obj )
    self.obj = 0
  def Height(self):
    return lib.GridHeight(self.obj)
  def Width(self):
    return lib.GridWidth(self.obj)
  def Size(self):
    return lib.GridSize(self.obj)
  def Row(self):
    return lib.GridRow(self.obj)
  def Col(self):
    return lib.GridCol(self.obj)
  def Rank(self):
    return lib.GridRank(self.obj)

class Mat(object):
  def __init__(self):
    self.obj = lib.CreateMat()
  def Resize(self,height,width):
    lib.ResizeMat(self.obj,height,width)
  def Attach(self,height,width,buf,ldim):
    lib.AttachToMat(self.obj,height,width,buf,ldim)
  def Free(self):
    lib.FreeMat( self.obj )
    self.obj = 0
  def Print(self,string):
    print string
    lib.PrintMat(self.obj)
  def Height(self):
    return lib.MatHeight(self.obj)
  def Width(self):
    return lib.MatWidth(self.obj)
  def LDim(self):
    return lib.MatLDim(self.obj)
  def Get(self,i,j):
    return lib.GetMatEntry(self.obj,i,j)
  def Set(self,i,j,alpha):
    lib.SetMatEntry(self.obj,i,j,alpha)

class CpxMat(object):
  def __init__(self):
    self.obj = lib.CreateCpxMat()
  def Resize(self,height,width):
    lib.ResizeCpxMat(self.obj,height,width)
  def Attach(self,height,width,buf,ldim):
    lib.AttachToCpxMat(self.obj,height,width,buf,ldim)
  def Free(self):
    lib.FreeCpxMat( self.obj )
    self.obj = 0
  def Print(self,string):
    print string
    lib.PrintCpxMat(self.obj)
  def Height(self):
    return lib.CpxMatHeight(self.obj)
  def Width(self):
    return lib.CpxMatWidth(self.obj)
  def LDim(self):
    return lib.CpxMatLDim(self.obj)
  # TODO: Figure out how to pass complex values

class DistMat(object):
  def __init__(self,grid):
    self.obj = lib.CreateDistMat(grid.obj)
    self.grid = grid
  def Grid(self):
    return self.grid
  def Resize(self,height,width):
    lib.ResizeDistMat(self.obj,height,width)
  def Attach(self,height,width,colAlign,rowAlign,buf,ldim,grid):
    lib.AttachToDistMat(self.obj,height,width,colAlign,rowAlign,buf,ldim,grid)
    self.grid = grid
  def Free(self):
    lib.FreeDistMat( self.obj )
    self.obj = 0
  def Print(self,string):
    if self.grid.Rank() == 0:
      print string
    lib.PrintDistMat(self.obj)
  def Height(self):
    return lib.DistMatHeight(self.obj)
  def Width(self):
    return lib.DistMatWidth(self.obj)
  def LocalHeight(self):
    return lib.DistMatLocalHeight(self.obj)
  def LocalWidth(self):
    return lib.DistMatLocalWidth(self.obj)
  def LDim(self):
    return lib.DistMatLDim(self.obj)
  def ColShift(self):
    return lib.DistMatColShift(self.obj)
  def RowShift(self):
    return lib.DistMatRowShift(self.obj)
  def ColStride(self):
    return lib.DistMatColStride(self.obj)
  def RowStride(self):
    return lib.DistMatRowStride(self.obj) 
  def Get(self,i,j):
    return lib.GetDistMatEntry(self.obj,i,j)
  def Set(self,i,j,alpha):
    lib.SetDistMatEntry(self.obj,i,j,alpha)
  def GetLocal(self,iLocal,jLocal):
    return lib.GetLocalDistMatEntry(self.obj,iLocal,jLocal)
  def SetLocal(self,iLocal,jLocal,alpha):
    lib.SetLocalDistMatEntry(self.obj,iLocal,jLocal,alpha)

class CpxDistMat(object):
  def __init__(self,grid):
    self.obj = lib.CreateCpxDistMat(grid.obj)
    self.grid = grid
  def Grid(self):
    return self.grid
  def Resize(self,height,width):
    lib.ResizeCpxDistMat(self.obj,height,width)
  def Attach(self,height,width,colAlign,rowAlign,buf,ldim,grid):
    lib.AttachToCpxDistMat(self.obj,height,width,colAlign,rowAlign,buf,ldim,grid)
    self.grid = grid
  def Free(self):
    lib.FreeCpxDistMat( self.obj )
    self.obj = 0
  def Print(self,string):
    if self.grid.Rank() == 0:
      print string
    lib.PrintCpxDistMat(self.obj)
  def Height(self):
    return lib.CpxDistMatHeight(self.obj)
  def Width(self):
    return lib.CpxDistMatWidth(self.obj)
  def LocalHeight(self):
    return lib.CpxDistMatLocalHeight(self.obj)
  def LocalWidth(self):
    return lib.CpxDistMatLocalWidth(self.obj)
  def LDim(self):
    return lib.CpxDistMatLDim(self.obj)
  def ColShift(self):
    return lib.CpxDistMatColShift(self.obj)
  def RowShift(self):
    return lib.CpxDistMatRowShift(self.obj)
  def ColStride(self):
    return lib.CpxDistMatColStride(self.obj)
  def RowStride(self):
    return lib.CpxDistMatRowStride(self.obj) 

class DistMat_VC_STAR(object):
  def __init__(self,grid):
    self.obj = lib.CreateDistMat_VC_STAR(grid.obj)
    self.grid = grid
  def Grid(self):
    return self.grid
  def Resize(self,height,width):
    lib.ResizeDistMat_VC_STAR(self.obj,height,width)
  def Free(self):
    lib.FreeDistMat_VC_STAR( self.obj )
    self.obj = 0
  def Print(self,string):
    if self.grid.Rank() == 0:
      print string
    lib.PrintDistMat_VC_STAR(self.obj)
  def Height(self):
    return lib.DistMatHeight_VC_STAR(self.obj)
  def Width(self):
    return lib.DistMatWidth_VC_STAR(self.obj)
  def LocalHeight(self):
    return lib.DistMatLocalHeight_VC_STAR(self.obj)
  def LocalWidth(self):
    return lib.DistMatLocalWidth_VC_STAR(self.obj)
  def LDim(self):
    return lib.DistMatLDim_VC_STAR(self.obj)
  def ColShift(self):
    return lib.DistMatColShift_VC_STAR(self.obj)
  def RowShift(self):
    return lib.DistMatRowShift_VC_STAR(self.obj)
  def ColStride(self):
    return lib.DistMatColStride_VC_STAR(self.obj)
  def RowStride(self):
    return lib.DistMatRowStride_VC_STAR(self.obj) 
  def Get(self,i,j):
    return lib.GetDistMatEntry_VC_STAR(self.obj,i,j)
  def Set(self,i,j,alpha):
    lib.SetDistMatEntry_VC_STAR(self.obj,i,j,alpha)
  def GetLocal(self,iLocal,jLocal):
    return lib.GetLocalDistMatEntry_VC_STAR(self.obj,iLocal,jLocal)
  def SetLocal(self,iLocal,jLocal,alpha):
    lib.SetLocalDistMatEntry_VC_STAR(self.obj,iLocal,jLocal,alpha)

class CpxDistMat_VC_STAR(object):
  def __init__(self,grid):
    self.obj = lib.CreateCpxDistMat_VC_STAR(grid.obj)
    self.grid = grid
  def Grid(self):
    return self.grid
  def Resize(self,height,width):
    lib.ResizeCpxDistMat_VC_STAR(self.obj,height,width)
  def Free(self):
    lib.FreeCpxDistMat_VC_STAR( self.obj )
    self.obj = 0
  def Print(self,string):
    if self.grid.Rank() == 0:
      print string
    lib.PrintCpxDistMat_VC_STAR(self.obj)
  def Height(self):
    return lib.CpxDistMatHeight_VC_STAR(self.obj)
  def Width(self):
    return lib.CpxDistMatWidth_VC_STAR(self.obj)
  def LocalHeight(self):
    return lib.CpxDistMatLocalHeight_VC_STAR(self.obj)
  def LocalWidth(self):
    return lib.CpxDistMatLocalWidth_VC_STAR(self.obj)
  def LDim(self):
    return lib.CpxDistMatLDim_VC_STAR(self.obj)
  def ColShift(self):
    return lib.CpxDistMatColShift_VC_STAR(self.obj)
  def RowShift(self):
    return lib.CpxDistMatRowShift_VC_STAR(self.obj)
  def ColStride(self):
    return lib.CpxDistMatColStride_VC_STAR(self.obj)
  def RowStride(self):
    return lib.CpxDistMatRowStride_VC_STAR(self.obj) 
  # TODO: Figure out how to pass complex data

class DistMat_VR_STAR(object):
  def __init__(self,grid):
    self.obj = lib.CreateDistMat_VR_STAR(grid.obj)
    self.grid = grid
  def Grid(self):
    return self.grid
  def Resize(self,height,width):
    lib.ResizeDistMat_VR_STAR(self.obj,height,width)
  def Free(self):
    lib.FreeDistMat_VR_STAR( self.obj )
    self.obj = 0
  def Print(self,string):
    if self.grid.Rank() == 0:
      print string
    lib.PrintDistMat_VR_STAR(self.obj)
  def Height(self):
    return lib.DistMatHeight_VR_STAR(self.obj)
  def Width(self):
    return lib.DistMatWidth_VR_STAR(self.obj)
  def LocalHeight(self):
    return lib.DistMatLocalHeight_VR_STAR(self.obj)
  def LocalWidth(self):
    return lib.DistMatLocalWidth_VR_STAR(self.obj)
  def LDim(self):
    return lib.DistMatLDim_VR_STAR(self.obj)
  def ColShift(self):
    return lib.DistMatColShift_VR_STAR(self.obj)
  def RowShift(self):
    return lib.DistMatRowShift_VR_STAR(self.obj)
  def ColStride(self):
    return lib.DistMatColStride_VR_STAR(self.obj)
  def RowStride(self):
    return lib.DistMatRowStride_VR_STAR(self.obj) 
  def Get(self,i,j):
    return lib.GetDistMatEntry_VR_STAR(self.obj,i,j)
  def Set(self,i,j,alpha):
    lib.SetDistMatEntry_VR_STAR(self.obj,i,j,alpha)
  def GetLocal(self,iLocal,jLocal):
    return lib.GetLocalDistMatEntry_VR_STAR(self.obj,iLocal,jLocal)
  def SetLocal(self,iLocal,jLocal,alpha):
    lib.SetLocalDistMatEntry_VR_STAR(self.obj,iLocal,jLocal,alpha)

class CpxDistMat_VR_STAR(object):
  def __init__(self,grid):
    self.obj = lib.CreateCpxDistMat_VR_STAR(grid.obj)
    self.grid = grid
  def Grid(self):
    return self.grid
  def Resize(self,height,width):
    lib.ResizeCpxDistMat_VR_STAR(self.obj,height,width)
  def Free(self):
    lib.FreeCpxDistMat_VR_STAR( self.obj )
    self.obj = 0
  def Print(self,string):
    if self.grid.Rank() == 0:
      print string
    lib.PrintCpxDistMat_VR_STAR(self.obj)
  def Height(self):
    return lib.CpxDistMatHeight_VR_STAR(self.obj)
  def Width(self):
    return lib.CpxDistMatWidth_VR_STAR(self.obj)
  def LocalHeight(self):
    return lib.CpxDistMatLocalHeight_VR_STAR(self.obj)
  def LocalWidth(self):
    return lib.CpxDistMatLocalWidth_VR_STAR(self.obj)
  def LDim(self):
    return lib.CpxDistMatLDim_VR_STAR(self.obj)
  def ColShift(self):
    return lib.CpxDistMatColShift_VR_STAR(self.obj)
  def RowShift(self):
    return lib.CpxDistMatRowShift_VR_STAR(self.obj)
  def ColStride(self):
    return lib.CpxDistMatColStride_VR_STAR(self.obj)
  def RowStride(self):
    return lib.CpxDistMatRowStride_VR_STAR(self.obj) 
  # TODO: Figure out how to pass complex data

#
# LAPACK-like routines
#

def SVD(A):
  s = DistMat_VR_STAR(A.grid)
  V = DistMat(A.grid) 
  lib.SVD( A.obj, s.obj, V.obj )
  return s, V
def CpxSVD(A):
  s = DistMat_VR_STAR(A.grid)
  V = CpxDistMat(A.grid) 
  lib.CpxSVD( A.obj, s.obj, V.obj )
  return s, V

def SingularValues(A):
  s = DistMat_VR_STAR(A.grid)
  lib.SingularValues( A.obj, s.obj )
  return s
def CpxSingularValues(A):
  s = DistMat_VR_STAR(A.grid)
  lib.SingularValues( A.obj, s.obj )
  return s

def ExplicitQR(A):
  R = DistMat(A.grid)
  lib.ExplicitQR( A.obj, R.obj )
  return R
def CpxExplicitQR(A):
  R = CpxDistMat(A.grid)
  lib.CpxExplicitQR( A.obj, R.obj )
  return R

#
# Special matrices
#

def UniformMat(A,height,width):
  lib.UniformMat(A.obj,height,width)
def UniformCpxMat(A,height,width):
  lib.UniformCpxMat(A.obj,height,width)
def UniformDistMat(A,height,width):
  lib.UniformDistMat(A.obj,height,width)
def UniformCpxDistMat(A,height,width):
  lib.UniformCpxDistMat(A.obj,height,width)

