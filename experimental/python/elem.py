import ctypes, sys
from ctypes import cdll
lib = cdll.LoadLibrary('./libElem.so')

# TODO: Figure out how to wrap MPI_Comm

class Grid(object):
  def __init__(self):
    self.obj = lib.DefaultGrid()
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

class DistMat(object):
  def __init__(self,grid):
    self.obj = lib.CreateDistMat(grid.obj)
    self.grid = grid
  def Attach(self,height,width,colAlign,rowAlign,buf,ldim,grid):
    lib.AttachToDistMat(self.obj,height,width,colAlign,rowAlign,buf,ldim,grid)
    self.grid = grid
  def Print(self):
    lib.PrintDistMat(self.obj)
  def MakeUniform(self,height,width):
    lib.UniformDistMat(self.obj,height,width)
  def Free(self):
    lib.FreeDistMat( self.obj )
    self.obj = 0
  def Grid(self):
    return self.grid
  def SVD(self):
    s = DistMat_VR_STAR(self.grid)
    V = DistMat(self.grid) 
    lib.SVD( self.obj, s.obj, V.obj )
    return s, V

class DistMat_VR_STAR(object):
  def __init__(self,grid):
    self.obj = lib.CreateDistMat_VR_STAR(grid.obj)
    self.grid = grid
  def Print(self):
    lib.PrintDistMat_VR_STAR(self.obj)
  def Free(self):
    lib.FreeDistMat_VR_STAR( self.obj )
    self.obj = 0
  def Grid(self):
    return self.grid

def Initialize():
  argc = ctypes.c_int(len(sys.argv))
  _argv = ""
  for arg in sys.argv:
    _argv += arg + ' '
  argv = ctypes.pointer(ctypes.c_char_p(_argv))
  lib.Initialize( ctypes.byref(argc), ctypes.byref(argv) )

def Finalize():
  lib.Finalize()
