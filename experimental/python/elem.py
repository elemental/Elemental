import ctypes, sys
from ctypes import cdll
elem = cdll.LoadLibrary('./libElem.so')

# TODO: Figure out how to wrap MPI_Comm

class Grid(object):
  def __init__(self):
    self.obj = elem.DefaultGrid()
  @classmethod
  def FromPointer(cls,gridPtr):
    grid = Grid()
    grid.obj = gridPtr
    return grid
  def Height(self):
    return elem.GridHeight(self.obj)
  def Width(self):
    return elem.GridWidth(self.obj)
  def Size(self):
    return elem.GridSize(self.obj)
  def Row(self):
    return elem.GridRow(self.obj)
  def Col(self):
    return elem.GridCol(self.obj)
  def Rank(self):
    return elem.GridRank(self.obj)

class DistMat(object):
  def __init__(self,grid):
    self.obj = elem.CreateDistMat(grid.obj)
    self.grid = grid
  def Attach(self,height,width,colAlign,rowAlign,buf,ldim,grid):
    elem.AttachToDistMat(self.obj,height,width,colAlign,rowAlign,buf,ldim,grid)
    self.grid = grid
  def Print(self):
    elem.PrintDistMat(self.obj)
  def MakeUniform(self,height,width):
    elem.UniformDistMat(self.obj,height,width)
  def Free(self):
    elem.FreeDistMat( ctypes.byref(self.obj) )
    elem.FreeGrid( ctypes.byref(self.grid) )
  def Grid(self):
    return self.grid
  def SVD(self):
    s = DistMat_VR_STAR(self.grid)
    V = DistMat(self.grid) 
    elem.SVD( self.obj, s.obj, V.obj )
    return s, V

class DistMat_VR_STAR(object):
  def __init__(self,grid):
    self.obj = elem.CreateDistMat_VR_STAR(grid.obj)
    self.grid = grid
  def Print(self):
    elem.PrintDistMat_VR_STAR(self.obj)
  def Free(self):
    elem.FreeDistMat_VR_STAR( ctypes.byref(self.obj) )
    elem.FreeGrid( ctypes.byref(self.grid) )
  def Grid(self):
    return self.grid

argc = ctypes.c_int(len(sys.argv))
_argv = ""
for arg in sys.argv:
  _argv += arg + ' '
argv = ctypes.pointer(ctypes.c_char_p(_argv))
elem.Initialize( ctypes.byref(argc), ctypes.byref(argv) )
grid = Grid()
rank = grid.Rank()

if rank == 0:
  print "Constructing DistMat"
A = DistMat( grid )

if rank == 0:
  print "Making uniform"
A.MakeUniform( 5, 5 )

if rank == 0:
  print "Printing A"
A.Print()

if rank == 0:
  print "Running SVD"
[s,V] = A.SVD()

if rank == 0:
  print "Printing SVD"
A.Print()
s.Print()
V.Print()

elem.Finalize()
