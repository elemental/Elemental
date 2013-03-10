import ctypes,elem

elem.Initialize()

grid = elem.Grid()
A = elem.DistMat( grid )
A.Resize(8,8)
localHeight = A.LocalHeight()
localWidth = A.LocalWidth()
colShift = A.ColShift()
rowShift = A.RowShift()
colStride = A.ColStride()
rowStride = A.RowStride()
for jLocal in xrange(0,localWidth):
  j = rowShift + jLocal*rowStride
  for iLocal in xrange(0,localHeight):
    i = colShift + iLocal*colStride
    A.SetLocal(iLocal,jLocal,ctypes.c_double(i-j))
A.Print("Original A")

rank = grid.Rank()
if rank == 0:
  print "Running SVD..."
[s,V] = elem.SVD(A)

A.Print("U")
s.Print("s")
V.Print("V")

elem.UniformDistMat( A, 8, 8 )
A.Print("New A")

if rank == 0:
  print "Running QR..."
R = elem.ExplicitQR(A)

A.Print("Q")
R.Print("R")

elem.Finalize()
