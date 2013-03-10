import elem

elem.Initialize()

grid = elem.Grid()
A = elem.DistMat( grid )
A.MakeUniform( 8, 8 )
A.Print("Original A")

rank = grid.Rank()
if rank == 0:
  print "Running SVD..."
[s,V] = A.SVD()

A.Print("U")
s.Print("s")
V.Print("V")

A.MakeUniform( 8, 8 )
A.Print("New A")

if rank == 0:
  print "Running QR..."
R = A.ExplicitQR()

A.Print("Q")
R.Print("R")

elem.Finalize()
