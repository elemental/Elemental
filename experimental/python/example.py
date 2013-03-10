import elem

elem.Initialize()

grid = elem.Grid()
rank = grid.Rank()

if rank == 0:
  print "Constructing DistMat"
A = elem.DistMat( grid )

if rank == 0:
  print "Making uniform"
A.MakeUniform( 8, 8 )

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
