#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

n0 = n1 = 20

def Laplacian(xSize,ySize):
  A = El.DistSparseMatrix(El.dTag)
  A.Resize(xSize*ySize,xSize*ySize)
  firstLocalRow = A.FirstLocalRow()
  localHeight = A.LocalHeight()
  A.Reserve(5*localHeight)
  hxInvSq = (1.*(xSize+1))**2
  hyInvSq = (1.*(ySize+1))**2
  for sLoc in xrange(localHeight):
    s = firstLocalRow + sLoc
    x = s % xSize
    y = s / xSize
    A.QueueLocalUpdate( sLoc, s, 2*(hxInvSq+hyInvSq) )
    if x != 0:       A.QueueLocalUpdate( sLoc, s-1,     -hxInvSq )
    if x != xSize-1: A.QueueLocalUpdate( sLoc, s+1,     -hxInvSq )
    if y != 0:       A.QueueLocalUpdate( sLoc, s-xSize, -hyInvSq )
    if y != ySize-1: A.QueueLocalUpdate( sLoc, s+xSize, -hyInvSq )

  A.MakeConsistent()
  return A

A = Laplacian(n0,n1)
x = El.DistMultiVec()
y = El.DistMultiVec()
El.Uniform( x, n0*n1, 1 )
El.Copy( x, y )

yNrm = El.Nrm2(y)
rank = El.mpi.WorldRank()
if rank == 0:
  print "|| y ||_2 =", yNrm

El.Display( A, "Laplacian" )
El.Display( A.DistGraph(), "Laplacian graph" )
El.Display( y, "y" )

El.SymmetricSolve(A,x)
El.Display( x, "x" )

xNrm = El.Nrm2(x)
if rank == 0:
  print "|| x ||_2 =", xNrm

El.SparseMultiply(El.NORMAL,-1.,A,x,1.,y)
El.Display( y, "A x - y" )
eNrm = El.Nrm2(y)
if rank == 0:
  print "|| A x - y ||_2 / || y ||_2 =", eNrm/yNrm

# Require the user to press a button before the figures are closed
commSize = El.mpi.Size( El.mpi.COMM_WORLD() )
El.Finalize()
if commSize == 1:
  raw_input('Press Enter to exit')
