#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

n0 = n1 = 20

def AsymmSparse(xSize,ySize):
  A = El.DistSparseMatrix()
  A.Resize(xSize*ySize,xSize*ySize)
  firstLocalRow = A.FirstLocalRow()
  localHeight = A.LocalHeight()
  A.Reserve(5*localHeight)
  for sLoc in xrange(localHeight):
    s = firstLocalRow + sLoc
    x = s % xSize
    y = s / xSize
    A.QueueLocalUpdate( sLoc, s, 11 )
    if x != 0:       A.QueueLocalUpdate( sLoc, s-1,     1 )
    if x != xSize-1: A.QueueLocalUpdate( sLoc, s+1,     2 )
    if y != 0:       A.QueueLocalUpdate( sLoc, s-xSize, 3 )
    if y != ySize-1: A.QueueLocalUpdate( sLoc, s+xSize, 4 )

  A.MakeConsistent()
  return A

A = AsymmSparse(n0,n1)

b = El.DistMultiVec()
c = El.DistMultiVec()
x = El.DistMultiVec()
l = El.DistMultiVec()
s = El.DistMultiVec()
El.Uniform(b,n0*n1,1,0.5,0.4999)
El.Uniform(c,n0*n1,1,0.5,0.4999)
El.Uniform(x,n0*n1,1,0.5,0.4999)
El.Uniform(l,n0*n1,1,0.5,0.4999)
El.Uniform(s,n0*n1,1,0.5,0.4999)

mu = El.Dot(x,s) / (1.*n0*n1)
sigma = 0.9
tau = mu*sigma
J, y = El.LinearProgramFormNormalSystem(A,b,c,x,l,s,tau)

El.Display( b, "b" )
El.Display( c, "c" )
El.Display( x, "x" )
El.Display( l, "l" )
El.Display( s, "s" )
El.Display( A, "A" )

El.Display( J, "J" )
El.Display( y, "y" )

dx, dl, ds = El.LinearProgramSolveNormalSystem(A,b,c,x,l,s,tau,J,y)

El.Display( dx, "dx" )
El.Display( dl, "dl" )
El.Display( ds, "ds" )

gamma = 0.5
beta = 1.1
psi = 100
progress = True
alpha = El.LinearProgramIPFLineSearch(A,b,c,x,l,s,dx,dl,ds, \
                                      gamma,beta,psi,progress)
if El.mpi.WorldRank() == 0:
  print "alpha =", alpha

# Require the user to press a button before the figures are closed
commSize = El.mpi.Size( El.mpi.COMM_WORLD() )
El.Finalize()
if commSize == 1:
  raw_input('Press Enter to exit')
