#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

n0 = n1 = 100
display = False
worldRank = El.mpi.WorldRank()
worldSize = El.mpi.WorldSize()

# Stack two 2D finite-difference matrices on top of each other
# and make the last column dense
#
# NOTE: Increasing the magnitudes of the off-diagonal entries by an order of
#       magnitude makes the condition number vastly higher.
def StackedFD2D(N0,N1):
  A = El.DistSparseMatrix()
  height = 2*N0*N1
  width = N0*N1
  A.Resize(height,width)
  localHeight = A.LocalHeight()
  A.Reserve(6*localHeight)
  for sLoc in xrange(localHeight):
    s = A.GlobalRow(sLoc)
    if s < N0*N1:
      x0 = s % N0
      x1 = s / N0
      A.QueueLocalUpdate( sLoc, s, 11 )
      if x0 > 0:
        A.QueueLocalUpdate( sLoc, s-1, -1 )
      if x0+1 < N0:
        A.QueueLocalUpdate( sLoc, s+1, 2 )
      if x1 > 0:
        A.QueueLocalUpdate( sLoc, s-N0, -3 )
      if x1+1 < N1:
        A.QueueLocalUpdate( sLoc, s+N0, 4 )
    else:
      sRel = s-N0*N1
      x0 = sRel % N0
      x1 = sRel / N0
      A.QueueLocalUpdate( sLoc, sRel, -20 )
      if x0 > 0:
        A.QueueLocalUpdate( sLoc, sRel-1, -1.7 )
      if x0+1 < N0:
        A.QueueLocalUpdate( sLoc, sRel+1, -2 )
      if x1 > 0:
        A.QueueLocalUpdate( sLoc, sRel-N0, -3 )
      if x1+1 < N1:
        A.QueueLocalUpdate( sLoc, sRel+N0, 3 )

    # The dense last column
    A.QueueLocalUpdate( sLoc, width-1, -10/height );

  A.ProcessQueues()
  return A

# Stack two 2D finite-difference-like matrices on top of each other
A = StackedFD2D(n0,n1)
if display:
  El.Display( A, "A" )
  El.Display( A[0:n0*n1,0:n0*n1], "AT" )
  El.Display( A[n0*n1:2*n0*n1,0:n0*n1], "AB" )

# Regularize with 3 I
G = El.DistSparseMatrix()
El.Identity( G, n0*n1, n0*n1 )
El.Scale( 3., G )

y = El.DistMultiVec()
El.Uniform( y, 2*n0*n1, 1 )
if display:
  El.Display( y, "y" )
yNrm = El.Nrm2(y)
if worldRank == 0:
  print "|| y ||_2 =", yNrm

startLS = El.mpi.Time()
x = El.Tikhonov(A,y,G)
endLS = El.mpi.Time()
if worldRank == 0:
  print "Tikhonov time:", endLS-startLS, "seconds"
xNrm = El.Nrm2(x)
if display:
  El.Display( x, "x" )
if worldRank == 0:
  print "|| x ||_2 =", xNrm
El.Multiply(El.NORMAL,-1.,A,x,1.,y)
if display:
  El.Display( y, "A x - y" )
eNrm = El.Nrm2(y)
if worldRank == 0:
  print "|| A x - y ||_2 / || y ||_2 =", eNrm/yNrm

# Require the user to press a button before the figures are closed

El.Finalize()
if worldSize == 1:
  raw_input('Press Enter to exit')
