#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El, time

n0 = n1 = 100
display = False
output = False
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
      A.QueueLocalUpdate( sLoc, s, 11.1 )
      if x0 > 0:
        A.QueueLocalUpdate( sLoc, s-1, -1.2 )
      if x0+1 < N0:
        A.QueueLocalUpdate( sLoc, s+1, 2.3 )
      if x1 > 0:
        A.QueueLocalUpdate( sLoc, s-N0, -3.4 )
      if x1+1 < N1:
        A.QueueLocalUpdate( sLoc, s+N0, 4.5 )
    else:
      sRel = s-N0*N1
      x0 = sRel % N0
      x1 = sRel / N0
      A.QueueLocalUpdate( sLoc, sRel, -20.1 )
      if x0 > 0:
        A.QueueLocalUpdate( sLoc, sRel-1, -1.7 )
      if x0+1 < N0:
        A.QueueLocalUpdate( sLoc, sRel+1, -2.2 )
      if x1 > 0:
        A.QueueLocalUpdate( sLoc, sRel-N0, -3.3 )
      if x1+1 < N1:
        A.QueueLocalUpdate( sLoc, sRel+N0, 3.4 )

    # The dense last column
    A.QueueLocalUpdate( sLoc, width-1, -10/height );

  A.ProcessQueues()
  return A

A = StackedFD2D(n0,n1)
if display:
  El.Display( A, "A" )
if output:
  El.Print( A, "A" )

y = El.DistMultiVec()
El.Uniform( y, 2*n0*n1, 1 )
if display:
  El.Display( y, "y" )
if output:
  El.Print( y, "y" )

AAdj = El.DistSparseMatrix()
El.Adjoint( A, AAdj )
z = El.DistMultiVec()
El.Uniform( z, A.Width(), 1 )
El.Multiply( El.ADJOINT, 1., A, y, 0., z )
zNrm2 = El.FrobeniusNorm( z )
El.Multiply( El.NORMAL, -1., AAdj, y, 1., z )
eNrm2 = El.FrobeniusNorm( z )
if worldRank == 0:
  print "|| A^H y ||_2 =", zNrm2
  print "|| error ||_2 =", eNrm2

# Require the user to press a button before the figures are closed
El.Finalize()
if worldSize == 1:
  raw_input('Press Enter to exit')
