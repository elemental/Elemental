#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

n0=100
n1=100
output = False
display = False
worldSize = El.mpi.WorldSize()
worldRank = El.mpi.WorldRank()

# Place two 2D finite-difference matrices next to each other
# and make the last column dense
def ConcatFD2D(N0,N1):
  A = El.DistSparseMatrix(El.zTag)
  height = N0*N1
  width = 2*N0*N1
  A.Resize(height,width)
  localHeight = A.LocalHeight()
  A.Reserve(11*localHeight)
  for sLoc in xrange(localHeight):
    s = A.GlobalRow(sLoc)
    x0 = s % N0
    x1 = s / N0
    sRel = s + N0*N1

    A.QueueUpdate( s, s,    El.ComplexDouble(1,1), passive=True )
    A.QueueUpdate( s, sRel, El.ComplexDouble(20,2), passive=True )
    if x0 > 0:
      A.QueueUpdate( s, s-1,    El.ComplexDouble(-1,3), passive=True )
      A.QueueUpdate( s, sRel-1, El.ComplexDouble(-17,4), passive=True )
    if x0+1 < N0:
      A.QueueUpdate( s, s+1,    El.ComplexDouble(2,5), passive=True )
      A.QueueUpdate( s, sRel+1, El.ComplexDouble(-20,6), passive=True )
    if x1 > 0:
      A.QueueUpdate( s, s-N0,    El.ComplexDouble(-30,7), passive=True )
      A.QueueUpdate( s, sRel-N0, El.ComplexDouble(-3,8), passive=True )
    if x1+1 < N1:
      A.QueueUpdate( s, s+N0,    El.ComplexDouble(4,9), passive=True )
      A.QueueUpdate( s, sRel+N0, El.ComplexDouble(3,10), passive=True )

  A.ProcessLocalQueues()
  return A

A = ConcatFD2D(n0,n1)
b = El.DistMultiVec(El.zTag)
#El.Gaussian( b, n0*n1, 1 )
El.Ones( b, n0*n1, 1 )

ctrl = El.BPCtrl_z()
ctrl.ipmCtrl.mehrotraCtrl.minTol = 1e-4
ctrl.ipmCtrl.mehrotraCtrl.targetTol = 1e-8
ctrl.ipmCtrl.mehrotraCtrl.time = True
ctrl.ipmCtrl.mehrotraCtrl.progress = True
ctrl.ipmCtrl.mehrotraCtrl.solveCtrl.progress = True
startBP = El.mpi.Time()
x = El.BP( A, b, ctrl )
endBP = El.mpi.Time()
if worldRank == 0:
  print "BP time:", endBP-startBP, "seconds"
if display:
  El.Display( x, "x" )
if output:
  El.Print( x, "x" )

xOneNorm = El.EntrywiseNorm( x, 1 )
e = El.DistMultiVec(El.zTag)
El.Copy( b, e )
El.Multiply \
( El.NORMAL, El.ComplexDouble(-1), A, x, El.ComplexDouble(1), e )
eTwoNorm = El.Nrm2( e )
if worldRank == 0:
  print "|| x ||_1       =", xOneNorm
  print "|| A x - b ||_2 =", eTwoNorm

# Require the user to press a button before the figures are closed
El.Finalize()
if worldSize == 1:
  raw_input('Press Enter to exit')
