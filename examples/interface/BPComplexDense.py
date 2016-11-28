#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

n0=15
n1=15
output = False
display = True
worldSize = El.mpi.WorldSize()
worldRank = El.mpi.WorldRank()

# Place two 2D finite-difference matrices next to each other
# and make the last column dense
def ConcatFD2D(N0,N1):
  A = El.DistMatrix(El.zTag)
  height = N0*N1
  width = 2*N0*N1
  El.Zeros(A,height,width)
  localHeight = A.LocalHeight()
  A.Reserve(11*localHeight)
  for iLoc in xrange(localHeight):
    i = A.GlobalRow(iLoc)
    x0 = i % N0
    x1 = i / N0
    iRel = i + N0*N1

    A.Update( i, i,    El.ComplexDouble(1,1) )
    A.Update( i, iRel, El.ComplexDouble(20,2) )
    if x0 > 0:
      A.Update( i, i-1,    El.ComplexDouble(-1,3) )
      A.Update( i, iRel-1, El.ComplexDouble(-17,4) )
    if x0+1 < N0:
      A.Update( i, i+1,    El.ComplexDouble(2,5) )
      A.Update( i, iRel+1, El.ComplexDouble(-20,6) )
    if x1 > 0:
      A.Update( i, i-N0,    El.ComplexDouble(-30,7) )
      A.Update( i, iRel-N0, El.ComplexDouble(-3,8) )
    if x1+1 < N1:
      A.Update( i, i+N0,    El.ComplexDouble(4,9) )
      A.Update( i, iRel+N0, El.ComplexDouble(3,10) )

    # The dense last column
    A.Update( i, width-1, El.ComplexDouble(-10/height) );

  return A

A = ConcatFD2D(n0,n1)
b = El.DistMatrix(El.zTag)
#El.Gaussian( b, n0*n1, 1 )
El.Ones( b, n0*n1, 1 )
if display:
  El.Display( A, "A" )
  El.Display( b, "b" )
if output:
  El.Print( A, "A" )
  El.Print( b, "b" )

ctrl = El.BPCtrl_z()
ctrl.ipmCtrl.mehrotraCtrl.minTol = 1e-5
ctrl.ipmCtrl.mehrotraCtrl.targetTol = 1e-8
ctrl.ipmCtrl.mehrotraCtrl.time = True
ctrl.ipmCtrl.mehrotraCtrl.progress = True
ctrl.ipmCtrl.mehrotraCtrl.solveCtrl.progress = True
startBP = El.mpi.Time()
x = El.BP( A, b, ctrl )
endBP = El.mpi.Time()
if worldRank == 0:
  print('BP time: {} seconds'.format(endBP-startBP))
if display:
  El.Display( x, "x" )
if output:
  El.Print( x, "x" )

xOneNorm = El.EntrywiseNorm( x, 1 )
e = El.DistMatrix(El.zTag)
El.Copy( b, e )
El.Gemv( El.NORMAL, El.ComplexDouble(-1), A, x, El.ComplexDouble(1), e )
if display:
  El.Display( e, "e" )
if output:
  El.Print( e, "e" )
eTwoNorm = El.Nrm2( e )
if worldRank == 0:
  print('|| x ||_1       = {}'.format(xOneNorm))
  print('|| A x - b ||_2 = {}'.format(eTwoNorm))

El.Finalize()
