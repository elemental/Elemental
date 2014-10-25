#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

n0 = n1 = 50

def Laplacian(xSize,ySize):
  A = El.DistSparseMatrix(El.dTag)
  A.Resize(xSize*ySize,xSize*ySize)
  firstLocalRow = A.FirstLocalRow()
  localHeight = A.LocalHeight()
  A.Reserve(5*localHeight)
  hxInvSq = (1.*(xSize+1))**2
  hyInvSq = (1.*(ySize+1))**2
  for iLoc in xrange(localHeight):
    s = firstLocalRow + iLoc
    x = s % xSize
    y = s / xSize
    A.QueueUpdate( s, s, 2*(hxInvSq+hyInvSq) )
    if x != 0:       A.QueueUpdate( s, s-1,     -hxInvSq )
    if x != xSize-1: A.QueueUpdate( s, s+1,     -hxInvSq )
    if y != 0:       A.QueueUpdate( s, s-xSize, -hyInvSq )
    if y != ySize-1: A.QueueUpdate( s, s+xSize, -hyInvSq )

  A.MakeConsistent()
  return A

A = Laplacian(n0,n1)
x = El.DistMultiVec()
y = El.DistMultiVec()
El.Uniform( x, n0*n1, 1 )
El.Copy( x, y )

El.Display( A, "Laplacian" )
El.Display( y, "y" )

El.SymmetricSolveSparse(A,x)
El.Display( x, "x" )

# Require the user to press a button before the figures are closed
El.Finalize()
raw_input('Press Enter to exit')
