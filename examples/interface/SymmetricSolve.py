#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

xSize = ySize = 50

def Laplacian(nx,ny):
  A = El.DistSparseMatrix(El.dTag)
  A.Resize(nx*ny,nx*ny)
  firstLocalRow = A.FirstLocalRow()
  localHeight = A.LocalHeight()
  A.Reserve(5*localHeight)
  hxInvSq = (1.*(nx+1))**2
  hyInvSq = (1.*(ny+1))**2
  for iLoc in xrange(localHeight):
    s = firstLocalRow + iLoc
    x = s % nx
    y = s / nx
    A.QueueUpdate( s, s, 2*(hxInvSq+hyInvSq) )
    if x != 0:    A.QueueUpdate( s, s-1, -hxInvSq )
    if x != nx-1: A.QueueUpdate( s, s+1, -hxInvSq )
    if y != 0:    A.QueueUpdate( s, s-nx, -hyInvSq )
    if y != ny-1: A.QueueUpdate( s, s+nx, -hyInvSq )

  A.MakeConsistent()
  return A

A = Laplacian(xSize,ySize)
x = El.DistMultiVec()
y = El.DistMultiVec()
x.Resize(xSize*ySize,1)
y.Resize(xSize*ySize,1)
xFirst = x.FirstLocalRow()
xLocalHeight = x.LocalHeight()
for iLoc in xrange(xLocalHeight):
  x.SetLocal(iLoc,0,-1.*iLoc)
  y.SetLocal(iLoc,0,-1.*iLoc)

El.Print( A, "Laplacian" )
El.Display( A, "Laplacian" )

El.Print( y, "y" )
El.Display( y, "y" )

El.SymmetricSolve(A,x)

El.Print( x, "x" )
El.Display( x, "x" )

# Require the user to press a button before the figures are closed
El.Finalize()
raw_input('Press Enter to exit')
