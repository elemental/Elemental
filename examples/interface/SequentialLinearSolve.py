#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

n0 = n1 = 200
display = False

def Square(xSize,ySize):
  A = El.SparseMatrix(El.dTag)
  n = xSize*ySize
  A.Resize(n,n)
  A.Reserve(5*n)
  hxInvSq = (1.*(xSize+1))**2
  hyInvSq = (1.*(ySize+1))**2
  for s in xrange(n):
    x = s % xSize
    y = s / xSize
    A.QueueUpdate( s, s, 8*(hxInvSq+hyInvSq) )
    if x != 0:       A.QueueUpdate( s, s-1,     -1*hxInvSq )
    if x != xSize-1: A.QueueUpdate( s, s+1,     -2*hxInvSq )
    if y != 0:       A.QueueUpdate( s, s-xSize, -4*hyInvSq )
    if y != ySize-1: A.QueueUpdate( s, s+xSize, -3*hyInvSq )

  A.ProcessQueues()
  return A

def Rectang(height,width):
  A = El.SparseMatrix()
  A.Resize(height,width)
  A.Reserve(5*height)
  for s in xrange(height):
    if s < width:
      A.QueueUpdate( s, s,        11 )
    if s >= 1 and s-1 < width:
      A.QueueUpdate( s, s-1,      -1 )
    if s+1 < width:
      A.QueueUpdate( s, s+1,       2 )
    if s >= height and s-height < width:
      A.QueueUpdate( s, s-height, -3 )
    if s+height < width:
      A.QueueUpdate( s, s+height,  4 )
    # The dense last column
    A.QueueUpdate( s, width-1, -5/height );

  A.ProcessQueues()
  return A

#A = Square(n0,n1)
A = Rectang(n0*n1,n0*n1)
x = El.Matrix()
y = El.Matrix()
El.Uniform( x, n0*n1, 1 )
El.Copy( x, y )

yNrm = El.Nrm2(y)
print('|| y ||_2 = {}'.format(yNrm))

if display:
  El.Display( y, "y" )

El.LinearSolve(A,x)
if display:
  El.Display( x, "x" )

xNrm = El.Nrm2(x)
print('|| x ||_2 = {}'.format(xNrm))

El.Multiply(El.NORMAL,-1.,A,x,1.,y)
if display:
  El.Display( y, "A x - y" )
eNrm = El.Nrm2(y)
print('|| A x - y ||_2 / || y ||_2 = {}'.format(eNrm/yNrm))

El.Finalize()
