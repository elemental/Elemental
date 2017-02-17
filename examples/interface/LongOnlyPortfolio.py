#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El, math

n = 1000
r = 30
gamma = 0.1
display = False
worldRank = El.mpi.WorldRank()
worldSize = El.mpi.WorldSize()

# Create a positive diagonal for the covariance
def CreateDiag(height):
  d = El.DistMultiVec()
  El.Uniform( d, height, 1, 0.45, 0.45 )
  return d

# Create the factor model (for now, make it dense)
def CreateFactor(height,width):
  F = El.DistSparseMatrix()
  El.Zeros( F, height, width )
  localHeight = F.LocalHeight()
  F.Reserve(localHeight*width)
  for iLoc in xrange(localHeight):
    i = F.GlobalRow(iLoc)
    for j in xrange(width):
      F.QueueLocalUpdate( iLoc, j, math.log(i+j+1.) )

  F.ProcessQueues()
  # NOTE: Without this rescaling, the problem is much harder, as the variables
  #       s becomes quite large 
  #       (|| F^T x ||_2 <= u, u^2 <= t, would imply u and t are large).
  FFrob = El.FrobeniusNorm( F )
  El.Scale( 1./FFrob, F )
  return F

# Create the vector of expected returns
def CreateExpected(height):
  c = El.DistMultiVec()
  #Zeros( c, height, 1 )
  #localHeight = c.LocalHeight()
  #for iLoc in xrange(localHeight):
  #  i = c.GlobalRow(iLoc)
  #  c.SetLocal(iLoc,0,1.+1./i)
  El.Gaussian( c, height, 1 )
  El.EntrywiseMap( c, lambda alpha : abs(alpha) )

  return c

d = CreateDiag(n)
F = CreateFactor(n,r)
c = CreateExpected(n)
if display:
  El.Display( d, "d" )
  El.Display( F, "F" )
  El.Display( c, "c" )

ctrl = El.SOCPAffineCtrl_d()
ctrl.mehrotraCtrl.progress = True
ctrl.mehrotraCtrl.time = False
ctrl.mehrotraCtrl.solveCtrl.progress = False
ctrl.mehrotraCtrl.solveCtrl.time = False

# Solve *with* resolving the regularization
ctrl.mehrotraCtrl.resolveReg = True
startLOP = El.mpi.Time()
x = El.LongOnlyPortfolio(d,F,c,gamma,ctrl)
endLOP = El.mpi.Time()
if worldRank == 0:
  print('LOP time (resolve reg. w/ equil): {} seconds'.format(endLOP-startLOP))
if display:
  El.Display( x, "x" )

# Solve without resolving the regularization
ctrl.mehrotraCtrl.resolveReg = False
startLOP = El.mpi.Time()
x = El.LongOnlyPortfolio(d,F,c,gamma,ctrl)
endLOP = El.mpi.Time()
if worldRank == 0:
  print('LOP time (no resolve reg. w/ equil): {} seconds'.format( \
    endLOP-startLOP))
if display:
  El.Display( x, "x" )

# Compute the risk-adjusted return
# ================================
e = El.DistMultiVec()
f = El.DistMultiVec()
El.Copy( x, e )
El.DiagonalScale( El.LEFT, El.NORMAL, d, e )
El.Zeros( f, r, 1 )
El.Multiply( El.TRANSPOSE, 1., F, x, 0., f )
El.Multiply( El.NORMAL, 1., F, f, 1., e )
rar = El.Dot(c,x) - gamma*El.Dot(x,e)
if worldRank == 0:
  print('c^T x - gamma x^T (D + F F^T) x = {}'.format(rar))

xOneNorm = El.EntrywiseNorm( x, 1 )
xTwoNorm = El.Nrm2( x )
if worldRank == 0:
  print('|| x ||_1 = {}'.format(xOneNorm))
  print('|| x ||_2 = {}'.format(xTwoNorm))

El.Finalize()
