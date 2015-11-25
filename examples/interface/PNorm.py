# A small SOCP formulation of minimizing the p-norm of a vector subject to
# constraints (in this trivial case, maximize the two-norm subject to the 
# entire vector being fixed). This example originally caused the ninth SOCP
# IPM KKT system to be excessively ill-conditioned and was reported by 
# Steven Diamond via CVXPY. The original script was:
#
#   from cvxpy import *
#   x = Variable(3)
#   prob = Problem(Maximize(pnorm(x,.5)),[x==[1.1,2,.1]])
#   prob.solve(verbose=True, solver=ELEMENTAL)
# 
# which should have a solution of roughly 7.72423. The problem data used
# in the following script can be found by then running
#
#   data = prob.get_problem_data(ELEMENTAL)
#   print( data['A'] )
#   print( data['G'] )
#   print( data['b'] )
#   print( data['c'] )
#   print( data['h'] )
#
# which has three SOCs of size three.

import El

orders = El.DistMultiVec(El.iTag)
firstInds = El.DistMultiVec(El.iTag)
orders.Resize(9,1)
firstInds.Resize(9,1)
for c in xrange(0,3):
  for i in xrange(0,3):
    orders.Set(3*c+i,0,3)
    firstInds.Set(3*c+i,0,3*c)

A = El.DistSparseMatrix(El.dTag)
A.Resize(4,7)
A.Reserve(7)
A.QueueUpdate( 1, 0,  1, passive=True )
A.QueueUpdate( 2, 1,  1, passive=True )
A.QueueUpdate( 3, 2,  1, passive=True )
A.QueueUpdate( 0, 3, -1, passive=True )
A.QueueUpdate( 0, 4,  1, passive=True )
A.QueueUpdate( 0, 5,  1, passive=True )
A.QueueUpdate( 0, 6,  1, passive=True )
A.ProcessLocalQueues()

G = El.DistSparseMatrix(El.dTag)
G.Resize(9,7)
G.Reserve(15)
G.QueueUpdate( 0, 0, -1, passive=True )
G.QueueUpdate( 1, 0, -1, passive=True )
G.QueueUpdate( 3, 1, -1, passive=True )
G.QueueUpdate( 4, 1, -1, passive=True )
G.QueueUpdate( 6, 2, -1, passive=True )
G.QueueUpdate( 7, 2, -1, passive=True )
G.QueueUpdate( 0, 3, -1, passive=True )
G.QueueUpdate( 1, 3,  1, passive=True )
G.QueueUpdate( 3, 3, -1, passive=True )
G.QueueUpdate( 4, 3,  1, passive=True )
G.QueueUpdate( 6, 3, -1, passive=True )
G.QueueUpdate( 7, 3,  1, passive=True )
G.QueueUpdate( 2, 4, -2, passive=True )
G.QueueUpdate( 5, 5, -2, passive=True )
G.QueueUpdate( 8, 6, -2, passive=True )
G.ProcessLocalQueues()

b = El.DistMultiVec(El.dTag)
b.Resize(4,1)
b.Set(0,0,-0.0)
b.Set(1,0, 1.1)
b.Set(2,0, 2.0)
b.Set(3,0, 0.1)

c = El.DistMultiVec(El.dTag)
c.Resize(7,1)
c.Set(0,0, 0)
c.Set(1,0, 0)
c.Set(2,0, 0)
c.Set(3,0,-1)
c.Set(4,0, 0)
c.Set(5,0, 0)
c.Set(6,0, 0)

h = El.DistMultiVec(El.dTag)
El.Zeros( h, 9, 1 )
  
x = El.DistMultiVec()
y = El.DistMultiVec()
z = El.DistMultiVec()
s = El.DistMultiVec()
ctrl = El.SOCPAffineCtrl_d()
ctrl.mehrotraCtrl.solveCtrl.progress = True
ctrl.mehrotraCtrl.progress = True
ctrl.mehrotraCtrl.outerEquil = True
ctrl.mehrotraCtrl.time = True
El.SOCPAffine(A,G,b,c,h,orders,firstInds,x,y,z,s,ctrl)

El.Finalize()
