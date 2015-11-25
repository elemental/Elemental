# A small SOCP formulation of minimizing the inverse of the top-left entry
# of the positive 2x2 matrix [1,2; 3,4]. This example originally caused the 
# 30th SOCP IPM iteration to fail within the symmetric tridiagonal eigensolver
# and was submitted by Steven Diamond. The original script was:
#
#   from cvxpy import *
#   x = Variable(2,2)
#   expr = inv_pos(x)
#   prob = Problem(Minimize(expr[0,0]), [x == [[1,2],[3,4]] ])

# which should have a solution of 1. The problem data used
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
orders.Resize(12,1)
firstInds.Resize(12,1)
for c in xrange(0,4):
  for i in xrange(0,3):
    orders.Set(3*c+i,0,3)
    firstInds.Set(3*c+i,0,3*c)

A = El.DistSparseMatrix(El.dTag)
A.Resize(4,8)
A.Reserve(4)
A.QueueUpdate( 0, 0, 1, passive=True )
A.QueueUpdate( 1, 1, 1, passive=True )
A.QueueUpdate( 2, 2, 1, passive=True )
A.QueueUpdate( 3, 3, 1, passive=True )
A.ProcessLocalQueues()

G = El.DistSparseMatrix(El.dTag)
G.Resize(12,8)
G.Reserve(16)
G.QueueUpdate( 0, 0, -1, passive=True )
G.QueueUpdate( 1, 0, -1, passive=True )
G.QueueUpdate( 3, 1, -1, passive=True )
G.QueueUpdate( 4, 1, -1, passive=True )
G.QueueUpdate( 6, 2, -1, passive=True )
G.QueueUpdate( 7, 2, -1, passive=True )
G.QueueUpdate( 9, 3, -1, passive=True )
G.QueueUpdate( 10, 3, -1, passive=True )
G.QueueUpdate( 0, 4, -1, passive=True )
G.QueueUpdate( 1, 4,  1, passive=True )
G.QueueUpdate( 3, 5, -1, passive=True )
G.QueueUpdate( 4, 5,  1, passive=True )
G.QueueUpdate( 6, 6, -1, passive=True )
G.QueueUpdate( 7, 6,  1, passive=True )
G.QueueUpdate( 9, 7, -1, passive=True )
G.QueueUpdate( 10, 7, 1, passive=True )
G.ProcessLocalQueues()

b = El.DistMultiVec(El.dTag)
b.Resize(4,1)
b.Set(0,0, 1.)
b.Set(1,0, 2.)
b.Set(2,0, 3.)
b.Set(3,0, 4.)

c = El.DistMultiVec(El.dTag)
c.Resize(8,1)
c.Set(0,0, 0)
c.Set(1,0, 0)
c.Set(2,0, 0)
c.Set(3,0, 0)
c.Set(4,0, 1)
c.Set(5,0, 0)
c.Set(6,0, 0)
c.Set(7,0, 0)

h = El.DistMultiVec(El.dTag)
El.Zeros( h, 12, 1 )
h.Set(2,0,2.)
h.Set(5,0,2.)
h.Set(8,0,2.)
h.Set(11,0,2.)

output = True
display = False
if output:
  El.Print( A, "A" ) 
  El.Print( G, "G" )
  El.Print( b, "b" )
  El.Print( c, "c" )
  El.Print( h, "h" )

if display:
  El.Display( A, "A" ) 
  El.Display( G, "G" )
  El.Display( b, "b" )
  El.Display( c, "c" )
  El.Display( h, "h" )
 
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
