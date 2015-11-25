# A problematic SOC problem from CVXPY's tests.

import El

orders = El.DistMultiVec(El.iTag)
firstInds = El.DistMultiVec(El.iTag)
orders.Resize(21,1)
firstInds.Resize(21,1)
for c in xrange(7):
  for i in xrange(3):
    orders.Set(3*c+i,0,3)
    firstInds.Set(3*c+i,0,3*c)

h = El.DistMultiVec(El.dTag)
El.Zeros( h, 21, 1 )
h.Set(  3, 0,  1 )
h.Set(  4, 0,  1 )
h.Set(  6, 0,  1 )
h.Set(  7, 0, -1 )
h.Set( 11, 0,  2 )
h.Set( 12, 0,  1 )
h.Set( 13, 0, -1 )

c = El.DistMultiVec(El.dTag)
El.Zeros( c, 8, 1 )
c.Set( 1, 0, 1 )

G = El.DistSparseMatrix(El.dTag)
El.Zeros( G, 21, 8 )
G.Reserve( 28 )
G.QueueUpdate( 0, 3, -1, passive=True )
G.QueueUpdate( 0, 4, -1, passive=True )
G.QueueUpdate( 1, 3, -1, passive=True )
G.QueueUpdate( 1, 4,  1, passive=True )
G.QueueUpdate( 2, 2, -2, passive=True )
G.QueueUpdate( 3, 5, -1, passive=True )
G.QueueUpdate( 4, 5,  1, passive=True )
G.QueueUpdate( 5, 4, -2, passive=True )
G.QueueUpdate( 6, 1, -1, passive=True )
G.QueueUpdate( 7, 1, -1, passive=True )
G.QueueUpdate( 8, 6, -2, passive=True )
G.QueueUpdate( 9, 6, -1, passive=True )
G.QueueUpdate( 9, 7, -1, passive=True )
G.QueueUpdate( 10, 6, -1, passive=True )
G.QueueUpdate( 10, 7, 1, passive=True )
G.QueueUpdate( 12, 0, -1, passive=True )
G.QueueUpdate( 13, 0, -1, passive=True )
G.QueueUpdate( 14, 5, -2, passive=True )
G.QueueUpdate( 15, 2, -1, passive=True )
G.QueueUpdate( 15, 5, -1, passive=True )
G.QueueUpdate( 16, 2,  1, passive=True )
G.QueueUpdate( 16, 5, -1, passive=True )
G.QueueUpdate( 17, 7, -2, passive=True )
G.QueueUpdate( 18, 0, -1, passive=True )
G.QueueUpdate( 18, 1, -1, passive=True )
G.QueueUpdate( 19, 0, -1, passive=True )
G.QueueUpdate( 19, 1, 1, passive=True )
G.QueueUpdate( 20, 3, -2, passive=True )
G.ProcessQueues()

b = El.DistMultiVec(El.dTag)
b.Resize( 1, 1 )
b.Set( 0, 0, 7.45 )

A = El.DistSparseMatrix(El.dTag)
El.Zeros( A, 1, 8 )
A.Reserve( 1 )
A.QueueUpdate( 0, 0, 1, passive=True )
A.ProcessQueues()

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
