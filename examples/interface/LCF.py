#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

def LCF(lcf):
  n = len(lcf)
  G = El.Graph()
  G.Resize(n,n)

  G.Reserve( 4*n )
  for s in xrange(n):
    # Build the links to this node in the Hamiltonian cycle
    tL = (s-1) % n
    tR = (s+1) % n
    G.QueueConnection(s,tL)
    G.QueueConnection(s,tR)
    # Build the LCF links
    t = (s + lcf[s]) % n
    G.QueueConnection(s,t)
    G.QueueConnection(t,s)

  G.ProcessQueues()
  return G

levi = [-13,-9,7,-7,9,13]*5
dodec = [10,7,4,-4,-7,10,-4,7,-7,4]*2
truncOct = [3,-7,7,-3]*6

if El.mpi.WorldRank() == 0:
  El.Display( LCF(levi), "Levi graph" )
  El.Display( LCF(dodec), "Dodecahedral graph" )
  El.Display( LCF(truncOct), "Trunacted octahedral graph" )

# Require the user to press a button before the figures are closed
worldSize = El.mpi.WorldSize()
El.Finalize()
if worldSize == 1:
  raw_input('Press Enter to exit')
