#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El
n = 20

# Generate a simple one-dimensional graph
G = El.Graph()
G.Resize(n,n)
for j in xrange(0,n):
  G.Insert(j,(j-1)%n)
  G.Insert(j,j)
  G.Insert(j,(j+1)%n)
G.MakeConsistent()

# Display the graph (with networkx if possible)
El.Display(G)

# Require the user to press a button before the figures are closed
El.Finalize()
raw_input('Press Enter to exit')
