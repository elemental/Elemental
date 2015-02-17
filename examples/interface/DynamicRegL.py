#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import math, El
n = 50                  # matrix size
realRes = imagRes = 100 # grid resolution

# Display an instance of the dynamic regularization problematic L
L = El.DistMatrix(El.dTag)
El.DynamicRegL(L,n)
El.Display(L,"Dynamic regularization L")

# Display its spectral portrait
portrait, box = El.SpectralPortrait(L,realRes,imagRes)
El.DisplayPortrait(portrait,box,"spectral portrait of L")

# Display its singular values
s = El.SVD(L)
El.EntrywiseMap(s,math.log10)
El.Display(s,"log10(svd(L))")

# Require the user to press a button before the figures are closed
commSize = El.mpi.Size( El.mpi.COMM_WORLD() )
El.Finalize()
if commSize == 1:
  raw_input('Press Enter to exit')
