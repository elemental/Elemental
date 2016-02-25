#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import math, El
n = 100                 # matrix size
realRes = imagRes = 100 # grid resolution

# Display an instance of the Fox-Li/Landau matrix
A = El.DistMatrix(El.zTag)
El.FoxLi(A,n)
El.Display(A,"Fox-Li/Landau matrix")

# Display its spectral portrait
portrait, box = El.SpectralPortrait(A,realRes,imagRes)
El.DisplayPortrait(portrait,box,"spectral portrait of Fox-Li/Landau matrix")

# Display its singular values
s = El.SingularValues(A)
El.EntrywiseMap(s,math.log10)
El.Display(s,"log10(svd(A))")

# Require the user to press a button before the figures are closed
worldSize = El.mpi.WorldSize()
El.Finalize()
if worldSize == 1:
  raw_input('Press Enter to exit')
