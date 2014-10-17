#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import math, El
n = 200                 # matrix size
realRes = imagRes = 100 # grid resolution

# Display an instance of the Fourier matrix
A = El.DistMatrix(El.zTag)
El.Fourier(A,n)
El.Display(A,"Fourier matrix")

# Display a submatrix and subvector of the Fourier matrix
El.Display(A[(n/4):(3*n/4),(n/4):(3*n/4)],"Middle submatrix")
El.Display(A[1,0:n],"Second row")

# Display the spectral portrait
portrait = El.SpectralPortrait(A,realRes,imagRes)
El.EntrywiseMap(portrait,math.log10)
El.Display(portrait,"spectral portrait of Fourier matrix")

# Require the user to press a button before the figures are closed
El.Finalize()
raw_input('Press Enter to exit')
