#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El
n = 200                 # matrix size
realRes = imagRes = 100 # grid resolution

# Display an instance of the Fourier matrix
A = El.DistMatrix(El.zTag)
El.Fourier(A,n)
El.Display(A,"Fourier matrix")

# Display a submatrix and subvector of the Fourier matrix
El.Display(A[(n/4):(3*n/4),(n/4):(3*n/4)],"Middle submatrix")
El.Display(A[1,:],"Second row")

# Display the spectral portrait
portrait, box = El.SpectralPortrait(A,realRes,imagRes)
El.DisplayPortrait(portrait,box,"spectral portrait of Fourier matrix")

El.Finalize()
