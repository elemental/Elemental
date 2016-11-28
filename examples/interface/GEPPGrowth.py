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

# Display an instance of the pathological example
A = El.DistMatrix()
El.GEPPGrowth(A,n)
El.Display(A,"GEPP growth matrix")

# Display the spectral portrait
portrait, box = El.SpectralPortrait(A,realRes,imagRes)
El.DisplayPortrait(portrait,box,"spectral portrait of GEPP growth matrix")

# Display the relevant pieces of pivoted LU factorization
p = El.LU(A)
El.Display(p,"LU permutation")
El.EntrywiseMap(A,lambda x:math.log10(max(abs(x),1)))
El.Display(A,"Logarithmically-scaled LU factors")
El.Display(A[0:n,n-1],"Last column of logarithmic U")

El.Finalize()
