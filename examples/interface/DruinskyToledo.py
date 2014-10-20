#
#  Copyright (c) 2009-2014, Jack Poulson
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
El.DruinskyToledo(A,n)
El.Display(A,"Bunch-Kaufman growth matrix")

# Display the spectral portrait
portrait = El.SpectralPortrait(A,realRes,imagRes)
El.EntrywiseMap(portrait,math.log10)
El.Display(portrait,"spectral portrait of Bunch-Kaufman growth matrix")

# Display the relevant pieces of pivoted LU factorization
dSub, p = El.LDL(A,False,El.BUNCH_KAUFMAN_A)
El.MakeTriangular(El.LOWER,A)
El.Display(dSub,"Subdiagonal of D from LDL")
El.Display(p,"LDL permutation")
El.Display(A[n-2:n,n-2:n],"Bottom-right 2x2 of L")
El.EntrywiseMap(A,lambda x:math.log10(max(abs(x),1)))
El.Display(A,"Logarithmically-scaled triangular factor")

# Require the user to press a button before the figures are closed
El.Finalize()
raw_input('Press Enter to exit')
