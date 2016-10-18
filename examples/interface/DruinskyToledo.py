#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import math, El
k = 140                 # matrix size
realRes = imagRes = 100 # grid resolution

# Display an instance of the pathological example
A = El.DistMatrix()
El.DruinskyToledo(A,k)
El.Display(A,"Bunch-Kaufman growth matrix")

# Display the spectral portrait
portrait, box = El.SpectralPortrait(A,realRes,imagRes)
El.DisplayPortrait(portrait,box,"spectral portrait of BK growth matrix")

# Make a copy before overwriting with LDL factorization
A_LU = El.DistMatrix()
El.Copy( A, A_LU )

# Display the relevant pieces of pivoted LDL factorization
dSub, p = El.LDL(A,False,El.BUNCH_KAUFMAN_A)
El.MakeTrapezoidal(El.LOWER,A)
El.Display(dSub,"Subdiagonal of D from LDL")
#P = El.DistMatrix(iTag,MC,MR,A.Grid())
# TODO: Construct P from p
#El.Display(P,"P")
El.EntrywiseMap(A,lambda x:math.log10(max(abs(x),1)))
El.Display(A,"Logarithmically-scaled LDL triangular factor")

# Display the relevant pieces of a pivoted LU factorization
p_LU = El.LU(A_LU)
El.Display(p_LU,"LU permutation")
El.EntrywiseMap(A_LU,lambda x:math.log10(max(abs(x),1)))
El.Display(A_LU,"Logarithmically-scaled LU factors")

# Require the user to press a button before the figures are closed
worldSize = El.mpi.WorldSize()
El.Finalize()
if worldSize == 1:
  raw_input('Press Enter to exit')
