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
A = El.DistMatrix(El.zTag)
El.FoxLi(A,n)
El.Display(A,"Fox-Li/Landau matrix")
portrait = El.SpectralPortrait(A,realRes,imagRes)
El.EntrywiseMap(portrait,math.log)
El.Display(portrait,"spectral portrait of Fox-Li/Landau matrix")
El.Finalize()

# Require the user to press a button before the figures are closed
raw_input('Press Enter to exit')
