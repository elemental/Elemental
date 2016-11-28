#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El, time, math

n=100
realRes = imagRes = 50 # grid resolution

A = El.DistMatrix()
El.JordanCholesky( A, n )
El.Display( A, "A" )

El.Cholesky( El.UPPER, A )
El.MakeTrapezoidal( El.UPPER, A )
El.Display( A, "U" )

portrait = El.SpectralWindow(A,0,6,4,realRes,imagRes)
box = El.SpectralBox_d()
box.center = El.TagToType(El.Complexify(A.tag))(0)
box.realWidth = 6
box.imagWidth = 4
El.DisplayPortrait(portrait,box,"spectral portrait of 2*J_{1/2}(n)")

El.Finalize()
