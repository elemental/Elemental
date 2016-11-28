#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El
import time

n=30
A = El.DistSparseMatrix()
El.DynamicRegCounter( A, n )
El.Display( A, "A" )

El.Finalize()
