/*
   Copyright (c) 2009-2013, Jack Poulson
                      2013, Michael C. Grant
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

%module elem_convex

%include "common.swg"
%import "elem.i"

/*
 * CONVEX OPTIMIZATION
 */
 
%include "elemental/convex/LogBarrier.hpp"
%include "elemental/convex/LogDetDivergence.hpp"
%include "elemental/convex/SVT.hpp"
%include "elemental/convex/SoftThreshold.hpp"
%include "elemental/convex/UnitaryCoherence.hpp" 

namespace elem {
OVERLOAD0(LogBarrier)
OVERLOAD0(LogDetDivergence)
OVERLOAD0(SVT)
OVERLOAD0(SoftThreshold)
OVERLOAD0(UnitaryCoherence)
}
