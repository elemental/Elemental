/*
   Copyright (c) 2009-2014, Jack Poulson
                      2013, Michael C. Grant
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

%module elem_control

%include "common.swg"
%import "elem.i"

// Control theory
// ==============
 
%include "elemental/control/Lyapunov.hpp"
%include "elemental/control/Ricatti.hpp"
%include "elemental/control/Sylvester.hpp"

namespace elem {
OVERLOAD0(Lyapunov)
OVERLOAD0(Ricatti)
OVERLOAD0(Sylvester)
}
