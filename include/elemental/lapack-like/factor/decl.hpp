/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_FACTOR_DECL_HPP
#define ELEM_FACTOR_DECL_HPP

namespace elem {

namespace ldl_pivot_type_wrapper {
enum LDLPivotType
{
    BUNCH_KAUFMAN_A=0,
    BUNCH_KAUFMAN_C=1,
    BUNCH_KAUFMAN_D=2,
    BUNCH_KAUFMAN_BOUNDED=3,
    BUNCH_PARLETT=4
};
}
using namespace ldl_pivot_type_wrapper;

struct LDLPivot
{
    Int nb;
    Int from[2];
};

} // namespace elem

#endif // ifndef ELEM_FACTOR_DECL_HPP
