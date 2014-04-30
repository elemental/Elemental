/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_FUNCS_DECL_HPP
#define ELEM_FUNCS_DECL_HPP

#include ELEM_DECOMP_DECL_INC

namespace elem {

namespace SignScalingNS {
enum SignScaling {
    SIGN_SCALE_NONE,
    SIGN_SCALE_DET,
    SIGN_SCALE_FROB
};
} 
using namespace SignScalingNS;

} // namespace elem

#endif // ifndef ELEM_FUNCS_DECL_HPP
