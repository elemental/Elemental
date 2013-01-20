/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CORE_RANDOM_DECL_HPP
#define CORE_RANDOM_DECL_HPP

namespace elem {

const double Pi = 3.141592653589793;

// Generate a sample from a uniform PDF over the (closed) unit ball about the 
// origin of the ring implied by the type T using the most natural metric.
template<typename T> T SampleUnitBall();

} // namespace elem

#endif // ifndef CORE_RANDOM_DECL_HPP
