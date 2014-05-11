/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_RANDOM_DECL_HPP
#define ELEM_RANDOM_DECL_HPP

namespace elem {

const double Pi = 3.141592653589793;

bool BooleanCoinFlip();
Int CoinFlip();

template<typename T>
T UnitCell();

template<typename T=double>
T SampleUniform( T a=0, T b=UnitCell<T>() );

template<>
Int SampleUniform<Int>( Int a, Int b );

// The complex extension of the normal distribution can actually be quite
// technical, and so we will use the simplest case, where both the real and
// imaginary components are independently drawn with the same standard 
// deviation, but different means.
template<typename T=double>
T SampleNormal( T mean=0, Base<T> stddev=1 );

// Generate a sample from a uniform PDF over the (closed) unit ball about the 
// origin of the ring implied by the type T using the most natural metric.
template<typename T> 
T SampleBall( T center=0, Base<T> radius=1 );

} // namespace elem

#endif // ifndef ELEM_RANDOM_DECL_HPP
