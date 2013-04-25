/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CONVEX_SOFTTHRESHOLD_HPP
#define CONVEX_SOFTTHRESHOLD_HPP

namespace elem {

template<typename F>
inline F
SoftThreshold( F alpha, BASE(F) tau )
{
#ifndef RELEASE
    CallStackEntry entry("SoftThreshold");
    if( tau < 0 )
        throw std::logic_error("Negative threshold does not make sense");
#endif
    typedef BASE(F) R;
    const R scale = Abs(alpha);
    return ( scale <= tau ? F(0) : alpha-(alpha/scale)*tau );
}

template<typename F>
inline void
SoftThreshold( Matrix<F>& A, BASE(F) tau )
{
#ifndef RELEASE
    CallStackEntry entry("SoftThreshold");
#endif
    const int height = A.Height();
    const int width = A.Width();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            A.Set( i, j, SoftThreshold(A.Get(i,j),tau) );
}

template<typename F,Distribution U,Distribution V>
inline void
SoftThreshold( DistMatrix<F,U,V>& A, BASE(F) tau )
{
#ifndef RELEASE
    CallStackEntry entry("SoftThreshold");
#endif
    SoftThreshold( A.Matrix(), tau );
}

} // namespace elem

#endif // ifndef CONVEX_SOFTTHRESHOLD_HPP
