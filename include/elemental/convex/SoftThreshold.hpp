/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CONVEX_SOFTTHRESHOLD_HPP
#define ELEM_CONVEX_SOFTTHRESHOLD_HPP

#include "elemental/lapack-like/Norm/Max.hpp"

namespace elem {

template<typename F>
inline F
SoftThreshold( F alpha, BASE(F) tau )
{
#ifndef RELEASE
    CallStackEntry cse("SoftThreshold");
    if( tau < 0 )
        LogicError("Negative threshold does not make sense");
#endif
    typedef Base<F> R;
    const R scale = Abs(alpha);
    return ( scale <= tau ? F(0) : alpha-(alpha/scale)*tau );
}

template<typename F>
inline void
SoftThreshold( Matrix<F>& A, BASE(F) tau, bool relative=false )
{
#ifndef RELEASE
    CallStackEntry cse("SoftThreshold");
#endif
    typedef Base<F> Real;
    if( relative )
    {
        const Real maxNorm = MaxNorm( A );
        tau *= maxNorm;
    }
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            A.Set( i, j, SoftThreshold(A.Get(i,j),tau) );
}

template<typename F,Distribution U,Distribution V>
inline void
SoftThreshold( DistMatrix<F,U,V>& A, BASE(F) tau, bool relative=false )
{
#ifndef RELEASE
    CallStackEntry cse("SoftThreshold");
#endif
    typedef Base<F> Real;
    if( relative )
    {
        const Real maxNorm = MaxNorm( A );
        tau *= maxNorm;
    }
    SoftThreshold( A.Matrix(), tau, false );
}

} // namespace elem

#endif // ifndef ELEM_CONVEX_SOFTTHRESHOLD_HPP
