/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SOFTTHRESHOLD_HPP
#define ELEM_SOFTTHRESHOLD_HPP

#include ELEM_MAXNORM_INC

namespace elem {

template<typename F>
inline F
SoftThreshold( F alpha, Base<F> tau )
{
    DEBUG_ONLY(
        CallStackEntry cse("SoftThreshold");
        if( tau < 0 )
            LogicError("Negative threshold does not make sense");
    )
    const Base<F> scale = Abs(alpha);
    return ( scale <= tau ? F(0) : alpha-(alpha/scale)*tau );
}

template<typename F>
inline void
SoftThreshold( Matrix<F>& A, Base<F> tau, bool relative=false )
{
    DEBUG_ONLY(CallStackEntry cse("SoftThreshold"))
    if( relative )
        tau *= MaxNorm(A);
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            A.Set( i, j, SoftThreshold(A.Get(i,j),tau) );
}

template<typename F,Dist U,Dist V>
inline void
SoftThreshold( DistMatrix<F,U,V>& A, Base<F> tau, bool relative=false )
{
    DEBUG_ONLY(CallStackEntry cse("SoftThreshold"))
    if( relative )
        tau *= MaxNorm(A);
    SoftThreshold( A.Matrix(), tau, false );
}

} // namespace elem

#endif // ifndef ELEM_SOFTTHRESHOLD_HPP
