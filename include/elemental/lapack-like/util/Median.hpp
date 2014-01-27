/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MEDIAN_HPP
#define ELEM_MEDIAN_HPP

namespace elem {

template<typename Real>
inline ValueInt<Real>
Median( const Matrix<Real>& x )
{
    DEBUG_ONLY(CallStackEntry cse("Median"))
    if( IsComplex<Real>::val )
        LogicError("Complex numbers do not have a natural ordering");
    const Int m = x.Height();
    const Int n = x.Width();
    if( m != 1 && n != 1 )
        LogicError("Median is meant for a single vector");

    const Int k = ( n==1 ? m : n );
    const Int stride = ( n==1 ? 1 : x.LDim() );
    const Real* xBuffer = x.LockedBuffer();

    std::vector<ValueInt<Real>> pairs( k );
    for( Int i=0; i<k; ++i )
    {
        pairs[i].value = xBuffer[i*stride];
        pairs[i].index = i;
    }
    std::sort( pairs.begin(), pairs.end(), ValueInt<Real>::Lesser );

    return pairs[k/2];
}

template<typename Real,Dist U,Dist V>
inline ValueInt<Real>
Median( const DistMatrix<Real,U,V>& x )
{
    DEBUG_ONLY(CallStackEntry cse("Median"))
    if( U==STAR && V==STAR )
    {
        return Median( x.LockedMatrix() );
    }
    else
    {
        DistMatrix<Real,STAR,STAR> x_STAR_STAR( x );
        return Median( x_STAR_STAR.LockedMatrix() );
    }
}

} // namespace elem

#endif // ifndef ELEM_MEDIAN_HPP
