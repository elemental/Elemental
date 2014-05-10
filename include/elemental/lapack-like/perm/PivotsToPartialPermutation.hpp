/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_PIVOTSTOPARTIALPERMUTATION_HPP
#define ELEM_LAPACK_PIVOTSTOPARTIALPERMUTATION_HPP

namespace elem {

inline void
PivotsToPartialPermutation
( const Matrix<Int>& pivots, Matrix<Int>& perm, Matrix<Int>& invPerm, 
  Int offset=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("PivotsToPartialPermutation");
        if( pivots.Width() != 1 )
            LogicError("pivots must be a column vector");
    )

    const Int b = pivots.Height();
    perm.Resize( b, 1 );
    invPerm.Resize( b, 1 );
 
    // Assume that an O(1) number of pivots is supplied and run an algorithm
    // which is quadratic in the number of pivots, but with a low coefficient.
    // A log-linear algorithm with a higher coefficient is also possible.

    for( Int i=0; i<b; ++i ) 
    {
        Int k = pivots.Get(i,0) - offset;
        for( Int j=i-1; j>=0; --j )
        {
            const Int relSwap = pivots.Get(j,0)-offset;
            if( k == relSwap )
                k = j;
            else if( k == j )
                k = relSwap;
        }
        perm.Set( i, 0, k );
    }

    for( Int i=0; i<b; ++i )
    {
        Int k = i;
        for( Int j=0; j<Min(k+1,b); ++j )
        {
            const Int relSwap = pivots.Get(j,0)-offset;
            if( k == relSwap )
                k = j; 
            else if( k == j )
                k = relSwap;
        }
        invPerm.Set( i, 0, k );
    }
}

template<Dist U,Dist UPerm>
inline void
PivotsToPartialPermutation
( const DistMatrix<Int,U,    STAR>& pivots, 
        DistMatrix<Int,UPerm,STAR>& perm, 
        DistMatrix<Int,UPerm,STAR>& invPerm, Int offset=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("PivotsToPartialPermutation");
        if( pivots.Width() != 1 )
            LogicError("pivots must be a column vector");
    )

    const Int b = pivots.Height();
    perm.SetGrid( pivots.Grid() );
    invPerm.SetGrid( pivots.Grid() );
    invPerm.AlignWith( perm );
    invPerm.Resize( b, 1 );
    perm.Resize( b, 1 );
 
    // Assume that an O(1) number of pivots is supplied and run an algorithm
    // which is quadratic in the number of pivots, but with a low coefficient.
    // A log-linear algorithm with a higher coefficient is also possible.

    for( Int i=0; i<b; ++i ) 
    {
        Int k = pivots.Get(i,0) - offset;
        for( Int j=i-1; j>=0; --j )
        {
            const Int relSwap = pivots.Get(j,0)-offset;
            if( k == relSwap )
                k = j;
            else if( k == j )
                k = relSwap;
        }
        perm.Set( i, 0, k );
    }

    // Construct the image using a similar algorithm.
    for( Int i=0; i<b; ++i )
    {
        Int k = i;
        for( Int j=0; j<Min(k+1,b); ++j )
        {
            const Int relSwap = pivots.Get(j,0)-offset;
            if( k == relSwap )
                k = j; 
            else if( k == j )
                k = relSwap;
        }
        invPerm.Set( i, 0, k );
    }
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_PIVOTSTOPARTIALPERMUTATION_HPP
