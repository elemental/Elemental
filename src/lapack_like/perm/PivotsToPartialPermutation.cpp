/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

void PivotsToPartialPermutation
( const Matrix<Int>& pivots, Matrix<Int>& p, Matrix<Int>& pInv, 
  Int offset )
{
    DEBUG_ONLY(
        CallStackEntry cse("PivotsToPartialPermutation");
        if( pivots.Width() != 1 )
            LogicError("pivots must be a column vector");
    )

    const Int b = pivots.Height();
    p.Resize( b, 1 );
    pInv.Resize( b, 1 );
 
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
        p.Set( i, 0, k );
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
        pInv.Set( i, 0, k );
    }
}

void PivotsToPartialPermutation
( const AbstractDistMatrix<Int>& pivots, 
        AbstractDistMatrix<Int>& p, 
        AbstractDistMatrix<Int>& pInv, Int offset )
{
    DEBUG_ONLY(
        CallStackEntry cse("PivotsToPartialPermutation");
        if( pivots.Width() != 1 )
            LogicError("pivots must be a column vector");
    )

    const Int b = pivots.Height();
    p.SetGrid( pivots.Grid() );
    pInv.SetGrid( pivots.Grid() );
    pInv.Resize( b, 1 );
    p.Resize( b, 1 );
 
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
        p.Set( i, 0, k );
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
        pInv.Set( i, 0, k );
    }
}

} // namespace El
