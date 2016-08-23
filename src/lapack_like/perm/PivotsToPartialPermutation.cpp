/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

void PivotsToPartialPermutation
( const Matrix<Int>& pivots,
        Matrix<Int>& p,
        Matrix<Int>& pInv, 
  Int offset )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( pivots.Width() != 1 )
          LogicError("pivots must be a column vector");
    )

    const Int b = pivots.Height();
    const Int* pivotsBuf = pivots.LockedBuffer();
    p.Resize( b, 1 );
    pInv.Resize( b, 1 );
    Int* pBuf = p.Buffer();
    Int* pInvBuf = pInv.Buffer();
 
    // Assume that an O(1) number of pivots is supplied and run an algorithm
    // which is quadratic in the number of pivots, but with a low coefficient.
    // A log-linear algorithm with a higher coefficient is also possible.

    for( Int i=0; i<b; ++i ) 
    {
        Int k = pivotsBuf[i] - offset;
        for( Int j=i-1; j>=0; --j )
        {
            const Int relSwap = pivotsBuf[j] - offset;
            if( k == relSwap )
                k = j;
            else if( k == j )
                k = relSwap;
        }
        pBuf[i] = k;
    }

    for( Int i=0; i<b; ++i )
    {
        Int k = i;
        // TODO: Double-check that the upper-bound should change
        for( Int j=0; j<Min(k+1,b); ++j )
        {
            const Int relSwap = pivotsBuf[j] - offset;
            if( k == relSwap )
                k = j; 
            else if( k == j )
                k = relSwap;
        }
        pInvBuf[i] = k;
    }
}

void PivotsToPartialPermutation
( const DistMatrix<Int,STAR,STAR>& pivots, 
        AbstractDistMatrix<Int>& p, 
        AbstractDistMatrix<Int>& pInv,
  Int offset )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( pivots.Width() != 1 )
          LogicError("pivots must be a column vector");
    )

    const Int b = pivots.Height();
    p.SetGrid( pivots.Grid() );
    pInv.SetGrid( pivots.Grid() );
    pInv.Resize( b, 1 );
    p.Resize( b, 1 );

    const Int* pivotsBuf = pivots.LockedBuffer();
 
    // Assume that an O(1) number of pivots is supplied and run an algorithm
    // which is quadratic in the number of pivots, but with a low coefficient.
    // A log-linear algorithm with a higher coefficient is also possible.

    for( Int i=0; i<b; ++i ) 
    {
        Int k = pivotsBuf[i] - offset;
        for( Int j=i-1; j>=0; --j )
        {
            const Int relSwap = pivotsBuf[j] - offset;
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
            const Int relSwap = pivotsBuf[j] - offset;
            if( k == relSwap )
                k = j; 
            else if( k == j )
                k = relSwap;
        }
        pInv.Set( i, 0, k );
    }
}

void PivotsToPartialPermutation
( const AbstractDistMatrix<Int>& pivotsPre,
        AbstractDistMatrix<Int>& p, 
        AbstractDistMatrix<Int>& pInv,
  Int offset )
{
    DEBUG_CSE
    DistMatrixReadProxy<Int,Int,STAR,STAR> pivotsProx( pivotsPre );
    auto& pivots = pivotsProx.GetLocked();
    PivotsToPartialPermutation( pivots, p, pInv, offset );
}

} // namespace El
