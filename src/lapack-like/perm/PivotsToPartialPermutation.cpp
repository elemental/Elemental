/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

void PivotsToPartialPermutation
( const Matrix<Int>& pivots, Matrix<Int>& perm, Matrix<Int>& invPerm, 
  Int offset )
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
void PivotsToPartialPermutation
( const DistMatrix<Int,U,    STAR>& pivots, 
        DistMatrix<Int,UPerm,STAR>& perm, 
        DistMatrix<Int,UPerm,STAR>& invPerm, Int offset )
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

#define PROTO_DIST_INTERNAL(U,UPERM) \
  template void PivotsToPartialPermutation \
  ( const DistMatrix<Int,U,    STAR>& pivots, \
          DistMatrix<Int,UPERM,STAR>& perm, \
          DistMatrix<Int,UPERM,STAR>& invPerm, Int offset );

#define PROTO_DIST(U) \
  PROTO_DIST_INTERNAL(U,MC  ) \
  PROTO_DIST_INTERNAL(U,MD  ) \
  PROTO_DIST_INTERNAL(U,MR  ) \
  PROTO_DIST_INTERNAL(U,STAR) \
  PROTO_DIST_INTERNAL(U,VC  ) \
  PROTO_DIST_INTERNAL(U,VR  ) \

PROTO_DIST(MC  )
PROTO_DIST(MD  )
PROTO_DIST(MR  )
PROTO_DIST(STAR)
PROTO_DIST(VC  )
PROTO_DIST(VR  )

} // namespace El
