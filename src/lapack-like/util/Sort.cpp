/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

// Sort each column of the real matrix X

template<typename Real>
void Sort( Matrix<Real>& X, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("Sort"))
    if( IsComplex<Real>::val )
        LogicError("Complex numbers do not have a natural ordering");
    if( sort == UNSORTED )
        return;
    const Int m = X.Height();
    const Int n = X.Width();
    for( Int j=0; j<n; ++j )
    {
        Real* XCol = X.Buffer(0,j);
        if( sort == ASCENDING )
            std::sort( XCol, XCol+m );
        else
            std::sort( XCol, XCol+m, std::greater<Real>() );
    }
}

template<typename Real,Dist U,Dist V>
void Sort( DistMatrix<Real,U,V>& X, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("Sort"))
    if( sort == UNSORTED )
        return;

    if( (U==STAR && V==STAR) || (U==CIRC && V==CIRC) )
    {
        if( X.Participating() )
            Sort( X.Matrix(), sort );
    }
    else
    {
        // Get a copy on a single process, sort, and then redistribute
        DistMatrix<Real,CIRC,CIRC> X_CIRC_CIRC( X );
        if( X_CIRC_CIRC.Participating() )
            Sort( X_CIRC_CIRC.Matrix(), sort );

        // Refill the distributed X with the sorted values
        X = X_CIRC_CIRC;
    }
}

// Tagged sort

template<typename Real>
std::vector<ValueInt<Real>> TaggedSort
( const Matrix<Real>& x, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("TaggedSort"))
    if( IsComplex<Real>::val )
        LogicError("Complex numbers do not have a natural ordering");
    const Int m = x.Height();
    const Int n = x.Width();
    if( m != 1 && n != 1 )
        LogicError("TaggedSort is meant for a single vector");

    const Int k = ( n==1 ? m : n );
    const Int stride = ( n==1 ? 1 : x.LDim() );
    const Real* xBuffer = x.LockedBuffer();

    std::vector<ValueInt<Real>> pairs( k );
    for( Int i=0; i<k; ++i )
    {
        pairs[i].value = xBuffer[i*stride];
        pairs[i].index = i;
    }

    if( sort == ASCENDING )
        std::sort( pairs.begin(), pairs.end(), ValueInt<Real>::Lesser );
    else if( sort == DESCENDING )
        std::sort( pairs.begin(), pairs.end(), ValueInt<Real>::Greater );

    return pairs;
}

template<typename Real,Dist U,Dist V>
std::vector<ValueInt<Real>> TaggedSort
( const DistMatrix<Real,U,V>& x, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("TaggedSort"))
    if( U==STAR && V==STAR )
    {
        return TaggedSort( x.LockedMatrix(), sort );
    }
    else
    {
        DistMatrix<Real,STAR,STAR> x_STAR_STAR( x );
        return TaggedSort( x_STAR_STAR.LockedMatrix(), sort );
    }
}

#define PROTO_DIST(Real,U,V) \
  template void Sort( DistMatrix<Real,U,V>& x, SortType sort ); \
  template std::vector<ValueInt<Real>> TaggedSort \
  ( const DistMatrix<Real,U,V>& x, SortType sort );

#define PROTO(Real) \
  template void Sort( Matrix<Real>& x, SortType sort ); \
  template std::vector<ValueInt<Real>> TaggedSort \
  ( const Matrix<Real>& x, SortType sort ); \
  PROTO_DIST(Real,CIRC,CIRC) \
  PROTO_DIST(Real,MC,  MR  ) \
  PROTO_DIST(Real,MC,  STAR) \
  PROTO_DIST(Real,MD,  STAR) \
  PROTO_DIST(Real,MR,  MC  ) \
  PROTO_DIST(Real,MR,  STAR) \
  PROTO_DIST(Real,STAR,MC  ) \
  PROTO_DIST(Real,STAR,MD  ) \
  PROTO_DIST(Real,STAR,MR  ) \
  PROTO_DIST(Real,STAR,STAR) \
  PROTO_DIST(Real,STAR,VC  ) \
  PROTO_DIST(Real,STAR,VR  ) \
  PROTO_DIST(Real,VC,  STAR) \
  PROTO_DIST(Real,VR,  STAR)

PROTO(Int)
PROTO(float)
PROTO(double)

} // namespace El
