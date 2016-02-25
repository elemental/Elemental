/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include <algorithm>

namespace El {

// Sort each column of the real matrix X

template<typename Real,typename>
void Sort( Matrix<Real>& X, SortType sort )
{
    DEBUG_ONLY(CSE cse("Sort"))
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

template<typename Real,typename>
void Sort( ElementalMatrix<Real>& X, SortType sort )
{
    DEBUG_ONLY(CSE cse("Sort"))
    if( sort == UNSORTED )
        return;

    if( (X.ColDist()==STAR && X.RowDist()==STAR) || 
        (X.ColDist()==CIRC && X.RowDist()==CIRC) )
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
        Copy( X_CIRC_CIRC, X );
    }
}

// Tagged sort

template<typename Real,typename>
vector<ValueInt<Real>> TaggedSort( const Matrix<Real>& x, SortType sort )
{
    DEBUG_ONLY(CSE cse("TaggedSort"))
    const Int m = x.Height();
    const Int n = x.Width();
    if( m != 1 && n != 1 )
        LogicError("TaggedSort is meant for a single vector");

    const Int k = ( n==1 ? m : n );
    const Int stride = ( n==1 ? 1 : x.LDim() );
    const Real* xBuffer = x.LockedBuffer();

    vector<ValueInt<Real>> pairs( k );
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

template<typename Real,typename>
vector<ValueInt<Real>>
TaggedSort( const ElementalMatrix<Real>& x, SortType sort )
{
    DEBUG_ONLY(CSE cse("TaggedSort"))
    if( x.ColDist()==STAR && x.RowDist()==STAR )
    {
        return TaggedSort( x.LockedMatrix(), sort );
    }
    else
    {
        DistMatrix<Real,STAR,STAR> x_STAR_STAR( x );
        return TaggedSort( x_STAR_STAR.LockedMatrix(), sort );
    }
}

#define PROTO(Real) \
  template void Sort( Matrix<Real>& x, SortType sort ); \
  template void Sort( ElementalMatrix<Real>& x, SortType sort ); \
  template vector<ValueInt<Real>> TaggedSort \
  ( const Matrix<Real>& x, SortType sort ); \
  template vector<ValueInt<Real>> TaggedSort \
  ( const ElementalMatrix<Real>& x, SortType sort );

#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
