/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include <algorithm>

namespace El {

// Sort each column of the real matrix X

template<typename Real,typename>
void Sort( Matrix<Real>& X, SortType sort, bool stable )
{
    DEBUG_CSE
    if( sort == UNSORTED )
        return;
    const Int m = X.Height();
    const Int n = X.Width();
    for( Int j=0; j<n; ++j )
    {
        Real* XCol = X.Buffer(0,j);
        if( sort == ASCENDING )
        {
            if( stable )
                std::stable_sort( XCol, XCol+m );
            else
                std::sort( XCol, XCol+m );
        }
        else
        {
            if( stable )
                std::stable_sort( XCol, XCol+m, std::greater<Real>() );
            else
                std::sort( XCol, XCol+m, std::greater<Real>() );
        }
    }
}

template<typename Real,typename>
void Sort( AbstractDistMatrix<Real>& X, SortType sort, bool stable )
{
    DEBUG_CSE
    if( sort == UNSORTED )
        return;

    if( (X.ColDist()==STAR && X.RowDist()==STAR) || 
        (X.ColDist()==CIRC && X.RowDist()==CIRC) )
    {
        if( X.Participating() )
            Sort( X.Matrix(), sort, stable );
    }
    else
    {
        // TODO(poulson): Distributed sort

        // Get a copy on a single process, sort, and then redistribute
        DistMatrix<Real,CIRC,CIRC> X_CIRC_CIRC( X );
        if( X_CIRC_CIRC.Participating() )
            Sort( X_CIRC_CIRC.Matrix(), sort, stable );

        // Refill the distributed X with the sorted values
        Copy( X_CIRC_CIRC, X );
    }
}

// Tagged sort

template<typename Real,typename>
vector<ValueInt<Real>>
TaggedSort( const Matrix<Real>& x, SortType sort, bool stable )
{
    DEBUG_CSE
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
    {
        if( stable )
            std::stable_sort
            ( pairs.begin(), pairs.end(), ValueInt<Real>::Lesser );
        else
            std::sort( pairs.begin(), pairs.end(), ValueInt<Real>::Lesser );
    }
    else if( sort == DESCENDING )
    {
        if( stable )
            std::stable_sort
            ( pairs.begin(), pairs.end(), ValueInt<Real>::Greater );
        else
            std::sort( pairs.begin(), pairs.end(), ValueInt<Real>::Greater );
    }

    return pairs;
}

template<typename Real,typename>
vector<ValueInt<Real>>
TaggedSort( const AbstractDistMatrix<Real>& x, SortType sort, bool stable )
{
    DEBUG_CSE
    if( x.ColDist()==STAR && x.RowDist()==STAR )
    {
        return TaggedSort( x.LockedMatrix(), sort, stable );
    }
    else
    {
        DistMatrix<Real,STAR,STAR> x_STAR_STAR( x );
        return TaggedSort( x_STAR_STAR.LockedMatrix(), sort, stable );
    }
}

template<typename Real,typename F>
void ApplyTaggedSortToEachRow
( const vector<ValueInt<Real>>& sortPairs,
        Matrix<F>& Z )
{
    DEBUG_CSE
    const Int m = Z.Height();
    const Int n = Z.Width();
    Matrix<F> ZPerm( m, n );
    for( Int j=0; j<n; ++j )
    {
        const Int source = sortPairs[j].index;
        MemCopy( ZPerm.Buffer(0,j), Z.LockedBuffer(0,source), m );
    }
    Z = ZPerm;
}

template<typename Real,typename F>
void ApplyTaggedSortToEachColumn
( const vector<ValueInt<Real>>& sortPairs,
        Matrix<F>& Z )
{
    DEBUG_CSE
    const Int m = Z.Height();
    const Int n = Z.Width();
    Matrix<F> ZPerm( m, n );
    for( Int i=0; i<m; ++i )
    {
        const Int source = sortPairs[i].index;
        for( Int j=0; j<n; ++j )
            ZPerm(i,j) = Z(source,j);
    }
    Z = ZPerm;
}

template<typename Real,typename F>
void ApplyTaggedSortToEachRow
( const vector<ValueInt<Real>>& sortPairs,
        AbstractDistMatrix<F>& Z )
{
    DEBUG_CSE
    const Int m = Z.Height();
    const Int n = Z.Width();
    DistMatrix<F,VC,STAR> Z_VC_STAR( Z );
    DistMatrix<F,VC,STAR> ZPerm_VC_STAR(Z.Grid());
    ZPerm_VC_STAR.AlignWith( Z_VC_STAR );
    ZPerm_VC_STAR.Resize( m, n );
    const Int mLocal = Z_VC_STAR.LocalHeight();
    for( Int j=0; j<n; ++j )
    {
        const Int source = sortPairs[j].index;
        MemCopy
        ( ZPerm_VC_STAR.Buffer(0,j),
          Z_VC_STAR.LockedBuffer(0,source), mLocal );
    }
    Z_VC_STAR.Empty();
    Copy( ZPerm_VC_STAR, Z );
}

template<typename Real,typename F>
void ApplyTaggedSortToEachColumn
( const vector<ValueInt<Real>>& sortPairs,
        AbstractDistMatrix<F>& Z )
{
    DEBUG_CSE
    const Int m = Z.Height();
    const Int n = Z.Width();
    DistMatrix<F,STAR,VR> Z_STAR_VR( Z );
    DistMatrix<F,STAR,VR> ZPerm_STAR_VR(Z.Grid());
    ZPerm_STAR_VR.AlignWith( Z_STAR_VR );
    ZPerm_STAR_VR.Resize( m, n );
    const Int nLocal = Z_STAR_VR.LocalWidth();
    for( Int i=0; i<m; ++i )
    {
        const Int source = sortPairs[i].index;
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            ZPerm_STAR_VR.SetLocal( i, jLoc, Z_STAR_VR.GetLocal(source,jLoc) );
    }
    Z_STAR_VR.Empty();
    Copy( ZPerm_STAR_VR, Z );
}

template<typename Real,typename>
void SortingPermutation
( const Matrix<Real>& x, Permutation& sortPerm, SortType sort, bool stable )
{
    DEBUG_CSE
    auto sortPairs = TaggedSort( x, sort, stable );
    const Int n = ( x.Width()==1 ? x.Height() : x.Width() );
    sortPerm.MakeIdentity( n );
    for( Int i=0; i<n; ++i )
        sortPerm.SetImage( sortPairs[i].index, i );
}

template<typename Real,typename>
void MergeSortingPermutation
( Int n0, Int n1, const Matrix<Real>& x, Permutation& sortPerm, SortType sort )
{
    DEBUG_CSE
    const Int m = x.Height();
    const Int n = x.Width();
    if( m != 1 && n != 1 )
        LogicError("MergeSortingPermutation meant for a vector");

    const Int k = ( n==1 ? m : n );
    if( k != n0+n1 )
        LogicError("Dimensions did not match");

    const Int stride = ( n==1 ? 1 : x.LDim() );
    const Real* xBuffer = x.LockedBuffer();

    vector<ValueInt<Real>> pairs( k );
    for( Int i=0; i<k; ++i )
    {
        pairs[i].value = xBuffer[i*stride];
        pairs[i].index = i;
    }

    if( sort == ASCENDING )
        std::inplace_merge
        ( pairs.begin(), pairs.begin()+n0, pairs.end(),
          ValueInt<Real>::Lesser );
    else if( sort == DESCENDING )
        std::inplace_merge
        ( pairs.begin(), pairs.begin()+n0, pairs.end(),
          ValueInt<Real>::Greater );

    sortPerm.MakeIdentity( k );
    for( Int i=0; i<k; ++i )
        sortPerm.SetImage( pairs[i].index, i );
}

#define PROTO_COMPLEX(F) \
  template void ApplyTaggedSortToEachRow \
  ( const vector<ValueInt<Base<F>>>& sortPairs, \
          Matrix<F>& Z ); \
  template void ApplyTaggedSortToEachColumn \
  ( const vector<ValueInt<Base<F>>>& sortPairs, \
          Matrix<F>& Z ); \
  template void ApplyTaggedSortToEachRow \
  ( const vector<ValueInt<Base<F>>>& sortPairs, \
          AbstractDistMatrix<F>& Z ); \
  template void ApplyTaggedSortToEachColumn \
  ( const vector<ValueInt<Base<F>>>& sortPairs, \
          AbstractDistMatrix<F>& Z );

#define PROTO(Real) \
  PROTO_COMPLEX(Real) \
  template void Sort( Matrix<Real>& x, SortType sort, bool stable ); \
  template void Sort \
  ( AbstractDistMatrix<Real>& x, SortType sort, bool stable ); \
  template vector<ValueInt<Real>> TaggedSort \
  ( const Matrix<Real>& x, SortType sort, bool stable ); \
  template vector<ValueInt<Real>> TaggedSort \
  ( const AbstractDistMatrix<Real>& x, SortType sort, bool stable ); \
  template void SortingPermutation \
  ( const Matrix<Real>& x, Permutation& sortPerm, SortType sort, bool stable );

// For support for double-precision MRRR with float eigenvectors

#define PROTO_FLOAT \
  PROTO(float) \
  template void ApplyTaggedSortToEachRow \
  ( const vector<ValueInt<float>>& sortPairs, \
    AbstractDistMatrix<double>& Z ); \
  template void ApplyTaggedSortToEachColumn \
  ( const vector<ValueInt<float>>& sortPairs, \
    AbstractDistMatrix<double>& Z );

#define PROTO_COMPLEX_FLOAT \
  PROTO_COMPLEX(Complex<float>) \
  template void ApplyTaggedSortToEachRow \
  ( const vector<ValueInt<float>>& sortPairs, \
    AbstractDistMatrix<Complex<double>>& Z ); \
  template void ApplyTaggedSortToEachColumn \
  ( const vector<ValueInt<float>>& sortPairs, \
    AbstractDistMatrix<Complex<double>>& Z );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
