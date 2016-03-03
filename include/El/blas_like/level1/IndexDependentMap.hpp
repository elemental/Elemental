/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_INDEXDEPENDENTMAP_HPP
#define EL_BLAS_INDEXDEPENDENTMAP_HPP

namespace El {

template<typename T>
void IndexDependentMap( Matrix<T>& A, function<T(Int,Int,T)> func )
{
    DEBUG_ONLY(CSE cse("IndexDependentMap"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, func(i,j,A.Get(i,j)) );
}

template<typename T>
void IndexDependentMap
( AbstractDistMatrix<T>& A, function<T(Int,Int,T)> func )
{
    DEBUG_ONLY(CSE cse("IndexDependentMap"))
    const Int mLoc = A.LocalHeight();
    const Int nLoc = A.LocalWidth();
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<mLoc; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            A.SetLocal( iLoc, jLoc, func(i,j,A.GetLocal(iLoc,jLoc)) );
        }
    }
}

template<typename S,typename T>
void IndexDependentMap
( const Matrix<S>& A, Matrix<T>& B, function<T(Int,Int,S)> func )
{
    DEBUG_ONLY(CSE cse("IndexDependentMap"))
    const Int m = A.Height();
    const Int n = A.Width();
    B.Resize( m, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            B.Set( i, j, func(i,j,A.Get(i,j)) );
}

template<typename S,typename T>
void IndexDependentMap
( const ElementalMatrix<S>& A, ElementalMatrix<T>& B, 
  function<T(Int,Int,S)> func )
{
    DEBUG_ONLY(CSE cse("IndexDependentMap"))
    const Int mLoc = A.LocalHeight();
    const Int nLoc = A.LocalWidth();
    B.AlignWith( A.DistData() );
    B.Resize( A.Height(), A.Width() );
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<mLoc; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            B.SetLocal( iLoc, jLoc, func(i,j,A.GetLocal(iLoc,jLoc)) );
        }
    }
}

template<typename S,typename T>
void IndexDependentMap
( const BlockMatrix<S>& A, BlockMatrix<T>& B, 
  function<T(Int,Int,S)> func )
{
    DEBUG_ONLY(CSE cse("IndexDependentMap"))
    const Int mLoc = A.LocalHeight();
    const Int nLoc = A.LocalWidth();
    B.AlignWith( A.DistData() );
    B.Resize( A.Height(), A.Width() );
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<mLoc; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            B.SetLocal( iLoc, jLoc, func(i,j,A.GetLocal(iLoc,jLoc)) );
        }
    }
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void IndexDependentMap \
  ( Matrix<T>& A, \
    function<T(Int,Int,T)> func ); \
  EL_EXTERN template void IndexDependentMap \
  ( AbstractDistMatrix<T>& A, \
    function<T(Int,Int,T)> func ); \
  EL_EXTERN template void IndexDependentMap \
  ( const Matrix<T>& A, \
          Matrix<T>& B, \
          function<T(Int,Int,T)> func ); \
  EL_EXTERN template void IndexDependentMap \
  ( const ElementalMatrix<T>& A, \
          ElementalMatrix<T>& B, \
          function<T(Int,Int,T)> func ); \
  EL_EXTERN template void IndexDependentMap \
  ( const BlockMatrix<T>& A, \
          BlockMatrix<T>& B, \
          function<T(Int,Int,T)> func );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_INDEXDEPENDENTMAP_HPP
