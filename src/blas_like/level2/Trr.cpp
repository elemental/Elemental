/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// A := A + alpha x y'
template<typename T>
void Trr
( UpperOrLower uplo,
  T alpha,
  const Matrix<T>& x,
  const Matrix<T>& y,
        Matrix<T>& A, 
  bool conjugate )
{
    DEBUG_ONLY(
      CSE cse("Trr");
      if( x.Width() != 1 || y.Width() != 1 )
          LogicError("x and y must be of width 1");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    DEBUG_ONLY(
      if( x.Height() != m || y.Height() != n )
          LogicError("x and y must conform with A");
    )
    const T* xCol = x.LockedBuffer();
    const T* yCol = y.LockedBuffer();
    if( uplo == LOWER )
    {
        for( Int j=0; j<n; ++j )
        {
            const T eta = alpha*( conjugate ? Conj(yCol[j]) : yCol[j] );
            T* ACol = A.Buffer(0,j);
            for( Int i=j; i<m; ++i )
                ACol[i] += xCol[i]*eta;
            if( conjugate )
                A.MakeReal( j, j );
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            const T eta = alpha*( conjugate ? Conj(yCol[j]) : yCol[j] );
            T* ACol = A.Buffer(0,j);
            for( Int i=0; i<=j; ++i )
                ACol[i] += xCol[i]*eta;
            if( conjugate )
                A.MakeReal( j, j );
        }
    }
}

template<typename T>
void Trr
( UpperOrLower uplo,
  T alpha,
  const ElementalMatrix<T>& x, 
  const ElementalMatrix<T>& y,
        ElementalMatrix<T>& APre,
  bool conjugate )
{
    DEBUG_ONLY(
      CSE cse("Trr");
      if( x.Width() != 1 || y.Width() != 1 )
          LogicError("x and y must be of width 1");
    )

    DistMatrixReadWriteProxy<T,T,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Grid& g = A.Grid();
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    DEBUG_ONLY(
      if( x.Height() != A.Height() || y.Height() != A.Width() )
          LogicError("x and y must conform with A");
    )

    DistMatrix<T,MC,STAR> x_MC_STAR(g);
    DistMatrix<T,MR,STAR> y_MR_STAR(g);
    x_MC_STAR.AlignWith( A );
    y_MR_STAR.AlignWith( A );
    x_MC_STAR = x;
    y_MR_STAR = y;

    const T* xLocCol = x_MC_STAR.LockedBuffer();
    const T* yLocCol = y_MR_STAR.LockedBuffer();

    if( uplo == LOWER )
    {
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            const Int mLocBefore = A.LocalRowOffset(j);

            const T eta = 
                alpha*( conjugate ? Conj(yLocCol[jLoc]) : yLocCol[jLoc] );
            T* ALocCol = A.Buffer(0,jLoc);
            for( Int iLoc=mLocBefore; iLoc<mLocal; ++iLoc )
                ALocCol[iLoc] += xLocCol[iLoc]*eta;
            if( conjugate )
                A.MakeReal( j, j );
        }
    }
    else
    {
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            const Int mLocBefore = A.LocalRowOffset(j+1);

            const T eta = 
                alpha*( conjugate ? Conj(yLocCol[jLoc]) : yLocCol[jLoc] );
            T* ALocCol = A.Buffer(0,jLoc);
            for( Int iLoc=0; iLoc<mLocBefore; ++iLoc )
                ALocCol[iLoc] += xLocCol[iLoc]*eta;
            if( conjugate )
                A.MakeReal( j, j );
        }
    }
}

#define PROTO(T) \
  template void Trr \
  ( UpperOrLower uplo, \
    T alpha, const Matrix<T>& x, const Matrix<T>& y, \
    Matrix<T>& A, bool conjugate ); \
  template void Trr \
  ( UpperOrLower uplo, \
    T alpha, const ElementalMatrix<T>& x, const ElementalMatrix<T>& y, \
    ElementalMatrix<T>& A, bool conjugate );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
