/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level2.hpp>

namespace El {

// A := A + alpha X Y'
template<typename T>
void Trr2
( UpperOrLower uplo,
  T alpha,
  const Matrix<T>& X,
  const Matrix<T>& Y,
        Matrix<T>& A, 
  bool conjugate )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( X.Width() != 2 || Y.Width() != 2 )
          LogicError("X and Y must be of width 2");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    EL_DEBUG_ONLY(
      if( X.Height() != m || Y.Height() != n )
          LogicError("X and Y must conform with A");
    )
    const T* XCol0 = X.LockedBuffer(0,0);
    const T* XCol1 = X.LockedBuffer(0,1);
    const T* YCol0 = Y.LockedBuffer(0,0);
    const T* YCol1 = Y.LockedBuffer(0,1);
    if( uplo == LOWER )
    {
        for( Int j=0; j<n; ++j )
        {
            const T eta0 = alpha*( conjugate ? Conj(YCol0[j]) : YCol0[j] );
            const T eta1 = alpha*( conjugate ? Conj(YCol1[j]) : YCol1[j] );
            T* ACol = A.Buffer( 0, j );
            for( Int i=j; i<m; ++i )
                ACol[i] += XCol0[i]*eta0 + XCol1[i]*eta1;
            if( conjugate )
                A.MakeReal( j, j ); 
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            const T eta0 = alpha*( conjugate ? Conj(YCol0[j]) : YCol0[j] );
            const T eta1 = alpha*( conjugate ? Conj(YCol1[j]) : YCol1[j] );
            T* ACol = A.Buffer( 0, j );
            for( Int i=0; i<=j; ++i )
                ACol[i] += XCol0[i]*eta0 + XCol1[i]*eta1;
            if( conjugate )
                A.MakeReal( j, j ); 
        }
    }
}

template<typename T>
void Trr2
( UpperOrLower uplo,
  T alpha,
  const AbstractDistMatrix<T>& XPre,
  const AbstractDistMatrix<T>& YPre,
        AbstractDistMatrix<T>& APre,
  bool conjugate )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( XPre.Width() != 2 || YPre.Width() != 2 )
          LogicError("X and Y must be of width 2");
    )

    DistMatrixReadProxy<T,T,MC,MR> XProx( XPre ), YProx( YPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> AProx( APre );
    auto& X = XProx.GetLocked();
    auto& Y = YProx.GetLocked();
    auto& A = AProx.Get();

    const Grid& g = A.Grid();
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    EL_DEBUG_ONLY(
        if( X.Height() != A.Height() || Y.Height() != A.Width() )
            LogicError("X and Y must conform with A");
    )
    DistMatrix<T,MC,STAR> X_MC_STAR(g);
    DistMatrix<T,MR,STAR> Y_MR_STAR(g);
    X_MC_STAR.AlignWith( A );
    X_MC_STAR = X;
    Y_MR_STAR.AlignWith( A );
    Y_MR_STAR = Y;

    const T* XLocCol0 = X_MC_STAR.LockedBuffer(0,0);
    const T* XLocCol1 = X_MC_STAR.LockedBuffer(0,1);
    const T* YLocCol0 = Y_MR_STAR.LockedBuffer(0,0);
    const T* YLocCol1 = Y_MR_STAR.LockedBuffer(0,1);

    if( uplo == LOWER )
    {
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            const Int mLocBefore = A.LocalRowOffset(j);

            const T value0 = YLocCol0[jLoc]; 
            const T value1 = YLocCol1[jLoc]; 
            const T eta0 = alpha*( conjugate ? Conj(value0) : value0 );
            const T eta1 = alpha*( conjugate ? Conj(value1) : value1 );
            T* ALocCol = A.Buffer(0,jLoc);
            for( Int iLoc=mLocBefore; iLoc<mLocal; ++iLoc )
                ALocCol[iLoc] += XLocCol0[iLoc]*eta0 + XLocCol1[iLoc]*eta1;
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

            const T value0 = YLocCol0[jLoc];
            const T value1 = YLocCol1[jLoc];
            const T eta0 = alpha*( conjugate ? Conj(value0) : value0 );
            const T eta1 = alpha*( conjugate ? Conj(value1) : value1 );
            T* ALocCol = A.Buffer(0,jLoc);
            for( Int iLoc=0; iLoc<mLocBefore; ++iLoc )
                ALocCol[iLoc] += XLocCol0[iLoc]*eta0 + XLocCol1[iLoc]*eta1;
            if( conjugate )
                A.MakeReal( j, j );
        }
    }
}

#define PROTO(T) \
  template void Trr2 \
  ( UpperOrLower uplo, \
    T alpha, const Matrix<T>& X, const Matrix<T>& Y, \
    Matrix<T>& A, bool conjugate ); \
  template void Trr2 \
  ( UpperOrLower uplo, \
    T alpha, const AbstractDistMatrix<T>& X, const AbstractDistMatrix<T>& Y, \
    AbstractDistMatrix<T>& A, bool conjugate );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
