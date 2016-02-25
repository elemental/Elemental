/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// A := Z Omega Z^H, where Omega is diagonal and real-valued

template<typename F>
void HermitianFromEVD
( UpperOrLower uplo,
  Matrix<F>& A,
  const Matrix<Base<F>>& w,
  const Matrix<F>& Z )
{
    DEBUG_ONLY(CSE cse("HermitianFromEVD"))
    Matrix<F> Z1Copy, Y1;

    const Int m = Z.Height();
    const Int n = Z.Width();
    A.Resize( m, m );
    if( uplo == LOWER )
        MakeTrapezoidal( UPPER, A, 1 );
    else
        MakeTrapezoidal( LOWER, A, -1 );
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        auto Z1 = Z( ALL,        IR(k,k+nb) );
        auto w1 = w( IR(k,k+nb), ALL        );

        Y1 = Z1Copy = Z1;
        DiagonalScale( RIGHT, NORMAL, w1, Y1 );
        Trrk( uplo, NORMAL, ADJOINT, F(1), Z1Copy, Y1, F(1), A );
    }
}

template<typename F>
void HermitianFromEVD
( UpperOrLower uplo, 
        ElementalMatrix<F>& APre,
  const ElementalMatrix<Base<F>>& wPre,
  const ElementalMatrix<F>& ZPre )
{
    DEBUG_ONLY(CSE cse("HermitianFromEVD"))
    typedef Base<F> Real;

    DistMatrixWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixReadProxy<Real,Real,VR,STAR> wProx( wPre );
    DistMatrixReadProxy<F,F,MC,MR> ZProx( ZPre );
    auto& A = AProx.Get();
    auto& w = wProx.GetLocked();
    auto& Z = ZProx.GetLocked();

    const Grid& g = A.Grid();
    DistMatrix<F,MC,  STAR> Z1_MC_STAR(g);
    DistMatrix<F,VR,  STAR> Z1_VR_STAR(g);
    DistMatrix<F,STAR,MR  > Z1Adj_STAR_MR(g);

    const Int m = Z.Height();
    const Int n = Z.Width();
    A.Resize( m, m );
    if( uplo == LOWER )
        MakeTrapezoidal( UPPER, A, 1 );
    else
        MakeTrapezoidal( LOWER, A, -1 );
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        auto Z1 = Z( ALL,        IR(k,k+nb) );
        auto w1 = w( IR(k,k+nb), ALL        );

        Z1_MC_STAR.AlignWith( A );
        Z1_MC_STAR = Z1;
        Z1_VR_STAR.AlignWith( A );
        Z1_VR_STAR = Z1_MC_STAR;

        DiagonalScale( RIGHT, NORMAL, w1, Z1_VR_STAR );

        Z1Adj_STAR_MR.AlignWith( A );
        Adjoint( Z1_VR_STAR, Z1Adj_STAR_MR );
        LocalTrrk( uplo, F(1), Z1_MC_STAR, Z1Adj_STAR_MR, F(1), A );
    }
}

#define PROTO(F) \
  template void HermitianFromEVD \
  ( UpperOrLower uplo, \
          Matrix<F>& A, \
    const Matrix<Base<F>>& w, \
    const Matrix<F>& Z ); \
  template void HermitianFromEVD \
  ( UpperOrLower uplo, \
          ElementalMatrix<F>& A, \
    const ElementalMatrix<Base<F>>& w, \
    const ElementalMatrix<F>& Z );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
