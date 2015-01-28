/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// A :=  Z Omega Z^T, where Omega is diagonal and real-valued

template<typename F>
void HermitianFromEVD
( UpperOrLower uplo, Matrix<F>& A,
  const Matrix<Base<F>>& w, const Matrix<F>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianFromEVD"))
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
        auto Z1 = Z( IR(0,m),    IR(k,k+nb) );
        auto w1 = w( IR(k,k+nb), IR(0,1)    );

        Y1 = Z1Copy = Z1;
        DiagonalScale( RIGHT, NORMAL, w1, Y1 );
        Trrk( uplo, NORMAL, ADJOINT, F(1), Z1Copy, Y1, F(1), A );
    }
}

template<typename F>
void HermitianFromEVD
( UpperOrLower uplo, AbstractDistMatrix<F>& APre,
  const AbstractDistMatrix<Base<F>>& wPre, const AbstractDistMatrix<F>& ZPre )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianFromEVD"))
    typedef Base<F> Real;

    auto APtr = WriteProxy<F,MC,MR>( &APre );     auto& A = *APtr;
    auto wPtr = ReadProxy<Real,VR,STAR>( &wPre ); auto& w = *wPtr;
    auto ZPtr = ReadProxy<F,MC,MR>( &ZPre );      auto& Z = *ZPtr;

    const Grid& g = A.Grid();
    DistMatrix<F,MC,  STAR> Z1_MC_STAR(g);
    DistMatrix<F,VR,  STAR> Z1_VR_STAR(g);
    DistMatrix<F,STAR,MR  > Z1Adj_STAR_MR(g);
    DistMatrix<Real,STAR,STAR> w1_STAR_STAR(g);

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
        auto Z1 = Z( IR(0,m),    IR(k,k+nb) );
        auto w1 = w( IR(k,k+nb), IR(0,1)    );

        Z1_MC_STAR.AlignWith( A );
        Z1_MC_STAR = Z1;
        Z1_VR_STAR.AlignWith( A );
        Z1_VR_STAR = Z1_MC_STAR;
        w1_STAR_STAR = w1;

        DiagonalScale( RIGHT, NORMAL, w1_STAR_STAR, Z1_VR_STAR );

        Z1Adj_STAR_MR.AlignWith( A );
        Adjoint( Z1_VR_STAR, Z1Adj_STAR_MR );
        LocalTrrk( uplo, F(1), Z1_MC_STAR, Z1Adj_STAR_MR, F(1), A );
    }
}

#define PROTO(F) \
  template void HermitianFromEVD \
  ( UpperOrLower uplo, Matrix<F>& A, \
    const Matrix<Base<F>>& w, const Matrix<F>& Z ); \
  template void HermitianFromEVD \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<Base<F>>& w, const AbstractDistMatrix<F>& Z );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
