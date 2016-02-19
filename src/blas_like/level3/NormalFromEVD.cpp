/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// A := Z Omega Z^H, where Omega is complex-valued and diagonal

template<typename Real>
void NormalFromEVD
(       Matrix<Complex<Real>>& A,
  const Matrix<Complex<Real>>& w,
  const Matrix<Complex<Real>>& Z )
{
    DEBUG_ONLY(CSE cse("NormalFromEVD"))
    typedef Complex<Real> C;

    Matrix<C> Y1, Z1Copy;

    const Int m = Z.Height();
    const Int n = Z.Width();
    Zeros( A, m, m );
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        auto Z1 = Z( ALL,        IR(k,k+nb) );
        auto w1 = w( IR(k,k+nb), ALL        );

        Y1 = Z1Copy = Z1;
        DiagonalScale( RIGHT, ADJOINT, w1, Y1 );
        Gemm( NORMAL, NORMAL, C(1), Z1Copy, Y1, C(1), A );
    }
}

template<typename Real>
void NormalFromEVD
(       ElementalMatrix<Complex<Real>>& APre,
  const ElementalMatrix<Complex<Real>>& wPre, 
  const ElementalMatrix<Complex<Real>>& ZPre )
{
    DEBUG_ONLY(CSE cse("NormalFromEVD"))
    typedef Complex<Real> C;

    DistMatrixWriteProxy<C,C,MC,MR> AProx( APre );
    DistMatrixReadProxy<C,C,VR,STAR> wProx( wPre );
    DistMatrixReadProxy<C,C,MC,MR> ZProx( ZPre );
    auto& A = AProx.Get();
    auto& w = wProx.GetLocked();
    auto& Z = ZProx.GetLocked();

    const Grid& g = A.Grid();
    DistMatrix<C,MC,  STAR> Z1_MC_STAR(g);
    DistMatrix<C,VR,  STAR> Z1_VR_STAR(g);
    DistMatrix<C,STAR,MR  > Z1Adj_STAR_MR(g);

    const Int m = Z.Height();
    const Int n = Z.Width();
    Zeros( A, m, m );
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

        DiagonalScale( RIGHT, ADJOINT, w1, Z1_VR_STAR );

        Z1Adj_STAR_MR.AlignWith( A );
        Adjoint( Z1_VR_STAR, Z1Adj_STAR_MR );
        LocalGemm( NORMAL, NORMAL, C(1), Z1_MC_STAR, Z1Adj_STAR_MR, C(1), A );
    }
}

#define PROTO(Real) \
  template void NormalFromEVD \
  (       Matrix<Complex<Real>>& A, \
    const Matrix<Complex<Real>>& w, \
    const Matrix<Complex<Real>>& Z ); \
  template void NormalFromEVD \
  (       ElementalMatrix<Complex<Real>>& A, \
    const ElementalMatrix<Complex<Real>>& w, \
    const ElementalMatrix<Complex<Real>>& Z );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
