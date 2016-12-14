/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>
#include <El/matrices.hpp>

namespace El {

// See Subsection 3.4 of Nguyen and Stehle's "LLL on the Average"

template<typename T>
void AjtaiTypeBasis( Matrix<T>& A, Int n, Base<T> alpha )
{
    EL_DEBUG_CSE
    typedef Base<T> Real;

    Zeros( A, n, n );

    Matrix<Real> d;
    d.Resize( n, 1 );
    for( Int j=0; j<n; ++j )
    {
        const Real exponent = Pow(Real(2)*n-j+1,alpha);
        const Real beta = Round(Pow(Real(2),exponent));
        d(j) = beta;
        A(j,j) = beta;

        for( Int i=0; i<j; ++i )
            A(i,j) = SampleUniform(T(0),T(d(j)/Real(2)));
    }
}

template<typename T>
void AjtaiTypeBasis( AbstractDistMatrix<T>& APre, Int n, Base<T> alpha )
{
    EL_DEBUG_CSE
    typedef Base<T> Real;
    DistMatrixWriteProxy<T,T,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Grid& g = A.Grid();
    Zeros( A, n, n );
    const Int ALocHeight = A.LocalHeight();
    const Int ALocWidth = A.LocalWidth();

    DistMatrix<Real,MR,STAR> d(g);
    d.AlignWith( A );
    d.Resize( n, 1 );
    for( Int jLoc=0; jLoc<ALocWidth; ++jLoc )
    {
        const Int j = A.GlobalCol( jLoc ); 
        const Real exponent = Pow(Real(2)*n-j+1,alpha);
        const Real beta = Round(Pow(Real(2),exponent));
        d.Set( j, 0, beta );
        A.Set( j, j, beta );
    }

    if( A.RedundantRank() == 0 )
    {
        auto& ALoc = A.Matrix();
        auto& dLoc = d.Matrix();
        for( Int jLoc=0; jLoc<ALocWidth; ++jLoc )
        {
            for( Int iLoc=0; iLoc<ALocHeight; ++iLoc )
            {
                ALoc(iLoc,jLoc) = SampleUniform(T(0),T(dLoc(jLoc,0)/Real(2)));
            }
        }
    }
    Broadcast( A, A.RedundantComm(), 0 );
}

#define PROTO(T) \
  template void AjtaiTypeBasis \
  ( Matrix<T>& A, Int n, Base<T> alpha ); \
  template void AjtaiTypeBasis \
  ( AbstractDistMatrix<T>& A, Int n, Base<T> alpha );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
