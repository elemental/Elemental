/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F>
void Covariance( const Matrix<F>& D, Matrix<F>& S )
{
    DEBUG_ONLY(CSE cse("Covariance"))
    const Int numObs = D.Height();
    const Int n = D.Width();

    // Compute the average column
    Matrix<F> ones, xMean;
    Ones( ones, numObs, 1 );
    Gemv( TRANSPOSE, F(1)/F(numObs), D, ones, xMean );

    // Subtract the mean from each column of D
    Matrix<F> DDev( D );
    for( Int i=0; i<numObs; ++i )
        blas::Axpy
        ( n, F(-1), xMean.LockedBuffer(), 1, DDev.Buffer(i,0), DDev.LDim() );

    // Form S := 1/(numObs-1) DDev DDev'
    Herk( LOWER, ADJOINT, Base<F>(1)/Base<F>(numObs-1), DDev, S );
    Conjugate( S );
    MakeHermitian( LOWER, S );
}

template<typename F>
void Covariance
( const ElementalMatrix<F>& DPre, ElementalMatrix<F>& SPre )
{
    DEBUG_ONLY(CSE cse("Covariance"))

    DistMatrixReadProxy<F,F,MC,MR>
      DProx( DPre );
    DistMatrixWriteProxy<F,F,MC,MR>
      SProx( SPre );
    auto& D = DProx.GetLocked();
    auto& S = SProx.Get();

    const Grid& g = D.Grid();
    const Int numObs = D.Height();

    // Compute the average column
    DistMatrix<F> ones(g), xMean(g);
    Ones( ones, numObs, 1 );
    Gemv( TRANSPOSE, F(1)/F(numObs), D, ones, xMean );
    DistMatrix<F,MR,STAR> xMean_MR_STAR(g);
    xMean_MR_STAR.AlignWith( D );
    xMean_MR_STAR = xMean;

    // Subtract the mean from each column of D
    DistMatrix<F> DDev( D );
    for( Int iLoc=0; iLoc<DDev.LocalHeight(); ++iLoc )
        blas::Axpy
        ( DDev.LocalWidth(), F(-1), 
          xMean_MR_STAR.LockedBuffer(), 1, 
          DDev.Buffer(iLoc,0),          DDev.LDim() );

    // Form S := 1/(numObs-1) DDev DDev'
    Herk( LOWER, ADJOINT, Base<F>(1)/Base<F>(numObs-1), DDev, S );
    Conjugate( S );
    MakeHermitian( LOWER, S );
}

#define PROTO(F) \
  template void Covariance( const Matrix<F>& D, Matrix<F>& S ); \
  template void Covariance \
  ( const ElementalMatrix<F>& D, ElementalMatrix<F>& S );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
