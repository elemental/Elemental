/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename Field>
void Covariance( const Matrix<Field>& D, Matrix<Field>& S )
{
    EL_DEBUG_CSE
    const Int numObs = D.Height();
    const Int n = D.Width();

    // Compute the average column
    Matrix<Field> ones, xMean;
    Ones( ones, numObs, 1 );
    Gemv( TRANSPOSE, Field(1)/Field(numObs), D, ones, xMean );

    // Subtract the mean from each column of D
    Matrix<Field> DDev( D );
    for( Int i=0; i<numObs; ++i )
        blas::Axpy
        ( n, Field(-1),
          xMean.LockedBuffer(), 1,
          DDev.Buffer(i,0), DDev.LDim() );

    // Form S := 1/(numObs-1) DDev DDev'
    Herk( LOWER, ADJOINT, Base<Field>(1)/Base<Field>(numObs-1), DDev, S );
    Conjugate( S );
    MakeHermitian( LOWER, S );
}

template<typename Field>
void Covariance
( const AbstractDistMatrix<Field>& DPre, AbstractDistMatrix<Field>& SPre )
{
    EL_DEBUG_CSE

    DistMatrixReadProxy<Field,Field,MC,MR>
      DProx( DPre );
    DistMatrixWriteProxy<Field,Field,MC,MR>
      SProx( SPre );
    auto& D = DProx.GetLocked();
    auto& S = SProx.Get();

    const Grid& g = D.Grid();
    const Int numObs = D.Height();

    // Compute the average column
    DistMatrix<Field> ones(g), xMean(g);
    Ones( ones, numObs, 1 );
    Gemv( TRANSPOSE, Field(1)/Field(numObs), D, ones, xMean );
    DistMatrix<Field,MR,STAR> xMean_MR(g);
    xMean_MR.AlignWith( D );
    xMean_MR = xMean;

    // Subtract the mean from each column of D
    DistMatrix<Field> DDev( D );
    for( Int iLoc=0; iLoc<DDev.LocalHeight(); ++iLoc )
        blas::Axpy
        ( DDev.LocalWidth(), Field(-1),
          xMean_MR.LockedBuffer(), 1,
          DDev.Buffer(iLoc,0),     DDev.LDim() );

    // Form S := 1/(numObs-1) DDev DDev'
    Herk( LOWER, ADJOINT, Base<Field>(1)/Base<Field>(numObs-1), DDev, S );
    Conjugate( S );
    MakeHermitian( LOWER, S );
}

#define PROTO(Field) \
  template void Covariance( const Matrix<Field>& D, Matrix<Field>& S ); \
  template void Covariance \
  ( const AbstractDistMatrix<Field>& D, AbstractDistMatrix<Field>& S );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
