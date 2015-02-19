/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Grab the full SVD of the general matrix A, A = U diag(s) V^H
// ============================================================
// NOTE: A is ovewritten with its packed reduction to tridiagonal form

template<typename F>
void HermitianSVD
( UpperOrLower uplo,
  Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& U, Matrix<F>& V )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSVD"))
#if 1
    // Grab an eigenvalue decomposition of A
    HermitianEig( uplo, A, s, V );

    // Copy V into U (flipping the sign as necessary)
    const Int n = A.Height();
    U.Resize( n, n );
    for( Int j=0; j<n; ++j )
    {
        const Base<F> sigma = s.Get( j, 0 );
        F* UCol = U.Buffer( 0, j );
        const F* VCol = V.LockedBuffer( 0, j );
        if( sigma >= 0 )
        {
            for( Int i=0; i<n; ++i )
                UCol[i] = VCol[i];
        }
        else
        {
            for( Int i=0; i<n; ++i )
                UCol[i] = -VCol[i];
            s.Set( j, 0, -sigma );
        }
    }

    // TODO: Descending sort of triplets
#else
    U = A;
    MakeHermitian( uplo, U );
    SVD( U, s, V );
#endif 
}

// TODO: Add support for HermitianEigSubset and HermitianEigCtrl
template<typename F>
void HermitianSVD
( UpperOrLower uplo, AbstractDistMatrix<F>& A, AbstractDistMatrix<Base<F>>& s, 
  AbstractDistMatrix<F>& U, AbstractDistMatrix<F>& V )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSVD"))

    // Grab an eigenvalue decomposition of A
    HermitianEig( uplo, A, s, V );

    // Copy V into U (flipping the sign as necessary)
    Copy( U, V );
    typedef Base<F> Real;
    DistMatrix<Real,VR,STAR> sSgn( s );
    auto sgnLambda = []( Real sigma ) { return Sgn(sigma,false); };
    EntrywiseMap( sSgn, function<Real(Real)>(sgnLambda) );
    DiagonalScale( RIGHT, NORMAL, sSgn, U );

    // Set the singular values to the absolute value of the eigenvalues
    auto absLambda = []( Real sigma ) { return Abs(sigma); };
    EntrywiseMap( s, function<Real(Real)>(absLambda) );

    // TODO: Descending sort of triplets
}

// Return the singular values
// ==========================
// NOTE: A is ovewritten with its packed reduction to tridiagonal form

template<typename F>
void HermitianSVD( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& s )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSVD"))
#if 1
    typedef Base<F> Real;
    // Grab the eigenvalues of A
    HermitianEig( uplo, A, s );

    // Set the singular values to the absolute value of the eigenvalues
    auto absLambda = []( Real sigma ) { return Abs(sigma); };
    EntrywiseMap( s, function<Real(Real)>(absLambda) );

    Sort( s, DESCENDING );
#else
    MakeHermitian( uplo, A );
    SVD( A, s );
#endif 
}

template<typename F>
void HermitianSVD
( UpperOrLower uplo, AbstractDistMatrix<F>& A, AbstractDistMatrix<Base<F>>& s )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSVD"))

    // Grab the eigenvalues of A
    HermitianEig( uplo, A, s );

    // Set the singular values to the absolute value of the eigenvalues
    typedef Base<F> Real;
    auto absLambda = []( Real sigma ) { return Abs(sigma); };
    EntrywiseMap( s, function<Real(Real)>(absLambda) );

    Sort( s, DESCENDING );
}

#define PROTO(F) \
  template void HermitianSVD \
  ( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& s ); \
  template void HermitianSVD \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<Base<F>>& s ); \
  template void HermitianSVD \
  ( UpperOrLower uplo, Matrix<F>& A, \
    Matrix<Base<F>>& s, Matrix<F>& U, Matrix<F>& V ); \
  template void HermitianSVD \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<Base<F>>& s, AbstractDistMatrix<F>& U, \
    AbstractDistMatrix<F>& V );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
