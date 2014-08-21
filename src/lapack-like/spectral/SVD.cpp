/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./SVD/Chan.hpp"
#include "./SVD/Thresholded.hpp"

namespace El {

// Grab the full SVD of the general matrix A, A = U diag(s) V^H
// ============================================================
// NOTE: On exit, A is overwritten with U

template<typename F>
void SVD( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V, bool useQR )
{
    DEBUG_ONLY(CallStackEntry cse("SVD"))
    if( useQR )
        svd::QRSVD( A, s, V );
    else
        svd::DivideAndConquerSVD( A, s, V );
}

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

template<typename F>
void SVD
( AbstractDistMatrix<F>& A, AbstractDistMatrix<Base<F>>& s, 
  AbstractDistMatrix<F>& V, double heightRatio )
{
    DEBUG_ONLY(CallStackEntry cse("SVD"))
    // TODO: Add more options
    svd::Chan( A, s, V, heightRatio );
}

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
    auto sgnLambda = []( Real sigma )
                     { return ( sigma >= 0 ? Real(1) : Real(-1) ); };
    EntrywiseMap( sSgn, std::function<Real(Real)>(sgnLambda) );
    DiagonalScale( RIGHT, NORMAL, sSgn, U );

    // Set the singular values to the absolute value of the eigenvalues
    auto absLambda = []( Real sigma ) { return Abs(sigma); };
    EntrywiseMap( s, std::function<Real(Real)>(absLambda) );

    // TODO: Descending sort of triplets
}

// Return the singular values
// ==========================

template<typename F>
void SVD( Matrix<F>& A, Matrix<Base<F>>& s )
{
    DEBUG_ONLY(CallStackEntry cse("SVD"))
    const Int m = A.Height();
    const Int n = A.Width();
    s.Resize( Min(m,n), 1 );
    lapack::SVD( m, n, A.Buffer(), A.LDim(), s.Buffer() );
}

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
    EntrywiseMap( s, std::function<Real(Real)>(absLambda) );

    Sort( s, DESCENDING );
#else
    MakeHermitian( uplo, A );
    SVD( A, s );
#endif 
}

template<typename F>
void SVD
( AbstractDistMatrix<F>& A, AbstractDistMatrix<Base<F>>& s, double heightRatio )
{
    DEBUG_ONLY(CallStackEntry cse("SVD"))
    // TODO: Add more options
    svd::Chan( A, s, heightRatio );
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
    EntrywiseMap( s, std::function<Real(Real)>(absLambda) );

    Sort( s, DESCENDING );
}

#define PROTO(F) \
  template void SVD \
  ( Matrix<F>& A, Matrix<Base<F>>& s ); \
  template void SVD \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<Base<F>>& s, \
    double heightRatio ); \
  template void SVD \
  ( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V, bool useQR ); \
  template void SVD \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<Base<F>>& s, \
    AbstractDistMatrix<F>& V, double heightRatio ); \
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
    AbstractDistMatrix<F>& V ); \
  template void svd::Thresholded \
  ( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V, \
    Base<F> tol, bool relative ); \
  template void svd::Thresholded \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<Base<F>>& s, \
    AbstractDistMatrix<F>& V, Base<F> tol, bool relative ); \
  template void svd::TallThresholded \
  ( DistMatrix<F,VC,STAR>& A, AbstractDistMatrix<Base<F>>& s, \
    AbstractDistMatrix<F>& V, Base<F> tol, bool relative );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
