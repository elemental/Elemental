/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SVD_HPP
#define ELEM_SVD_HPP

#include ELEM_MAKEHERMITIAN_INC
#include ELEM_HERMITIANEIG_INC

#include "./SVD/Chan.hpp"
#include "./SVD/Thresholded.hpp"

namespace elem {

// Grab the full SVD of the general matrix A, A = U diag(s) V^H
// ============================================================
// NOTE: On exit, A is overwritten with U

template<typename F>
inline void
SVD( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V, bool useQR=false )
{
    DEBUG_ONLY(CallStackEntry cse("SVD"))
    if( useQR )
        svd::QRSVD( A, s, V );
    else
        svd::DivideAndConquerSVD( A, s, V );
}

template<typename F>
inline void HermitianSVD
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
inline void
SVD
( DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& s, DistMatrix<F>& V,
  double heightRatio=1.5 )
{
    DEBUG_ONLY(CallStackEntry cse("SVD"))
    // TODO: Add more options
    svd::Chan( A, s, V, heightRatio );
}

template<typename F>
inline void HermitianSVD
( UpperOrLower uplo, DistMatrix<F>& A, 
  DistMatrix<Base<F>,VR,STAR>& s, DistMatrix<F>& U, DistMatrix<F>& V )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSVD"))
    // Grab an eigenvalue decomposition of A
    HermitianEig( uplo, A, s, V );

    // Redistribute the singular values into an [MR,* ] distribution
    typedef Base<F> Real;
    const Grid& grid = A.Grid();
    DistMatrix<Real,MR,STAR> s_MR_STAR( grid );
    s_MR_STAR.AlignWith( V.DistData() );
    s_MR_STAR = s;

    // Copy V into U (flipping the sign as necessary)
    U.AlignWith( V );
    U.Resize( V.Height(), V.Width() );
    const Int localHeight = V.LocalHeight();
    const Int localWidth = V.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Real sigma = s_MR_STAR.GetLocal( jLoc, 0 );
        F* UCol = U.Buffer( 0, jLoc );
        const F* VCol = V.LockedBuffer( 0, jLoc );
        if( sigma >= 0 )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                UCol[iLoc] = VCol[iLoc];
        else
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                UCol[iLoc] = -VCol[iLoc];
    }

    // Set the singular values to the absolute value of the eigenvalues
    const Int numLocalVals = s.LocalHeight();
    for( Int iLoc=0; iLoc<numLocalVals; ++iLoc )
    {
        const Real sigma = s.GetLocal(iLoc,0);
        s.SetLocal(iLoc,0,Abs(sigma));
    }

    // TODO: Descending sort of triplets
}

// Return the singular values
// ==========================

template<typename F>
inline void
SVD( Matrix<F>& A, Matrix<Base<F>>& s )
{
    DEBUG_ONLY(CallStackEntry cse("SVD"))
    const Int m = A.Height();
    const Int n = A.Width();
    s.Resize( Min(m,n), 1 );
    lapack::SVD( m, n, A.Buffer(), A.LDim(), s.Buffer() );
}

template<typename F>
inline void HermitianSVD
( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& s )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSVD"))
#if 1
    // Grab the eigenvalues of A
    HermitianEig( uplo, A, s );

    // Set the singular values to the absolute value of the eigenvalues
    for( Int i=0; i<s.Height(); ++i )
        s.Set(i,0,Abs(s.Get(i,0)));

    Sort( s, DESCENDING );
#else
    MakeHermitian( uplo, A );
    SVD( A, s );
#endif 
}

template<typename F>
inline void
SVD( DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& s, double heightRatio=1.2 )
{
    DEBUG_ONLY(CallStackEntry cse("SVD"))
    // TODO: Add more options
    svd::Chan( A, s, heightRatio );
}

template<typename F>
inline void HermitianSVD
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& s )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSVD"))
    // Grab the eigenvalues of A
    HermitianEig( uplo, A, s );

    // Replace the eigenvalues with their absolute values
    const Int numLocalVals = s.LocalHeight();
    for( Int iLoc=0; iLoc<numLocalVals; ++iLoc )
    {
        const Base<F> sigma = s.GetLocal(iLoc,0);
        s.SetLocal(iLoc,0,Abs(sigma));
    }

    Sort( s, DESCENDING );
}

} // namespace elem

#endif // ifndef ELEM_SVD_HPP
