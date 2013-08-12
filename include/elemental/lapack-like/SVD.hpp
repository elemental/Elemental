/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_SVD_HPP
#define ELEM_LAPACK_SVD_HPP

#include "elemental/blas-like/level1/MakeHermitian.hpp"
#include "elemental/lapack-like/HermitianEig.hpp"
#include "elemental/lapack-like/SVD/Chan.hpp"
#include "elemental/lapack-like/SVD/Thresholded.hpp"

namespace elem {

//----------------------------------------------------------------------------//
// Grab the full SVD of the general matrix A, A = U diag(s) V^H.              //
// On exit, A is overwritten with U.                                          //
//----------------------------------------------------------------------------//

template<typename F>
inline void
SVD( Matrix<F>& A, Matrix<BASE(F)>& s, Matrix<F>& V, bool useQR=false )
{
#ifndef RELEASE
    CallStackEntry entry("SVD");
#endif
    if( useQR )
        svd::QRSVD( A, s, V );
    else
        svd::DivideAndConquerSVD( A, s, V );
}

template<typename F>
inline void HermitianSVD
( UpperOrLower uplo,
  Matrix<F>& A, Matrix<BASE(F)>& s, Matrix<F>& U, Matrix<F>& V )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianSVD");
#endif
#if 1
    typedef BASE(F) R;

    // Grab an eigenvalue decomposition of A
    HermitianEig( uplo, A, s, V );

    // Set the singular values to the absolute value of the eigenvalues
    for( Int i=0; i<s.Height(); ++i )
        s.Set(i,0,Abs(s.Get(i,0)));

    // Copy V into U (flipping the sign as necessary)
    const Int n = A.Height();
    U.ResizeTo( n, n );
    for( Int j=0; j<n; ++j )
    {
        const R sigma = s.Get( j, 0 );
        F* UCol = U.Buffer( 0, j );
        const F* VCol = V.LockedBuffer( 0, j );
        if( sigma >= 0 )
            for( Int i=0; i<n; ++i )
                UCol[i] = VCol[i];
        else
            for( Int i=0; i<n; ++i )
                UCol[i] = -VCol[i];
    }
#else
    U = A;
    MakeHermitian( uplo, U );
    SVD( U, s, V );
#endif 
}

template<typename F>
inline void
SVD
( DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& s, DistMatrix<F>& V,
  double heightRatio=1.5 )
{
#ifndef RELEASE
    CallStackEntry entry("SVD");
#endif
    // TODO: Add more options
    svd::Chan( A, s, V, heightRatio );
}

template<typename F>
inline void HermitianSVD
( UpperOrLower uplo, DistMatrix<F>& A, 
  DistMatrix<BASE(F),VR,STAR>& s, DistMatrix<F>& U, DistMatrix<F>& V )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianSVD");
#endif
#ifdef HAVE_PMRRR
    typedef BASE(F) R;

    // Grab an eigenvalue decomposition of A
    HermitianEig( uplo, A, s, V );

    // Redistribute the singular values into an [MR,* ] distribution
    const Grid& grid = A.Grid();
    DistMatrix<R,MR,STAR> s_MR_STAR( grid );
    s_MR_STAR.AlignWith( V.DistData() );
    s_MR_STAR = s;

    // Set the singular values to the absolute value of the eigenvalues
    const Int numLocalVals = s.LocalHeight();
    for( Int iLoc=0; iLoc<numLocalVals; ++iLoc )
    {
        const R sigma = s.GetLocal(iLoc,0);
        s.SetLocal(iLoc,0,Abs(sigma));
    }

    // Copy V into U (flipping the sign as necessary)
    U.AlignWith( V );
    U.ResizeTo( V.Height(), V.Width() );
    const Int localHeight = V.LocalHeight();
    const Int localWidth = V.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const R sigma = s_MR_STAR.GetLocal( jLoc, 0 );
        F* UCol = U.Buffer( 0, jLoc );
        const F* VCol = V.LockedBuffer( 0, jLoc );
        if( sigma >= 0 )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                UCol[iLoc] = VCol[iLoc];
        else
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                UCol[iLoc] = -VCol[iLoc];
    }
#else
    U = A;
    MakeHermitian( uplo, U );
    SVD( U, s, V );
#endif // ifdef HAVE_PMRRR
}

//----------------------------------------------------------------------------//
// Grab the singular values of the general matrix A using the QR algorithm.   //
//----------------------------------------------------------------------------//

template<typename F>
inline void
SVD( Matrix<F>& A, Matrix<BASE(F)>& s )
{
#ifndef RELEASE
    CallStackEntry entry("SVD");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    s.ResizeTo( std::min(m,n), 1 );
    lapack::SVD( m, n, A.Buffer(), A.LDim(), s.Buffer() );
}

template<typename F>
inline void HermitianSVD
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& s )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianSVD");
#endif
#if 1
    // Grab the eigenvalues of A
    HermitianEig( uplo, A, s );

    // Set the singular values to the absolute value of the eigenvalues
    for( Int i=0; i<s.Height(); ++i )
        s.Set(i,0,Abs(s.Get(i,0)));
#else
    MakeHermitian( uplo, A );
    SVD( A, s );
#endif 
}

template<typename F>
inline void
SVD( DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& s, double heightRatio=1.2 )
{
#ifndef RELEASE
    CallStackEntry entry("SVD");
#endif
    // TODO: Add more options
    svd::Chan( A, s, heightRatio );
}

template<typename F>
inline void HermitianSVD
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& s )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianSVD");
#endif
#ifdef HAVE_PMRRR
    typedef BASE(F) R;

    // Grab the eigenvalues of A
    HermitianEig( uplo, A, s );

    // Replace the eigenvalues with their absolute values
    const Int numLocalVals = s.LocalHeight();
    for( Int iLoc=0; iLoc<numLocalVals; ++iLoc )
    {
        const R sigma = s.GetLocal(iLoc,0);
        s.SetLocal(iLoc,0,Abs(sigma));
    }
#else
    MakeHermitian( uplo, A );
    SVD( A, s );
#endif // ifdef HAVE_PMRRR
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_SVD_HPP
