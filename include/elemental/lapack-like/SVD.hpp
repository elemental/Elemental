/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_SVD_HPP
#define LAPACK_SVD_HPP

#include "elemental/lapack-like/SVD/Chan.hpp"
#include "elemental/lapack-like/SVD/Thresholded.hpp"

namespace elem {

//----------------------------------------------------------------------------//
// Grab the full SVD of the general matrix A, A = U diag(s) V^H.              //
// On exit, A is overwritten with U.                                          //
//----------------------------------------------------------------------------//

template<typename F>
inline void
SVD
( Matrix<F>& A, Matrix<typename Base<F>::type>& s, Matrix<F>& V, 
  bool useQR=false )
{
#ifndef RELEASE
    PushCallStack("SVD");
#endif
    if( useQR )
        svd::QRSVD( A, s, V );
    else
        svd::DivideAndConquerSVD( A, s, V );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void HermitianSVD
( UpperOrLower uplo,
  Matrix<F>& A, Matrix<typename Base<F>::type>& s,
  Matrix<F>& U, Matrix<F>& V )
{
#ifndef RELEASE
    PushCallStack("HermitianSVD");
#endif
    // TODO: Make use of a sequential Hermitian eigensolver
    U = A;
    MakeHermitian( uplo, U );
    SVD( U, s, V );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
SVD
( DistMatrix<F>& A,
  DistMatrix<typename Base<F>::type,VR,STAR>& s,
  DistMatrix<F>& V,
  double heightRatio=1.5 )
{
#ifndef RELEASE
    PushCallStack("SVD");
#endif
    // TODO: Add more options
    svd::Chan( A, s, V, heightRatio );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void HermitianSVD
( UpperOrLower uplo,
  DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& s,
  DistMatrix<F>& U, DistMatrix<F>& V )
{
#ifndef RELEASE
    PushCallStack("HermitianSVD");
#endif
#ifdef HAVE_PMRRR
    typedef typename Base<F>::type R;

    // Grab an eigenvalue decomposition of A
    HermitianEig( uplo, A, s, V );

    // Redistribute the singular values into an [MR,* ] distribution
    const Grid& grid = A.Grid();
    DistMatrix<R,MR,STAR> s_MR_STAR( grid );
    s_MR_STAR.AlignWith( V.DistData() );
    s_MR_STAR = s;

    // Set the singular values to the absolute value of the eigenvalues
    const int numLocalVals = s.LocalHeight();
    for( int iLocal=0; iLocal<numLocalVals; ++iLocal )
    {
        const R sigma = s.GetLocal(iLocal,0);
        s.SetLocal(iLocal,0,Abs(sigma));
    }

    // Copy V into U (flipping the sign as necessary)
    U.AlignWith( V );
    U.ResizeTo( V.Height(), V.Width() );
    const int localHeight = V.LocalHeight();
    const int localWidth = V.LocalWidth();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const R sigma = s_MR_STAR.GetLocal( jLocal, 0 );
        F* UCol = U.Buffer( 0, jLocal );
        const F* VCol = V.LockedBuffer( 0, jLocal );
        if( sigma >= 0 )
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                UCol[iLocal] = VCol[iLocal];
        else
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                UCol[iLocal] = -VCol[iLocal];
    }
#else
    U = A;
    MakeHermitian( uplo, U );
    SVD( U, s, V );
#endif // ifdef HAVE_PMRRR
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Grab the singular values of the general matrix A using the QR algorithm.   //
//----------------------------------------------------------------------------//

template<typename F>
inline void
SVD( Matrix<F>& A, Matrix<typename Base<F>::type>& s )
{
#ifndef RELEASE
    PushCallStack("SVD");
#endif
    typedef typename Base<F>::type R;

    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min(m,n);
    s.ResizeTo( k, 1 );
    lapack::SVD( m, n, A.Buffer(), A.LDim(), s.Buffer() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void HermitianSVD
( UpperOrLower uplo,
  Matrix<F>& A, Matrix<typename Base<F>::type>& s )
{
#ifndef RELEASE
    PushCallStack("HermitianSVD");
#endif
    // TODO: Make use of a sequential Hermitian eigensolver
    MakeHermitian( uplo, A );
    SVD( A, s );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
SVD
( DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& s, 
  double heightRatio=1.2 )
{
#ifndef RELEASE
    PushCallStack("SVD");
#endif
    // TODO: Add more options
    svd::Chan( A, s, heightRatio );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void HermitianSVD
( UpperOrLower uplo,
  DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& s )
{
#ifndef RELEASE
    PushCallStack("HermitianSVD");
#endif
#ifdef HAVE_PMRRR
    typedef typename Base<F>::type R;

    // Grab an eigenvalue decomposition of A
    HermitianEig( uplo, A, s );

    // Replace the eigenvalues with their absolute values
    const int numLocalVals = s.LocalHeight();
    for( int iLocal=0; iLocal<numLocalVals; ++iLocal )
    {
        const R sigma = s.GetLocal(iLocal,0);
        s.SetLocal(iLocal,0,Abs(sigma));
    }
#else
    MakeHermitian( uplo, A );
    SVD( A, s );
#endif // ifdef HAVE_PMRRR
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef LAPACK_SVD_HPP
