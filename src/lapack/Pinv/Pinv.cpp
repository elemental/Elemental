/*
   Copyright (c) 2009-2010, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental/blas_internal.hpp"
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;

// Assumes A is a tall, skinny m x n matrix that is small enough for all cores 
// to keep the upper n x n portion of A in main memory.
template<typename R>
void
elemental::lapack::Pinv
( DistMatrix<R,MC,MR>& A, DistMatrix<R,MC,MR>& PinvA )
{
#ifndef RELEASE
    PushCallStack("lapack::Pinv");
    if( A.Height() < A.Width() )
        throw std::runtime_error("A is assumed to be tall and skinny.");
#endif
    const unsigned m = A.Height();
    const unsigned n = A.Width();
    const Grid& g = A.GetGrid();

    // Fill the lower triangle of A with Householder reflectors defining Q^T
    // and the upper triangle with R
    lapack::QR( A );

    // Partition off the top n x n portion of A
    DistMatrix<R,MC,MR> ATop(g), ABottom(g);
    PartitionDown( A, ATop,
                      ABottom, n );

    // Give all processes a copy of R (with zeroes below the diagonal)
    DistMatrix<R,Star,Star> R_Star_Star(g);
    R_Star_Star = ATop;
    R_Star_Star.MakeTrapezoidal( Left, Upper );

    // Redundantly overwrite R with its pseudo-inverse
    {
        std::vector<R> SigmaDiag;
        Matrix<R> U, VT;
        lapack::SVD( R_Star_Star.LocalMatrix(), U, VT, SigmaDiag );

        // Threshold and invert Sigma with max(m,n) ||A||_2 eps, where epsilon 
        // is machine-precision
        R ANorm = 0;
        for( unsigned i=0; i<n; ++i )
            ANorm = std::max( ANorm, SigmaDiag[i] );
        R eps = std::numeric_limits<R>::epsilon();
        R threshold = m*ANorm*eps;
        for( unsigned i=0; i<n; ++i )
        {
            if( SigmaDiag[i] >= threshold )
                SigmaDiag[i] = 1 / SigmaDiag[i];
            else
                SigmaDiag[i] = 0;
        }

        // Scale U by our thresholded and inverted Sigma
        for( unsigned j=0; j<n; ++j )
        {
            R alpha = SigmaDiag[j];
            R* thisColumn = U.Buffer(0,j);
            for( unsigned i=0; i<m; ++i )
                thisColumn[i] *= alpha;
        }

        // Overwrite R with our scaled U times VT
        blas::Gemm
        ( Normal, Normal, (R)1, U, VT, (R)0, R_Star_Star.LocalMatrix() );
    }

    // Form the pseudo-inverse of A:
    //   We have the Householder vectors for Q^T stored in the lower portion of
    //   A, so we need to form    
    //     PinvA = H_{n-1} ... H_0 | I | PinvR
    //                             | 0 |
    //   We can efficiently do so by multiplying the two right-most matrices 
    //   first,    | I | PinvR = | PinvR |,
    //             | 0 |         |   0   |
    //   and then applying the Householder transforms to this matrix.
    PinvA.AlignWith( A );
    PinvA.ResizeTo( m, n );
    DistMatrix<R,MC,MR> PinvATop(g), PinvABottom(g);
    PartitionDown( PinvA, PinvATop,
                          PinvABottom, n );
    PinvATop = R_Star_Star;
    PinvABottom.SetToZero();
    lapack::UT( Left, Lower, Transpose, 0, A, PinvA );

#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
// Assumes A is a tall, skinny m x n matrix that is small enough for all cores 
// to keep the upper n x n portion of A in main memory.
template<typename R>
void
elemental::lapack::Pinv
( DistMatrix< std::complex<R>,MC,MR>& A, 
  DistMatrix< std::complex<R>,MC,MR>& PinvA )
{
#ifndef RELEASE
    PushCallStack("lapack::Pinv");
    if( A.Height() < A.Width() )
        throw std::runtime_error("A is assumed to be tall and skinny.");
#endif
    typedef std::complex<R> C;
    const unsigned m = A.Height();
    const unsigned n = A.Width();
    const Grid& g = A.GetGrid();

    // Fill the lower triangle of A with Householder reflectors defining Q^T
    // and the upper triangle with R
    DistMatrix<C,MD,Star> t(g);
    lapack::QR( A, t );

    // Partition off the top n x n portion of A
    DistMatrix<C,MC,MR> ATop(g), ABottom(g);
    PartitionDown( A, ATop,
                      ABottom, n );

    // Give all processes a copy of R (with zeroes below the diagonal)
    DistMatrix<C,Star,Star> R_Star_Star(g);
    R_Star_Star = ATop;
    R_Star_Star.MakeTrapezoidal( Left, Upper );

    // Redundantly overwrite R with its pseudo-inverse
    {
        std::vector<R> SigmaDiag;
        Matrix<C> U, VT;
        lapack::SVD( R_Star_Star.LocalMatrix(), U, VT, SigmaDiag );

        // Threshold and invert Sigma with max(m,n) ||A||_2 eps, where epsilon 
        // is machine-precision
        R ANorm = 0;
        for( unsigned i=0; i<n; ++i )
            ANorm = std::max( ANorm, SigmaDiag[i] );
        R eps = std::numeric_limits<R>::epsilon();
        R threshold = m*ANorm*eps;
        for( unsigned i=0; i<n; ++i )
        {
            if( SigmaDiag[i] >= threshold )
                SigmaDiag[i] = 1 / SigmaDiag[i];
            else
                SigmaDiag[i] = 0;
        }

        // Scale U by our thresholded and inverted Sigma
        for( unsigned j=0; j<n; ++j )
        {
            R alpha = SigmaDiag[j];
            C* thisColumn = U.Buffer(0,j);
            for( unsigned i=0; i<m; ++i )
                thisColumn[i] *= alpha;
        }

        // Overwrite R with our scaled U times VT
        blas::Gemm
        ( Normal, Normal, (C)1, U, VT, (C)0, R_Star_Star.LocalMatrix() );
    }

    // Form the pseudo-inverse of A:
    //   We have the Householder vectors for Q^T stored in the lower portion of
    //   A, so we need to form    
    //     PinvA = H_{n-1} ... H_0 | I | PinvR
    //                             | 0 |
    //   We can efficiently do so by multiplying the two right-most matrices 
    //   first,    | I | PinvR = | PinvR |,
    //             | 0 |         |   0   |
    //   and then applying the Householder transforms to this matrix.
    PinvA.AlignWith( A );
    PinvA.ResizeTo( m, n );
    DistMatrix<C,MC,MR> PinvATop(g), PinvABottom(g);
    PartitionDown( PinvA, PinvATop,
                          PinvABottom, n );
    PinvATop = R_Star_Star;
    PinvABottom.SetToZero();
    lapack::UT( Left, Lower, ConjugateTranspose, 0, A, t, PinvA );

#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

template void
elemental::lapack::Pinv
( DistMatrix<float,MC,MR>& A, DistMatrix<float,MC,MR>& PinvA );

template void
elemental::lapack::Pinv
( DistMatrix<double,MC,MR>& A, DistMatrix<double,MC,MR>& PinvA );

#ifndef WITHOUT_COMPLEX
template void
elemental::lapack::Pinv
( DistMatrix<scomplex,MC,MR>& A, DistMatrix<scomplex,MC,MR>& PinvA );

template void
elemental::lapack::Pinv
( DistMatrix<dcomplex,MC,MR>& A, DistMatrix<dcomplex,MC,MR>& PinvA );
#endif

