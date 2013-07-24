/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_SVD_CHAN_HPP
#define ELEM_LAPACK_SVD_CHAN_HPP

#include "elemental/blas-like/level1/Adjoint.hpp"
#include "elemental/blas-like/level1/MakeTriangular.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/lapack-like/QR.hpp"

#include "elemental/lapack-like/SVD/GolubReinsch.hpp"

namespace elem {
namespace svd {

template<typename F>
inline void
ChanUpper
( DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& s, DistMatrix<F>& V,
  double heightRatio=1.5 )
{
#ifndef RELEASE
    CallStackEntry entry("svd::ChanUpper");
    if( A.Height() < A.Width() )
        throw std::logic_error("A must be at least as tall as it is wide");
    if( heightRatio <= 1.0 )
        throw std::logic_error("Nonsensical switchpoint for SVD");
#endif
    const Grid& g = A.Grid();
    const int m = A.Height();
    const int n = A.Width();
    if( m > heightRatio*n )
    {
        DistMatrix<F> R(g);
        qr::Explicit( A, R );
        svd::GolubReinschUpper( R, s, V );
        // Unfortunately, extra memory is used in forming A := A R,
        // where A has been overwritten with the Q from the QR factorization
        // of the original state of A, and R has been overwritten with the U 
        // from the SVD of the R from the QR factorization of A
        //
        // Perhaps this should be broken into pieces.
        DistMatrix<F> ACopy( A );
        Gemm( NORMAL, NORMAL, F(1), ACopy, R, F(0), A );
    }
    else
    {
        svd::GolubReinschUpper( A, s, V );
    }
}

template<typename F>
inline void
ChanUpper
( DistMatrix<F>& A,
  DistMatrix<BASE(F),VR,STAR>& s,
  double heightRatio=1.2 )
{
#ifndef RELEASE
    CallStackEntry entry("svd::ChanUpper");    
    if( heightRatio <= 1.0 )
        throw std::logic_error("Nonsensical switchpoint");
#endif
    const Grid& g = A.Grid();
    const int m = A.Height();
    const int n = A.Width();
    if( m >= heightRatio*n )
    {
        QR( A );
        DistMatrix<F> AT(g),
                      AB(g);
        PartitionDown
        ( A, AT,
             AB, n );
        MakeTriangular( UPPER, AT );
        GolubReinschUpper( AT, s );
    }
    else
    {
        GolubReinschUpper( A, s );
    }
}

//----------------------------------------------------------------------------//
// Grab the full SVD of the general matrix A, A = U diag(s) V^H using Chan's  //
// algorithm. On exit, A is overwritten with U.                               //
//----------------------------------------------------------------------------//

template<typename F>
inline void
Chan
( DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& s, DistMatrix<F>& V,
  double heightRatio=1.5 )
{
#ifndef RELEASE
    CallStackEntry entry("svd::Chan");
    if( heightRatio <= 1.0 )
        throw std::logic_error("Nonsensical switchpoint for SVD");
#endif
    // Check if we need to rescale the matrix, and do so if necessary
    bool needRescaling;
    BASE(F) scale;
    svd::CheckScale( A, needRescaling, scale );
    if( needRescaling )
        Scale( scale, A );

    // TODO: Switch between different algorithms. For instance, starting 
    //       with a QR decomposition of tall-skinny matrices.
    if( A.Height() >= A.Width() )
    {
        svd::ChanUpper( A, s, V, heightRatio );
    }
    else
    {
        // Lower bidiagonalization is not yet supported, so we instead play a 
        // trick to get the SVD of A.
        Adjoint( A, V );
        svd::ChanUpper( V, s, A, heightRatio );
    }

    // Rescale the singular values if necessary
    if( needRescaling )
        Scale( 1/scale, s );
}

//----------------------------------------------------------------------------//
// Grab the singular values of the general matrix A.                          //
//----------------------------------------------------------------------------//

template<typename F>
inline void
Chan( DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& s, double heightRatio=1.2 )
{
#ifndef RELEASE
    CallStackEntry entry("svd::Chan");
#endif
    // Check if we need to rescale the matrix, and do so if necessary
    typedef BASE(F) R;
    bool needRescaling;
    R scale;
    svd::CheckScale( A, needRescaling, scale );
    if( needRescaling )
        Scale( scale, A );

    // TODO: Switch between different algorithms. For instance, starting 
    //       with a QR decomposition of tall-skinny matrices.
    if( A.Height() >= A.Width() )
    {
        svd::ChanUpper( A, s, heightRatio );
    }
    else
    {
        // Lower bidiagonalization is not yet supported, so we instead play a 
        // trick to get the SVD of A.
        DistMatrix<F> AAdj( A.Grid() );
        Adjoint( A, AAdj );
        svd::ChanUpper( AAdj, s, heightRatio );
    }

    // Rescale the singular values if necessary
    if( needRescaling )
        Scale( 1/scale, s );
}

} // namespace svd
} // namespace elem

#endif // ifndef ELEM_LAPACK_SVD_CHAN_HPP
