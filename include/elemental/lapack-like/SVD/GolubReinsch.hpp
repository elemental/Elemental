/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_SVD_GOLUBREINSCH_HPP
#define LAPACK_SVD_GOLUBREINSCH_HPP

#include "elemental/blas-like/level1/Adjoint.hpp"
#include "elemental/blas-like/level1/Transpose.hpp"
#include "elemental/lapack-like/ApplyPackedReflectors.hpp"
#include "elemental/lapack-like/Bidiag.hpp"
#include "elemental/matrices/Identity.hpp"
#include "elemental/matrices/Zeros.hpp"

#include "elemental/lapack-like/SVD/Util.hpp"

namespace elem {
namespace svd {

// TODO: Wrap below routines to produce svd::GolubReinsch?

template<typename F>
inline void
GolubReinschUpper
( DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& s, DistMatrix<F>& V )
{
#ifndef RELEASE
    CallStackEntry entry("svd::GolubReinschUpper");
#endif
    typedef BASE(F) Real;
    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min( m, n );
    const int offdiagonal = ( m>=n ? 1 : -1 );
    const char uplo = ( m>=n ? 'U' : 'L' );
    const Grid& g = A.Grid();

    // Bidiagonalize A
    DistMatrix<F,STAR,STAR> tP( g ), tQ( g );
    Bidiag( A, tP, tQ );

    // Grab copies of the diagonal and sub/super-diagonal of A
    DistMatrix<Real,MD,STAR> d_MD_STAR( g ),
                             e_MD_STAR( g );
    A.GetRealPartOfDiagonal( d_MD_STAR );
    A.GetRealPartOfDiagonal( e_MD_STAR, offdiagonal );

    // NOTE: lapack::BidiagQRAlg expects e to be of length k
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                               eHat_STAR_STAR( k, 1, g ),
                               e_STAR_STAR( g );
    View( e_STAR_STAR, eHat_STAR_STAR, 0, 0, k-1, 1 );
    e_STAR_STAR = e_MD_STAR;

    // Initialize U and VAdj to the appropriate identity matrices
    DistMatrix<F,VC,STAR> U_VC_STAR( g );
    DistMatrix<F,STAR,VC> VAdj_STAR_VC( g );
    U_VC_STAR.AlignWith( A );
    VAdj_STAR_VC.AlignWith( V );
    Identity( U_VC_STAR, m, k );
    Identity( VAdj_STAR_VC, k, n );

    // Compute the SVD of the bidiagonal matrix and accumulate the Givens
    // rotations into our local portion of U and VAdj
    Matrix<F>& ULoc = U_VC_STAR.Matrix();
    Matrix<F>& VAdjLoc = VAdj_STAR_VC.Matrix();
    lapack::BidiagQRAlg
    ( uplo, k, VAdjLoc.Width(), ULoc.Height(),
      d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(), 
      VAdjLoc.Buffer(), VAdjLoc.LDim(), 
      ULoc.Buffer(), ULoc.LDim() );

    // Make a copy of A (for the Householder vectors) and pull the necessary 
    // portions of U and VAdj into a standard matrix dist.
    DistMatrix<F> B( A );
    if( m >= n )
    {
        DistMatrix<F> AT( g ),
                      AB( g );
        DistMatrix<F,VC,STAR> UT_VC_STAR( g ),
                              UB_VC_STAR( g );
        PartitionDown( A, AT,
                          AB, n );
        PartitionDown( U_VC_STAR, UT_VC_STAR,
                                  UB_VC_STAR, n );
        AT = UT_VC_STAR;
        MakeZeros( AB );
        Adjoint( VAdj_STAR_VC, V );
    }
    else
    {
        DistMatrix<F> VT( g ), 
                      VB( g );
        DistMatrix<F,STAR,VC> VAdjL_STAR_VC( g ), VAdjR_STAR_VC( g );
        PartitionDown( V, VT, 
                          VB, m );
        PartitionRight( VAdj_STAR_VC, VAdjL_STAR_VC, VAdjR_STAR_VC, m );
        Adjoint( VAdjL_STAR_VC, VT );
        MakeZeros( VB );
    }

    // Backtransform U and V
    if( m >= n )
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, 0, B, tQ, A );
        ApplyPackedReflectors
        ( LEFT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 1, B, tP, V );
    }
    else
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, -1, B, tQ, A );
        ApplyPackedReflectors
        ( LEFT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 0, B, tP, V );
    }

    // Copy out the appropriate subset of the singular values
    s = d_STAR_STAR;
}

#ifdef HAVE_FLA_BSVD
template<typename F>
inline void
GolubReinschUpper_FLA
( DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& s, DistMatrix<F>& V )
{
#ifndef RELEASE
    CallStackEntry entry("svd::GolubReinschUpper_FLA");
#endif
    typedef double Real;
    typedef Complex<Real> C;

    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min( m, n );
    const int offdiagonal = ( m>=n ? 1 : -1 );
    const Grid& g = A.Grid();

    // Bidiagonalize A
    DistMatrix<C,STAR,STAR> tP( g ), tQ( g );
    Bidiag( A, tP, tQ );

    // Grab copies of the diagonal and sub/super-diagonal of A
    DistMatrix<Real,MD,STAR> d_MD_STAR( g ),
                             e_MD_STAR( g );
    A.GetRealPartOfDiagonal( d_MD_STAR );
    A.GetRealPartOfDiagonal( e_MD_STAR, offdiagonal );

    // In order to use serial QR kernels, we need the full bidiagonal matrix
    // on each process
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                               e_STAR_STAR( e_MD_STAR );

    // Initialize U and VAdj to the appropriate identity matrices
    DistMatrix<C,VC,STAR> U_VC_STAR( g );
    DistMatrix<C,VC,STAR> V_VC_STAR( g );
    U_VC_STAR.AlignWith( A );
    V_VC_STAR.AlignWith( V );
    Identity( U_VC_STAR, m, k );
    Identity( V_VC_STAR, n, k );

    FlaSVD
    ( k, U_VC_STAR.LocalHeight(), V_VC_STAR.LocalHeight(),
      d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
      U_VC_STAR.Buffer(), U_VC_STAR.LDim(),
      V_VC_STAR.Buffer(), V_VC_STAR.LDim() );

    // Make a copy of A (for the Householder vectors) and pull the necessary 
    // portions of U and V into a standard matrix dist.
    DistMatrix<C> B( A );
    if( m >= n )
    {
        DistMatrix<C> AT( g ),
                      AB( g );
        DistMatrix<C,VC,STAR> UT_VC_STAR( g ), 
                              UB_VC_STAR( g );
        PartitionDown( A, AT,
                          AB, n );
        PartitionDown( U_VC_STAR, UT_VC_STAR,
                                  UB_VC_STAR, n );
        AT = UT_VC_STAR;
        MakeZeros( AB );
        V = V_VC_STAR;
    }
    else
    {
        DistMatrix<C> VT( g ), 
                      VB( g );
        DistMatrix<C,VC,STAR> VT_VC_STAR( g ), 
                              VB_VC_STAR( g );
        PartitionDown( V, VT, 
                          VB, m );
        PartitionDown
        ( V_VC_STAR, VT_VC_STAR, 
                     VB_VC_STAR, m );
        VT = VT_VC_STAR;
        MakeZeros( VB );
    }

    // Backtransform U and V
    if( m >= n )
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, 0, B, tQ, A );
        ApplyPackedReflectors
        ( LEFT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 1, B, tP, V );
    }
    else
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, -1, B, tQ, A );
        ApplyPackedReflectors
        ( LEFT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 0, B, tP, V );
    }

    // Copy out the appropriate subset of the singular values
    s = d_STAR_STAR;
}

template<>
inline void
GolubReinschUpper
( DistMatrix<double>& A, DistMatrix<double,VR,STAR>& s, DistMatrix<double>& V )
{
#ifndef RELEASE
    CallStackEntry entry("svd::GolubReinschUpper");
#endif
    GolubReinschUpper_FLA( A, s, V );
}

template<>
inline void
GolubReinschUpper
( DistMatrix<Complex<double> >& A, 
  DistMatrix<double,VR,STAR>& s, DistMatrix<Complex<double> >& V )
{
#ifndef RELEASE
    CallStackEntry entry("svd::GolubReinschUpper");
#endif
    GolubReinschUpper_FLA( A, s, V );
}
#endif // HAVE_FLA_BSVD

template<typename F>
inline void
GolubReinschUpper
( DistMatrix<F>& A,
  DistMatrix<BASE(F),VR,STAR>& s )
{
#ifndef RELEASE
    CallStackEntry entry("svd::GolubReinschUpper");
#endif
    typedef BASE(F) Real;
    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min( m, n );
    const int offdiagonal = ( m>=n ? 1 : -1 );
    const Grid& g = A.Grid();

    // Bidiagonalize A
    DistMatrix<F,STAR,STAR> tP(g), tQ(g);
    Bidiag( A, tP, tQ );

    // Grab copies of the diagonal and sub/super-diagonal of A
    DistMatrix<Real,MD,STAR> d_MD_STAR( g ), 
                             e_MD_STAR( g );
    A.GetRealPartOfDiagonal( d_MD_STAR );
    A.GetRealPartOfDiagonal( e_MD_STAR, offdiagonal );

    // In order to use serial DQDS kernels, we need the full bidiagonal matrix
    // on each process
    //
    // NOTE: lapack::BidiagDQDS expects e to be of length k
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                               eHat_STAR_STAR( k, 1, g ),
                               e_STAR_STAR( g );
    View( e_STAR_STAR, eHat_STAR_STAR, 0, 0, k-1, 1 );
    e_STAR_STAR = e_MD_STAR;

    // Compute the singular values of the bidiagonal matrix via DQDS
    lapack::BidiagDQDS( k, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer() );

    // Copy out the appropriate subset of the singular values
    s = d_STAR_STAR;
}

} // namespace svd
} // namespace elem

#endif // ifndef LAPACK_SVD_GOLUBREINSCH_HPP
