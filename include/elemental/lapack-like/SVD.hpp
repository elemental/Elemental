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

#include "elemental/lapack-like/ApplyPackedReflectors.hpp"
#include "elemental/lapack-like/Bidiag.hpp"
#include "elemental/lapack-like/ExplicitQR.hpp"
#include "elemental/lapack-like/QR.hpp"

namespace elem {

namespace svd {

template<typename Real>
inline void
SimpleSVDUpper
( DistMatrix<Real>& A,
  DistMatrix<Real,VR,STAR>& s,
  DistMatrix<Real>& V )
{
#ifndef RELEASE
    PushCallStack("svd::SimpleSVDUpper");
#endif
    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min( m, n );
    const int offdiagonal = ( m>=n ? 1 : -1 );
    const char uplo = ( m>=n ? 'U' : 'L' );
    const Grid& g = A.Grid();

    // Bidiagonalize A
    Bidiag( A );

    // Grab copies of the diagonal and sub/super-diagonal of A
    DistMatrix<Real,MD,STAR> d_MD_STAR( g ), 
                             e_MD_STAR( g );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, offdiagonal );

    // In order to use serial QR kernels, we need the full bidiagonal matrix
    // on each process.
    //
    // NOTE: lapack::BidiagQRAlg expects e to be of length k
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d_MD_STAR );
    DistMatrix<Real,STAR,STAR> eHat_STAR_STAR( k, 1, g );
    DistMatrix<Real,STAR,STAR> e_STAR_STAR( g );
    View( e_STAR_STAR, eHat_STAR_STAR, 0, 0, k-1, 1 );
    e_STAR_STAR = e_MD_STAR;

    // Initialize U and VTrans to the appropriate identity matrices.
    DistMatrix<Real,VC,STAR> U_VC_STAR( g );
    DistMatrix<Real,STAR,VC> VTrans_STAR_VC( g );
    U_VC_STAR.AlignWith( A );
    VTrans_STAR_VC.AlignWith( V );
    Identity( m, k, U_VC_STAR );
    Identity( k, n, VTrans_STAR_VC );

    // Compute the SVD of the bidiagonal matrix and accumulate the Givens
    // rotations into our local portion of U and VTrans
    Matrix<Real>& ULocal = U_VC_STAR.LocalMatrix();
    Matrix<Real>& VTransLocal = VTrans_STAR_VC.LocalMatrix();
    lapack::BidiagQRAlg
    ( uplo, k, VTransLocal.Width(), ULocal.Height(),
      d_STAR_STAR.LocalBuffer(), e_STAR_STAR.LocalBuffer(), 
      VTransLocal.Buffer(), VTransLocal.LDim(), 
      ULocal.Buffer(), ULocal.LDim() );

    // Make a copy of A (for the Householder vectors) and pull the necessary 
    // portions of U and VTrans into a standard matrix dist.
    DistMatrix<Real> B( A );
    if( m >= n )
    {
        DistMatrix<Real> AT( g ),
                         AB( g );
        DistMatrix<Real,VC,STAR> UT_VC_STAR( g ), 
                                 UB_VC_STAR( g );
        PartitionDown( A, AT,
                          AB, n );
        PartitionDown( U_VC_STAR, UT_VC_STAR,
                                  UB_VC_STAR, n );
        AT = UT_VC_STAR;
        MakeZeros( AB );
        Transpose( VTrans_STAR_VC, V );
    }
    else
    {
        DistMatrix<Real> VT( g ), 
                         VB( g );
        DistMatrix<Real,STAR,VC> VTransL_STAR_VC( g ), VTransR_STAR_VC( g );
        PartitionDown( V, VT, 
                          VB, m );
        PartitionRight( VTrans_STAR_VC, VTransL_STAR_VC, VTransR_STAR_VC, m );
        Transpose( VTransL_STAR_VC, VT );
        MakeZeros( VB );
    }

    // Backtransform U and V
    if( m >= n )
    {
        ApplyPackedReflectors( LEFT, LOWER, VERTICAL, BACKWARD, 0, B, A );
        ApplyPackedReflectors( LEFT, UPPER, HORIZONTAL, BACKWARD, 1, B, V );
    }
    else
    {
        ApplyPackedReflectors( LEFT, LOWER, VERTICAL, BACKWARD, -1, B, A );
        ApplyPackedReflectors( LEFT, UPPER, HORIZONTAL, BACKWARD, 0, B, V );
    }

    // Copy out the appropriate subset of the singular values
    s = d_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Real>
inline void
SimpleSVDUpper
( DistMatrix<Complex<Real> >& A,
  DistMatrix<Real,VR,STAR>& s,
  DistMatrix<Complex<Real> >& V )
{
#ifndef RELEASE
    PushCallStack("svd::SimpleSVDUpper");
#endif
    typedef Complex<Real> C;
    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min( m, n );
    const int offdiagonal = ( m>=n ? 1 : -1 );
    const char uplo = ( m>=n ? 'U' : 'L' );
    const Grid& g = A.Grid();

    // Bidiagonalize A
    DistMatrix<C,STAR,STAR> tP( g ), tQ( g );
    Bidiag( A, tP, tQ );

    // Grab copies of the diagonal and sub/super-diagonal of A
    DistMatrix<Real,MD,STAR> d_MD_STAR( g ),
                             e_MD_STAR( g );
    A.GetRealPartOfDiagonal( d_MD_STAR );
    A.GetRealPartOfDiagonal( e_MD_STAR, offdiagonal );

    // NOTE: lapack::BidiagQRAlg expects e to be of length k
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d_MD_STAR );
    DistMatrix<Real,STAR,STAR> eHat_STAR_STAR( k, 1, g );
    DistMatrix<Real,STAR,STAR> e_STAR_STAR( g );
    View( e_STAR_STAR, eHat_STAR_STAR, 0, 0, k-1, 1 );
    e_STAR_STAR = e_MD_STAR;

    // Initialize U and VAdj to the appropriate identity matrices
    DistMatrix<C,VC,STAR> U_VC_STAR( g );
    DistMatrix<C,STAR,VC> VAdj_STAR_VC( g );
    U_VC_STAR.AlignWith( A );
    VAdj_STAR_VC.AlignWith( V );
    Identity( m, k, U_VC_STAR );
    Identity( k, n, VAdj_STAR_VC );

    // Compute the SVD of the bidiagonal matrix and accumulate the Givens
    // rotations into our local portion of U and VAdj
    Matrix<C>& ULocal = U_VC_STAR.LocalMatrix();
    Matrix<C>& VAdjLocal = VAdj_STAR_VC.LocalMatrix();
    lapack::BidiagQRAlg
    ( uplo, k, VAdjLocal.Width(), ULocal.Height(),
      d_STAR_STAR.LocalBuffer(), e_STAR_STAR.LocalBuffer(), 
      VAdjLocal.Buffer(), VAdjLocal.LDim(), 
      ULocal.Buffer(), ULocal.LDim() );

    // Make a copy of A (for the Householder vectors) and pull the necessary 
    // portions of U and VAdj into a standard matrix dist.
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
        Adjoint( VAdj_STAR_VC, V );
    }
    else
    {
        DistMatrix<C> VT( g ), 
                      VB( g );
        DistMatrix<C,STAR,VC> VAdjL_STAR_VC( g ), VAdjR_STAR_VC( g );
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
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifdef HAVE_FLA_BSVD
template<>
inline void
SimpleSVDUpper
( DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& s,
  DistMatrix<double>& V )
{
#ifndef RELEASE
    PushCallStack("svd::SimpleSVDUpper");
#endif
    typedef double Real;
    typedef Complex<Real> C;

    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min( m, n );
    const int offdiagonal = ( m>=n ? 1 : -1 );
    const char uplo = ( m>=n ? 'U' : 'L' );
    const Grid& g = A.Grid();

    // Bidiagonalize A
    Bidiag( A );

    // Grab copies of the diagonal and sub/super-diagonal of A
    DistMatrix<Real,MD,STAR> d_MD_STAR( g ), 
                             e_MD_STAR( g );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, offdiagonal );

    // In order to use serial QR kernels, we need the full bidiagonal matrix
    // on each process
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                               e_STAR_STAR( e_MD_STAR );

    // Initialize U and V to the appropriate identity matrices.
    DistMatrix<Real,VC,STAR> U_VC_STAR( g );
    DistMatrix<Real,VC,STAR> V_VC_STAR( g );
    U_VC_STAR.AlignWith( A );
    V_VC_STAR.AlignWith( V );
    Identity( m, k, U_VC_STAR );
    Identity( n, k, V_VC_STAR );

    // Compute the SVD of the bidiagonal matrix and accumulate the Givens
    // rotations into our local portion of U and V
    // NOTE: This _only_ works in the case where m >= n
    const int numAccum = 32;
    const int maxNumIts = 30;
    const int bAlg = 512;
    std::vector<C> GBuffer( (k-1)*numAccum ),
                   HBuffer( (k-1)*numAccum );
    FLA_Bsvd_v_opd_var1
    ( k, U_VC_STAR.LocalHeight(), V_VC_STAR.LocalHeight(), 
      numAccum, maxNumIts,
      d_STAR_STAR.LocalBuffer(), 1,
      e_STAR_STAR.LocalBuffer(), 1,
      &GBuffer[0], 1, k-1,
      &HBuffer[0], 1, k-1,
      U_VC_STAR.LocalBuffer(), 1, U_VC_STAR.LocalLDim(),
      V_VC_STAR.LocalBuffer(), 1, V_VC_STAR.LocalLDim(),
      bAlg );

    // Make a copy of A (for the Householder vectors) and pull the necessary 
    // portions of U and V into a standard matrix dist.
    DistMatrix<Real> B( A );
    if( m >= n )
    {
        DistMatrix<Real> AT( g ),
                         AB( g );
        DistMatrix<Real,VC,STAR> UT_VC_STAR( g ), 
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
        DistMatrix<Real> VT( g ), 
                         VB( g );
        DistMatrix<Real,VC,STAR> VT_VC_STAR( g ), 
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
        ApplyPackedReflectors( LEFT, LOWER, VERTICAL, BACKWARD, 0, B, A );
        ApplyPackedReflectors( LEFT, UPPER, HORIZONTAL, BACKWARD, 1, B, V );
    }
    else
    {
        ApplyPackedReflectors( LEFT, LOWER, VERTICAL, BACKWARD, -1, B, A );
        ApplyPackedReflectors( LEFT, UPPER, HORIZONTAL, BACKWARD, 0, B, V );
    }

    // Copy out the appropriate subset of the singular values
    s = d_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<>
inline void
SimpleSVDUpper
( DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& s,
  DistMatrix<Complex<double> >& V )
{
#ifndef RELEASE
    PushCallStack("svd::SimpleSVDUpper");
#endif
    typedef double Real;
    typedef Complex<Real> C;

    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min( m, n );
    const int offdiagonal = ( m>=n ? 1 : -1 );
    const char uplo = ( m>=n ? 'U' : 'L' );
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
    Identity( m, k, U_VC_STAR );
    Identity( n, k, V_VC_STAR );

    // Compute the SVD of the bidiagonal matrix and accumulate the Givens
    // rotations into our local portion of U and V
    // NOTE: This _only_ works in the case where m >= n
    const int numAccum = 32;
    const int maxNumIts = 30;
    const int bAlg = 512;
    std::vector<C> GBuffer( (k-1)*numAccum ),
                   HBuffer( (k-1)*numAccum );
    FLA_Bsvd_v_opz_var1
    ( k, U_VC_STAR.LocalHeight(), V_VC_STAR.LocalHeight(), 
      numAccum, maxNumIts,
      d_STAR_STAR.LocalBuffer(), 1,
      e_STAR_STAR.LocalBuffer(), 1,
      &GBuffer[0], 1, k-1,
      &HBuffer[0], 1, k-1,
      U_VC_STAR.LocalBuffer(), 1, U_VC_STAR.LocalLDim(),
      V_VC_STAR.LocalBuffer(), 1, V_VC_STAR.LocalLDim(),
      bAlg );

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
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // HAVE_FLA_BSVD

template<typename F>
inline void
SVDUpper
( DistMatrix<F>& A,
  DistMatrix<typename Base<F>::type,VR,STAR>& s,
  DistMatrix<F>& V,
  double heightRatio=1.5 )
{
#ifndef RELEASE
    PushCallStack("svd::SVDUpper");
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
        ExplicitQR( A, R );
        svd::SimpleSVDUpper( R, s, V );
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
        svd::SimpleSVDUpper( A, s, V );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Real>
inline void
SimpleSingularValuesUpper
( DistMatrix<Real>& A,
  DistMatrix<Real,VR,STAR>& s )
{
#ifndef RELEASE
    PushCallStack("svd::SimpleSingularValuesUpper");
#endif
    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min( m, n );
    const int offdiagonal = ( m>=n ? 1 : -1 );
    const char uplo = ( m>=n ? 'U' : 'L' );
    const Grid& g = A.Grid();

    // Bidiagonalize A
    Bidiag( A );

    // Grab copies of the diagonal and sub/super-diagonal of A
    DistMatrix<Real,MD,STAR> d_MD_STAR( g ), 
                             e_MD_STAR( g );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, offdiagonal );


    // In order to use serial QR kernels, we need the full bidiagonal matrix
    // on each process
    //
    // NOTE: lapack::BidiagQRAlg expects e to be of length k
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d_MD_STAR );
    DistMatrix<Real,STAR,STAR> eHat_STAR_STAR( k, 1, g );
    DistMatrix<Real,STAR,STAR> e_STAR_STAR( g );
    View( e_STAR_STAR, eHat_STAR_STAR, 0, 0, k-1, 1 );
    e_STAR_STAR = e_MD_STAR;

    // Compute the singular values of the bidiagonal matrix
    lapack::BidiagQRAlg
    ( uplo, k, 0, 0,
      d_STAR_STAR.LocalBuffer(), e_STAR_STAR.LocalBuffer(), 
      (Real*)0, 1, (Real*)0, 1 );

    // Copy out the appropriate subset of the singular values
    s = d_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Real>
inline void
SimpleSingularValuesUpper
( DistMatrix<Complex<Real> >& A,
  DistMatrix<Real,VR,STAR>& s )
{
#ifndef RELEASE
    PushCallStack("svd::SimpleSingularValuesUpper");
#endif
    typedef Complex<Real> C;
    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min( m, n );
    const int offdiagonal = ( m>=n ? 1 : -1 );
    const char uplo = ( m>=n ? 'U' : 'L' );
    const Grid& g = A.Grid();

    // Bidiagonalize A
    DistMatrix<C,STAR,STAR> tP(g), tQ(g);
    Bidiag( A, tP, tQ );

    // Grab copies of the diagonal and sub/super-diagonal of A
    DistMatrix<Real,MD,STAR> d_MD_STAR( g ), 
                             e_MD_STAR( g );
    A.GetRealPartOfDiagonal( d_MD_STAR );
    A.GetRealPartOfDiagonal( e_MD_STAR, offdiagonal );

    // In order to use serial QR kernels, we need the full bidiagonal matrix
    // on each process
    //
    // NOTE: lapack::BidiagQRAlg expects e to be of length k
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d_MD_STAR );
    DistMatrix<Real,STAR,STAR> eHat_STAR_STAR( k, 1, g );
    DistMatrix<Real,STAR,STAR> e_STAR_STAR( g );
    View( e_STAR_STAR, eHat_STAR_STAR, 0, 0, k-1, 1 );
    e_STAR_STAR = e_MD_STAR;

    // Compute the singular values of the bidiagonal matrix
    lapack::BidiagQRAlg
    ( uplo, k, 0, 0,
      d_STAR_STAR.LocalBuffer(), e_STAR_STAR.LocalBuffer(), 
      (C*)0, 1, (C*)0, 1 );

    // Copy out the appropriate subset of the singular values
    s = d_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Real>
inline void
SingularValuesUpper
( DistMatrix<Real>& A,
  DistMatrix<Real,VR,STAR>& s,
  double heightRatio=1.2 )
{
#ifndef RELEASE
    PushCallStack("svd::SingularValuesUpper");    
    if( heightRatio <= 1.0 )
        throw std::logic_error("Nonsensical switchpoint for SingularValues");
#endif
    const Grid& g = A.Grid();
    const int m = A.Height();
    const int n = A.Width();
    if( m >= heightRatio*n )
    {
        QR( A );
        DistMatrix<Real> AT(g),
                         AB(g);
        PartitionDown
        ( A, AT,
             AB, n );
        MakeTrapezoidal( LEFT, UPPER, 0, AT );
        SimpleSingularValuesUpper( AT, s );
    }
    else
    {
        SimpleSingularValuesUpper( A, s );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Real>
inline void
SingularValuesUpper
( DistMatrix<Complex<Real> >& A,
  DistMatrix<Real,VR,STAR>& s,
  double heightRatio=1.2 )
{
#ifndef RELEASE
    PushCallStack("svd::SingularValuesUpper");
    if( heightRatio <= 1.0 )
        throw std::logic_error("Nonsensical switchpoint for SingularValues");
#endif
    typedef Complex<Real> C;
    const Grid& g = A.Grid();
    const int m = A.Height();
    const int n = A.Width();
    if( m >= heightRatio*n )
    {
        DistMatrix<C,MD,STAR> t(g);
        QR( A, t );
        DistMatrix<C> AT(g),
                      AB(g);
        PartitionDown
        ( A, AT,
             AB, n );
        MakeTrapezoidal( LEFT, UPPER, 0, AT );
        SimpleSingularValuesUpper( AT, s );
    }
    else
    {
        SimpleSingularValuesUpper( A, s );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
CheckScale
( DistMatrix<F>& A, bool& needRescaling, typename Base<F>::type& scale )
{
    typedef typename Base<F>::type R;

    scale = 1;
    needRescaling = false;
    const R oneNormOfA = Norm( A, ONE_NORM );
    const R safeMin = lapack::MachineSafeMin<R>();
    const R precision = lapack::MachinePrecision<R>();
    const R smallNumber = safeMin/precision;
    const R bigNumber = 1/smallNumber;
    const R rhoMin = Sqrt(smallNumber);
    const R rhoMax = std::min( Sqrt(bigNumber), 1/Sqrt(Sqrt(safeMin)) );

    if( oneNormOfA > 0 && oneNormOfA < rhoMin )
    {
        needRescaling = true;
        scale = rhoMin/oneNormOfA;
    }
    else if( oneNormOfA > rhoMax )
    {
        needRescaling = true;
        scale = rhoMax/oneNormOfA;
    }
}

template<typename F>
inline void
DivideAndConquerSVD
( Matrix<F>& A, Matrix<typename Base<F>::type>& s, Matrix<F>& V )
{
#ifndef RELEASE
    PushCallStack("svd::DivideAndConquerSVD");
#endif
    typedef typename Base<F>::type R;

    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min(m,n);
    s.ResizeTo( k, 1 );
    Matrix<F> U( m, k );
    Matrix<F> VAdj( k, n );
    lapack::DivideAndConquerSVD
    ( m, n, A.Buffer(), A.LDim(), s.Buffer(), U.Buffer(), U.LDim(),
      VAdj.Buffer(), VAdj.LDim() );

    A = U;
    Adjoint( VAdj, V );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
QRSVD( Matrix<F>& A, Matrix<typename Base<F>::type>& s, Matrix<F>& V )
{
#ifndef RELEASE
    PushCallStack("svd::QRSVD");
#endif
    typedef typename Base<F>::type R;

    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min(m,n);
    s.ResizeTo( k, 1 );
    Matrix<F> U( m, k );
    Matrix<F> VAdj( k, n );
    lapack::QRSVD
    ( m, n, A.Buffer(), A.LDim(), s.Buffer(), U.Buffer(), U.LDim(),
      VAdj.Buffer(), VAdj.LDim() );

    A = U;
    Adjoint( VAdj, V );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace svd

//----------------------------------------------------------------------------//
// Grab the full SVD of the general matrix A, A = U diag(s) V^H.              //
// On exit, A is overwritten with U.                                          //
//----------------------------------------------------------------------------//

template<typename F>
inline void
SVD( Matrix<F>& A, Matrix<typename Base<F>::type>& s, Matrix<F>& V, bool useQR )
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
inline void
SVD
( DistMatrix<F>& A,
  DistMatrix<typename Base<F>::type,VR,STAR>& s,
  DistMatrix<F>& V,
  double heightRatio )
{
#ifndef RELEASE
    PushCallStack("SVD");
    if( heightRatio <= 1.0 )
        throw std::logic_error("Nonsensical switchpoint for SVD");
#endif
    typedef typename Base<F>::type Real;

    // Check if we need to rescale the matrix, and do so if necessary
    bool needRescaling;
    Real scale;
    svd::CheckScale( A, needRescaling, scale );
    if( needRescaling )
        Scale( scale, A );

    // TODO: Switch between different algorithms. For instance, starting 
    //       with a QR decomposition of tall-skinny matrices.
    if( A.Height() >= A.Width() )
    {
        svd::SVDUpper( A, s, V, heightRatio );
    }
    else
    {
        // Lower bidiagonalization is not yet supported, so we instead play a 
        // trick to get the SVD of A.
        Adjoint( A, V );
        svd::SVDUpper( V, s, A, heightRatio );
    }

    // Rescale the singular values if necessary
    if( needRescaling )
        Scal( 1/scale, s );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Grab the singular values of the general matrix A.                          //
//----------------------------------------------------------------------------//
template<typename F>
inline void
SingularValues( Matrix<F>& A, Matrix<typename Base<F>::type>& s )
{
#ifndef RELEASE
    PushCallStack("SingularValues");
#endif
    typedef typename Base<F>::type R;

    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min(m,n);
    s.ResizeTo( k, 1 );
    lapack::SingularValues( m, n, A.Buffer(), A.LDim(), s.Buffer() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
SingularValues
( DistMatrix<F>& A,
  DistMatrix<typename Base<F>::type,VR,STAR>& s,
  double heightRatio )
{
#ifndef RELEASE
    PushCallStack("SingularValues");
#endif
    typedef typename Base<F>::type R;

    // Check if we need to rescale the matrix, and do so if necessary
    bool needRescaling;
    R scale;
    svd::CheckScale( A, needRescaling, scale );
    if( needRescaling )
        Scale( scale, A );

    // TODO: Switch between different algorithms. For instance, starting 
    //       with a QR decomposition of tall-skinny matrices.
    if( A.Height() >= A.Width() )
    {
        svd::SingularValuesUpper( A, s, heightRatio );
    }
    else
    {
        // Lower bidiagonalization is not yet supported, so we instead play a 
        // trick to get the SVD of A.
        DistMatrix<F> AAdj( A.Grid() );
        Adjoint( A, AAdj );
        svd::SingularValuesUpper( AAdj, s, heightRatio );
    }

    // Rescale the singular values if necessary
    if( needRescaling )
        Scal( 1/scale, s );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef LAPACK_SVD_HPP
