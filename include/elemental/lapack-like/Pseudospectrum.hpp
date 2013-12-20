/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_PSEUDOSPECTRUM_HPP
#define ELEM_LAPACK_PSEUDOSPECTRUM_HPP

#include "elemental/lapack-like/Schur.hpp"

namespace elem {

namespace pspec {

template<typename F>
inline void
FrobNorms( const Matrix<F>& X, Matrix<BASE(F)>& norms )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::FrobNorms"))
    const Int m = X.Height();
    const Int n = X.Width();
    norms.ResizeTo( n, 1 );
    for( Int j=0; j<n; ++j )
        norms.Set( j, 0, blas::Nrm2(m,X.LockedBuffer(0,j),1) );
}

template<typename F>
inline void
FrobNorms( const DistMatrix<F>& X, DistMatrix<BASE(F),MR,STAR>& norms )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::FrobNorms"))
    const Int n = X.Width();
    const Int mLocal = X.LocalHeight();
    const Int nLocal = X.LocalWidth();
    const Grid& g = X.Grid();

    // TODO: Switch to more stable parallel norm computation using scaling
    norms.AlignWith( X );
    norms.ResizeTo( n, 1 ); 
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Base<F> localNorm = blas::Nrm2(mLocal,X.LockedBuffer(0,jLoc),1);
        norms.SetLocal( jLoc, 0, localNorm*localNorm );
    }

    mpi::AllReduce( norms.Buffer(), nLocal, mpi::SUM, X.ColComm() );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        norms.SetLocal( jLoc, 0, Sqrt(norms.GetLocal(jLoc,0)) );
}

template<typename F>
inline void
FixColumns( Matrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::FixColumns"))
    typedef Base<F> Real;
    Matrix<Real> norms;
    FrobNorms( X, norms );
    const Int m = X.Height();
    const Int n = X.Width();
    for( Int j=0; j<n; ++j )
    {
        auto x = View( X, 0, j, m, 1 );
        Real norm = norms.Get(j,0);
        if( norm == Real(0) )
        {
            MakeGaussian( x );
            norm = FrobeniusNorm( x );
        }
        Scale( Real(1)/norm, x );
    }
}

template<typename F>
inline void
FixColumns( DistMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::FixZeroColumns"))
    typedef Base<F> Real;
    DistMatrix<Real,MR,STAR> norms( X.Grid() );
    FrobNorms( X, norms );
    const Int m = X.Height();
    const Int nLocal = X.LocalWidth();
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Int j = X.RowShift() + jLoc*X.RowStride();
        auto x = View( X, 0, j, m, 1 );
        Real norm = norms.GetLocal(jLoc,0);
        if( norm == Real(0) )
        {
            MakeGaussian( x );
            norm = FrobeniusNorm( x );
        }
        Scale( Real(1)/norm, x );
    }
}

template<typename F>
inline void
ShiftedTrsmLUNUnb
( Matrix<F>& U, const Matrix<Complex<BASE(F)> >& shifts, Matrix<F>& X ) 
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::ShiftedTrsmLUNUnb");
        if( shifts.Height() != X.Width() )
            LogicError("Incompatible number of shifts");
    )
    auto diag = U.GetDiagonal();
    const Int n = U.Height();
    const Int ldim = U.LDim();
    const Int numShifts = shifts.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        UpdateDiagonal( U, -shifts.Get(j,0) );
        blas::Trsv
        ( 'U', 'N', 'N', n, U.LockedBuffer(), ldim, X.Buffer(0,j), 1 );
        U.SetDiagonal( diag );
    }
}

template<typename F>
inline void
ShiftedTrsmLUN
( Matrix<F>& U, const Matrix<Complex<BASE(F)> >& shifts, Matrix<F>& X ) 
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ShiftedTrsmLUN"))

    Matrix<F>
        UTL, UTR,  U00, U01, U02,
        UBL, UBR,  U10, U11, U12,
                   U20, U21, U22;
    Matrix<F> XT,  X0,
              XB,  X1,
                   X2;

    PartitionUpDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    PartitionUp
    ( X, XT,
         XB, 0 );
    while( XT.Height() > 0 )
    {
        RepartitionUpDiagonal
        ( UTL, /**/ UTR,   U00, U01, /**/ U02,
               /**/        U10, U11, /**/ U12,
         /*************/  /******************/
          UBL, /**/ UBR,   U20, U21, /**/ U22 );

        RepartitionUp
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );

        //--------------------------------------------------------------------//
        ShiftedTrsmLUNUnb( U11, shifts, X1 );
        Gemm( NORMAL, NORMAL, F(-1), U01, X1, F(1), X0 );
        //--------------------------------------------------------------------//

        SlidePartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        SlidePartitionUp
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );
    }
}

template<typename F>
inline void
ShiftedTrsmLUN
( const DistMatrix<F>& U, 
  const DistMatrix<Complex<BASE(F)>,VR,STAR>& shifts,
        DistMatrix<F>& X ) 
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ShiftedTrsmLUN"))
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<F>
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);
    DistMatrix<F> XT(g),  X0(g),
                  XB(g),  X1(g),
                          X2(g);

    // Temporary distributions
    DistMatrix<F,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    // Start the algorithm
    LockedPartitionUpDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    PartitionUp
    ( X, XT,
         XB, 0 );
    while( XT.Height() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( UTL, /**/ UTR,   U00, U01, /**/ U02,
               /**/        U10, U11, /**/ U12,
         /*************/  /******************/
          UBL, /**/ UBR,   U20, U21, /**/ U22 );

        RepartitionUp
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );

        U01_MC_STAR.AlignWith( X0 );
        X1_STAR_MR.AlignWith( X0 );
        X1_STAR_VR.AlignWith( shifts );
        //--------------------------------------------------------------------//
        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        X1_STAR_VR    = X1;  // X1[* ,VR] <- X1[MC,MR]

        // X1[* ,VR] := U11^-1[* ,* ] X1[* ,VR]
        ShiftedTrsmLUN
        ( U11_STAR_STAR.Matrix(), shifts.LockedMatrix(), X1_STAR_VR.Matrix() );

        X1_STAR_MR  = X1_STAR_VR; // X1[* ,MR]  <- X1[* ,VR]
        X1          = X1_STAR_MR; // X1[MC,MR] <- X1[* ,MR]
        U01_MC_STAR = U01;        // U01[MC,* ] <- U01[MC,MR]

        // X0[MC,MR] -= U01[MC,* ] X1[* ,MR]
        LocalGemm( NORMAL, NORMAL, F(-1), U01_MC_STAR, X1_STAR_MR, F(1), X0 );
        //--------------------------------------------------------------------//

        SlideLockedPartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        SlidePartitionUp
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );
    }
}

template<typename F>
inline void
ShiftedTrsmLUTUnb
( Matrix<F>& U, const Matrix<Complex<BASE(F)> >& shifts, Matrix<F>& X ) 
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::ShiftedTrsmLUTUnb");
        if( shifts.Height() != X.Width() )
            LogicError("Incompatible number of shifts");
    )
    auto diag = U.GetDiagonal();
    const Int n = U.Height();
    const Int ldim = U.LDim();
    const Int numShifts = shifts.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        UpdateDiagonal( U, -shifts.Get(j,0) );
        blas::Trsv
        ( 'U', 'C', 'N', n, U.LockedBuffer(), ldim, X.Buffer(0,j), 1 );
        U.SetDiagonal( diag );
    }
}

template<typename F>
inline void
ShiftedTrsmLUT
( Matrix<F>& U, const Matrix<Complex<BASE(F)> >& shifts, Matrix<F>& X ) 
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ShiftedTrsmLUT"))

    Matrix<F>
        UTL, UTR,  U00, U01, U02,
        UBL, UBR,  U10, U11, U12,
                   U20, U21, U22;

    Matrix<F> XT,  X0,
              XB,  X1,
                   X2;

    PartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XB.Height() > 0 )
    {
        RepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );

        //--------------------------------------------------------------------//
        ShiftedTrsmLUTUnb( U11, shifts, X1 );
        Gemm( ADJOINT, NORMAL, F(-1), U12, X1, F(1), X2 );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( UTL, /**/ UTR,   U00, U01, /**/ U02,
               /**/        U10, U11, /**/ U12,
         /*************/  /******************/
          UBL, /**/ UBR,   U20, U21, /**/ U22 );

        SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
}

template<typename F>
inline void
ShiftedTrsmLUT
( const DistMatrix<F>& U, const DistMatrix<Complex<BASE(F)>,VR,STAR>& shifts,
        DistMatrix<F>& X ) 
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ShiftedTrsmLUT"))
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<F>
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    DistMatrix<F> XT(g),  X0(g),
                  XB(g),  X1(g),
                          X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    // Start the algorithm
    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );

        U12_STAR_MC.AlignWith( X2 );
        X1_STAR_MR.AlignWith( X2 );
        X1_STAR_VR.AlignWith( shifts );
        //--------------------------------------------------------------------//
        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        X1_STAR_VR    = X1;  // X1[* ,VR] <- X1[MC,MR]

        // X1[* ,VR] := U11^-'[*,*] X1[* ,VR]
        ShiftedTrsmLUT
        ( U11_STAR_STAR.Matrix(), shifts.LockedMatrix(), X1_STAR_VR.Matrix() );

        X1_STAR_MR  = X1_STAR_VR; // X1[* ,MR]  <- X1[* ,VR]
        X1          = X1_STAR_MR; // X1[MC,MR]  <- X1[* ,MR]
        U12_STAR_MC = U12;        // U12[* ,MC] <- U12[MC,MR]

        // X2[MC,MR] -= (U12[* ,MC])' X1[* ,MR]
        //            = U12'[MC,*] X1[* ,MR]
        LocalGemm
        ( ADJOINT, NORMAL, F(-1), U12_STAR_MC, X1_STAR_MR, F(1), X2 );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,   U00, U01, /**/ U02,
               /**/        U10, U11, /**/ U12,
         /*************/  /******************/
          UBL, /**/ UBR,   U20, U21, /**/ U22 );

        SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
}

} // namespace pspec

template<typename F>
inline void
Pseudospectrum
( const Matrix<F>& A, const Matrix<Complex<BASE(F)> >& shifts, 
  Matrix<BASE(F)>& invNorms, Int maxIts=1000, BASE(F) tol=1e-6 )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;
    const Int n = A.Height();

    Matrix<C> T( n, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i )
            T.Set( i, j, A.Get(i,j) );

    Matrix<C> w;
    schur::QR( T, w );

    // Simultaneously run inverse iteration for various shifts
    const int numShifts = shifts.Height();
    Matrix<C> X;
    Gaussian( X, n, numShifts );
    Int numIts=0;
    Matrix<Real> estimates(numShifts,1),
                 lastEsts(numShifts,1), 
                 diffs(numShifts,1);
    Zeros( estimates, numShifts, 1 );
    while( true )
    {
        lastEsts = estimates;
        pspec::ShiftedTrsmLUN( T, shifts, X );
        pspec::FixColumns( X );
        pspec::ShiftedTrsmLUT( T, shifts, X );
        pspec::FrobNorms( X, estimates );

        ++numIts;
        if( numIts >= maxIts )
            break;

        diffs = estimates;
        Axpy( Real(-1), lastEsts, diffs );
        const Real maxDiff = MaxNorm( diffs );
        if( maxDiff <= tol*n )
            break;
    } 
    
    diffs = estimates;
    Axpy( Real(-1), lastEsts, diffs );
    const Real maxDiff = MaxNorm( diffs );
    if( maxDiff > tol*n )
        RuntimeError("Two-norm estimate did not converge in time");

    invNorms = estimates;
}

template<typename F>
inline void
Pseudospectrum
( const DistMatrix<F>& A, const DistMatrix<Complex<BASE(F)>,VR,STAR>& shifts,
  DistMatrix<BASE(F),VR,STAR>& invNorms, Int maxIts=1000, BASE(F) tol=1e-6 )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;
    const Int n = A.Height();
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    const Grid& g = A.Grid();

    DistMatrix<C> T( g );
    T.AlignWith( A );
    T.ResizeTo( n, n );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            T.SetLocal( iLoc, jLoc, A.GetLocal(iLoc,jLoc) );

    // We don't actually need the Schur vectors, but SDC requires their 
    // computation in order to form the full triangular factor
    DistMatrix<C> X( g );
    DistMatrix<C,VR,STAR> w( g );
    schur::SDC( T, w, X );

    // Simultaneously run inverse iteration for various shifts
    const int numShifts = shifts.Height();
    Gaussian( X, n, numShifts );
    Int numIts=0;
    DistMatrix<Real,MR,STAR> estimates(g), lastEsts(g), diffs(g);
    estimates.AlignWith( X );
    lastEsts.AlignWith( X );
    diffs.AlignWith( X );
    Zeros( estimates, numShifts, 1 );
    lastEsts.ResizeTo( numShifts, 1 );
    diffs.ResizeTo( numShifts, 1 );
    while( true )
    {
        lastEsts = estimates;
        pspec::ShiftedTrsmLUN( T, shifts, X );
        pspec::FixColumns( X );
        pspec::ShiftedTrsmLUT( T, shifts, X );
        pspec::FrobNorms( X, estimates );

        ++numIts;
        if( numIts >= maxIts )
            break;

        diffs = estimates;
        Axpy( Real(-1), lastEsts, diffs );
        const Real maxDiff = MaxNorm( diffs );
        if( maxDiff <= tol*n )
            break;
    } 
    
    diffs = estimates;
    Axpy( Real(-1), lastEsts, diffs );
    const Real maxDiff = MaxNorm( diffs );
    if( maxDiff > tol*n )
        RuntimeError("Two-norm estimate did not converge in time");

    invNorms = estimates;
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_PSEUDOSPECTRUM_HPP
