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
FrobNorms( const Matrix<Complex<BASE(F)> >& X, Matrix<BASE(F)>& norms )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::FrobNorms"));
    const Int m = X.Height();
    const Int n = X.Width();
    norms.ResizeTo( n, 1 );
    for( Int j=0; j<n; ++j )
        norms.Set( j, 0, Nrm2(m,X.LockedBuffer(0,j),1) );
}

template<typename F>
inline void
FrobNorms
( const DistMatrix<Complex<BASE(F)> >& X, DistMatrix<BASE(F),MR,STAR>& norms )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::FrobNorms"));
    const Int mLocal = X.LocalHeight();
    const Int nLocal = X.LocalWidth();
    const Grid& g = X.Grid();

    norms.AlignWith( X );
    norms.ResizeTo( nLocal, 1 ); 
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        norms.SetLocal( jLoc, 0, Nrms(mLocal,X.LockedBuffer(0,jLoc),1) );

    mpi::AllReduce( norms.Buffer(), nLocal, mpi::SUM, X.ColComm() );
}

template<typename F>
inline void
FixZeroColumns( Matrix<Complex<BASE(F)> >& X )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::FixZeroColumns"));
    Matrix<Complex<Base<F>>> norms;
    FrobNorms( X, norms );
    const Int m = X.Height();
    const Int n = X.Width();
    for( Int j=0; j<n; ++j )
    {
        if( norms.Get(j,0) == Base<F>(0) )
        {
            auto x = View( X, 0, j, m, 1 );
            MakeGaussian( x );
        }
    }
}

template<typename F>
inline void
FixZeroColumns( DistMatrix<Complex<BASE(F)> >& X )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::FixZeroColumns"));
    DistMatrix<Complex<Base<F>>,MR,STAR> norms( X.Grid() );
    FrobNorms( X, norms );
    const Int mLocal = X.LocalHeight();
    const Int nLocal = X.LocalWidth();
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        if( norms.Get(jLoc,0) == Base<F>(0) )
        {
            auto xLoc = View( X.Matrix(), 0, jLoc, mLocal, 1 );
            MakeGaussian( xLoc );
        }
    }
}

template<typename F>
inline void
ShiftedTrsmLUNUnb
( Matrix<F>& U, const Matrix<Complex<BASE(F)> >& shifts, Matrix<F>& X ) 
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ShiftedTrsmLUNUnb"));
    auto diag = U.GetDiagonal();
    const Int n = U.Height();
    const Int ldim = U.LDim();
    for( Int j=0; j<n; ++j )
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
    DEBUG_ONLY(CallStackEntry cse("pspec::ShiftedTrsmLUN"));

    Matrix<F>
        UTL, UTR,  U00, U01, U02,
        UBL, UBR,  U10, U11, U12,
                   U20, U21, U22;
    Matrix<F> XT,  X0,
              XB,  X1,
                   X2;

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

        //--------------------------------------------------------------------//
        ShiftedTrsmLUNUnb( U11, shifts, X1 );
        Gemm( NORMAL, NORMAL, F(-1), U01, X1, F(1), X0 );
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
ShiftedTrsmLUN
( const DistMatrix<F>& U, 
  const DistMatrix<Complex<BASE(F)>,VR,STAR>& shifts,
        DistMatrix<F>& X ) 
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ShiftedTrsmLUN"));
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
    DEBUG_ONLY(CallStackEntry cse("pspec::ShiftedTrsmLUTUnb"));
    auto diag = U.GetDiagonal();
    const Int n = U.Height();
    const Int ldim = U.LDim();
    for( Int j=0; j<n; ++j )
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
( const Matrix<F>& U, const Matrix<Complex<BASE(F)> >& shifts,
        Matrix<F>& X ) 
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ShiftedTrsmLUT"));

    Matrix<F>
        UTL, UTR,  U00, U01, U02,
        UBL, UBR,  U10, U11, U12,
                   U20, U21, U22;

    Matrix<F> XT,  X0,
              XB,  X1,
                   X2;

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

        //--------------------------------------------------------------------//
        ShiftedTrsmLUTUnb( U11, shifts, X1 );
        Gemm( ADJOINT, NORMAL, F(-1), U12, X1, F(1), X2 );
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

template<typename F>
inline void
ShiftedTrsmLUT
( const DistMatrix<F>& U, const DistMatrix<Complex<BASE(F)>,VR,STAR>& shifts,
        DistMatrix<F>& X ) 
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ShiftedTrsmLUT"));
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

        // X1[* ,VR] := U11^-[T/H][*,*] X1[* ,VR]
        ShiftedTrsmLUT
        ( U11_STAR_STAR.Matrix(), shifts.LockedMatrix(), X1_STAR_VR.Matrix() );

        X1_STAR_MR  = X1_STAR_VR; // X1[* ,MR]  <- X1[* ,VR]
        X1          = X1_STAR_MR; // X1[MC,MR]  <- X1[* ,MR]
        U12_STAR_MC = U12;        // U12[* ,MC] <- U12[MC,MR]

        // X2[MC,MR] -= (U12[* ,MC])^(T/H) X1[* ,MR]
        //            = U12^(T/H)[MC,*] X1[* ,MR]
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
( const Matrix<F>& A, Matrix<Complex<BASE(F)> >& shifts, 
  Int maxIts=1000, BASE(F) tol=1e-6 )
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
    while( true )
    {
        lastEsts = estimates;
        pspec::ShiftedTrsmLUN( T, shifts, X );
        pspec::FixZeroColumns( X );
        pspec::ShiftedTrsmLUT( T, shifts, X );

        ++numIts;
        if( numIts >= maxIts )
            break;
        
        pspec::FrobNorms( X, estimates );
        diffs = estimates;
        Axpy( Real(-1), lastEsts, diffs );
        const Real maxDiff = MaxNorm( diffs );
        if( maxDiff <= tol*n )
            break;
    } 
    
    Axpy( Real(-1), lastEsts, estimates );
    const Real maxDiff = MaxNorm( estimates );
    if( maxDiff > tol*n )
        RuntimeError("Two-norm estimate did not converge in time");
}

template<typename F>
inline void
Pseudospectrum
( const Matrix<F>& A, DistMatrix<Complex<BASE(F)>,VR,STAR>& shifts,
  Int maxIts=1000, BASE(F) tol=1e-6 )
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
    T.ResizeTo( A.Height(), A.Width() );
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
    estimates.ResizeTo( numShifts, 1 );
    lastEsts.ResizeTo( numShifts, 1 );
    diffs.ResizeTo( numShifts, 1 );
    while( true )
    {
        lastEsts = estimates;
        pspec::ShiftedTrsmLUN( T, shifts, X );
        pspec::FixZeroColumns( X );
        pspec::ShiftedTrsmLUT( T, shifts, X );

        ++numIts;
        if( numIts >= maxIts )
            break;
        
        pspec::FrobNorms( X, estimates );
        diffs = estimates;
        Axpy( Real(-1), lastEsts, diffs );
        const Real maxDiff = MaxNorm( diffs );
        if( maxDiff <= tol*n )
            break;
    } 
    
    Axpy( Real(-1), lastEsts, estimates );
    const Real maxDiff = MaxNorm( estimates );
    if( maxDiff > tol*n )
        RuntimeError("Two-norm estimate did not converge in time");
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_PSEUDOSPECTRUM_HPP
