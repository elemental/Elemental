/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_LDL_VAR3_HPP
#define LAPACK_LDL_VAR3_HPP

#include "elemental/blas-like/level1/DiagonalSolve.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"

namespace elem {
namespace ldl {

// Unblocked serial LDL _without_ partial pivoting
template<typename F> 
inline void
Var3Unb( Orientation orientation, Matrix<F>& A, Matrix<F>& d )
{
#ifndef RELEASE
    PushCallStack("ldl::Var3Unb");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( d.Viewing() && (d.Height() != A.Height() || d.Width() != 1) )
        throw std::logic_error
        ("d must be a column vector the same height as A");
    if( orientation == NORMAL )
        throw std::logic_error("Can only perform LDL^T or LDL^H");
#endif
    const int n = A.Height();
    if( !d.Viewing() )
        d.ResizeTo( n, 1 );

    F* ABuffer = A.Buffer();
    F* dBuffer = d.Buffer();
    const int ldim = A.LDim();
    for( int j=0; j<n; ++j )
    {
        const int a21Height = n - (j+1);

        // Extract and store the diagonal of D
        const F alpha11 = ABuffer[j+j*ldim];
        if( alpha11 == F(0) )
            throw SingularMatrixException();
        dBuffer[j] = alpha11; 

        F* RESTRICT a21 = &ABuffer[(j+1)+j*ldim];
        if( orientation == ADJOINT )
        {
            // A22 := A22 - a21 (a21 / alpha11)^H
            for( int k=0; k<a21Height; ++k )
            {
                const F beta = Conj(a21[k]/alpha11);
                F* RESTRICT A22Col = &ABuffer[(j+1)+(j+1+k)*ldim];
                for( int i=k; i<a21Height; ++i )
                    A22Col[i] -= a21[i]*beta;
            }
        }
        else
        {
            // A22 := A22 - a21 (a21 / alpha11)^T
            for( int k=0; k<a21Height; ++k )
            {
                const F beta = a21[k]/alpha11;
                F* RESTRICT A22Col = &ABuffer[(j+1)+(j+1+k)*ldim];
                for( int i=k; i<a21Height; ++i )
                    A22Col[i] -= a21[i]*beta;
            }
        }
        
        // a21 := a21 / alpha11
        for( int i=0; i<a21Height; ++i )
            a21[i] /= alpha11;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Blocked serial LDL _without_ partial pivoting
template<typename F>
inline void
Var3( Orientation orientation, Matrix<F>& A, Matrix<F>& d )
{
#ifndef RELEASE
    PushCallStack("ldl::Var3");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( d.Viewing() && (d.Height() != A.Height() || d.Width() != 1) )
        throw std::logic_error
        ("d must be a column vector the same height as A");
    if( orientation == NORMAL )
        throw std::logic_error("Can only perform LDL^T or LDL^H");
#endif
    const int n = A.Height();
    if( !d.Viewing() )
        d.ResizeTo( n, 1 );

    Matrix<F>
        ATL, ATR,  A00, A01, A02,
        ABL, ABR,  A10, A11, A12,
                   A20, A21, A22;
    Matrix<F> 
        dT,  d0,
        dB,  d1,
             d2;
    Matrix<F> S21;

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( d, dT,
         dB, 0 );
    while( ABR.Height() > 0 )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        RepartitionDown
        ( dT,  d0,
         /**/ /**/
               d1,
          dB,  d2 );

        //--------------------------------------------------------------------//
        ldl::Var3Unb( orientation, A11, d1 );
        Trsm( RIGHT, LOWER, orientation, UNIT, F(1), A11, A21 );
        S21 = A21;
        DiagonalSolve( RIGHT, NORMAL, d1, A21 );
        internal::TrrkNT( LOWER, orientation, F(-1), S21, A21, F(1), A22 );
        //--------------------------------------------------------------------//

        SlidePartitionDown
        ( dT,  d0,
               d1,
         /**/ /**/
          dB,  d2 );

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
Var3( Orientation orientation, DistMatrix<F>& A, DistMatrix<F,MC,STAR>& d )
{
#ifndef RELEASE
    PushCallStack("ldl::Var3");
    if( orientation == NORMAL )
        throw std::logic_error("Can only perform LDL^T and LDL^H");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Grid() != d.Grid() )
        throw std::logic_error("A and d must use the same grid");
    if( d.Viewing() && (d.Height() != A.Height() || d.Width() != 1) )
        throw std::logic_error
        ("d must be a column vector of the same height as A");
    if( d.Viewing() && d.ColAlignment() != A.ColAlignment() )
        throw std::logic_error("d must be aligned with A");
#endif
    const Grid& g = A.Grid();
    if( !d.Viewing() )
    {
        d.AlignWith( A );
        d.ResizeTo( A.Height(), 1 );
    }

    // Matrix views
    DistMatrix<F>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<F,MC,STAR>
        dT(g),  d0(g),
        dB(g),  d1(g),
                d2(g);

    // Temporary matrices
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> d1_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> A21_VR_STAR(g);
    DistMatrix<F,STAR,MC  > S21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A21AdjOrTrans_STAR_MR(g);

    const bool conjugate = ( orientation == ADJOINT );

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( d, dT,
         dB, 0 );
    while( ABR.Height() > 0 )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        RepartitionDown
        ( dT,  d0,
         /**/ /**/
               d1,
          dB,  d2 );

        A21_VC_STAR.AlignWith( A22 );
        A21_VR_STAR.AlignWith( A22 );
        S21Trans_STAR_MC.AlignWith( A22 );
        A21AdjOrTrans_STAR_MR.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_STAR_STAR = A11;
        LocalLDL( orientation, A11_STAR_STAR, d1_STAR_STAR );
        A11 = A11_STAR_STAR;
        d1 = d1_STAR_STAR;

        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, orientation, UNIT,
          F(1), A11_STAR_STAR, A21_VC_STAR );

        S21Trans_STAR_MC.TransposeFrom( A21_VC_STAR );
        DiagonalSolve( RIGHT, NORMAL, d1_STAR_STAR, A21_VC_STAR );
        A21_VR_STAR = A21_VC_STAR;
        A21AdjOrTrans_STAR_MR.TransposeFrom( A21_VR_STAR, conjugate );
        LocalTrrk
        ( LOWER, TRANSPOSE,
          F(-1), S21Trans_STAR_MC, A21AdjOrTrans_STAR_MR, F(1), A22 );

        A21 = A21_VC_STAR;
        //--------------------------------------------------------------------//
        A21_VC_STAR.FreeAlignments();
        A21_VR_STAR.FreeAlignments();
        S21Trans_STAR_MC.FreeAlignments();
        A21AdjOrTrans_STAR_MR.FreeAlignments();

        SlidePartitionDown
        ( dT,  d0,
               d1,
         /**/ /**/
          dB,  d2 );

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace ldl
} // namespace elem

#endif // ifndef LAPACK_LDL_VAR3_HPP
