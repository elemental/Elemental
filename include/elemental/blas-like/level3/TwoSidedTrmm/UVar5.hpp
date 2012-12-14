/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

// The only requirement that this is a field comes from the necessity for 
// the existence of 1/2, which is artifact of the algorithm...
template<typename F> 
inline void
TwoSidedTrmmUVar5( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& U )
{
#ifndef RELEASE
    PushCallStack("internal::TwoSidedTrmmUVar5");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( U.Height() != U.Width() )
        throw std::logic_error("Triangular matrices must be square");
    if( A.Height() != U.Height() )
        throw std::logic_error("A and U must be the same size");
#endif
    // Matrix views
    Matrix<F>
        ATL, ATR,  A00, A01, A02,
        ABL, ABR,  A10, A11, A12,
                   A20, A21, A22;
    Matrix<F>
        UTL, UTR,  U00, U01, U02,
        UBL, UBR,  U10, U11, U12,
                   U20, U21, U22;

    // Temporary products
    Matrix<F> Y01;

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        //--------------------------------------------------------------------//
        // Y01 := U01 A11
        Zeros( A01.Height(), A01.Width(), Y01 );
        Hemm( RIGHT, UPPER, F(1), A11, U01, F(0), Y01 );

        // A01 := U00 A01
        Trmm( LEFT, UPPER, NORMAL, diag, F(1), U00, A01 );

        // A01 := A01 + 1/2 Y01
        Axpy( F(1)/F(2), Y01, A01 );

        // A00 := A00 + (U01 A01' + A01 U01')
        Her2k( UPPER, NORMAL, F(1), U01, A01, F(1), A00 );

        // A01 := A01 + 1/2 Y01
        Axpy( F(1)/F(2), Y01, A01 );

        // A01 := A01 U11'
        Trmm( RIGHT, UPPER, ADJOINT, diag, F(1), U11, A01 );

        // A11 := U11 A11 U11'
        TwoSidedTrmmUUnb( diag, A11, U11 );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
TwoSidedTrmmUVar5
( UnitOrNonUnit diag, DistMatrix<F>& A, const DistMatrix<F>& U )
{
#ifndef RELEASE
    PushCallStack("internal::TwoSidedTrmmUVar5");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( U.Height() != U.Width() )
        throw std::logic_error("Triangular matrices must be square");
    if( A.Height() != U.Height() )
        throw std::logic_error("A and U must be the same size");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<F>
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> A01_MC_STAR(g);
    DistMatrix<F,MR,  STAR> A01_MR_STAR(g);
    DistMatrix<F,VC,  STAR> A01_VC_STAR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<F,MR,  STAR> U01_MR_STAR(g);
    DistMatrix<F,VC,  STAR> U01_VC_STAR(g);
    DistMatrix<F,VC,  STAR> Y01_VC_STAR(g);
    DistMatrix<F> Y01(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        A01_MC_STAR.AlignWith( A00 );
        A01_MR_STAR.AlignWith( A00 );
        A01_VC_STAR.AlignWith( A00 );
        U01_MC_STAR.AlignWith( A00 );
        U01_MR_STAR.AlignWith( A00 );
        U01_VC_STAR.AlignWith( A00 );
        Y01.AlignWith( A01 );
        Y01_VC_STAR.AlignWith( A01 );
        //--------------------------------------------------------------------//
        // Y01 := U01 A11
        A11_STAR_STAR = A11;
        U01_VC_STAR = U01;
        Y01_VC_STAR.ResizeTo( A01.Height(), A01.Width() );
        Hemm
        ( RIGHT, UPPER,
          F(1), A11_STAR_STAR.LocalMatrix(), U01_VC_STAR.LocalMatrix(),
          F(0), Y01_VC_STAR.LocalMatrix() );
        Y01 = Y01_VC_STAR;

        // A01 := U00 A01
        Trmm( LEFT, UPPER, NORMAL, diag, F(1), U00, A01 );

        // A01 := A01 + 1/2 Y01
        Axpy( F(1)/F(2), Y01, A01 );

        // A00 := A00 + (U01 A01' + A01 U01')
        A01_MC_STAR = A01;
        U01_MC_STAR = U01;
        A01_VC_STAR = A01_MC_STAR;
        A01_MR_STAR = A01_VC_STAR;
        U01_MR_STAR = U01_MC_STAR;
        LocalTrr2k
        ( UPPER, ADJOINT, ADJOINT,
          F(1), U01_MC_STAR, A01_MR_STAR, 
                A01_MC_STAR, U01_MR_STAR,
          F(1), A00 );

        // A01 := A01 + 1/2 Y01
        Axpy( F(1)/F(2), Y01_VC_STAR, A01_VC_STAR );

        // A01 := A01 U11'
        U11_STAR_STAR = U11;
        LocalTrmm
        ( RIGHT, UPPER, ADJOINT, diag, F(1), U11_STAR_STAR, A01_VC_STAR );
        A01 = A01_VC_STAR;

        // A11 := U11 A11 U11'
        LocalTwoSidedTrmm( UPPER, diag, A11_STAR_STAR, U11_STAR_STAR );
        A11 = A11_STAR_STAR;
        //--------------------------------------------------------------------//
        A01_MC_STAR.FreeAlignments();
        A01_MR_STAR.FreeAlignments();
        A01_VC_STAR.FreeAlignments();
        U01_MC_STAR.FreeAlignments();
        U01_MR_STAR.FreeAlignments();
        U01_VC_STAR.FreeAlignments();
        Y01.FreeAlignments();
        Y01_VC_STAR.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
