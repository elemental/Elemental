/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

// The only reason a field is required is for the existence of 1/2, which is 
// an artifact of the algorithm...
template<typename F> 
inline void
TwoSidedTrmmLVar5( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& L )
{
#ifndef RELEASE
    PushCallStack("internal::TwoSidedTrmmLVar5");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( L.Height() != L.Width() )
        throw std::logic_error("Triangular matrices must be square");
    if( A.Height() != L.Height() )
        throw std::logic_error("A and L must be the same size");
#endif
    // Matrix views
    Matrix<F>
        ATL, ATR,  A00, A01, A02,
        ABL, ABR,  A10, A11, A12,
                   A20, A21, A22;
    Matrix<F>
        LTL, LTR,  L00, L01, L02,
        LBL, LBR,  L10, L11, L12,
                   L20, L21, L22;

    // Temporary products
    Matrix<F> Y10;

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        //--------------------------------------------------------------------//
        // Y10 := A11 L10
        Zeros( A10.Height(), A10.Width(), Y10 );
        Hemm( LEFT, LOWER, F(1), A11, L10, F(0), Y10 );

        // A10 := A10 L00
        Trmm( RIGHT, LOWER, NORMAL, diag, F(1), L00, A10 );

        // A10 := A10 + 1/2 Y10
        Axpy( F(1)/F(2), Y10, A10 );

        // A00 := A00 + (L10' A10 + A10' L10)
        Her2k( LOWER, ADJOINT, F(1), L10, A10, F(1), A00 );

        // A10 := A10 + 1/2 Y10
        Axpy( F(1)/F(2), Y10, A10 );

        // A10 := L11' A10
        Trmm( LEFT, LOWER, ADJOINT, diag, F(1), L11, A10 );

        // A11 := L11' A11 L11
        TwoSidedTrmmLUnb( diag, A11, L11 );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /**********************************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
TwoSidedTrmmLVar5
( UnitOrNonUnit diag, DistMatrix<F>& A, const DistMatrix<F>& L )
{
#ifndef RELEASE
    PushCallStack("internal::TwoSidedTrmmLVar5");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( L.Height() != L.Width() )
        throw std::logic_error("Triangular matrices must be square");
    if( A.Height() != L.Height() )
        throw std::logic_error("A and L must be the same size");
#endif
    const Grid& g = A.Grid();
    
    // Matrix views
    DistMatrix<F>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<F>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,MC  > A10_STAR_MC(g);
    DistMatrix<F,MR,  STAR> A10Trans_MR_STAR(g);
    DistMatrix<F,STAR,VR  > A10_STAR_VR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<F,STAR,MR  > L10_STAR_MR(g);
    DistMatrix<F,STAR,VR  > L10_STAR_VR(g);
    DistMatrix<F,STAR,VR  > Y10_STAR_VR(g);
    DistMatrix<F> Y10(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        A10_STAR_MC.AlignWith( A00 );
        A10Trans_MR_STAR.AlignWith( A00 );
        A10_STAR_VR.AlignWith( A10 );
        L10_STAR_MC.AlignWith( A00 );
        L10_STAR_MR.AlignWith( A00 );
        L10_STAR_VR.AlignWith( A00 );
        Y10.AlignWith( A10 );
        Y10_STAR_VR.AlignWith( A10 );
        //--------------------------------------------------------------------//
        // Y10 := A11 L10
        A11_STAR_STAR = A11;
        L10_STAR_VR = L10;
        Y10_STAR_VR.ResizeTo( A10.Height(), A10.Width() );
        Hemm
        ( LEFT, LOWER, 
          F(1), A11_STAR_STAR.LocalMatrix(), L10_STAR_VR.LocalMatrix(),
          F(0), Y10_STAR_VR.LocalMatrix() );
        Y10 = Y10_STAR_VR;

        // A10 := A10 L00
        Trmm( RIGHT, LOWER, NORMAL, diag, F(1), L00, A10 );

        // A10 := A10 + 1/2 Y10
        Axpy( F(1)/F(2), Y10, A10 );

        // A00 := A00 + (L10' A10 + A10' L10)
        A10Trans_MR_STAR.TransposeFrom( A10 );
        L10_STAR_MR = L10_STAR_VR;
        A10_STAR_VR = A10;
        A10_STAR_MC = A10_STAR_VR;
        L10_STAR_MC = L10_STAR_VR;
        LocalTrr2k
        ( LOWER, ADJOINT, TRANSPOSE, ADJOINT,
          F(1), L10_STAR_MC, A10Trans_MR_STAR, 
                A10_STAR_MC, L10_STAR_MR, 
          F(1), A00 );

        // A10 := A10 + 1/2 Y10
        Axpy( F(1)/F(2), Y10_STAR_VR, A10_STAR_VR );

        // A10 := L11' A10
        L11_STAR_STAR = L11;
        LocalTrmm
        ( LEFT, LOWER, ADJOINT, diag, F(1), L11_STAR_STAR, A10_STAR_VR );
        A10 = A10_STAR_VR;

        // A11 := L11' A11 L11
        LocalTwoSidedTrmm( LOWER, diag, A11_STAR_STAR, L11_STAR_STAR );
        A11 = A11_STAR_STAR;
        //--------------------------------------------------------------------//
        A10_STAR_MC.FreeAlignments();
        A10Trans_MR_STAR.FreeAlignments();
        A10_STAR_VR.FreeAlignments();
        L10_STAR_MC.FreeAlignments();
        L10_STAR_MR.FreeAlignments();
        L10_STAR_VR.FreeAlignments();
        Y10.FreeAlignments();
        Y10_STAR_VR.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /**********************************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
