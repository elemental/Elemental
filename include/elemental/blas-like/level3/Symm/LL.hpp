/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

template<typename T>
inline void
SymmLLA
( T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::SymmLLA");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
#endif
    const Grid& g = A.Grid();

    DistMatrix<T> 
        BL(g), BR(g),
        B0(g), B1(g), B2(g);

    DistMatrix<T>
        CL(g), CR(g),
        C0(g), C1(g), C2(g);

    DistMatrix<T,MC,STAR> B1_MC_STAR(g);
    DistMatrix<T,VR,STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR> B1Trans_STAR_MR(g);
    DistMatrix<T> Z1(g);
    DistMatrix<T,MC,STAR> Z1_MC_STAR(g);
    DistMatrix<T,MR,STAR> Z1_MR_STAR(g);
    DistMatrix<T,MR,MC  > Z1_MR_MC(g);

    B1_MC_STAR.AlignWith( A );
    B1_VR_STAR.AlignWith( A );
    B1Trans_STAR_MR.AlignWith( A );
    Z1_MC_STAR.AlignWith( A );
    Z1_MR_STAR.AlignWith( A );

    Scale( beta, C );
    LockedPartitionRight
    ( B, BL, BR, 0 );
    PartitionRight
    ( C, CL, CR, 0 );
    while( CL.Width() < C.Width() )
    {
        LockedRepartitionRight 
        ( BL, /**/ BR,
          B0, /**/ B1, B2 );

        RepartitionRight
        ( CL, /**/ CR,
          C0, /**/ C1, C2 );

        Z1.AlignWith( C1 );
        Zeros( C1.Height(), C1.Width(), Z1_MC_STAR );
        Zeros( C1.Height(), C1.Width(), Z1_MR_STAR );
        //--------------------------------------------------------------------//
        B1_MC_STAR = B1;
        B1_VR_STAR = B1_MC_STAR;
        B1Trans_STAR_MR.TransposeFrom( B1_VR_STAR );
        LocalSymmetricAccumulateLL
        ( TRANSPOSE, 
          alpha, A, B1_MC_STAR, B1Trans_STAR_MR, Z1_MC_STAR, Z1_MR_STAR );

        Z1_MR_MC.SumScatterFrom( Z1_MR_STAR );
        Z1 = Z1_MR_MC;
        Z1.SumScatterUpdate( T(1), Z1_MC_STAR );
        Axpy( T(1), Z1, C1 );
        //--------------------------------------------------------------------//
        Z1.FreeAlignments();

        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );

        SlidePartitionRight
        ( CL,     /**/ CR,
          C0, C1, /**/ C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
SymmLLC
( T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::SymmLLC");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  AColPan(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),  ARowPan(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<T> 
        BT(g),  B0(g),
        BB(g),  B1(g),
                B2(g);
    DistMatrix<T> 
        CT(g),  C0(g),  CAbove(g),
        CB(g),  C1(g),  CBelow(g),
                C2(g);

    // Temporary distributions
    DistMatrix<T,MC,  STAR> AColPan_MC_STAR(g);
    DistMatrix<T,STAR,MC  > ARowPan_STAR_MC(g);
    DistMatrix<T,MR,  STAR> B1Trans_MR_STAR(g);

    B1Trans_MR_STAR.AlignWith( C );

    // Start the algorithm
    Scale( beta, C );
    LockedPartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionDown
    ( B, BT,
         BB, 0 );
    PartitionDown
    ( C, CT,
         CB, 0 );
    while( CB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionDown
        ( BT,  B0,
         /**/ /**/
               B1,
          BB,  B2 );

        RepartitionDown
        ( CT,  C0,
         /**/ /**/
               C1,
          CB,  C2 );

        LockedView1x2( ARowPan, A10, A11 );
        LockedView2x1
        ( AColPan, A11,
                   A21 );

        View2x1
        ( CAbove, C0,
                  C1 );
        View2x1
        ( CBelow, C1,
                  C2 );

        AColPan_MC_STAR.AlignWith( CBelow );
        ARowPan_STAR_MC.AlignWith( CAbove );
        //--------------------------------------------------------------------//
        AColPan_MC_STAR = AColPan;
        ARowPan_STAR_MC = ARowPan;
        MakeTrapezoidal( LEFT,  LOWER,  0, AColPan_MC_STAR );
        MakeTrapezoidal( RIGHT, LOWER, -1, ARowPan_STAR_MC );

        B1Trans_MR_STAR.TransposeFrom( B1 );

        LocalGemm
        ( NORMAL, TRANSPOSE, 
          alpha, AColPan_MC_STAR, B1Trans_MR_STAR, T(1), CBelow );

        LocalGemm
        ( TRANSPOSE, TRANSPOSE, 
          alpha, ARowPan_STAR_MC, B1Trans_MR_STAR, T(1), CAbove );
        //--------------------------------------------------------------------//
        AColPan_MC_STAR.FreeAlignments();
        ARowPan_STAR_MC.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDown
        ( BT,  B0,
               B1,
         /**/ /**/
          BB,  B2 );

        SlidePartitionDown
        ( CT,  C0,
               C1,
         /**/ /**/
          CB,  C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
SymmLL
( T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::SymmLL");
#endif
    // TODO: Come up with a better routing mechanism
    if( A.Height() > 5*B.Width() )
        SymmLLA( alpha, A, B, beta, C );
    else
        SymmLLC( alpha, A, B, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
