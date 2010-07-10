/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
*/
#include "elemental/blas_internal.hpp"
#include "elemental/lapack_internal.hpp"
using namespace elemental;
using namespace std;

namespace {

template<typename T>
void
SetDiagonalToOne( int offset, DistMatrix<T,MC,MR>& H )
{
#ifndef RELEASE
    PushCallStack("SetDiagonalToOne");
#endif
    const int localHeight = H.LocalHeight();
    const int localWidth = H.LocalWidth();
    const int r = H.GetGrid().Height();
    const int c = H.GetGrid().Width();
    const int colShift = H.ColShift();
    const int rowShift = H.RowShift();

    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*r + offset;
        if( i % c == rowShift )
        {
            const int jLoc = (i-rowShift) / c;
            if( jLoc < localWidth )
                H.LocalEntry(iLoc,jLoc) = (T)1;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void 
HalveMainDiagonal( DistMatrix<R,Star,Star>& SInv )
{
#ifndef RELEASE
    PushCallStack("HalveMainDiagonal");
#endif
    for( int j=0; j<SInv.Height(); ++j )
        SInv.LocalEntry(j,j) /= (R)2;
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void 
FixDiagonal
( const DistMatrix<complex<R>,Star,Star>& t,
        DistMatrix<complex<R>,Star,Star>& SInv )
{
#ifndef RELEASE
    PushCallStack("FixDiagonal");
#endif
    for( int j=0; j<SInv.Height(); ++j )
        SInv.LocalEntry(j,j) = complex<R>(1)/t.LocalEntry(j,0);
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

} // anonymous namespace

// This routine applies the accumulation of Householder transforms stored 
// in the portion of H below the diagonal marked by 'offset'. It is assumed 
// that the Householder transforms were accumulated left-to-right.

template<typename R>
void
elemental::lapack::internal::UTRLN
( int offset, 
  const DistMatrix<R,MC,MR>& H,
        DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::UTRLN");
    if( H.GetGrid() != A.GetGrid() )
        throw logic_error( "H and A must be distributed over the same grid." );
    if( offset > 0 )
        throw logic_error( "Transforms cannot extend above matrix." );
    if( offset < -H.Height() )
        throw logic_error( "Transforms cannot extend below matrix." );
    if( H.Height() != A.Width() )
        throw logic_error
              ( "Height of transforms must equal width of target matrix." );
#endif
    const Grid& g = H.GetGrid();

    // Matrix views    
    DistMatrix<R,MC,MR>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g), HPanCopy(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<R,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  ARight(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<R,VC,  Star> HPan_VC_Star(g);
    DistMatrix<R,MR,  Star> HPan_MR_Star(g);
    DistMatrix<R,Star,Star> SInv_Star_Star(g);
    DistMatrix<R,MC,  Star> Z_MC_Star(g);
    DistMatrix<R,VC,  Star> Z_VC_Star(g);

    LockedPartitionDownDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( HTL.Height() < H.Height() && HTL.Width() < H.Width() )
    {
        LockedRepartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );

        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        HPan.LockedView2x1
        ( H11,
          H21 );

        ARight.View2x1
        ( ATR,
          ABR );

        HPan_MR_Star.AlignWith( ARight );
        Z_MC_Star.AlignWith( ARight );
        Z_VC_Star.AlignWith( ARight );
        Z_MC_Star.ResizeTo( ARight.Height(), HPan.Width() );
        SInv_Star_Star.ResizeTo( HPan.Width(), HPan.Width() );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        HPanCopy.MakeTrapezoidal( Left, Lower, offset );
        SetDiagonalToOne( offset, HPanCopy );

        HPan_VC_Star = HPanCopy;
        blas::Herk
        ( Upper, ConjugateTranspose, 
          (R)1, HPan_VC_Star.LockedLocalMatrix(),
          (R)0, SInv_Star_Star.LocalMatrix() );     
        SInv_Star_Star.AllSum();
        HalveMainDiagonal( SInv_Star_Star );

        HPan_MR_Star = HPanCopy;
        blas::internal::LocalGemm
        ( Normal, Normal,
          (R)1, ARight, HPan_MR_Star, (R)0, Z_MC_Star );
        Z_VC_Star.SumScatterFrom( Z_MC_Star );
        
        blas::internal::LocalTrsm
        ( Right, Upper, ConjugateTranspose, NonUnit, 
          (R)1, SInv_Star_Star, Z_VC_Star );

        Z_MC_Star = Z_VC_Star;
        blas::internal::LocalGemm
        ( Normal, ConjugateTranspose,
          (R)-1, Z_MC_Star, HPan_MR_Star, (R)1, ARight );
        //--------------------------------------------------------------------//
        HPan_MR_Star.FreeAlignments();
        Z_MC_Star.FreeAlignments();
        Z_VC_Star.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
               /**/       A10, /**/ A11, A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, /**/ A21, A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::lapack::internal::UTRLN
( int offset, 
  const DistMatrix<complex<R>,MC,MR  >& H,
  const DistMatrix<complex<R>,MD,Star>& t,
        DistMatrix<complex<R>,MC,MR  >& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::UTRLN");
    if( H.GetGrid() != t.GetGrid() || t.GetGrid() != A.GetGrid() )
        throw logic_error
              ( "H, t, and A must be distributed over the same grid." );
    if( offset > 0 )
        throw logic_error( "Transforms cannot extend above matrix." );
    if( offset < -H.Height() )
        throw logic_error( "Transforms cannot extend below matrix." );
    if( H.Height() != A.Width() )
        throw logic_error
              ( "Height of transforms must equal width of target matrix." );
    if( t.Height() != H.DiagonalLength( offset ) )
        throw logic_error( "t must be the same length as H's 'offset' diag." );
    if( !t.AlignedWithDiag( H, offset ) )
        throw logic_error( "t must be aligned with H's 'offset' diagonal." );
#endif
    typedef complex<R> C;
    const Grid& g = H.GetGrid();

    // Matrix views    
    DistMatrix<C,MC,MR>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g), HPanCopy(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<C,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  ARight(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<C,MD,Star>
        tT(g),  t0(g),
        tB(g),  t1(g),
                t2(g);

    DistMatrix<C,VC,  Star> HPan_VC_Star(g);
    DistMatrix<C,MR,  Star> HPan_MR_Star(g);
    DistMatrix<C,Star,Star> t1_Star_Star(g);
    DistMatrix<C,Star,Star> SInv_Star_Star(g);
    DistMatrix<C,MC,  Star> Z_MC_Star(g);
    DistMatrix<C,VC,  Star> Z_VC_Star(g);

    LockedPartitionDownDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    LockedPartitionDown
    ( t, tT,
         tB, 0 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( HTL.Height() < H.Height() && HTL.Width() < H.Width() )
    {
        LockedRepartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );

        LockedRepartitionDown
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );

        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        HPan.LockedView2x1
        ( H11,
          H21 );

        ARight.View2x1
        ( ATR,
          ABR );

        HPan_MR_Star.AlignWith( ARight );
        Z_MC_Star.AlignWith( ARight );
        Z_VC_Star.AlignWith( ARight );
        Z_MC_Star.ResizeTo( ARight.Height(), HPan.Width() );
        SInv_Star_Star.ResizeTo( HPan.Width(), HPan.Width() );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        HPanCopy.MakeTrapezoidal( Left, Lower, offset );
        SetDiagonalToOne( offset, HPanCopy );

        HPan_VC_Star = HPanCopy;
        blas::Herk
        ( Upper, ConjugateTranspose, 
          (C)1, HPan_VC_Star.LockedLocalMatrix(),
          (C)0, SInv_Star_Star.LocalMatrix() );     
        SInv_Star_Star.AllSum();
        t1_Star_Star = t1;
        FixDiagonal( t1_Star_Star, SInv_Star_Star );

        HPan_MR_Star = HPanCopy;
        blas::internal::LocalGemm
        ( Normal, Normal,
          (C)1, ARight, HPan_MR_Star, (C)0, Z_MC_Star );
        Z_VC_Star.SumScatterFrom( Z_MC_Star );
        
        blas::internal::LocalTrsm
        ( Right, Upper, ConjugateTranspose, NonUnit, 
          (C)1, SInv_Star_Star, Z_VC_Star );

        Z_MC_Star = Z_VC_Star;
        blas::internal::LocalGemm
        ( Normal, ConjugateTranspose,
          (C)-1, Z_MC_Star, HPan_MR_Star, (C)1, ARight );
        //--------------------------------------------------------------------//
        HPan_MR_Star.FreeAlignments();
        Z_MC_Star.FreeAlignments();
        Z_VC_Star.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        SlideLockedPartitionDown
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2 ); 

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
               /**/       A10, /**/ A11, A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, /**/ A21, A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

template void elemental::lapack::internal::UTRLN
( int offset,
  const DistMatrix<float,MC,MR>& H,
        DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::UTRLN
( int offset,
  const DistMatrix<double,MC,MR>& H,
        DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::UTRLN
( int offset,
  const DistMatrix<scomplex,MC,MR  >& H,
  const DistMatrix<scomplex,MD,Star>& t,
        DistMatrix<scomplex,MC,MR  >& A );

template void elemental::lapack::internal::UTRLN
( int offset,
  const DistMatrix<dcomplex,MC,MR  >& H,
  const DistMatrix<dcomplex,MD,Star>& t,
        DistMatrix<dcomplex,MC,MR  >& A );
#endif

