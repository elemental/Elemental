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
SetDiagonalToOne( int offset, DistMatrix<T,Star,VR>& H )
{
#ifndef RELEASE
    PushCallStack("SetDiagonalToOne");
#endif
    const int height = H.Height();
    const int localWidth = H.LocalWidth();
    const int p = H.GetGrid().Size();
    const int rowShift = H.RowShift();

    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*p;
        if( (j-offset) >= 0 && (j-offset) < height )
            H.LocalEntry(j-offset,jLoc) = (T)1;
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

// This routine reverses the accumulation of Householder transforms stored 
// in the portion of H above the diagonal marked by 'offset'. It is assumed 
// that the Householder transforms were accumulated right-to-left.

template<typename R>
void
elemental::lapack::internal::UTLUH
( int offset, 
  const DistMatrix<R,MC,MR>& H,
        DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::UTLUH");
    if( H.GetGrid() != A.GetGrid() )
        throw logic_error( "H and A must be distributed over the same grid." );
    if( offset > H.Height() )
        throw logic_error( "Transforms cannot extend above matrix." );
    if( offset < 0 )
        throw logic_error( "Transforms cannot extend below matrix." );
    if( H.Width() != A.Height() )
        throw logic_error
              ( "Width of transforms must equal height of target matrix." );
#endif
    const Grid& g = H.GetGrid();

    // Matrix views    
    DistMatrix<R,MC,MR>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<R,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  ABottom(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<R,Star,VR  > HPan_Star_VR(g);
    DistMatrix<R,Star,MC  > HPan_Star_MC(g);
    DistMatrix<R,Star,Star> SInv_Star_Star(g);
    DistMatrix<R,Star,MR  > Z_Star_MR(g);
    DistMatrix<R,Star,VR  > Z_Star_VR(g);

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

        HPan.LockedView1x2( H11, H12 );

        ABottom.View1x2( ABL, ABR );

        HPan_Star_MC.AlignWith( ABottom );
        Z_Star_MR.AlignWith( ABottom );
        Z_Star_VR.AlignWith( ABottom );
        Z_Star_MR.ResizeTo( HPan.Width(), ABottom.Width() );
        SInv_Star_Star.ResizeTo( HPan.Width(), HPan.Width() );
        //--------------------------------------------------------------------//
        HPan_Star_VR = HPan;
        HPan_Star_VR.MakeTrapezoidal( Left, Upper, offset );
        SetDiagonalToOne( offset, HPan_Star_VR );
        blas::Herk
        ( Upper, Normal, 
          (R)1, HPan_Star_VR.LockedLocalMatrix(),
          (R)0, SInv_Star_Star.LocalMatrix() ); 
        SInv_Star_Star.AllSum();
        HalveMainDiagonal( SInv_Star_Star );

        HPan_Star_MC = HPan_Star_VR;
        blas::internal::LocalGemm
        ( Normal, Normal, (R)1, HPan_Star_MC, ABottom, (R)0, Z_Star_MR );
        Z_Star_VR.SumScatterFrom( Z_Star_MR );
        
        blas::internal::LocalTrsm
        ( Left, Upper, ConjugateTranspose, NonUnit, 
          (R)1, SInv_Star_Star, Z_Star_VR );

        Z_Star_MR = Z_Star_VR;
        blas::internal::LocalGemm
        ( ConjugateTranspose, Normal, 
          (R)-1, HPan_Star_MC, Z_Star_MR, (R)1, ABottom );
        //--------------------------------------------------------------------//
        HPan_Star_MC.FreeAlignments();
        Z_Star_MR.FreeAlignments();
        Z_Star_VR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
         /*************/ /******************/
               /**/       A10, A11, /**/ A12,
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::lapack::internal::UTLUH
( int offset, 
  const DistMatrix<complex<R>,MC,MR  >& H,
  const DistMatrix<complex<R>,MD,Star>& t,
        DistMatrix<complex<R>,MC,MR  >& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::UTLUH");
    if( H.GetGrid() != t.GetGrid() || t.GetGrid() != A.GetGrid() )
        throw logic_error
              ( "H, t, and A must be distributed over the same grid." );
    if( offset > H.Height() )
        throw logic_error( "Transforms cannot extend above matrix." );
    if( offset < 0 )
        throw logic_error( "Transforms cannot extend below matrix." );
    if( H.Width() != A.Height() )
        throw logic_error
              ( "Width of transforms must equal height of target matrix." );
    if( t.Height() != H.DiagonalLength( offset ) )
        throw logic_error( "t must be the same length as H's 'offset' diag." );
    if( !t.AlignedWithDiag( H, offset ) )
        throw logic_error( "t must be aligned with H's 'offset' diagonal." );
#endif
    typedef complex<R> C;
    const Grid& g = H.GetGrid();

    // Matrix views    
    DistMatrix<C,MC,MR>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<C,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  ABottom(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<C,MD,Star>
        tT(g),  t0(g),
        tB(g),  t1(g),
                t2(g);

    DistMatrix<C,Star,VR  > HPan_Star_VR(g);
    DistMatrix<C,Star,MC  > HPan_Star_MC(g);
    DistMatrix<C,Star,Star> t1_Star_Star(g);
    DistMatrix<C,Star,Star> SInv_Star_Star(g);
    DistMatrix<C,Star,MR  > Z_Star_MR(g);
    DistMatrix<C,Star,VR  > Z_Star_VR(g);

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

        HPan.LockedView1x2( H11, H12 );

        ABottom.View1x2( ABL, ABR );

        HPan_Star_MC.AlignWith( ABottom );
        Z_Star_MR.AlignWith( ABottom );
        Z_Star_VR.AlignWith( ABottom );
        Z_Star_MR.ResizeTo( HPan.Width(), ABottom.Width() );
        SInv_Star_Star.ResizeTo( HPan.Width(), HPan.Width() );
        //--------------------------------------------------------------------//
        HPan_Star_VR = HPan;
        HPan_Star_VR.MakeTrapezoidal( Left, Upper, offset );
        SetDiagonalToOne( offset, HPan_Star_VR );
        blas::Herk
        ( Upper, Normal, 
          (C)1, HPan_Star_VR.LockedLocalMatrix(),
          (C)0, SInv_Star_Star.LocalMatrix() ); 
        SInv_Star_Star.AllSum();
        t1_Star_Star = t1;
        FixDiagonal( t1_Star_Star, SInv_Star_Star );

        HPan_Star_MC = HPan_Star_VR;
        blas::internal::LocalGemm
        ( Normal, Normal, (C)1, HPan_Star_MC, ABottom, (C)0, Z_Star_MR );
        Z_Star_VR.SumScatterFrom( Z_Star_MR );
        
        blas::internal::LocalTrsm
        ( Left, Upper, ConjugateTranspose, NonUnit, 
          (C)1, SInv_Star_Star, Z_Star_VR );

        Z_Star_MR = Z_Star_VR;
        blas::internal::LocalGemm
        ( ConjugateTranspose, Normal, 
          (C)-1, HPan_Star_MC, Z_Star_MR, (C)1, ABottom );
        //--------------------------------------------------------------------//
        HPan_Star_MC.FreeAlignments();
        Z_Star_MR.FreeAlignments();
        Z_Star_VR.FreeAlignments();

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

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
         /*************/ /******************/
               /**/       A10, A11, /**/ A12,
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

template void elemental::lapack::internal::UTLUH
( int offset,
  const DistMatrix<float,MC,MR>& H,
        DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::UTLUH
( int offset,
  const DistMatrix<double,MC,MR>& H,
        DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::UTLUH
( int offset,
  const DistMatrix<scomplex,MC,MR  >& H,
  const DistMatrix<scomplex,MD,Star>& t,
        DistMatrix<scomplex,MC,MR  >& A );

template void elemental::lapack::internal::UTLUH
( int offset,
  const DistMatrix<dcomplex,MC,MR  >& H,
  const DistMatrix<dcomplex,MD,Star>& t,
        DistMatrix<dcomplex,MC,MR  >& A );
#endif

