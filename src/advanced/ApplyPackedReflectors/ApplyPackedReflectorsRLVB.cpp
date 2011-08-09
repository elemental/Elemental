/*
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental/basic_internal.hpp"
#include "elemental/advanced_internal.hpp"
using namespace elemental;
using namespace std;

#include "./UTUtil.hpp"

template<typename R> // representation of a real number
void
elemental::advanced::internal::ApplyPackedReflectorsRLVB
( int offset, 
  const DistMatrix<R,MC,MR>& H,
        DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::ApplyPackedReflectorsRLVB");
    if( H.Grid() != A.Grid() )
        throw logic_error( "H and A must be distributed over the same grid." );
    if( offset > 0 )
        throw logic_error( "Transforms cannot extend above matrix." );
    if( offset < -H.Height() )
        throw logic_error( "Transforms cannot extend below matrix." );
    if( H.Height() != A.Width() )
        throw logic_error
              ( "Height of transforms must equal width of target matrix." );
#endif
    const Grid& g = H.Grid();

    // Matrix views    
    DistMatrix<R,MC,MR>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g), HPanCopy(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<R,MC,MR> ARight(g);

    DistMatrix<R,VC,  Star> HPan_VC_Star(g);
    DistMatrix<R,MR,  Star> HPan_MR_Star(g);
    DistMatrix<R,Star,Star> SInv_Star_Star(g);
    DistMatrix<R,MC,  Star> Z_MC_Star(g);
    DistMatrix<R,VC,  Star> Z_VC_Star(g);

    LockedPartitionUpDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    while( HBR.Height() < H.Height() && HBR.Width() < H.Width() )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        int HPanHeight = H11.Height() + H21.Height();
        int HPanWidth = min( H11.Width(), max(HPanHeight+offset,0) );
        int leftover = A.Width()-HPanHeight;
        HPan.LockedView( H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        ARight.View( A, 0, leftover, A.Height(), HPanHeight );

        HPan_MR_Star.AlignWith( ARight );
        Z_MC_Star.AlignWith( ARight );
        Z_VC_Star.AlignWith( ARight );
        Z_MC_Star.ResizeTo( ARight.Height(), HPanWidth );
        SInv_Star_Star.ResizeTo( HPanWidth, HPanWidth );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        HPanCopy.MakeTrapezoidal( Left, Lower, offset );
        SetDiagonalToOne( Left, offset, HPanCopy );

        HPan_VC_Star = HPanCopy;
        basic::Syrk
        ( Upper, Transpose, 
          (R)1, HPan_VC_Star.LockedLocalMatrix(),
          (R)0, SInv_Star_Star.LocalMatrix() );     
        SInv_Star_Star.SumOverGrid();
        HalveMainDiagonal( SInv_Star_Star );

        HPan_MR_Star = HPan_VC_Star;
        basic::internal::LocalGemm
        ( Normal, Normal, 
          (R)1, ARight, HPan_MR_Star, (R)0, Z_MC_Star );
        Z_VC_Star.SumScatterFrom( Z_MC_Star );
 
        basic::internal::LocalTrsm
        ( Right, Upper, Transpose, NonUnit,
          (R)1, SInv_Star_Star, Z_VC_Star );

        Z_MC_Star = Z_VC_Star;
        basic::internal::LocalGemm
        ( Normal, Transpose, 
          (R)-1, Z_MC_Star, HPan_MR_Star, (R)1, ARight );
        //--------------------------------------------------------------------//
        HPan_MR_Star.FreeAlignments();
        Z_MC_Star.FreeAlignments();
        Z_VC_Star.FreeAlignments();

        SlideLockedPartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
void
elemental::advanced::internal::ApplyPackedReflectorsRLVB
( Conjugation conjugation, int offset, 
  const DistMatrix<complex<R>,MC,MR  >& H,
  const DistMatrix<complex<R>,MD,Star>& t,
        DistMatrix<complex<R>,MC,MR  >& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::ApplyPackedReflectorsRLVB");
    if( H.Grid() != t.Grid() || t.Grid() != A.Grid() )
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
    const Grid& g = H.Grid();

    // Matrix views    
    DistMatrix<C,MC,MR>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g), HPanCopy(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<C,MC,MR> ARight(g);
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

    LockedPartitionUpDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    LockedPartitionUp
    ( t, tT,
         tB, 0 );
    while( HBR.Height() < H.Height() && HBR.Width() < H.Width() )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        int HPanHeight = H11.Height() + H21.Height();
        int HPanWidth = min( H11.Width(), max(HPanHeight+offset,0) );
        int leftover = A.Width()-HPanHeight;
        HPan.LockedView( H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        LockedRepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2, HPanWidth );

        ARight.View( A, 0, leftover, A.Height(), HPanHeight );

        HPan_MR_Star.AlignWith( ARight );
        Z_MC_Star.AlignWith( ARight );
        Z_VC_Star.AlignWith( ARight );
        Z_MC_Star.ResizeTo( ARight.Height(), HPan.Width() );
        SInv_Star_Star.ResizeTo( HPan.Width(), HPan.Width() );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        HPanCopy.MakeTrapezoidal( Left, Lower, offset );
        SetDiagonalToOne( Left, offset, HPanCopy );

        HPan_VC_Star = HPanCopy;
        basic::Herk
        ( Upper, ConjugateTranspose, 
          (C)1, HPan_VC_Star.LockedLocalMatrix(),
          (C)0, SInv_Star_Star.LocalMatrix() );     
        SInv_Star_Star.SumOverGrid();
        t1_Star_Star = t1;
        FixDiagonal( conjugation, t1_Star_Star, SInv_Star_Star );

        HPan_MR_Star = HPan_VC_Star;
        basic::internal::LocalGemm
        ( Normal, Normal, 
          (C)1, ARight, HPan_MR_Star, (C)0, Z_MC_Star );
        Z_VC_Star.SumScatterFrom( Z_MC_Star );
 
        basic::internal::LocalTrsm
        ( Right, Upper, ConjugateTranspose, NonUnit, 
          (C)1, SInv_Star_Star, Z_VC_Star );

        Z_MC_Star = Z_VC_Star;
        basic::internal::LocalGemm
        ( Normal, ConjugateTranspose, 
          (C)-1, Z_MC_Star, HPan_MR_Star, (C)1, ARight );
        //--------------------------------------------------------------------//
        HPan_MR_Star.FreeAlignments();
        Z_MC_Star.FreeAlignments();
        Z_VC_Star.FreeAlignments();

        SlideLockedPartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );

        SlideLockedPartitionUp
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

template void elemental::advanced::internal::ApplyPackedReflectorsRLVB
( int offset,
  const DistMatrix<float,MC,MR>& H,
        DistMatrix<float,MC,MR>& A );

template void elemental::advanced::internal::ApplyPackedReflectorsRLVB
( int offset,
  const DistMatrix<double,MC,MR>& H,
        DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::advanced::internal::ApplyPackedReflectorsRLVB
( Conjugation conjugation, int offset,
  const DistMatrix<scomplex,MC,MR  >& H,
  const DistMatrix<scomplex,MD,Star>& t,
        DistMatrix<scomplex,MC,MR  >& A );

template void elemental::advanced::internal::ApplyPackedReflectorsRLVB
( Conjugation conjugation, int offset,
  const DistMatrix<dcomplex,MC,MR  >& H,
  const DistMatrix<dcomplex,MD,Star>& t,
        DistMatrix<dcomplex,MC,MR  >& A );
#endif

