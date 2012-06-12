/*
   Copyright (c) 2009-2012, Jack Poulson
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

namespace elem {
namespace internal {

template<typename T>
inline void
AddOneToDiagonal( LeftOrRight side, int offset, Matrix<T>& H )
{
#ifndef RELEASE
    PushCallStack("AddOneToDiagonal");
#endif
    const int height = H.Height();
    const int width = H.Width();

    if( side == LEFT )
    {
        for( int j=0; j<width; ++j )
        {
            const int i = j-offset;
            if( i >= 0 && i < height )
                H.Update(i,j,1);
        }
    }
    else
    {
        for( int j=0; j<width; ++j )
        {
            const int i = j-offset+height-width;
            if( i >= 0 && i < height )
                H.Update(i,j,1);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
AddOneToDiagonal( LeftOrRight side, int offset, DistMatrix<T>& H )
{
#ifndef RELEASE
    PushCallStack("AddOneToDiagonal");
#endif
    const int height = H.Height();
    const int width = H.Width();
    const int localWidth = H.LocalWidth();
    const int r = H.Grid().Height();
    const int c = H.Grid().Width();
    const int colShift = H.ColShift();
    const int rowShift = H.RowShift();

    if( side == LEFT )
    {
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*c;
            const int i = j-offset;
            if( i >= 0 && i < height && (i-colShift) % r == 0 )
            {
                const int iLoc = (i-colShift)/r;
                H.UpdateLocal(iLoc,jLoc,1);
            }
        }
    }
    else
    {
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*c;
            const int i = j-offset+height-width;
            if( i >= 0 && i < height && (i-colShift) % r == 0 )
            {
                const int iLoc = (i-colShift)/r;
                H.UpdateLocal(iLoc,jLoc,1);
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Since applying Householder transforms from vectors stored right-to-left
// implies that we will be forming a generalization of
//
//   (I - tau_0 u_0 u_0^H) (I - tau_1 u_1 u_1^H) = 
//   I - tau_0 u_0 u_0^H - tau_1 u_1 u_1^H + (tau_0 tau_1 u_0^H u_1) u_0 u_1^H =
//   I - [ u_0, u_1 ] [ tau_0, -tau_0 tau_1 u_0^H u_1 ] [ u_0^H ]
//                    [ 0,      tau_1                 ] [ u_1^H ],
//
// which has an upper-triangular center matrix, say S, we will form S as 
// the inverse of a matrix T, which can easily be formed as
// 
//   triu(T) = triu( U^H U ),  diag(T) = 1/t or 1/conj(t),
//
// where U is the matrix of Householder vectors and t is the vector of scalars.
//

template<typename R> 
inline void
ExpandPackedReflectorsLV( int offset, Matrix<R>& H )
{
#ifndef RELEASE
    PushCallStack("internal::ExpandPackedReflectorsLV");
    if( offset > 0 || offset < -H.Height() )
        throw std::logic_error("Transforms out of bounds");
#endif
    // Start by zeroing everything above the offset and setting that diagonal
    // to all ones. We can also ensure that H is not wider than it is tall.
    if( H.Width() > H.Height() )
        H.ResizeTo( H.Height(), H.Height() );
    MakeTrapezoidal( LEFT, LOWER, offset, H );
    SetDiagonalToOne( LEFT, offset, H );
    const int dimDiff = H.Height() - H.Width();

    Matrix<R>
        HTL, HTR,  H00, H01, H02,  HPan, HPanCopy, HPanT,
        HBL, HBR,  H10, H11, H12,                  HPanB,
                   H20, H21, H22;
    Matrix<R> HEffectedNew, HEffectedOld, HEffectedOldB;

    Matrix<R> SInv, Z;

    LockedPartitionUpDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    int oldEffectedHeight=dimDiff;
    while( HBR.Height() < H.Height() && HBR.Width() < H.Width() )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        const int HPanHeight = H11.Height() + H21.Height();
        const int effectedHeight = std::max(HPanHeight+offset,0);
        const int HPanWidth = std::min( H11.Width(), effectedHeight );

        const int oldEffectedWidth = oldEffectedHeight - dimDiff;
        const int effectedWidth = effectedHeight - dimDiff;

        HPan.LockedView
        ( H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );
        LockedPartitionDown
        ( HPan, HPanT,
                HPanB, effectedHeight-oldEffectedHeight );

        HEffectedOld.View
        ( H, H.Height()-effectedHeight, H.Width()-oldEffectedWidth, 
          effectedHeight, oldEffectedWidth );
        HEffectedOldB.View
        ( H, H.Height()-oldEffectedHeight, H.Width()-oldEffectedWidth,
          oldEffectedHeight, oldEffectedWidth );

        HEffectedNew.View
        ( H, H.Height()-effectedHeight, H.Width()-effectedWidth, 
          effectedHeight, effectedWidth-oldEffectedWidth );

        Zeros( HPanWidth, oldEffectedWidth, Z );
        Zeros( HPanWidth, HPanWidth, SInv );
        //--------------------------------------------------------------------//
        Syrk( UPPER, TRANSPOSE, (R)1, HPan, (R)0, SInv );
        HalveMainDiagonal( SInv );

        // Update the already effected portion of the matrix
        Gemm( TRANSPOSE, NORMAL, (R)1, HPanB, HEffectedOldB, (R)0, Z );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (R)1, SInv, Z );
        Gemm( NORMAL, NORMAL, (R)-1, HPan, Z, (R)1, HEffectedOld );

        // Update the newly effected portion of the matrix
        Transpose( HPanT, Z );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (R)1, SInv, Z );
        HPanCopy = HPan;
        Gemm( NORMAL, NORMAL, (R)-1, HPanCopy, Z, (R)0, HEffectedNew );
        AddOneToDiagonal( LEFT, 0, HEffectedNew );
        //--------------------------------------------------------------------//

        oldEffectedHeight = effectedHeight;

        SlideLockedPartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );
    }

    // Take care of any untouched columns on the left side of H
    const int oldEffectedWidth = oldEffectedHeight - dimDiff;
    if( oldEffectedWidth < H.Width() )
    {
        HEffectedNew.View( H, 0, 0, H.Height(), H.Width()-oldEffectedWidth );
        MakeZeros( HEffectedNew );
        SetDiagonalToOne( LEFT, 0, HEffectedNew );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
ExpandPackedReflectorsLV
( Conjugation conjugation, int offset,
  Matrix<Complex<R> >& H, const Matrix<Complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("internal::ExpandPackedReflectorsLLVB");
    if( offset > 0 || offset < -H.Height() )
        throw std::logic_error("Transforms out of bounds");
    if( t.Height() != H.DiagonalLength( offset ) )
        throw std::logic_error("t must be the same length as H's offset diag");
#endif
    typedef Complex<R> C;

    // Start by zeroing everything above the offset and setting that diagonal
    // to all ones. We can also ensure that H is not wider than it is tall.
    if( H.Width() > H.Height() )
        H.ResizeTo( H.Height(), H.Height() );
    MakeTrapezoidal( LEFT, LOWER, offset, H );
    SetDiagonalToOne( LEFT, offset, H );
    const int dimDiff = H.Height() - H.Width();

    Matrix<C>
        HTL, HTR,  H00, H01, H02,  HPan, HPanCopy, HPanT,
        HBL, HBR,  H10, H11, H12,                  HPanB,
                   H20, H21, H22;
    Matrix<C> HEffectedNew, HEffectedOld, HEffectedOldB;
    Matrix<C>
        tT,  t0,
        tB,  t1,
             t2;

    Matrix<C> SInv, Z;

    LockedPartitionUpDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    LockedPartitionUp
    ( t, tT,
         tB, 0 );
    int oldEffectedHeight=dimDiff;
    while( HBR.Height() < H.Height() && HBR.Width() < H.Width() )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        const int HPanHeight = H11.Height() + H21.Height();
        const int effectedHeight = std::max(HPanHeight+offset,0);
        const int HPanWidth = std::min( H11.Width(), effectedHeight );

        const int oldEffectedWidth = oldEffectedHeight - dimDiff;
        const int effectedWidth = effectedHeight - dimDiff;

        HPan.LockedView
        ( H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );
        LockedPartitionDown
        ( HPan, HPanT,
                HPanB, effectedHeight-oldEffectedHeight );

        HEffectedOld.View
        ( H, H.Height()-effectedHeight, H.Width()-oldEffectedWidth,
          effectedHeight, oldEffectedWidth );
        HEffectedOldB.View
        ( H, H.Height()-oldEffectedHeight, H.Width()-oldEffectedWidth,
          oldEffectedHeight, oldEffectedWidth );

        HEffectedNew.View
        ( H, H.Height()-effectedHeight, H.Width()-effectedWidth,
          effectedHeight, effectedWidth-oldEffectedWidth );

        LockedRepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2, HPanWidth );

        Zeros( HPanWidth, oldEffectedWidth, Z );
        Zeros( HPanWidth, HPanWidth, SInv );
        //--------------------------------------------------------------------//
        Herk( UPPER, ADJOINT, (C)1, HPan, (C)0, SInv );
        FixDiagonal( conjugation, t1, SInv );

        // Update the already effected portion of the matrix
        Gemm( ADJOINT, NORMAL, (C)1, HPanB, HEffectedOldB, (C)0, Z );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (C)1, SInv, Z );
        Gemm( NORMAL, NORMAL, (C)-1, HPan, Z, (C)1, HEffectedOld );

        // Update the newly effected portion of the matrix
        Adjoint( HPanT, Z );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (C)1, SInv, Z );
        HPanCopy = HPan;
        Gemm( NORMAL, NORMAL, (C)-1, HPanCopy, Z, (C)0, HEffectedNew );
        AddOneToDiagonal( LEFT, 0, HEffectedNew );
        //--------------------------------------------------------------------//

        oldEffectedHeight = effectedHeight;

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

    // Take care of any untouched columns on the left side of H
    const int oldEffectedWidth = oldEffectedHeight - dimDiff;
    if( oldEffectedWidth < H.Width() )
    {
        HEffectedNew.View( H, 0, 0, H.Height(), H.Width()-oldEffectedWidth );
        MakeZeros( HEffectedNew );
        SetDiagonalToOne( LEFT, 0, HEffectedNew );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
