/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef LAPACK_HERMITIANTRIDIAG_PANELL_HPP
#define LAPACK_HERMITIANTRIDIAG_PANELL_HPP

#include "elemental/blas-like/level1/Zero.hpp"
#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/blas-like/level2/Symv.hpp"
#include "elemental/lapack-like/Reflector/Col.hpp"

namespace elem {
namespace hermitian_tridiag {

template<typename F>
void PanelL
( DistMatrix<F>& A,
  DistMatrix<F>& W,
  DistMatrix<F,MD,STAR>& t,
  DistMatrix<F,MC,STAR>& APan_MC_STAR, 
  DistMatrix<F,MR,STAR>& APan_MR_STAR,
  DistMatrix<F,MC,STAR>& W_MC_STAR,
  DistMatrix<F,MR,STAR>& W_MR_STAR )
{
    const Int panelSize = W.Width();
    const Int bottomSize = W.Height()-panelSize;
#ifndef RELEASE
    PushCallStack("hermitian_tridiag::PanelL");
    if( A.Grid() != W.Grid() || W.Grid() != t.Grid() )
        LogicError("A, W, and t must be distributed over the same grid.");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( A.Height() != W.Height() )
        LogicError("A and W must be the same height");
    if( W.Height() < panelSize )
        LogicError("W must be a column panel");
    if( W.ColAlignment() != A.ColAlignment() || 
        W.RowAlignment() != A.RowAlignment() )
        LogicError("W and A must be aligned");
    if( t.Height() != W.Width() || t.Width() != 1 )
        LogicError("t must be a column vector of the same length as W's width");
    if( !t.AlignedWithDiagonal(A,-1) )
        LogicError("t is not aligned with A's subdiagonal.");
#endif
    typedef BASE(F) R;

    const Grid& g = A.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();

    // Create a distributed matrix for storing the subdiagonal
    DistMatrix<R,MD,STAR> e(g);
    e.AlignWithDiagonal( A.DistData(), -1 );
    e.ResizeTo( panelSize, 1 );

    // Matrix views 
    DistMatrix<F> 
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  ACol(g),  alpha21T(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),            a21B(g),
                         A20(g), a21(g),     A22(g),  A20B(g);
    DistMatrix<F> 
        WTL(g), WTR(g),  W00(g), w01(g),     W02(g),  WCol(g),
        WBL(g), WBR(g),  w10(g), omega11(g), w12(g),
                         W20(g), w21(g),     W22(g),  W20B(g), w21Last(g);
    DistMatrix<R,MD,STAR> eT(g),  e0(g),
                          eB(g),  epsilon1(g),
                                  e2(g);
    DistMatrix<F,MD,STAR>
        tT(g),  t0(g),
        tB(g),  tau1(g),
                t2(g);

    // Temporary distributions
    std::vector<F> w21LastBuffer(A.Height()/r+1);
    DistMatrix<F,MC,STAR> a21_MC_STAR(g), a21B_MC_STAR(g), a21Last_MC_STAR(g);
    DistMatrix<F,MR,STAR> a21_MR_STAR(g), a21Last_MR_STAR(g);
    DistMatrix<F,MC,STAR> p21_MC_STAR(g), p21B_MC_STAR(g);
    DistMatrix<F,MR,STAR> q21_MR_STAR(g);
    DistMatrix<F,MR,STAR> x01_MR_STAR(g);
    DistMatrix<F,MR,STAR> y01_MR_STAR(g);
    DistMatrix<F,MC,STAR> w21Last_MC_STAR(g);
    DistMatrix<F,MR,STAR> w21Last_MR_STAR(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, 0 );
    PartitionDown
    ( e, eT,
         eB, 0 );
    PartitionDown
    ( t, tT,
         tB, 0 );
    bool firstIteration = true;
    F tau = 0;
    F w21LastFirstEntry = 0;
    while( WTL.Width() < panelSize )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12, 
          ABL, /**/ ABR,  A20, /**/ a21,     A22, 1 );
        
        RepartitionDownDiagonal
        ( WTL, /**/ WTR,  W00, /**/ w01,     W02,
         /*************/ /**********************/
               /**/       w10, /**/ omega11, w12,
          WBL, /**/ WBR,  W20, /**/ w21,     W22, 1 );

        RepartitionDown
        ( eT,  e0,
         /**/ /********/
               epsilon1,
          eB,  e2, 1 );

        RepartitionDown
        ( tT,  t0,
         /**/ /****/
               tau1,
          tB,  t2, 1 );

        View2x1
        ( ACol, alpha11,
                a21 );

        View2x1
        ( WCol, omega11,
                w21 );

        // View the portions of A20 and W20 outside of this panel's square
        View( A20B, A, panelSize, 0, bottomSize, A20.Width() );
        View( W20B, W, panelSize, 0, bottomSize, W20.Width() );

        if( !firstIteration )
        {
            View
            ( a21Last_MC_STAR,
              APan_MC_STAR, WTL.Height(), WTL.Width()-1, WBL.Height(), 1 );
            View
            ( a21Last_MR_STAR,
              APan_MR_STAR, WTL.Height(), WTL.Width()-1, WBL.Height(), 1 );
            View( w21Last, W, WTL.Height(), WTL.Width()-1, WBL.Height(), 1 );
        }
            
        PartitionDown
        ( a21, alpha21T,
               a21B,     1 );

        a21_MC_STAR.AlignWith( A22 );
        a21_MR_STAR.AlignWith( A22 );
        p21_MC_STAR.AlignWith( A22 );
        q21_MR_STAR.AlignWith( A22 );
        x01_MR_STAR.AlignWith( W20B );
        y01_MR_STAR.AlignWith( W20B );
        
        a21_MC_STAR.ResizeTo( a21.Height(), 1 );
        a21_MR_STAR.ResizeTo( a21.Height(), 1 );
        p21_MC_STAR.ResizeTo( a21.Height(), 1 );

        // View the portions of a21[MC,* ] and p21[MC,* ] below the current
        // panel's square
        View
        ( a21B_MC_STAR,
          a21_MC_STAR, a21_MC_STAR.Height()-bottomSize, 0, bottomSize, 1 );
        View
        ( p21B_MC_STAR,
          p21_MC_STAR, p21_MC_STAR.Height()-bottomSize, 0, bottomSize, 1 );
        //--------------------------------------------------------------------//
        const bool thisIsMyCol = ( g.Col() == alpha11.RowAlignment() );
        if( thisIsMyCol )
        {
            if( !firstIteration )
            {
                // Finish updating the current column with two axpy's
                const Int AColLocalHeight = ACol.LocalHeight();
                F* AColBuffer = ACol.Buffer();
                const F* a21Last_MC_STAR_Buffer = a21Last_MC_STAR.Buffer();
                for( Int i=0; i<AColLocalHeight; ++i )
                    AColBuffer[i] -=
                        w21LastBuffer[i] + 
                        a21Last_MC_STAR_Buffer[i]*Conj(w21LastFirstEntry);
            }
            // Compute the Householder reflector
            tau = reflector::Col( alpha21T, a21B );
            if( g.Row() == alpha21T.ColAlignment() )
                tau1.SetLocal(0,0,tau);
        }

        // Store the subdiagonal value and turn a21 into a proper scaled 
        // reflector by explicitly placing the implicit one in its first entry.
        alpha21T.GetRealPartOfDiagonal( epsilon1 );
        alpha21T.Set( 0, 0, F(1) );

        // If this is the first iteration, have each member of the owning 
        // process column broadcast tau and a21 within its process row. 
        // Otherwise, also add w21 into the broadcast.
        if( firstIteration )
        {
            const Int a21LocalHeight = a21.LocalHeight();
            std::vector<F> rowBroadcastBuffer(a21LocalHeight+1);
            if( thisIsMyCol )
            {
                // Pack the broadcast buffer with a21 and tau
                MemCopy( &rowBroadcastBuffer[0], a21.Buffer(), a21LocalHeight );
                rowBroadcastBuffer[a21LocalHeight] = tau;
            }
            // Broadcast a21 and tau across the process row
            mpi::Broadcast
            ( &rowBroadcastBuffer[0], 
              a21LocalHeight+1, a21.RowAlignment(), g.RowComm() );
            // Store a21[MC,* ] into its DistMatrix class and also store a copy
            // for the next iteration
            MemCopy
            ( a21_MC_STAR.Buffer(), &rowBroadcastBuffer[0], a21LocalHeight );
            // Store a21[MC,* ] into APan[MC,* ]
            const Int APan_MC_STAR_Offset = 
                APan_MC_STAR.LocalHeight()-a21LocalHeight;
            MemCopy
            ( APan_MC_STAR.Buffer(APan_MC_STAR_Offset,0), 
              &rowBroadcastBuffer[0],
              APan_MC_STAR.LocalHeight()-APan_MC_STAR_Offset );
            // Store tau
            tau = rowBroadcastBuffer[a21LocalHeight];
            
            a21_MR_STAR = a21_MC_STAR;
            // Store a21[MR,* ]
            const Int APan_MR_STAR_Offset = 
                APan_MR_STAR.LocalHeight()-a21_MR_STAR.LocalHeight();
            MemCopy
            ( APan_MR_STAR.Buffer(APan_MR_STAR_Offset,A00.Width()),
              a21_MR_STAR.Buffer(),
              APan_MR_STAR.LocalHeight()-APan_MR_STAR_Offset );
        }
        else
        {
            const Int a21LocalHeight = a21.LocalHeight();
            const Int w21LastLocalHeight = ACol.LocalHeight();
            std::vector<F> 
                rowBroadcastBuffer(a21LocalHeight+w21LastLocalHeight+1);
            if( thisIsMyCol ) 
            {
                // Pack the broadcast buffer with a21, w21Last, and tau
                MemCopy( &rowBroadcastBuffer[0], a21.Buffer(), a21LocalHeight );
                MemCopy
                ( &rowBroadcastBuffer[a21LocalHeight], 
                  &w21LastBuffer[0], w21LastLocalHeight );
                rowBroadcastBuffer[a21LocalHeight+w21LastLocalHeight] = tau;
            }
            // Broadcast a21, w21Last, and tau across the process row
            mpi::Broadcast
            ( &rowBroadcastBuffer[0], 
              a21LocalHeight+w21LastLocalHeight+1, 
              a21.RowAlignment(), g.RowComm() );
            // Store a21[MC,* ] into its DistMatrix class 
            MemCopy
            ( a21_MC_STAR.Buffer(), 
              &rowBroadcastBuffer[0], a21LocalHeight );
            // Store a21[MC,* ] into APan[MC,* ]
            const Int APan_MC_STAR_Offset = 
                APan_MC_STAR.LocalHeight()-a21LocalHeight;
            MemCopy
            ( APan_MC_STAR.Buffer(APan_MC_STAR_Offset,A00.Width()), 
              &rowBroadcastBuffer[0],
              APan_MC_STAR.LocalHeight()-APan_MC_STAR_Offset );
            // Store w21Last[MC,* ] into its DistMatrix class
            w21Last_MC_STAR.AlignWith( alpha11 );
            w21Last_MC_STAR.ResizeTo( a21.Height()+1, 1 );
            MemCopy
            ( w21Last_MC_STAR.Buffer(), 
              &rowBroadcastBuffer[a21LocalHeight], w21LastLocalHeight );
            // Store the bottom part of w21Last[MC,* ] into WB[MC,* ] and, 
            // if necessary, w21.
            const Int W_MC_STAR_Offset = 
                W_MC_STAR.LocalHeight()-w21LastLocalHeight;
            MemCopy
            ( W_MC_STAR.Buffer(W_MC_STAR_Offset,A00.Width()-1),
              &rowBroadcastBuffer[a21LocalHeight],
              W_MC_STAR.LocalHeight()-W_MC_STAR_Offset );
            if( g.Col() == w21Last.RowAlignment() )
            {
                MemCopy
                ( w21Last.Buffer(),
                  &rowBroadcastBuffer[a21LocalHeight], w21LastLocalHeight );
            }
            // Store tau
            tau = rowBroadcastBuffer[a21LocalHeight+w21LastLocalHeight];

            // Form a21[MR,* ] and w21Last[MR,* ] by combining the 
            // communications needed for taking a vector from 
            // [MC,* ] -> [MR,* ]: 
            //   local copy to [VC,* ], 
            //   Send/Recv to [VR,* ], 
            //   AllGather to [MR,* ]
            // We can combine the two by treating a21 as [0; a21] so that its 
            // alignments are equivalent to w21Last.

            // alpha11 is the top-left entry of the A22 from the last iteration,
            // so its alignments are identical
            const Int colAlignSource = alpha11.ColAlignment();
            const Int colAlignDest = alpha11.RowAlignment();
            const Int colShiftSource = alpha11.ColShift();
            const Int colShiftDest = alpha11.RowShift();

            const Int height = a21.Height()+1;
            const Int portionSize = mpi::Pad( 2*MaxLength(height,p) );

            const Int colShiftVRDest = Shift(g.VRRank(),colAlignDest,p);
            const Int colShiftVCSource = Shift(g.VCRank(),colAlignSource,p);
            const Int sendRankRM = 
                (g.VRRank()+(p+colShiftVCSource-colShiftVRDest))%p;
            const Int recvRankCM = 
                (g.VCRank()+(p+colShiftVRDest-colShiftVCSource))%p;
            const Int recvRankRM = 
                (recvRankCM/r)+c*(recvRankCM%r);

            std::vector<F> transposeBuffer( (r+1)*portionSize );
            F* sendBuf = &transposeBuffer[0];
            F* recvBuf = &transposeBuffer[r*portionSize];

            // (w21Last[VC,* ] <- w21Last[MC,* ]) and
            // ([0; a21][VC,* ] <- [0; a21][MC,* ])
            {
                // Pack the necessary portion of w21Last[MC,* ]
                const Int w21Shift = Shift(g.VCRank(),colAlignSource,p);
                const Int w21Offset = (w21Shift-colShiftSource)/r;
                const Int w21VCLocalHeight = Length(height,w21Shift,p);
                const F* w21Buffer = w21Last_MC_STAR.Buffer(w21Offset,0);
                for( Int i=0; i<w21VCLocalHeight; ++i )
                    sendBuf[i] = w21Buffer[i*c];
                
                // Pack the necessary portion of a21[MC,* ]
                const Int a21Shift = (w21Shift+p-1) % p;
                const Int a21Offset = (a21Shift-((colShiftSource+r-1)%r))/r;
                const Int a21VCLocalHeight = Length(height-1,a21Shift,p);
                const F* a21Buffer = a21_MC_STAR.Buffer(a21Offset,0);
                for( Int i=0; i<a21VCLocalHeight; ++i )
                    sendBuf[w21VCLocalHeight+i] = a21Buffer[i*c];
            }

            // [VR,* ] <- [VC,* ]
            mpi::SendRecv
            ( sendBuf, portionSize, sendRankRM, 
              recvBuf, portionSize, recvRankRM, g.VRComm() );

            // [MR,* ] <- [VR,* ]
            mpi::AllGather
            ( recvBuf, portionSize,
              sendBuf, portionSize, g.ColComm() );

            // Unpack
            w21Last_MR_STAR.AlignWith( alpha11 );
            w21Last_MR_STAR.ResizeTo( a21.Height()+1, 1 );
            for( Int k=0; k<r; ++k )
            {
                // Unpack into w21Last[MR,* ]
                const F* w21Data = &sendBuf[k*portionSize];
                const Int w21Shift = Shift(g.Col()+c*k,colAlignDest,p);
                const Int w21Offset = (w21Shift-colShiftDest) / c;
                const Int w21VCLocalHeight = Length(height,w21Shift,p);
                F* w21Buffer = w21Last_MR_STAR.Buffer(w21Offset,0);
                for( Int i=0; i<w21VCLocalHeight; ++i )
                    w21Buffer[i*r] = w21Data[i];

                // Unpack into a21[MR,* ]
                const F* a21Data = &sendBuf[k*portionSize+w21VCLocalHeight];
                const Int a21Shift = (w21Shift+p-1) % p;
                const Int a21Offset = (a21Shift-((colShiftDest+c-1)%c))/c;
                const Int a21VCLocalHeight = Length(height-1,a21Shift,p);
                F* a21Buffer = a21_MR_STAR.Buffer(a21Offset,0);
                for( Int i=0; i<a21VCLocalHeight; ++i )
                    a21Buffer[i*r] = a21Data[i];
            }
            // Store w21Last[MR,* ]
            const Int W_MR_STAR_Offset = 
                W_MR_STAR.LocalHeight()-w21Last_MR_STAR.LocalHeight();
            MemCopy
            ( W_MR_STAR.Buffer(W_MR_STAR_Offset,A00.Width()-1),
              w21Last_MR_STAR.Buffer(),
              W_MR_STAR.LocalHeight()-W_MR_STAR_Offset );
            // Store a21[MR,* ]
            const Int APan_MR_STAR_Offset = 
                APan_MR_STAR.LocalHeight()-a21_MR_STAR.LocalHeight();
            MemCopy
            ( APan_MR_STAR.Buffer(APan_MR_STAR_Offset,A00.Width()),
              a21_MR_STAR.Buffer(),
              APan_MR_STAR.LocalHeight()-APan_MR_STAR_Offset );

            // Update the portion of A22 that is in our current panel with 
            // w21Last and a21Last using two gers. We do not need their top 
            // entries. We trash the upper triangle of our panel of A since we 
            // are only doing slightly more work and we can replace it
            // afterwards.
            DistMatrix<F,MC,STAR> a21Last_MC_STAR_Bottom(g),
                                  w21Last_MC_STAR_Bottom(g);
            DistMatrix<F,MR,STAR> a21Last_MR_STAR_Bottom(g),
                                  w21Last_MR_STAR_Bottom(g);
            View
            ( a21Last_MC_STAR_Bottom,
              a21Last_MC_STAR, 1, 0, a21Last_MC_STAR.Height()-1, 1 );
            View
            ( w21Last_MC_STAR_Bottom,
              w21Last_MC_STAR, 1, 0, w21Last_MC_STAR.Height()-1, 1 );
            View
            ( a21Last_MR_STAR_Bottom,
              a21Last_MR_STAR, 1, 0, a21Last_MR_STAR.Height()-1, 1 );
            View
            ( w21Last_MR_STAR_Bottom,
              w21Last_MR_STAR, 1, 0, w21Last_MR_STAR.Height()-1, 1 );
            const F* a21_MC_STAR_Buffer = a21Last_MC_STAR_Bottom.Buffer();
            const F* w21_MC_STAR_Buffer = w21Last_MC_STAR_Bottom.Buffer();
            const F* a21_MR_STAR_Buffer = a21Last_MR_STAR_Bottom.Buffer();
            const F* w21_MR_STAR_Buffer = w21Last_MR_STAR_Bottom.Buffer();
            F* A22Buffer = A22.Buffer();
            const Int localHeight = W22.LocalHeight();
            const Int localWidth = W22.LocalWidth();
            const Int lDim = A22.LDim();
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    A22Buffer[iLocal+jLocal*lDim] -=
                        w21_MC_STAR_Buffer[iLocal]*
                        Conj(a21_MR_STAR_Buffer[jLocal]) +
                        a21_MC_STAR_Buffer[iLocal]*
                        Conj(w21_MR_STAR_Buffer[jLocal]);
        }

        // Form the local portions of (A22 a21) into p21[MC,* ] and q21[MR,* ]:
        //   p21[MC,* ] := tril(A22)[MC,MR] a21[MR,* ]
        //   q21[MR,* ] := tril(A22,-1)'[MR,MC] a21[MC,* ]
        Zero( p21_MC_STAR );
        Zeros( q21_MR_STAR, a21.Height(), 1 );
        internal::LocalSymvColAccumulateL
        ( F(1), A22, a21_MC_STAR, a21_MR_STAR, p21_MC_STAR, q21_MR_STAR, true );

        Zeros( x01_MR_STAR, W20B.Width(), 1 );
        Zeros( y01_MR_STAR, W20B.Width(), 1 );
        LocalGemv( ADJOINT, F(1), W20B, a21B_MC_STAR, F(0), x01_MR_STAR );
        LocalGemv( ADJOINT, F(1), A20B, a21B_MC_STAR, F(0), y01_MR_STAR );

        // Combine the AllReduce column summations of x01[MR,* ], y01[MR,* ],
        // and q21[MR,* ]
        {
            const Int x01LocalHeight = x01_MR_STAR.LocalHeight();
            const Int q21LocalHeight = q21_MR_STAR.LocalHeight();
            std::vector<F> colSumSendBuffer(2*x01LocalHeight+q21LocalHeight),
                           colSumRecvBuffer(2*x01LocalHeight+q21LocalHeight);
            MemCopy
            ( &colSumSendBuffer[0], x01_MR_STAR.Buffer(), x01LocalHeight );
            MemCopy
            ( &colSumSendBuffer[x01LocalHeight],
              y01_MR_STAR.Buffer(), x01LocalHeight );
            MemCopy
            ( &colSumSendBuffer[2*x01LocalHeight],
              q21_MR_STAR.Buffer(), q21LocalHeight );
            mpi::AllReduce
            ( &colSumSendBuffer[0], &colSumRecvBuffer[0],
              2*x01LocalHeight+q21LocalHeight, g.ColComm() );
            MemCopy
            ( x01_MR_STAR.Buffer(), 
              &colSumRecvBuffer[0], x01LocalHeight );
            MemCopy
            ( y01_MR_STAR.Buffer(), 
              &colSumRecvBuffer[x01LocalHeight], x01LocalHeight );
            MemCopy
            ( q21_MR_STAR.Buffer(), 
              &colSumRecvBuffer[2*x01LocalHeight], q21LocalHeight );
        }

        LocalGemv( NORMAL, F(-1), A20B, x01_MR_STAR, F(1), p21B_MC_STAR );
        LocalGemv( NORMAL, F(-1), W20B, y01_MR_STAR, F(1), p21B_MC_STAR );

        if( W22.Width() > 0 )
        {
            // This is not the last iteration of the panel factorization, 
            // combine the Reduce to one of p21[MC,* ] with the redistribution 
            // of q21[MR,* ] -> q21[MC,MR] to the next process column.
            const Int localHeight = p21_MC_STAR.LocalHeight();
            std::vector<F> reduceToOneSendBuffer(2*localHeight),
                           reduceToOneRecvBuffer(2*localHeight);

            // Pack p21[MC,* ]
            MemCopy
            ( &reduceToOneSendBuffer[0], p21_MC_STAR.Buffer(), localHeight );

            // Fill in contributions to q21[MC,MR] from q21[MR,* ]
            const bool contributing = 
                ( q21_MR_STAR.ColShift() % g.GCD() ==
                  p21_MC_STAR.ColShift() % g.GCD() );
            if( contributing )
            {
                if( r == c )
                {
                    MemCopy
                    ( &reduceToOneSendBuffer[localHeight],
                      q21_MR_STAR.Buffer(), localHeight );
                }
                else
                {
                    // Zero the entire buffer first
                    MemZero( &reduceToOneSendBuffer[localHeight], localHeight );
                    // Fill in the entries that we contribute to.
                    // We seek to find the minimum s in N such that
                    //   s*c = a0-b0 (mod r)
                    // where a0 is the column shift of MC, b0 is the row shift
                    // of MR, and s is our first local entry of MR that will 
                    // contribute to MC. I cannot think of an O(1) method, so
                    // I will instead use the worst-case O(lcm(c,r)/c) method.
                    const Int sourcePeriod = g.LCM() / c;
                    const Int targetPeriod = g.LCM() / r;
                    const Int a0 = p21_MC_STAR.ColShift();
                    const Int b0 = q21_MR_STAR.ColShift();

                    Int sourceStart = 0;
                    const Int f = (r+a0-b0) % r;
                    for( Int s=0; s<sourcePeriod; ++s )
                    {
                        if( (s*c) % r == f )
                        {
                            sourceStart = s;
                            break;
                        }
                    }

                    const Int globalShift = b0+sourceStart*c;
                    const Int targetStart = (globalShift-a0)/r;
                    const Int localLength =
                        Length(localHeight,targetStart,targetPeriod);
                    const F* q21_MR_STAR_Buffer = q21_MR_STAR.Buffer();
                    const Int offset = localHeight + targetStart;
                    for( Int i=0; i<localLength; ++i )                        
                        reduceToOneSendBuffer[offset+i*targetPeriod] = 
                            q21_MR_STAR_Buffer[sourceStart+i*sourcePeriod];
                }
            }
            else
                MemZero( &reduceToOneSendBuffer[localHeight], localHeight );

            const Int nextProcessRow = (alpha11.ColAlignment()+1) % r;
            const Int nextProcessCol = (alpha11.RowAlignment()+1) % c;
            mpi::Reduce
            ( &reduceToOneSendBuffer[0], &reduceToOneRecvBuffer[0],
              2*localHeight, nextProcessCol, g.RowComm() );
            if( g.Col() == nextProcessCol )
            {
                // Combine the second half into the first half        
                for( Int i=0; i<localHeight; ++i )
                    reduceToOneRecvBuffer[i] +=
                        reduceToOneRecvBuffer[i+localHeight];

                // Finish computing w21. During its computation, ensure that 
                // every process has a copy of the first element of the w21.
                // We know a priori that the first element of a21 is one.
                const F* a21_MC_STAR_Buffer = a21_MC_STAR.Buffer();
                F myDotProduct = blas::Dot
                    ( localHeight, &reduceToOneRecvBuffer[0],   1, 
                                   &a21_MC_STAR_Buffer[0], 1 );
                F sendBuffer[2], recvBuffer[2];
                sendBuffer[0] = myDotProduct;
                sendBuffer[1] = ( g.Row()==nextProcessRow ? 
                                  reduceToOneRecvBuffer[0] : 0 );
                mpi::AllReduce( sendBuffer, recvBuffer, 2, g.ColComm() );
                F dotProduct = recvBuffer[0];

                // Set up for the next iteration by filling in the values for:
                // - w21LastBuffer
                // - w21LastFirstEntry
                F scale = dotProduct*Conj(tau) / F(2);
                for( Int i=0; i<localHeight; ++i )
                    w21LastBuffer[i] = tau*
                        ( reduceToOneRecvBuffer[i]-
                          scale*a21_MC_STAR_Buffer[i] );
                w21LastFirstEntry = tau*( recvBuffer[1]-scale );
            }
        }
        else
        {
            // This is the last iteration, our last task is to finish forming
            // w21[MC,* ] and w21[MR,* ] so that we may place them into W[MC,* ]
            // and W[MR,* ]
            const Int localHeight = p21_MC_STAR.LocalHeight();
            std::vector<F> allReduceSendBuffer(2*localHeight),
                           allReduceRecvBuffer(2*localHeight);

            // Pack p21[MC,* ]
            MemCopy
            ( &allReduceSendBuffer[0], p21_MC_STAR.Buffer(), localHeight );

            // Fill in contributions to q21[MC,* ] from q21[MR,* ]
            const bool contributing = 
                ( q21_MR_STAR.ColShift() % g.GCD() ==
                  p21_MC_STAR.ColShift() % g.GCD() );
            if( contributing )
            {
                if( r == c )
                {
                    MemCopy
                    ( &allReduceSendBuffer[localHeight],
                      q21_MR_STAR.Buffer(), localHeight );
                }
                else
                {
                    // Zero the entire buffer first
                    MemZero( &allReduceSendBuffer[localHeight], localHeight );
                    // Fill in the entries that we contribute to.
                    // We seek to find the minimum s in N such that
                    //   s*c = a0-b0 (mod r)
                    // where a0 is the column shift of MC, b0 is the row shift
                    // of MR, and s is our first local entry of MR that will 
                    // contribute to MC. I cannot think of an O(1) method, so
                    // I will instead use the worst-case O(lcm(c,r)/c) method.
                    const Int sourcePeriod = g.LCM() / c;
                    const Int targetPeriod = g.LCM() / r;
                    const Int a0 = p21_MC_STAR.ColShift();
                    const Int b0 = q21_MR_STAR.ColShift();

                    Int sourceStart = 0;
                    const Int f = (r+a0-b0) % r;
                    for( Int s=0; s<sourcePeriod; ++s )
                    {
                        if( (s*c) % r == f )
                        {
                            sourceStart = s;
                            break;
                        }
                    }

                    const Int globalShift = b0+sourceStart*c;
                    const Int targetStart = (globalShift-a0)/r;
                    const Int localLength = 
                        Length(localHeight,targetStart,targetPeriod);
                    const F* q21_MR_STAR_Buffer = q21_MR_STAR.Buffer();
                    const Int offset = localHeight + targetStart;
                    for( Int i=0; i<localLength; ++i )
                        allReduceSendBuffer[offset+i*targetPeriod] = 
                            q21_MR_STAR_Buffer[sourceStart+i*sourcePeriod];
                }
            }
            else
                MemZero( &allReduceSendBuffer[localHeight], localHeight );

            mpi::AllReduce
            ( &allReduceSendBuffer[0], &allReduceRecvBuffer[0],
              2*localHeight, g.RowComm() );

            // Combine the second half into the first half        
            for( Int i=0; i<localHeight; ++i )
                allReduceRecvBuffer[i] += allReduceRecvBuffer[i+localHeight];
 
            // Finish computing w21.
            const F* a21_MC_STAR_Buffer = a21_MC_STAR.Buffer();
            F myDotProduct = blas::Dot
                ( localHeight, &allReduceRecvBuffer[0], 1, 
                               a21_MC_STAR_Buffer, 1 );
            const F dotProduct = mpi::AllReduce( myDotProduct, g.ColComm() );

            // Grab views into W[MC,* ] and W[MR,* ]
            DistMatrix<F,MC,STAR> w21_MC_STAR(g);
            DistMatrix<F,MR,STAR> w21_MR_STAR(g);
            View
            ( w21_MC_STAR, 
              W_MC_STAR, W00.Height()+1, W00.Width(), w21.Height(), 1 );
            View
            ( w21_MR_STAR,
              W_MR_STAR, W00.Height()+1, W00.Width(), w21.Height(), 1 );

            // Store w21[MC,* ]
            F scale = dotProduct*Conj(tau) / F(2);
            F* w21_MC_STAR_Buffer = w21_MC_STAR.Buffer();
            for( Int i=0; i<localHeight; ++i )
                w21_MC_STAR_Buffer[i] = 
                    tau*( allReduceRecvBuffer[i]-scale*a21_MC_STAR_Buffer[i] );

            // Form w21[MR,* ]
            w21_MR_STAR = w21_MC_STAR;
        }
        //--------------------------------------------------------------------//

        SlidePartitionDown
        ( tT,  t0,
               tau1,
         /**/ /****/
          tB,  t2 );
       
        SlidePartitionDown
        ( eT,  e0,
               epsilon1,
         /**/ /********/
          eB,  e2 );

        SlidePartitionDownDiagonal
        ( WTL, /**/ WTR,  W00, w01,     /**/ W02,
               /**/       w10, omega11, /**/ w12,
         /*************/ /**********************/
          WBL, /**/ WBR,  W20, w21,     /**/ W22 );
 
        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );

        firstIteration = false;
    }

    // View the portion of A that e is the subdiagonal of, then place e into it
    DistMatrix<F> expandedATL(g);
    View( expandedATL, A, 0, 0, panelSize+1, panelSize+1 );
    expandedATL.SetRealPartOfDiagonal( e, -1 );
}

} // namespace hermitian_tridiag
} // namespace elem

#endif // ifndef LAPACK_HERMITIANTRIDIAG_PANELL_HPP
