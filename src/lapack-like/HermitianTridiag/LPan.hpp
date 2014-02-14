/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef ELEM_HERMITIANTRIDIAG_LPAN_HPP
#define ELEM_HERMITIANTRIDIAG_LPAN_HPP

#include ELEM_ZERO_INC
#include ELEM_GEMV_INC
#include ELEM_SYMV_INC
#include ELEM_REFLECTOR_INC

namespace elem {
namespace herm_tridiag {

template<typename F>
void LPan
( DistMatrix<F>& A,
  DistMatrix<F>& W,
  DistMatrix<F,MD,STAR>& t,
  DistMatrix<F,MC,STAR>& APan_MC_STAR, 
  DistMatrix<F,MR,STAR>& APan_MR_STAR,
  DistMatrix<F,MC,STAR>& W_MC_STAR,
  DistMatrix<F,MR,STAR>& W_MR_STAR )
{
    const Int n = A.Height();
    const Int nW = W.Width();
    DEBUG_ONLY(
        CallStackEntry cse("herm_tridiag::LPan");
        if( A.Grid() != W.Grid() || W.Grid() != t.Grid() )
            LogicError("A, W, and t must be distributed over the same grid.");
        if( n != A.Width() )
            LogicError("A must be square");
        if( n != W.Height() )
            LogicError("A and W must be the same height");
        if( n <= nW )
            LogicError("W must be a column panel");
        if( W.ColAlign() != A.ColAlign() || 
            W.RowAlign() != A.RowAlign() )
            LogicError("W and A must be aligned");
        if( t.Height() != nW || t.Width() != 1 )
            LogicError
            ("t must be a column vector of the same length as W's width");
        if( !A.DiagonalAlignedWith(t,-1) )
            LogicError("t is not aligned with A's subdiagonal.");
    )
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();

    // Create a distributed matrix for storing the subdiagonal
    DistMatrix<Real,MD,STAR> e(g);
    e.SetRoot( A.DiagonalRoot(-1) );
    e.AlignCols( A.DiagonalAlign(-1) );
    e.Resize( nW, 1 );

    std::vector<F> w21LastBuffer(n/r+1);
    DistMatrix<F> w21Last(g);
    DistMatrix<F,MC,STAR> a21_MC_STAR(g), a21B_MC_STAR(g),
                          p21_MC_STAR(g), p21B_MC_STAR(g),
                          a21Last_MC_STAR(g), w21Last_MC_STAR(g);
    DistMatrix<F,MR,STAR> a21_MR_STAR(g), q21_MR_STAR(g),
                          x01_MR_STAR(g), y01_MR_STAR(g),
                          a21Last_MR_STAR(g), w21Last_MR_STAR(g);

    F tau = 0;
    F w21LastFirstEntry = 0;
    for( Int k=0; k<nW; ++k )
    {
        if( k > 0 )
        {
            a21Last_MC_STAR = ViewRange( APan_MC_STAR, k, k-1, n, k );
            a21Last_MR_STAR = ViewRange( APan_MR_STAR, k, k-1, n, k );
            w21Last         = ViewRange( W,            k, k-1, n, k );
        }

        auto A00      = ViewRange( A, 0,   0,   k,   k   );
        auto W00      = ViewRange( W, 0,   0,   k,   k   );
        auto alpha11  = ViewRange( A, k,   k,   k+1, k+1 );
        auto a21      = ViewRange( A, k+1, k,   n,   k+1 );
        auto alpha21T = ViewRange( A, k+1, k,   k+2, k+1 );
        auto a21B     = ViewRange( A, k+2, k,   n,   k+1 );
        auto A22      = ViewRange( A, k+1, k+1, n,   n   );
        auto W22      = ViewRange( W, k+1, k+1, n,   nW  );
        auto ACol     = ViewRange( A, k,   k,   n,   k+1 );
        auto WCol     = ViewRange( W, k,   k,   n,   k+1 );
        auto A20B     = ViewRange( A, nW,  0,   n,   k   );
        auto W20B     = ViewRange( W, nW,  0,   n,   k   );
        auto tau1     = View( t, k, 0, 1, 1 );
        auto epsilon1 = View( e, k, 0, 1, 1 );
           
        a21_MC_STAR.AlignWith( A22 );
        a21_MR_STAR.AlignWith( A22 );
        p21_MC_STAR.AlignWith( A22 );
        a21_MC_STAR.Resize( n-(k+1), 1 );
        a21_MR_STAR.Resize( n-(k+1), 1 );
        p21_MC_STAR.Resize( n-(k+1), 1 );

        // View the portions of a21[MC,* ] and p21[MC,* ] below the current
        // panel's square
        auto a21B_MC_STAR = View( a21_MC_STAR, nW-(k+1), 0, n-nW, 1 );
        auto p21B_MC_STAR = View( p21_MC_STAR, nW-(k+1), 0, n-nW, 1 );

        const bool thisIsMyCol = ( g.Col() == alpha11.RowAlign() );
        if( thisIsMyCol )
        {
            if( k > 0 )
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
            if( g.Row() == alpha21T.ColAlign() )
                tau1.SetLocal(0,0,tau);
        }

        // Store the subdiagonal value and turn a21 into a proper scaled 
        // reflector by explicitly placing the implicit one in its first entry.
        alpha21T.GetRealPartOfDiagonal( epsilon1 );
        alpha21T.Set( 0, 0, F(1) );

        // If this is the first iteration, have each member of the owning 
        // process column broadcast tau and a21 within its process row. 
        // Otherwise, also add w21 into the broadcast.
        if( k == 0 )
        {
            const Int a21LocalHeight = a21.LocalHeight();
            std::vector<F> rowBroadcastBuffer(a21LocalHeight+1);
            if( thisIsMyCol )
            {
                // Pack the broadcast buffer with a21 and tau
                MemCopy
                ( rowBroadcastBuffer.data(), a21.Buffer(), a21LocalHeight );
                rowBroadcastBuffer[a21LocalHeight] = tau;
            }
            // Broadcast a21 and tau across the process row
            mpi::Broadcast
            ( rowBroadcastBuffer.data(), 
              a21LocalHeight+1, a21.RowAlign(), g.RowComm() );
            // Store a21[MC,* ] into its DistMatrix class and also store a copy
            // for the next iteration
            MemCopy
            ( a21_MC_STAR.Buffer(), rowBroadcastBuffer.data(), a21LocalHeight );
            // Store a21[MC,* ] into APan[MC,* ]
            const Int APan_MC_STAR_Offset = 
                APan_MC_STAR.LocalHeight()-a21LocalHeight;
            MemCopy
            ( APan_MC_STAR.Buffer(APan_MC_STAR_Offset,0), 
              rowBroadcastBuffer.data(),
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
                MemCopy
                ( rowBroadcastBuffer.data(), a21.Buffer(), a21LocalHeight );
                MemCopy
                ( &rowBroadcastBuffer[a21LocalHeight], 
                  w21LastBuffer.data(), w21LastLocalHeight );
                rowBroadcastBuffer[a21LocalHeight+w21LastLocalHeight] = tau;
            }
            // Broadcast a21, w21Last, and tau across the process row
            mpi::Broadcast
            ( rowBroadcastBuffer.data(), 
              a21LocalHeight+w21LastLocalHeight+1, 
              a21.RowAlign(), g.RowComm() );
            // Store a21[MC,* ] into its DistMatrix class 
            MemCopy
            ( a21_MC_STAR.Buffer(), 
              rowBroadcastBuffer.data(), a21LocalHeight );
            // Store a21[MC,* ] into APan[MC,* ]
            const Int APan_MC_STAR_Offset = 
                APan_MC_STAR.LocalHeight()-a21LocalHeight;
            MemCopy
            ( APan_MC_STAR.Buffer(APan_MC_STAR_Offset,A00.Width()), 
              rowBroadcastBuffer.data(),
              APan_MC_STAR.LocalHeight()-APan_MC_STAR_Offset );
            // Store w21Last[MC,* ] into its DistMatrix class
            w21Last_MC_STAR.AlignWith( alpha11 );
            w21Last_MC_STAR.Resize( n-k, 1 );
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
            if( g.Col() == w21Last.RowAlign() )
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
            const Int colAlignSource = alpha11.ColAlign();
            const Int colAlignDest = alpha11.RowAlign();
            const Int colShiftSource = alpha11.ColShift();
            const Int colShiftDest = alpha11.RowShift();

            const Int height = n-k;
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
            w21Last_MR_STAR.Resize( n-k, 1 );
            for( Int row=0; row<r; ++row )
            {
                // Unpack into w21Last[MR,* ]
                const F* w21Data = &sendBuf[row*portionSize];
                const Int w21Shift = Shift(g.Col()+c*row,colAlignDest,p);
                const Int w21Offset = (w21Shift-colShiftDest) / c;
                const Int w21VCLocalHeight = Length(height,w21Shift,p);
                F* w21Buffer = w21Last_MR_STAR.Buffer(w21Offset,0);
                for( Int i=0; i<w21VCLocalHeight; ++i )
                    w21Buffer[i*r] = w21Data[i];

                // Unpack into a21[MR,* ]
                const F* a21Data = &sendBuf[row*portionSize+w21VCLocalHeight];
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
            auto a21Last_MC_STAR_Bottom = 
                ViewRange( a21Last_MC_STAR, 1, 0, n-k, 1 );
            auto w21Last_MC_STAR_Bottom =  
                ViewRange( w21Last_MC_STAR, 1, 0, n-k, 1 );
            auto a21Last_MR_STAR_Bottom = 
                ViewRange( a21Last_MR_STAR, 1, 0, n-k, 1 );
            auto w21Last_MR_STAR_Bottom =
                ViewRange( w21Last_MR_STAR, 1, 0, n-k, 1 );
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
        q21_MR_STAR.AlignWith( A22 );
        Zeros( q21_MR_STAR, a21.Height(), 1 );
        internal::LocalSymvColAccumulateL
        ( F(1), A22, a21_MC_STAR, a21_MR_STAR, p21_MC_STAR, q21_MR_STAR, true );

        x01_MR_STAR.AlignWith( W20B );
        y01_MR_STAR.AlignWith( W20B );
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
            ( colSumSendBuffer.data(), x01_MR_STAR.Buffer(), x01LocalHeight );
            MemCopy
            ( &colSumSendBuffer[x01LocalHeight],
              y01_MR_STAR.Buffer(), x01LocalHeight );
            MemCopy
            ( &colSumSendBuffer[2*x01LocalHeight],
              q21_MR_STAR.Buffer(), q21LocalHeight );
            mpi::AllReduce
            ( colSumSendBuffer.data(), colSumRecvBuffer.data(),
              2*x01LocalHeight+q21LocalHeight, g.ColComm() );
            MemCopy
            ( x01_MR_STAR.Buffer(), 
              colSumRecvBuffer.data(), x01LocalHeight );
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
            ( reduceToOneSendBuffer.data(), p21_MC_STAR.Buffer(), localHeight );

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

            const Int nextProcessRow = (alpha11.ColAlign()+1) % r;
            const Int nextProcessCol = (alpha11.RowAlign()+1) % c;
            mpi::Reduce
            ( reduceToOneSendBuffer.data(), reduceToOneRecvBuffer.data(),
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
                    ( localHeight, reduceToOneRecvBuffer.data(), 1, 
                                   a21_MC_STAR_Buffer,           1 );
                F sendBuffer[2], recvBuffer[2];
                sendBuffer[0] = myDotProduct;
                sendBuffer[1] = ( g.Row()==nextProcessRow ? 
                                  reduceToOneRecvBuffer[0] : 0 );
                mpi::AllReduce( sendBuffer, recvBuffer, 2, g.ColComm() );
                F dotProduct = recvBuffer[0];

                // Set up for the next iteration by filling in the values for:
                // - w21LastBuffer
                // - w21LastFirstEntry
                F scale = dotProduct*tau / F(2);
                for( Int i=0; i<localHeight; ++i )
                    w21LastBuffer[i] = Conj(tau)*
                        ( reduceToOneRecvBuffer[i]-
                          scale*a21_MC_STAR_Buffer[i] );
                w21LastFirstEntry = Conj(tau)*( recvBuffer[1]-scale );
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
            ( allReduceSendBuffer.data(), p21_MC_STAR.Buffer(), localHeight );

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
            ( allReduceSendBuffer.data(), allReduceRecvBuffer.data(),
              2*localHeight, g.RowComm() );

            // Combine the second half into the first half        
            for( Int i=0; i<localHeight; ++i )
                allReduceRecvBuffer[i] += allReduceRecvBuffer[i+localHeight];
 
            // Finish computing w21.
            const F* a21_MC_STAR_Buffer = a21_MC_STAR.Buffer();
            F myDotProduct = blas::Dot
                ( localHeight, allReduceRecvBuffer.data(), 1, 
                               a21_MC_STAR_Buffer,         1 );
            const F dotProduct = mpi::AllReduce( myDotProduct, g.ColComm() );

            // Grab views into W[MC,* ] and W[MR,* ]
            auto w21_MC_STAR = ViewRange( W_MC_STAR, k+1, k, n, k+1 );
            auto w21_MR_STAR = ViewRange( W_MR_STAR, k+1, k, n, k+1 );

            // Store w21[MC,* ]
            F scale = dotProduct*tau / F(2);
            F* w21_MC_STAR_Buffer = w21_MC_STAR.Buffer();
            for( Int i=0; i<localHeight; ++i )
                w21_MC_STAR_Buffer[i] = Conj(tau)*
                    ( allReduceRecvBuffer[i]-scale*a21_MC_STAR_Buffer[i] );

            // Form w21[MR,* ]
            w21_MR_STAR = w21_MC_STAR;
        }
    }

    // View the portion of A that e is the subdiagonal of, then place e into it
    auto expandedATL = View( A, 0, 0, nW+1, nW+1 );
    expandedATL.SetRealPartOfDiagonal( e, -1 );
}

} // namespace herm_tridiag
} // namespace elem

#endif // ifndef ELEM_HERMITIANTRIDIAG_LPAN_HPP
