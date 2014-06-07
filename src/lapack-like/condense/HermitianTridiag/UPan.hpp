/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HERMITIANTRIDIAG_UPAN_HPP
#define EL_HERMITIANTRIDIAG_UPAN_HPP

namespace El {
namespace herm_tridiag {

template<typename F>
void UPan
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
        CallStackEntry cse("herm_tridiag::UPan");
        if( A.Grid() != W.Grid() || W.Grid() != t.Grid() )
            LogicError("A, W, and t must be distributed over the same grid.");
        if( n != A.Width() )
            LogicError("A must be square.");
        if( n != W.Height() )
            LogicError( "A and W must be the same height.");
        if( n <= nW )
            LogicError("W must be a column panel.");
        if( t.Height() != nW || t.Width() != 1 )
            LogicError
            ("t must be a column vector of the same length as W's width.");
    )
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int off = n-nW;

    // Create a distributed matrix for storing the superdiagonal
    auto expandedABR = ViewRange( A, off-1, off-1, n, n );
    DistMatrix<Real,MD,STAR> e(g);
    e.SetRoot( expandedABR.DiagonalRoot(1) );
    e.AlignCols( expandedABR.DiagonalAlign(1) );
    e.Resize( nW, 1 );

    std::vector<F> w01LastBuffer(n/r+1);
    DistMatrix<F> w01Last(g);
    DistMatrix<F,MC,STAR> a01_MC_STAR(g), p01_MC_STAR(g),
                          a01Last_MC_STAR(g), w01Last_MC_STAR(g);
    DistMatrix<F,MR,STAR> a01_MR_STAR(g), q01_MR_STAR(g),
                          x21_MR_STAR(g), y21_MR_STAR(g),
                          a01Last_MR_STAR(g), w01Last_MR_STAR(g);

    F tau = 0;
    F w01LastBottomEntry = 0;
    for( Int k=nW-1; k>=0; --k )
    {
        const Int kA = k+off;
        const bool firstIteration = ( k == nW-1 );
        if( !firstIteration ) 
        {
            a01Last_MC_STAR = View( APan_MC_STAR, 0, k+1, kA+1, 1 );
            a01Last_MR_STAR = View( APan_MR_STAR, 0, k+1, kA+1, 1 );
            w01Last         = View( W,            0, k+1, kA+1, 1 );
        }

        auto A00      = ViewRange( A, 0,    0,    kA,   kA   );
        auto W00      = ViewRange( W, 0,    0,    kA,   k    );
        auto A00Pan   = ViewRange( A, 0,    off,  kA,   kA   );
        auto a01      = ViewRange( A, 0,    kA,   kA,   kA+1 );
        auto a01T     = ViewRange( A, 0,    kA,   kA-1, kA+1 );
        auto alpha01B = ViewRange( A, kA-1, kA,   kA,   kA+1 );
        auto alpha11  = ViewRange( A, kA,   kA,   kA+1, kA+1 );
        auto ACol     = ViewRange( A, 0,    kA,   kA+1, kA+1 );
        auto WCol     = ViewRange( W, 0,    k,    kA+1, k+1  );
        auto A02      = ViewRange( A, 0,    kA+1, kA,   n    );
        auto A02T     = ViewRange( A, 0,    kA+1, off,  n    );
        auto W02T     = ViewRange( W, 0,    k+1,  off,  nW   );
        auto tau1     = View( t, k, 0, 1, 1 );
        auto epsilon1 = View( e, k, 0, 1, 1 );

        a01_MC_STAR.AlignWith( A00 );
        a01_MR_STAR.AlignWith( A00 );
        p01_MC_STAR.AlignWith( A00 );
        a01_MC_STAR.Resize( kA, 1 );
        a01_MR_STAR.Resize( kA, 1 );
        p01_MC_STAR.Resize( kA, 1 );

        // View the portions of a01[MC,* ] and p01[MC,* ] above the current
        // panel's square
        auto a01T_MC_STAR = View( a01_MC_STAR, 0, 0, off, 1 );
        auto p01T_MC_STAR = View( p01_MC_STAR, 0, 0, off, 1 );

        const bool thisIsMyCol = ( g.Col() == alpha11.RowAlign() );
        if( thisIsMyCol )
        {
            if( !firstIteration )
            {
                // Finish updating the current column with two axpy's
                const Int AColLocalHeight = ACol.LocalHeight();
                F* AColBuffer = ACol.Buffer();
                const F* a01Last_MC_STAR_Buffer = a01Last_MC_STAR.Buffer();
                for( Int i=0; i<AColLocalHeight; ++i )
                    AColBuffer[i] -=
                        w01LastBuffer[i] + 
                        a01Last_MC_STAR_Buffer[i]*Conj(w01LastBottomEntry);
            }
            // Compute the Householder reflector
            tau = reflector::Col( alpha01B, a01T );
            if( g.Row() == alpha01B.ColAlign() )
                tau1.SetLocal(0,0,tau);
        }

        // Store the subdiagonal value and turn a01 into a proper scaled 
        // reflector by explicitly placing the implicit one in its first entry.
        alpha01B.GetRealPartOfDiagonal( epsilon1 );
        alpha01B.Set( 0, 0, F(1) );

        // If this is the first iteration, have each member of the owning 
        // process column broadcast tau and a01 within its process row. 
        // Otherwise, also add w01 into the broadcast.
        if( firstIteration )
        {
            const Int a01LocalHeight = a01.LocalHeight();
            std::vector<F> rowBroadcastBuffer(a01LocalHeight+1);
            if( thisIsMyCol )
            {
                // Pack the broadcast buffer with a01 and tau
                MemCopy
                ( rowBroadcastBuffer.data(), a01.Buffer(), a01LocalHeight );
                rowBroadcastBuffer[a01LocalHeight] = tau;
            }
            // Broadcast a01 and tau across the process row
            mpi::Broadcast
            ( rowBroadcastBuffer.data(), 
              a01LocalHeight+1, a01.RowAlign(), g.RowComm() );
            // Store a01[MC,* ] into its DistMatrix class and also store a copy
            // for the next iteration
            MemCopy
            ( a01_MC_STAR.Buffer(), rowBroadcastBuffer.data(), a01LocalHeight );
            // Store a01[MC,* ] into APan[MC,* ]
            MemCopy
            ( APan_MC_STAR.Buffer(0,W00.Width()), 
              rowBroadcastBuffer.data(), a01LocalHeight );
            // Store tau
            tau = rowBroadcastBuffer[a01LocalHeight];
            
            a01_MR_STAR = a01_MC_STAR;
            // Store a01[MR,* ]
            MemCopy
            ( APan_MR_STAR.Buffer(0,W00.Width()),
              a01_MR_STAR.Buffer(),
              a01_MR_STAR.LocalHeight() );
        }
        else
        {
            const Int a01LocalHeight = a01.LocalHeight();
            const Int w01LastLocalHeight = ACol.LocalHeight();
            std::vector<F> 
                rowBroadcastBuffer(a01LocalHeight+w01LastLocalHeight+1);
            if( thisIsMyCol ) 
            {
                // Pack the broadcast buffer with a01, w01Last, and tau
                MemCopy
                ( rowBroadcastBuffer.data(), a01.Buffer(), a01LocalHeight );
                MemCopy
                ( &rowBroadcastBuffer[a01LocalHeight], 
                  w01LastBuffer.data(), w01LastLocalHeight );
                rowBroadcastBuffer[a01LocalHeight+w01LastLocalHeight] = tau;
            }
            // Broadcast a01, w01Last, and tau across the process row
            mpi::Broadcast
            ( rowBroadcastBuffer.data(), 
              a01LocalHeight+w01LastLocalHeight+1, 
              a01.RowAlign(), g.RowComm() );
            // Store a01[MC,* ] into its DistMatrix class 
            MemCopy
            ( a01_MC_STAR.Buffer(), rowBroadcastBuffer.data(), a01LocalHeight );
            // Store a01[MC,* ] into APan[MC,* ]
            MemCopy
            ( APan_MC_STAR.Buffer(0,W00.Width()), 
              rowBroadcastBuffer.data(), a01LocalHeight );
            // Store w01Last[MC,* ] into its DistMatrix class
            w01Last_MC_STAR.AlignWith( A00 );
            w01Last_MC_STAR.Resize( a01.Height()+1, 1 );
            MemCopy
            ( w01Last_MC_STAR.Buffer(), 
              &rowBroadcastBuffer[a01LocalHeight], w01LastLocalHeight );
            // Store the bottom part of w01Last[MC,* ] into WB[MC,* ] and, 
            // if necessary, w01.
            MemCopy
            ( W_MC_STAR.Buffer(0,W00.Width()+1),
              &rowBroadcastBuffer[a01LocalHeight], w01LastLocalHeight );
            if( g.Col() == w01Last.RowAlign() )
            {
                MemCopy
                ( w01Last.Buffer(),
                  &rowBroadcastBuffer[a01LocalHeight], w01LastLocalHeight );
            }
            // Store tau
            tau = rowBroadcastBuffer[a01LocalHeight+w01LastLocalHeight];

            // Form a01[MR,* ] and w01Last[MR,* ] by combining the 
            // communications needed for taking a vector from 
            // [MC,* ] -> [MR,* ]: 
            //   local copy to [VC,* ], 
            //   Send/Recv to [VR,* ], 
            //   AllGather to [MR,* ]
            // We can combine the two by treating a01 as [a01; 0] 

            const Int colAlignSource = A00.ColAlign();
            const Int colAlignDest = A00.RowAlign();
            const Int colShiftSource = A00.ColShift();
            const Int colShiftDest = A00.RowShift();

            const Int height = a01.Height()+1;
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

            // (w01Last[VC,* ] <- w01Last[MC,* ]) and
            // ([a01; 0][VC,* ] <- [a01; 0][MC,* ])
            {
                // Pack the necessary portion of w01Last[MC,* ]
                const Int shift = Shift(g.VCRank(),colAlignSource,p);
                const Int offset = (shift-colShiftSource)/r;
                const Int w01VCLocalHeight = Length(height,shift,p);
                const F* w01Buffer = w01Last_MC_STAR.Buffer(offset,0);
                for( Int i=0; i<w01VCLocalHeight; ++i )
                    sendBuf[i] = w01Buffer[i*c];
                
                // Pack the necessary portion of a01[MC,* ]
                const Int a01VCLocalHeight = Length(height-1,shift,p);
                const F* a01Buffer = a01_MC_STAR.Buffer(offset,0);
                for( Int i=0; i<a01VCLocalHeight; ++i )
                    sendBuf[w01VCLocalHeight+i] = a01Buffer[i*c];
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
            w01Last_MR_STAR.AlignWith( A00 );
            w01Last_MR_STAR.Resize( a01.Height()+1, 1 );
            for( Int row=0; row<r; ++row )
            {
                // Unpack into w01Last[MR,* ]
                const F* w01Data = &sendBuf[row*portionSize];
                const Int shift = Shift(g.Col()+c*row,colAlignDest,p);
                const Int offset = (shift-colShiftDest) / c;
                const Int w01VCLocalHeight = Length(height,shift,p);
                F* w01Buffer = w01Last_MR_STAR.Buffer(offset,0);
                for( Int i=0; i<w01VCLocalHeight; ++i )
                    w01Buffer[i*r] = w01Data[i];

                // Unpack into a01[MR,* ]
                const F* a01Data = &sendBuf[row*portionSize+w01VCLocalHeight];
                const Int a01VCLocalHeight = Length(height-1,shift,p);
                F* a01Buffer = a01_MR_STAR.Buffer(offset,0);
                for( Int i=0; i<a01VCLocalHeight; ++i )
                    a01Buffer[i*r] = a01Data[i];
            }
            // Store w01Last[MR,* ]
            MemCopy
            ( W_MR_STAR.Buffer(0,W00.Width()+1),
              w01Last_MR_STAR.Buffer(),
              w01Last_MR_STAR.LocalHeight() );
            // Store a01[MR,* ]
            MemCopy
            ( APan_MR_STAR.Buffer(0,W00.Width()),
              a01_MR_STAR.Buffer(),
              a01_MR_STAR.LocalHeight() );

            // Update the portion of A00 that is in our current panel with 
            // w01Last and a01Last using two gers. We do not need their bottom
            // entries. We trash the lower triangle of our panel of A since we 
            // are only doing slightly more work and we can replace it
            // afterwards.
            auto a01Last_MC_STAR_Top = View( a01Last_MC_STAR, 0, 0, kA, 1 );
            auto w01Last_MC_STAR_Top = View( w01Last_MC_STAR, 0, 0, kA, 1 );
            auto a01Last_MR_STAR_TopPan = View( a01Last_MR_STAR, off, 0, k, 1 );
            auto w01Last_MR_STAR_TopPan = View( w01Last_MR_STAR, off, 0, k, 1 );
            const F* a01_MC_STAR_Buffer = a01Last_MC_STAR_Top.Buffer();
            const F* w01_MC_STAR_Buffer = w01Last_MC_STAR_Top.Buffer();
            const F* a01_MR_STAR_Buffer = a01Last_MR_STAR_TopPan.Buffer();
            const F* w01_MR_STAR_Buffer = w01Last_MR_STAR_TopPan.Buffer();
            F* A00PanBuffer = A00Pan.Buffer();
            const Int localHeight = A00Pan.LocalHeight();
            const Int localWidth = A00Pan.LocalWidth();
            const Int lDim = A00Pan.LDim();
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    A00PanBuffer[iLocal+jLocal*lDim] -=
                        w01_MC_STAR_Buffer[iLocal]*
                        Conj(a01_MR_STAR_Buffer[jLocal]) +
                        a01_MC_STAR_Buffer[iLocal]*
                        Conj(w01_MR_STAR_Buffer[jLocal]);
        }

        // Form the local portions of (A00 a01) into p01[MC,* ] and q01[MR,* ]:
        //   p01[MC,* ] := triu(A00)[MC,MR] a01[MR,* ]
        //   q01[MR,* ] := triu(A00,+1)'[MR,MC] a01[MC,* ]
        Zero( p01_MC_STAR );
        q01_MR_STAR.AlignWith( A00 );
        Zeros( q01_MR_STAR, a01.Height(), 1 );
        symv::LocalColAccumulate
        ( UPPER, F(1), 
          A00, a01_MC_STAR, a01_MR_STAR, p01_MC_STAR, q01_MR_STAR, true );

        x21_MR_STAR.AlignWith( A02T );
        y21_MR_STAR.AlignWith( A02T );
        Zeros( x21_MR_STAR, A02.Width(), 1 );
        Zeros( y21_MR_STAR, A02.Width(), 1 );
        LocalGemv( ADJOINT, F(1), W02T, a01T_MC_STAR, F(0), x21_MR_STAR );
        LocalGemv( ADJOINT, F(1), A02T, a01T_MC_STAR, F(0), y21_MR_STAR );

        // Combine the AllReduce column summations of x21[MR,* ], y21[MR,* ],
        // and q01[MR,* ]
        {
            const Int x21LocalHeight = x21_MR_STAR.LocalHeight();
            const Int y21LocalHeight = y21_MR_STAR.LocalHeight();
            const Int q01LocalHeight = q01_MR_STAR.LocalHeight();
            const Int reduceSize = x21LocalHeight+y21LocalHeight+q01LocalHeight;
            std::vector<F> colSumSendBuffer(reduceSize),
                           colSumRecvBuffer(reduceSize);
            MemCopy
            ( colSumSendBuffer.data(), x21_MR_STAR.Buffer(), x21LocalHeight );
            MemCopy
            ( &colSumSendBuffer[x21LocalHeight],
              y21_MR_STAR.Buffer(), y21LocalHeight );
            MemCopy
            ( &colSumSendBuffer[x21LocalHeight+y21LocalHeight],
              q01_MR_STAR.Buffer(), q01LocalHeight );
            mpi::AllReduce
            ( colSumSendBuffer.data(), colSumRecvBuffer.data(),
              reduceSize, g.ColComm() );
            MemCopy
            ( x21_MR_STAR.Buffer(), colSumRecvBuffer.data(), x21LocalHeight );
            MemCopy
            ( y21_MR_STAR.Buffer(), 
              &colSumRecvBuffer[x21LocalHeight], y21LocalHeight );
            MemCopy
            ( q01_MR_STAR.Buffer(), 
              &colSumRecvBuffer[x21LocalHeight+y21LocalHeight], 
              q01LocalHeight );
        }

        LocalGemv( NORMAL, F(-1), A02T, x21_MR_STAR, F(1), p01T_MC_STAR );
        LocalGemv( NORMAL, F(-1), W02T, y21_MR_STAR, F(1), p01T_MC_STAR );

        if( W00.Width() > 0 )
        {
            // This is not the last iteration of the panel factorization, 
            // combine the Reduce to one of p01[MC,* ] with the redistribution 
            // of q01[MR,* ] -> q01[MC,MR] to the next process column.
            const Int localHeight = p01_MC_STAR.LocalHeight();
            std::vector<F> reduceToOneSendBuffer(2*localHeight),
                           reduceToOneRecvBuffer(2*localHeight);

            // Pack p01[MC,* ]
            MemCopy
            ( reduceToOneSendBuffer.data(), p01_MC_STAR.Buffer(), localHeight );

            // Fill in contributions to q01[MC,MR] from q01[MR,* ]
            const bool contributing = 
                ( q01_MR_STAR.ColShift() % g.GCD() ==
                  p01_MC_STAR.ColShift() % g.GCD() );
            if( contributing )
            {
                if( r == c )
                {
                    MemCopy
                    ( &reduceToOneSendBuffer[localHeight],
                      q01_MR_STAR.Buffer(), localHeight );
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
                    const Int a0 = p01_MC_STAR.ColShift();
                    const Int b0 = q01_MR_STAR.ColShift();

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
                    const F* q01_MR_STAR_Buffer = q01_MR_STAR.Buffer();
                    const Int offset = localHeight + targetStart;
                    for( Int i=0; i<localLength; ++i )                        
                        reduceToOneSendBuffer[offset+i*targetPeriod] = 
                            q01_MR_STAR_Buffer[sourceStart+i*sourcePeriod];
                }
            }
            else
                MemZero( &reduceToOneSendBuffer[localHeight], localHeight );

            const Int nextProcessRow = (alpha11.ColAlign()+r-1) % r;
            const Int nextProcessCol = (alpha11.RowAlign()+c-1) % c;
            mpi::Reduce
            ( reduceToOneSendBuffer.data(), reduceToOneRecvBuffer.data(),
              2*localHeight, nextProcessCol, g.RowComm() );
            if( g.Col() == nextProcessCol )
            {
                // Combine the second half into the first half        
                for( Int i=0; i<localHeight; ++i )
                    reduceToOneRecvBuffer[i] +=
                        reduceToOneRecvBuffer[i+localHeight];

                // Finish computing w01. During its computation, ensure that 
                // every process has a copy of the last element of the w01.
                // We know a priori that the last element of a01 is one.
                const F* a01_MC_STAR_Buffer = a01_MC_STAR.Buffer();
                F myDotProduct = blas::Dot
                    ( localHeight, reduceToOneRecvBuffer.data(), 1, 
                                   a01_MC_STAR_Buffer,           1 );
                F sendBuffer[2], recvBuffer[2];
                sendBuffer[0] = myDotProduct;
                sendBuffer[1] = ( g.Row()==nextProcessRow ? 
                                  reduceToOneRecvBuffer[localHeight-1] : 0 );
                mpi::AllReduce( sendBuffer, recvBuffer, 2, g.ColComm() );
                F dotProduct = recvBuffer[0];

                // Set up for the next iteration by filling in the values for:
                // - w01LastBuffer
                // - w01LastBottomEntry
                F scale = dotProduct*tau / F(2);
                for( Int i=0; i<localHeight; ++i )
                    w01LastBuffer[i] = Conj(tau)*
                        ( reduceToOneRecvBuffer[i]-
                          scale*a01_MC_STAR_Buffer[i] );
                w01LastBottomEntry = Conj(tau)*( recvBuffer[1]-scale );
            }
        }
        else
        {
            // This is the last iteration, our last task is to finish forming
            // w01[MC,* ] and w01[MR,* ] so that we may place them into W[MC,* ]
            // and W[MR,* ]
            const Int localHeight = p01_MC_STAR.LocalHeight();
            std::vector<F> allReduceSendBuffer(2*localHeight),
                           allReduceRecvBuffer(2*localHeight);

            // Pack p01[MC,* ]
            MemCopy
            ( allReduceSendBuffer.data(), p01_MC_STAR.Buffer(), localHeight );

            // Fill in contributions to q01[MC,* ] from q01[MR,* ]
            const bool contributing = 
                ( q01_MR_STAR.ColShift() % g.GCD() ==
                  p01_MC_STAR.ColShift() % g.GCD() );
            if( contributing )
            {
                if( r == c )
                {
                    MemCopy
                    ( &allReduceSendBuffer[localHeight],
                      q01_MR_STAR.Buffer(), localHeight );
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
                    const Int a0 = p01_MC_STAR.ColShift();
                    const Int b0 = q01_MR_STAR.ColShift();

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
                    const F* q01_MR_STAR_Buffer = q01_MR_STAR.Buffer();
                    const Int offset = localHeight + targetStart;
                    for( Int i=0; i<localLength; ++i )
                        allReduceSendBuffer[offset+i*targetPeriod] = 
                            q01_MR_STAR_Buffer[sourceStart+i*sourcePeriod];
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
 
            // Finish computing w01. During its computation, ensure that 
            // every process has a copy of the last element of the w01.
            // We know a priori that the last element of a01 is one.
            const F* a01_MC_STAR_Buffer = a01_MC_STAR.Buffer();
            F myDotProduct = blas::Dot
                ( localHeight, allReduceRecvBuffer.data(), 1, 
                               a01_MC_STAR_Buffer,         1 );
            const F dotProduct = mpi::AllReduce( myDotProduct, g.ColComm() );

            // Grab views into W[MC,* ] and W[MR,* ]
            auto w01_MC_STAR = View( W_MC_STAR, 0, k, kA, 1 );
            auto w01_MR_STAR = View( W_MR_STAR, 0, k, kA, 1 );

            // Store w01[MC,* ]
            F scale = dotProduct*tau / F(2);
            F* w01_MC_STAR_Buffer = w01_MC_STAR.Buffer();
            for( Int i=0; i<localHeight; ++i )
                w01_MC_STAR_Buffer[i] = Conj(tau)*
                    ( allReduceRecvBuffer[i]-scale*a01_MC_STAR_Buffer[i] );

            // Form w01[MR,* ]
            w01_MR_STAR = w01_MC_STAR;
        }
    }

    expandedABR.SetRealPartOfDiagonal( e, 1 );
}

} // namespace herm_tridiag
} // namespace El

#endif // ifndef EL_HERMITIANTRIDIAG_UPAN_HPP
