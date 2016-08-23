/*
   Copyright (c) 2009-2016, Jack Poulson
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
  DistMatrix<F,MC,STAR>& B_MC_STAR, 
  DistMatrix<F,MR,STAR>& B_MR_STAR,
  DistMatrix<F,MC,STAR>& W_MC_STAR,
  DistMatrix<F,MR,STAR>& W_MR_STAR,
  const SymvCtrl<F>& ctrl )
{
    DEBUG_CSE
    const Int n = A.Height();
    const Int nW = W.Width();
    DEBUG_ONLY(
      AssertSameGrids( A, W, t );
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
    auto expandedABR = A( IR(off-1,n), IR(off-1,n) ); 
    DistMatrix<Real,MD,STAR> e(g);
    e.SetRoot( expandedABR.DiagonalRoot(1) );
    e.AlignCols( expandedABR.DiagonalAlign(1) );
    e.Resize( nW, 1 );

    vector<F> w01LastBuf;
    FastResize( w01LastBuf, n/r+1 );

    DistMatrix<F> w01Last(g);
    DistMatrix<F,MC,STAR> a01_MC(g), p01_MC(g),
                          a01Last_MC(g), w01Last_MC(g);
    DistMatrix<F,MR,STAR> a01_MR(g), q01_MR(g),
                          x21_MR(g), y21_MR(g),
                          a01Last_MR(g), w01Last_MR(g);

    F tau = 0;
    F w01LastBottomEntry = 0;
    for( Int k=nW-1; k>=0; --k )
    {
        const Int kA = k+off;
        const bool firstIteration = ( k == nW-1 );
        if( !firstIteration ) 
        {
            // TODO: Move these and make them auto
            View( a01Last_MC, B_MC_STAR, IR(0,kA+1), IR(k+1) );
            View( a01Last_MR, B_MR_STAR, IR(0,kA+1), IR(k+1) );
            View( w01Last,         W,         IR(0,kA+1), IR(k+1) );
        }

        const Range<Int> ind0( 0, kA   ),
                         indT( 0, kA+1 ),
                         ind1( kA ),
                         ind2( kA+1, n );

        auto A00     = A( ind0, ind0 );
        auto a01     = A( ind0, ind1 );
        auto A02     = A( ind0, ind2 );
        auto aT1     = A( indT, ind1 );
        auto alpha11 = A( ind1, ind1 );

        auto a01T     = A( IR(0,kA-1), ind1       );
        auto alpha01B = A( IR(kA-1),   ind1       );
        auto A02T     = A( IR(0,off),  ind2       );
        auto A00Pan   = A( ind0,       IR(off,kA) );

        auto W02T     = W( IR(0,off), IR(k+1,nW) );

        auto tau1     = t( ind1-off, ALL );
        auto epsilon1 = e( ind1-off, ALL );

        a01_MC.AlignWith( A00 );
        a01_MR.AlignWith( A00 );
        p01_MC.AlignWith( A00 );
        a01_MC.Resize( kA, 1 );
        a01_MR.Resize( kA, 1 );
        p01_MC.Resize( kA, 1 );

        // View the portions of a01[MC] and p01[MC] above the current
        // panel's square
        auto a01T_MC = a01_MC( IR(0,off), ALL );
        auto p01T_MC = p01_MC( IR(0,off), ALL ); 

        const bool thisIsMyCol = ( g.Col() == alpha11.RowAlign() );
        if( thisIsMyCol )
        {
            if( !firstIteration )
            {
                // Finish updating the current column with two axpy's
                const Int aT1LocalHeight = aT1.LocalHeight();
                F* aT1Buf = aT1.Buffer();
                const F* a01Last_MC_Buf = a01Last_MC.Buffer();
                for( Int i=0; i<aT1LocalHeight; ++i )
                    aT1Buf[i] -=
                      w01LastBuf[i] + 
                      a01Last_MC_Buf[i]*Conj(w01LastBottomEntry);
            }
            // Compute the Householder reflector
            tau = reflector::Col( alpha01B, a01T );
            tau1.Set(0,0,tau);
        }

        // Store the subdiagonal value and turn a01 into a proper scaled 
        // reflector by explicitly placing the implicit one in its first entry.
        GetRealPartOfDiagonal( alpha01B, epsilon1 );
        alpha01B.Set( 0, 0, F(1) );

        // If this is the first iteration, have each member of the owning 
        // process column broadcast tau and a01 within its process row. 
        // Otherwise, also add w01 into the broadcast.
        if( firstIteration )
        {
            const Int a01LocalHeight = a01.LocalHeight();
            vector<F> rowBcastBuf;
            FastResize( rowBcastBuf, a01LocalHeight+1 );

            if( thisIsMyCol )
            {
                // Pack the broadcast buffer with a01 and tau
                MemCopy
                ( rowBcastBuf.data(), a01.Buffer(), a01LocalHeight );
                rowBcastBuf[a01LocalHeight] = tau;
            }
            // Broadcast a01 and tau across the process row
            mpi::Broadcast
            ( rowBcastBuf.data(), 
              a01LocalHeight+1, a01.RowAlign(), g.RowComm() );
            // Store a01[MC] into its DistMatrix class and also store a copy
            // for the next iteration
            MemCopy
            ( a01_MC.Buffer(), rowBcastBuf.data(), a01LocalHeight );
            // Store a01[MC] into B[MC,* ]
            MemCopy
            ( B_MC_STAR.Buffer(0,k), 
              rowBcastBuf.data(), a01LocalHeight );
            // Store tau
            tau = rowBcastBuf[a01LocalHeight];
            
            a01_MR = a01_MC;
            // Store a01[MR]
            MemCopy
            ( B_MR_STAR.Buffer(0,k),
              a01_MR.Buffer(),
              a01_MR.LocalHeight() );
        }
        else
        {
            const Int a01LocalHeight = a01.LocalHeight();
            const Int w01LastLocalHeight = aT1.LocalHeight();
            vector<F> rowBcastBuf;
            FastResize( rowBcastBuf, a01LocalHeight+w01LastLocalHeight+1 );

            if( thisIsMyCol ) 
            {
                // Pack the broadcast buffer with a01, w01Last, and tau
                MemCopy
                ( rowBcastBuf.data(), a01.Buffer(), a01LocalHeight );
                MemCopy
                ( &rowBcastBuf[a01LocalHeight], 
                  w01LastBuf.data(), w01LastLocalHeight );
                rowBcastBuf[a01LocalHeight+w01LastLocalHeight] = tau;
            }
            // Broadcast a01, w01Last, and tau across the process row
            mpi::Broadcast
            ( rowBcastBuf.data(), 
              a01LocalHeight+w01LastLocalHeight+1, 
              a01.RowAlign(), g.RowComm() );
            // Store a01[MC] into its DistMatrix class 
            MemCopy
            ( a01_MC.Buffer(), rowBcastBuf.data(), a01LocalHeight );
            // Store a01[MC] into B[MC,* ]
            MemCopy
            ( B_MC_STAR.Buffer(0,k), 
              rowBcastBuf.data(), a01LocalHeight );
            // Store w01Last[MC] into its DistMatrix class
            w01Last_MC.AlignWith( A00 );
            w01Last_MC.Resize( kA+1, 1 );
            MemCopy
            ( w01Last_MC.Buffer(), 
              &rowBcastBuf[a01LocalHeight], w01LastLocalHeight );
            // Store the bottom part of w01Last[MC] into WB[MC,* ] and, 
            // if necessary, w01.
            MemCopy
            ( W_MC_STAR.Buffer(0,k+1),
              &rowBcastBuf[a01LocalHeight], w01LastLocalHeight );
            if( g.Col() == w01Last.RowAlign() )
            {
                MemCopy
                ( w01Last.Buffer(),
                  &rowBcastBuf[a01LocalHeight], w01LastLocalHeight );
            }
            // Store tau
            tau = rowBcastBuf[a01LocalHeight+w01LastLocalHeight];

            // Form a01[MR] and w01Last[MR] by combining the 
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

            const Int height = kA+1;
            const Int portionSize = mpi::Pad( 2*MaxLength(height,p) );

            const Int colShiftVRDest = Shift(g.VRRank(),colAlignDest,p);
            const Int colShiftVCSource = Shift(g.VCRank(),colAlignSource,p);
            const Int sendRankRM = 
                (g.VRRank()+(p+colShiftVCSource-colShiftVRDest))%p;
            const Int recvRankCM = 
                (g.VCRank()+(p+colShiftVRDest-colShiftVCSource))%p;
            const Int recvRankRM = 
                (recvRankCM/r)+c*(recvRankCM%r);

            vector<F> transBuf;
            FastResize( transBuf, (r+1)*portionSize );
            F* sendBuf = &transBuf[0];
            F* recvBuf = &transBuf[r*portionSize];

            // (w01Last[VC] <- w01Last[MC]) and
            // ([a01; 0][VC] <- [a01; 0][MC])
            {
                // Pack the necessary portion of w01Last[MC]
                const Int shift = Shift(g.VCRank(),colAlignSource,p);
                const Int offset = (shift-colShiftSource)/r;
                const Int w01VCLocalHeight = Length(height,shift,p);
                StridedMemCopy
                ( sendBuf,                           1,
                  w01Last_MC.LockedBuffer(offset,0), c, w01VCLocalHeight );
                
                // Pack the necessary portion of a01[MC]
                const Int a01VCLocalHeight = Length(height-1,shift,p);
                StridedMemCopy
                ( &sendBuf[w01VCLocalHeight],    1,
                  a01_MC.LockedBuffer(offset,0), c, a01VCLocalHeight );
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
            w01Last_MR.AlignWith( A00 );
            w01Last_MR.Resize( kA+1, 1 );
            for( Int row=0; row<r; ++row )
            {
                // Unpack into w01Last[MR]
                const Int shift = Shift(g.Col()+c*row,colAlignDest,p);
                const Int offset = (shift-colShiftDest) / c;
                const Int w01VCLocalHeight = Length(height,shift,p);
                StridedMemCopy
                ( w01Last_MR.Buffer(offset,0), r,
                  &sendBuf[row*portionSize],   1, w01VCLocalHeight );

                // Unpack into a01[MR]
                const F* a01Data = &sendBuf[row*portionSize+w01VCLocalHeight];
                const Int a01VCLocalHeight = Length(height-1,shift,p);
                StridedMemCopy
                ( a01_MR.Buffer(offset,0), r,
                  a01Data,                 1, a01VCLocalHeight );
            }
            // Store w01Last[MR]
            MemCopy
            ( W_MR_STAR.Buffer(0,k+1),
              w01Last_MR.Buffer(),
              w01Last_MR.LocalHeight() );
            // Store a01[MR]
            MemCopy
            ( B_MR_STAR.Buffer(0,k),
              a01_MR.Buffer(),
              a01_MR.LocalHeight() );

            // Update the portion of A00 that is in our current panel with 
            // w01Last and a01Last using two gers. We do not need their bottom
            // entries. We trash the lower triangle of our panel of A since we 
            // are only doing slightly more work and we can replace it
            // afterwards.
            auto a01Last_MC_Top = a01Last_MC( ind0, ALL );
            auto w01Last_MC_Top = w01Last_MC( ind0, ALL );
            auto a01Last_MR_TopPan = a01Last_MR( IR(off,off+k), ALL );
            auto w01Last_MR_TopPan = w01Last_MR( IR(off,off+k), ALL );
            Ger2Sub
            ( a01Last_MC_Top,
              w01Last_MC_Top,
              a01Last_MR_TopPan,
              w01Last_MR_TopPan,
              A00Pan );
        }

        // Form the local portions of (A00 a01) into p01[MC] and q01[MR]:
        //   p01[MC] := triu(A00)[MC,MR] a01[MR]
        //   q01[MR] := triu(A00,+1)'[MR,MC] a01[MC]
        Zero( p01_MC );
        q01_MR.AlignWith( A00 );
        Zeros( q01_MR, kA, 1 );
        symv::LocalColAccumulate
        ( UPPER, F(1), A00, a01_MC, a01_MR, p01_MC, q01_MR, true, ctrl );

        x21_MR.AlignWith( A02T );
        y21_MR.AlignWith( A02T );
        Zeros( x21_MR, A02.Width(), 1 );
        Zeros( y21_MR, A02.Width(), 1 );
        LocalGemv( ADJOINT, F(1), W02T, a01T_MC, F(0), x21_MR );
        LocalGemv( ADJOINT, F(1), A02T, a01T_MC, F(0), y21_MR );

        // Combine the AllReduce column summations of x21[MR], y21[MR],
        // and q01[MR]
        {
            const Int x21LocalHeight = x21_MR.LocalHeight();
            const Int y21LocalHeight = y21_MR.LocalHeight();
            const Int q01LocalHeight = q01_MR.LocalHeight();
            const Int reduceSize = x21LocalHeight+y21LocalHeight+q01LocalHeight;
            vector<F> colSumSendBuf, colSumRecvBuf;
            FastResize( colSumSendBuf, reduceSize );
            FastResize( colSumRecvBuf, reduceSize );

            MemCopy
            ( colSumSendBuf.data(), x21_MR.Buffer(), x21LocalHeight );
            MemCopy
            ( &colSumSendBuf[x21LocalHeight],
              y21_MR.Buffer(), y21LocalHeight );
            MemCopy
            ( &colSumSendBuf[x21LocalHeight+y21LocalHeight],
              q01_MR.Buffer(), q01LocalHeight );
            mpi::AllReduce
            ( colSumSendBuf.data(), colSumRecvBuf.data(),
              reduceSize, g.ColComm() );
            MemCopy
            ( x21_MR.Buffer(), colSumRecvBuf.data(), x21LocalHeight );
            MemCopy
            ( y21_MR.Buffer(), 
              &colSumRecvBuf[x21LocalHeight], y21LocalHeight );
            MemCopy
            ( q01_MR.Buffer(), 
              &colSumRecvBuf[x21LocalHeight+y21LocalHeight], q01LocalHeight );
        }

        LocalGemv( NORMAL, F(-1), A02T, x21_MR, F(1), p01T_MC );
        LocalGemv( NORMAL, F(-1), W02T, y21_MR, F(1), p01T_MC );

        if( k > 0 )
        {
            // This is not the last iteration of the panel factorization, 
            // combine the Reduce to one of p01[MC] with the redistribution 
            // of q01[MR,* ] -> q01[MC,MR] to the next process column.
            const Int localHeight = p01_MC.LocalHeight();
            vector<F> reduceToOneSendBuf, reduceToOneRecvBuf;
            FastResize( reduceToOneSendBuf, 2*localHeight );
            FastResize( reduceToOneRecvBuf, 2*localHeight );

            // Pack p01[MC]
            MemCopy
            ( reduceToOneSendBuf.data(), p01_MC.Buffer(), localHeight );

            // Fill in contributions to q01[MC,MR] from q01[MR,* ]
            const bool contributing = 
                ( q01_MR.ColShift() % g.GCD() ==
                  p01_MC.ColShift() % g.GCD() );
            if( contributing )
            {
                if( r == c )
                {
                    MemCopy
                    ( &reduceToOneSendBuf[localHeight],
                      q01_MR.Buffer(), localHeight );
                }
                else
                {
                    // Zero the entire buffer first
                    MemZero( &reduceToOneSendBuf[localHeight], localHeight );
                    // Fill in the entries that we contribute to.
                    // We seek to find the minimum s in N such that
                    //   s*c = a0-b0 (mod r)
                    // where a0 is the column shift of MC, b0 is the row shift
                    // of MR, and s is our first local entry of MR that will 
                    // contribute to MC. I cannot think of an O(1) method, so
                    // I will instead use the worst-case O(lcm(c,r)/c) method.
                    const Int sourcePeriod = g.LCM() / c;
                    const Int targetPeriod = g.LCM() / r;
                    const Int a0 = p01_MC.ColShift();
                    const Int b0 = q01_MR.ColShift();

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
                    const Int offset = localHeight + targetStart;
                    StridedMemCopy
                    ( &reduceToOneSendBuf[offset],       targetPeriod,
                      q01_MR.Buffer(sourceStart,0), sourcePeriod,
                      localLength );
                }
            }
            else
                MemZero( &reduceToOneSendBuf[localHeight], localHeight );

            const Int nextProcessRow = (alpha11.ColAlign()+r-1) % r;
            const Int nextProcessCol = (alpha11.RowAlign()+c-1) % c;
            mpi::Reduce
            ( reduceToOneSendBuf.data(), reduceToOneRecvBuf.data(),
              2*localHeight, nextProcessCol, g.RowComm() );
            if( g.Col() == nextProcessCol )
            {
                // Combine the second half into the first half        
                blas::Axpy
                ( localHeight, F(1),
                  &reduceToOneRecvBuf[localHeight], 1,
                  &reduceToOneRecvBuf[0],           1 );

                // Finish computing w01. During its computation, ensure that 
                // every process has a copy of the last element of the w01.
                // We know a priori that the last element of a01 is one.
                const F* a01_MC_Buf = a01_MC.Buffer();
                F myDotProduct = blas::Dot
                    ( localHeight, reduceToOneRecvBuf.data(), 1, 
                                   a01_MC_Buf,                1 );
                F sendBuf[2], recvBuf[2];
                sendBuf[0] = myDotProduct;
                sendBuf[1] = ( g.Row()==nextProcessRow ? 
                               reduceToOneRecvBuf[localHeight-1] : 0 );
                mpi::AllReduce( sendBuf, recvBuf, 2, g.ColComm() );
                F dotProduct = recvBuf[0];

                // Set up for the next iteration by filling in the values for:
                // - w01LastBuf
                // - w01LastBottomEntry
                F scale = dotProduct*tau / F(2);
                for( Int i=0; i<localHeight; ++i )
                    w01LastBuf[i] = Conj(tau)*
                      ( reduceToOneRecvBuf[i]-scale*a01_MC_Buf[i] );
                w01LastBottomEntry = Conj(tau)*( recvBuf[1]-scale );
            }
        }
        else
        {
            // This is the last iteration, our last task is to finish forming
            // w01[MC] and w01[MR] so that we may place them into W[MC,* ]
            // and W[MR,* ]
            const Int localHeight = p01_MC.LocalHeight();
            vector<F> allReduceSendBuf, allReduceRecvBuf;
            FastResize( allReduceSendBuf, 2*localHeight );
            FastResize( allReduceRecvBuf, 2*localHeight );

            // Pack p01[MC]
            MemCopy
            ( allReduceSendBuf.data(), p01_MC.Buffer(), localHeight );

            // Fill in contributions to q01[MC] from q01[MR]
            const bool contributing = 
                ( q01_MR.ColShift() % g.GCD() ==
                  p01_MC.ColShift() % g.GCD() );
            if( contributing )
            {
                if( r == c )
                {
                    MemCopy
                    ( &allReduceSendBuf[localHeight],
                      q01_MR.Buffer(), localHeight );
                }
                else
                {
                    // Zero the entire buffer first
                    MemZero( &allReduceSendBuf[localHeight], localHeight );
                    // Fill in the entries that we contribute to.
                    // We seek to find the minimum s in N such that
                    //   s*c = a0-b0 (mod r)
                    // where a0 is the column shift of MC, b0 is the row shift
                    // of MR, and s is our first local entry of MR that will 
                    // contribute to MC. I cannot think of an O(1) method, so
                    // I will instead use the worst-case O(lcm(c,r)/c) method.
                    const Int sourcePeriod = g.LCM() / c;
                    const Int targetPeriod = g.LCM() / r;
                    const Int a0 = p01_MC.ColShift();
                    const Int b0 = q01_MR.ColShift();

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
                    const Int offset = localHeight + targetStart;
                    StridedMemCopy
                    ( &allReduceSendBuf[offset],    targetPeriod,
                      q01_MR.Buffer(sourceStart,0), sourcePeriod,
                      localLength );
                }
            }
            else
                MemZero( &allReduceSendBuf[localHeight], localHeight );

            mpi::AllReduce
            ( allReduceSendBuf.data(), allReduceRecvBuf.data(),
              2*localHeight, g.RowComm() );

            // Combine the second half into the first half        
            blas::Axpy
            ( localHeight, F(1), 
              &allReduceRecvBuf[localHeight], 1, 
              &allReduceRecvBuf[0],           1 );
 
            // Finish computing w01. During its computation, ensure that 
            // every process has a copy of the last element of the w01.
            // We know a priori that the last element of a01 is one.
            const F* a01_MC_Buf = a01_MC.Buffer();
            F myDotProduct = blas::Dot
                ( localHeight, allReduceRecvBuf.data(), 1, 
                               a01_MC_Buf,              1 );
            const F dotProduct = mpi::AllReduce( myDotProduct, g.ColComm() );

            // Grab views into W[MC,* ] and W[MR,* ]
            auto w01_MC = W_MC_STAR( ind0, ALL );
            auto w01_MR = W_MR_STAR( ind0, ALL );

            // Store w01[MC]
            F scale = dotProduct*tau / F(2);
            F* w01_MC_Buf = w01_MC.Buffer();
            for( Int i=0; i<localHeight; ++i )
                w01_MC_Buf[i] = Conj(tau)*
                  ( allReduceRecvBuf[i]-scale*a01_MC_Buf[i] );

            // Form w01[MR]
            w01_MR = w01_MC;
        }
    }

    SetRealPartOfDiagonal( expandedABR, e, 1 );
}

} // namespace herm_tridiag
} // namespace El

#endif // ifndef EL_HERMITIANTRIDIAG_UPAN_HPP
