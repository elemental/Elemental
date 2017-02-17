/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HERMITIANTRIDIAG_UPPER_PANEL_SQUARE_HPP
#define EL_HERMITIANTRIDIAG_UPPER_PANEL_SQUARE_HPP

namespace El {
namespace herm_tridiag {

template<typename F>
void UpperPanelSquare
( DistMatrix<F>& A,
  DistMatrix<F>& W,
  DistMatrix<F,MD,STAR>& t,
  DistMatrix<F,MC,STAR>& B_MC_STAR, 
  DistMatrix<F,MR,STAR>& B_MR_STAR,
  DistMatrix<F,MC,STAR>& W_MC_STAR,
  DistMatrix<F,MR,STAR>& W_MR_STAR,
  const SymvCtrl<F>& ctrl )
{
    EL_DEBUG_CSE
    const Int n = A.Height();
    const Int nW = W.Width();
    EL_DEBUG_ONLY(
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
    const Int off = n-nW;

    // Find the process holding our transposed data
    Int transposeRank;
    {
        const Int colAlign = A.ColAlign();
        const Int rowAlign = A.RowAlign();
        const Int colShift = A.ColShift();
        const Int rowShift = A.RowShift();

        const Int transposeRow = (colAlign+rowShift) % r;
        const Int transposeCol = (rowAlign+colShift) % r;
        transposeRank = transposeRow + r*transposeCol;
    }
    const bool onDiagonal = ( transposeRank == g.VCRank() );

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
                          a01Last_MR(g), w01Last_MR(g),
                          x21_MR(g), y21_MR(g);

    F tau = 0;
    F w01LastBottomEntry = 0;
    for( Int k=nW-1; k>=0; --k )
    {
        const Int kA = k+off;
        const bool firstIteration = ( k == nW-1 );

        const Range<Int> ind0( 0,    kA   ),
                         indT( 0,    kA+1 ), indL( 0, kA+1 ),
                         ind1( kA,   kA+1 ),
                         ind2( kA+1, n    );

        auto A00     = A( ind0, ind0 );
        auto a01     = A( ind0, ind1 );
        auto ATL     = A( indT, indL );
        auto aT1     = A( indT, ind1 );
        auto A02     = A( ind0, ind2 );
        auto alpha11 = A( ind1, ind1 );

        auto A00Pan   = A( ind0, IR(off,kA) );

        auto a01T     = A( IR(0,kA-1), ind1 );
        auto alpha01B = A( IR(kA-1),   ind1 );
        auto A02T     = A( IR(0,off),  ind2 );

        auto W02T = W( IR(0,off), IR(k+1,END) );

        auto tau1     = t( ind1-off, ALL );
        auto epsilon1 = e( ind1-off, ALL );

        a01_MC.AlignWith( A00 );
        a01_MR.AlignWith( A00 );
        p01_MC.AlignWith( A00 );
        q01_MR.AlignWith( A00 );
        a01_MC.Resize( kA, 1 );
        a01_MR.Resize( kA, 1 );
        Zeros( p01_MC, kA, 1 );
        Zeros( q01_MR, kA, 1 );

        // View the portions of A02 and W0T outside of this panel's square
        auto a01T_MC = a01_MC( IR(0,off), ALL );
        auto p01T_MC = p01_MC( IR(0,off), ALL ); 

        if( !firstIteration )
        {
            // TODO: Move these and make them auto
            View( a01Last_MC, B_MC_STAR, indT, IR(k+1) );
            View( a01Last_MR, B_MR_STAR, indT, IR(k+1) );
            View( w01Last,    W,         indT, IR(k+1) );
        }

        const bool thisIsMyCol = ( g.Col() == alpha11.RowAlign() );
        if( thisIsMyCol )
        {
            if( !firstIteration )
            {
                // Finish updating the current column with two axpy's
                const Int aT1LocalHeight = aT1.LocalHeight();
                F* aT1Buffer = aT1.Buffer();
                const F* a01Last_MC_Buf = a01Last_MC.Buffer();
                for( Int i=0; i<aT1LocalHeight; ++i )
                    aT1Buffer[i] -=
                        w01LastBuf[i] + 
                        a01Last_MC_Buf[i]*Conj(w01LastBottomEntry);
            }
            // Compute the Householder reflector
            tau = reflector::Col( alpha01B, a01T );
            if( g.Row() == alpha01B.ColAlign() )
                tau1.SetLocal(0,0,tau);
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
            
            // Take advantage of the square process grid in order to form
            // a01[MR,* ] from a01[MC,* ]
            if( onDiagonal )
            {
                MemCopy
                ( a01_MR.Buffer(),
                  a01_MC.Buffer(), a01LocalHeight );
            }
            else
            {
                // Pairwise exchange
                const Int sendSize = A00.LocalHeight();
                const Int recvSize = A00.LocalWidth();
                mpi::SendRecv
                ( a01_MC.Buffer(), sendSize, transposeRank,
                  a01_MR.Buffer(), recvSize, transposeRank, g.VCComm() );
            }
            // Store a01[MR]
            MemCopy
            ( B_MR_STAR.Buffer(0,k),
              a01_MR.Buffer(), a01_MR.LocalHeight() );
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
                MemCopy( rowBcastBuf.data(), a01.Buffer(), a01LocalHeight );
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
            w01Last_MC.Resize( a01.Height()+1, 1 );
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

            // Take advantage of the square process grid in order to quickly
            // form a01[MR,* ] and w01Last[MR,* ] from their [MC,* ]
            // counterparts
            w01Last_MR.AlignWith( A00 );
            w01Last_MR.Resize( w01Last.Height(), 1 );
            if( onDiagonal )
            {
                MemCopy
                ( a01_MR.Buffer(),
                  a01_MC.Buffer(), a01LocalHeight );
                MemCopy
                ( w01Last_MR.Buffer(),
                  w01Last_MC.Buffer(), w01LastLocalHeight );
            }
            else
            {
                const Int sendSize = A00.LocalHeight()+ATL.LocalHeight();
                const Int recvSize = A00.LocalWidth()+ATL.LocalWidth();
                vector<F> sendBuf, recvBuf;
                FastResize( sendBuf, sendSize );
                FastResize( recvBuf, recvSize );

                // Pack the send buffer
                MemCopy
                ( sendBuf.data(), a01_MC.Buffer(), A00.LocalHeight() );
                MemCopy
                ( &sendBuf[A00.LocalHeight()],
                  w01Last_MC.Buffer(), ATL.LocalHeight() );

                // Pairwise exchange
                mpi::SendRecv
                ( sendBuf.data(), sendSize, transposeRank,
                  recvBuf.data(), recvSize, transposeRank, g.VCComm() );

                // Unpack the recv buffer
                MemCopy
                ( a01_MR.Buffer(), recvBuf.data(), A00.LocalWidth() );
                MemCopy
                ( w01Last_MR.Buffer(),
                  &recvBuf[A00.LocalWidth()], ATL.LocalWidth() );
            }

            // Store w01Last[MR]
            MemCopy
            ( W_MR_STAR.Buffer(0,k+1),
              w01Last_MR.Buffer(),
              w01Last_MR.LocalHeight() );
            // Store a01[MR]
            MemCopy
            ( B_MR_STAR.Buffer(0,k),
              a01_MR.Buffer(), a01_MR.LocalHeight() );

            // Update the portion of A00 that is in our current panel with 
            // w01Last and a01Last using two gers. We do not need their bottom
            // entries. We trash the lower triangle of our panel of A since we 
            // are only doing slightly more work and we can replace it
            // afterwards.
            auto a01Last_MC_Top = a01Last_MC( ind0, ALL );
            auto w01Last_MC_Top = w01Last_MC( ind0, ALL );
            auto a01Last_MR_TopPan = a01Last_MR( IR(off,kA), ALL );
            auto w01Last_MR_TopPan = w01Last_MR( IR(off,kA), ALL );
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
        symv::LocalColAccumulate
        ( UPPER, F(1), A00, a01_MC, a01_MR, p01_MC, q01_MR, true, ctrl );

        x21_MR.AlignWith( A02T );
        y21_MR.AlignWith( A02T );
        Zeros( x21_MR, A02.Width(), 1 );
        Zeros( y21_MR, A02.Width(), 1 );
        LocalGemv( ADJOINT, F(1), W02T, a01T_MC, F(0), x21_MR );
        LocalGemv( ADJOINT, F(1), A02T, a01T_MC, F(0), y21_MR );

        // Combine the AllReduce column summations of x21[MR] and y21[MR]
        {
            const Int x21LocalHeight = x21_MR.LocalHeight();
            const Int y21LocalHeight = y21_MR.LocalHeight();
            const Int reduceSize = x21LocalHeight+y21LocalHeight;
            vector<F> colSumSendBuf, colSumRecvBuf;
            FastResize( colSumSendBuf, reduceSize );
            FastResize( colSumRecvBuf, reduceSize );

            MemCopy
            ( colSumSendBuf.data(), x21_MR.Buffer(), x21LocalHeight );
            MemCopy
            ( &colSumSendBuf[x21LocalHeight],
              y21_MR.Buffer(), y21LocalHeight );
            mpi::AllReduce
            ( colSumSendBuf.data(), colSumRecvBuf.data(),
              reduceSize, g.ColComm() );
            MemCopy
            ( x21_MR.Buffer(), colSumRecvBuf.data(), x21LocalHeight );
            MemCopy
            ( y21_MR.Buffer(), 
              &colSumRecvBuf[x21LocalHeight], y21LocalHeight );
        }

        LocalGemv( NORMAL, F(-1), A02T, x21_MR, F(1), p01T_MC );
        LocalGemv( NORMAL, F(-1), W02T, y21_MR, F(1), p01T_MC );

        // Fast transpose the unsummed q01[MR] -> q01[MC], so that
        // it needs to be summed over process rows instead of process 
        // columns. We immediately add it onto p01[MC], which also needs
        // to be summed within process rows.
        if( onDiagonal )
        {
            const Int a01LocalHeight = a01.LocalHeight();
            F* p01_MC_Buf = p01_MC.Buffer();
            const F* q01_MR_Buf = q01_MR.Buffer();
            for( Int i=0; i<a01LocalHeight; ++i )
                p01_MC_Buf[i] += q01_MR_Buf[i];
        }
        else
        {
            // Pairwise exchange with the transpose process
            const Int sendSize = A00.LocalWidth();
            const Int recvSize = A00.LocalHeight();
            vector<F> recvBuf;
            FastResize( recvBuf, recvSize );

            mpi::SendRecv
            ( q01_MR.Buffer(), sendSize, transposeRank,
              recvBuf.data(),  recvSize, transposeRank, g.VCComm() );

            // Unpack the recv buffer directly onto p01[MC]
            F* p01_MC_Buf = p01_MC.Buffer();
            for( Int i=0; i<recvSize; ++i )
                p01_MC_Buf[i] += recvBuf[i];
        }

        if( k > 0 )
        {
            // This is not the last iteration of the panel factorization, 
            // Reduce to one p01[MC] to the next process column.
            const Int a01LocalHeight = a01.LocalHeight();

            const Int nextProcRow = Mod( alpha11.ColAlign()-1, r );
            const Int nextProcCol = Mod( alpha11.RowAlign()-1, r );

            vector<F> reduceToOneRecvBuf;
            FastResize( reduceToOneRecvBuf, a01LocalHeight );

            mpi::Reduce
            ( p01_MC.Buffer(), reduceToOneRecvBuf.data(),
              a01LocalHeight, nextProcCol, g.RowComm() );
            if( g.Col() == nextProcCol )
            {
                // Finish computing w01. During its computation, ensure that 
                // every process has a copy of the last element of the w01.
                // We know a priori that the last element of a01 is one.
                const F* a01_MC_Buf = a01_MC.Buffer();
                F myDotProduct = blas::Dot
                    ( a01LocalHeight, reduceToOneRecvBuf.data(), 1, 
                                      a01_MC_Buf,                1 );
                F sendBuf[2], recvBuf[2];
                sendBuf[0] = myDotProduct;
                sendBuf[1] = ( g.Row()==nextProcRow ? 
                               reduceToOneRecvBuf[a01LocalHeight-1] : 0 );
                mpi::AllReduce( sendBuf, recvBuf, 2, g.ColComm() );
                F dotProduct = recvBuf[0];

                // Set up for the next iteration by filling in the values for:
                // - w01LastBuf
                // - w01LastBottomEntry
                F scale = dotProduct*tau / F(2);
                for( Int i=0; i<a01LocalHeight; ++i )
                    w01LastBuf[i] = Conj(tau)*
                        ( reduceToOneRecvBuf[i]-scale*a01_MC_Buf[i] );
                w01LastBottomEntry = Conj(tau)*( recvBuf[1]-scale );
            }
        }
        else
        {
            // This is the last iteration, our last task is to finish forming
            // w01[MC] and w01[MR]
            const Int a01LocalHeight = a01.LocalHeight();

            // AllReduce sum p01[MC] over process rows
            vector<F> allReduceRecvBuf;
            FastResize( allReduceRecvBuf, a01LocalHeight );

            mpi::AllReduce
            ( p01_MC.Buffer(), allReduceRecvBuf.data(), 
              a01LocalHeight, g.RowComm() );

            // Finish computing w01. During its computation, ensure that 
            // every process has a copy of the last element of the w01.
            // We know a priori that the last element of a01 is one.
            const F* a01_MC_Buf = a01_MC.Buffer();
            F myDotProduct = blas::Dot
                ( a01LocalHeight, allReduceRecvBuf.data(), 1, 
                                  a01_MC_Buf,              1 );
            const F dotProduct = mpi::AllReduce( myDotProduct, g.ColComm() );

            // Grab views into W[MC,* ] and W[MR,* ]
            auto w01_MC = W_MC_STAR( ind0, ind1-off );
            auto w01_MR = W_MR_STAR( ind0, ind1-off );

            // Store w01[MC]
            F scale = dotProduct*tau / F(2);
            F* w01_MC_Buf = w01_MC.Buffer();
            for( Int i=0; i<a01LocalHeight; ++i )
                w01_MC_Buf[i] = Conj(tau)*
                    ( allReduceRecvBuf[i]-scale*a01_MC_Buf[i] );

            // Fast transpose w01[MC] -> w01[MR]
            if( onDiagonal )
            {
                MemCopy
                ( w01_MR.Buffer(), w01_MC.Buffer(), a01LocalHeight );
            }
            else
            {
                // Pairwise exchange with the transpose process
                const Int sendSize = A00.LocalHeight();
                const Int recvSize = A00.LocalWidth();
                mpi::SendRecv
                ( w01_MC.Buffer(), sendSize, transposeRank,
                  w01_MR.Buffer(), recvSize, transposeRank, g.VCComm() );
            }
        }
    }

    SetRealPartOfDiagonal( expandedABR, e, 1 );
}

} // namespace herm_tridiag
} // namespace El

#endif // ifndef EL_HERMITIANTRIDIAG_UPPER_PANEL_SQUARE_HPP
