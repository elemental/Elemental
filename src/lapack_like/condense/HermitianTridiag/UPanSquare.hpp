/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HERMITIANTRIDIAG_UPANSQUARE_HPP
#define EL_HERMITIANTRIDIAG_UPANSQUARE_HPP

namespace El {
namespace herm_tridiag {

template<typename F>
void UPanSquare
( DistMatrix<F>& A,
  DistMatrix<F>& W,
  DistMatrix<F,MD,STAR>& t,
  DistMatrix<F,MC,STAR>& B_MC_STAR, 
  DistMatrix<F,MR,STAR>& B_MR_STAR,
  DistMatrix<F,MC,STAR>& W_MC_STAR,
  DistMatrix<F,MR,STAR>& W_MR_STAR,
  const SymvCtrl<F>& ctrl )
{
    const Int n = A.Height();
    const Int nW = W.Width();
    DEBUG_ONLY(
        CallStackEntry cse("herm_tridiag::UPanSquare");
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

    vector<F> w01LastBuffer(n/r+1);
    DistMatrix<F> w01Last(g);
    DistMatrix<F,MC,STAR> a01_MC_STAR(g), p01_MC_STAR(g),
                          a01Last_MC_STAR(g), w01Last_MC_STAR(g);
    DistMatrix<F,MR,STAR> a01_MR_STAR(g), q01_MR_STAR(g),
                          a01Last_MR_STAR(g), w01Last_MR_STAR(g),
                          x21_MR_STAR(g), y21_MR_STAR(g);

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

        auto a01T     = A( IR(0,kA-1),  ind1 );
        auto alpha01B = A( IR(kA-1,kA), ind1 );
        auto A02T     = A( IR(0,off),   ind2 );

        auto W02T = W( IR(0,off), IR(k+1,nW) );

        auto tau1     = t( ind1-off, IR(0,1) );
        auto epsilon1 = e( ind1-off, IR(0,1) );

        a01_MC_STAR.AlignWith( A00 );
        a01_MR_STAR.AlignWith( A00 );
        p01_MC_STAR.AlignWith( A00 );
        q01_MR_STAR.AlignWith( A00 );
        a01_MC_STAR.Resize( kA, 1 );
        a01_MR_STAR.Resize( kA, 1 );
        Zeros( p01_MC_STAR, kA, 1 );
        Zeros( q01_MR_STAR, kA, 1 );

        // View the portions of A02 and W0T outside of this panel's square
        auto a01T_MC_STAR = a01_MC_STAR( IR(0,off), IR(0,1) );
        auto p01T_MC_STAR = p01_MC_STAR( IR(0,off), IR(0,1) ); 

        if( !firstIteration )
        {
            // TODO: Move these and make them auto
            View( a01Last_MC_STAR, B_MC_STAR, indT, IR(k+1,k+2) );
            View( a01Last_MR_STAR, B_MR_STAR, indT, IR(k+1,k+2) );
            View( w01Last,         W,         indT, IR(k+1,k+2) );
        }

        const bool thisIsMyCol = ( g.Col() == alpha11.RowAlign() );
        if( thisIsMyCol )
        {
            if( !firstIteration )
            {
                // Finish updating the current column with two axpy's
                const Int aT1LocalHeight = aT1.LocalHeight();
                F* aT1Buffer = aT1.Buffer();
                const F* a01Last_MC_STAR_Buffer = a01Last_MC_STAR.Buffer();
                for( Int i=0; i<aT1LocalHeight; ++i )
                    aT1Buffer[i] -=
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
        GetRealPartOfDiagonal( alpha01B, epsilon1 );
        alpha01B.Set( 0, 0, F(1) );

        // If this is the first iteration, have each member of the owning 
        // process column broadcast tau and a01 within its process row. 
        // Otherwise, also add w01 into the broadcast.
        if( firstIteration )
        {
            const Int a01LocalHeight = a01.LocalHeight();
            vector<F> rowBroadcastBuffer(a01LocalHeight+1);
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
            // Store a01[MC,* ] into B[MC,* ]
            MemCopy
            ( B_MC_STAR.Buffer(0,k), 
              rowBroadcastBuffer.data(), a01LocalHeight );
            // Store tau
            tau = rowBroadcastBuffer[a01LocalHeight];
            
            // Take advantage of the square process grid in order to form
            // a01[MR,* ] from a01[MC,* ]
            if( onDiagonal )
            {
                MemCopy
                ( a01_MR_STAR.Buffer(),
                  a01_MC_STAR.Buffer(), a01LocalHeight );
            }
            else
            {
                // Pairwise exchange
                const Int sendSize = A00.LocalHeight();
                const Int recvSize = A00.LocalWidth();
                mpi::SendRecv
                ( a01_MC_STAR.Buffer(), sendSize, transposeRank,
                  a01_MR_STAR.Buffer(), recvSize, transposeRank, g.VCComm() );
            }
            // Store a01[MR,* ]
            MemCopy
            ( B_MR_STAR.Buffer(0,k),
              a01_MR_STAR.Buffer(), a01_MR_STAR.LocalHeight() );
        }
        else
        {
            const Int a01LocalHeight = a01.LocalHeight();
            const Int w01LastLocalHeight = aT1.LocalHeight();
            vector<F> 
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
            // Store a01[MC,* ] into B[MC,* ]
            MemCopy
            ( B_MC_STAR.Buffer(0,k), 
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
            ( W_MC_STAR.Buffer(0,k+1),
              &rowBroadcastBuffer[a01LocalHeight], w01LastLocalHeight );
            if( g.Col() == w01Last.RowAlign() )
            {
                MemCopy
                ( w01Last.Buffer(),
                  &rowBroadcastBuffer[a01LocalHeight], w01LastLocalHeight );
            }
            // Store tau
            tau = rowBroadcastBuffer[a01LocalHeight+w01LastLocalHeight];

            // Take advantage of the square process grid in order to quickly
            // form a01[MR,* ] and w01Last[MR,* ] from their [MC,* ]
            // counterparts
            w01Last_MR_STAR.AlignWith( A00 );
            w01Last_MR_STAR.Resize( w01Last.Height(), 1 );
            if( onDiagonal )
            {
                MemCopy
                ( a01_MR_STAR.Buffer(),
                  a01_MC_STAR.Buffer(), a01LocalHeight );
                MemCopy
                ( w01Last_MR_STAR.Buffer(),
                  w01Last_MC_STAR.Buffer(), w01LastLocalHeight );
            }
            else
            {
                const Int sendSize = A00.LocalHeight()+ATL.LocalHeight();
                const Int recvSize = A00.LocalWidth()+ATL.LocalWidth();
                vector<F> sendBuffer(sendSize), recvBuffer(recvSize);

                // Pack the send buffer
                MemCopy
                ( sendBuffer.data(), a01_MC_STAR.Buffer(), A00.LocalHeight() );
                MemCopy
                ( &sendBuffer[A00.LocalHeight()],
                  w01Last_MC_STAR.Buffer(), ATL.LocalHeight() );

                // Pairwise exchange
                mpi::SendRecv
                ( sendBuffer.data(), sendSize, transposeRank,
                  recvBuffer.data(), recvSize, transposeRank, g.VCComm() );

                // Unpack the recv buffer
                MemCopy
                ( a01_MR_STAR.Buffer(), recvBuffer.data(), A00.LocalWidth() );
                MemCopy
                ( w01Last_MR_STAR.Buffer(),
                  &recvBuffer[A00.LocalWidth()], ATL.LocalWidth() );
            }

            // Store w01Last[MR,* ]
            MemCopy
            ( W_MR_STAR.Buffer(0,k+1),
              w01Last_MR_STAR.Buffer(),
              w01Last_MR_STAR.LocalHeight() );
            // Store a01[MR,* ]
            MemCopy
            ( B_MR_STAR.Buffer(0,k),
              a01_MR_STAR.Buffer(), a01_MR_STAR.LocalHeight() );

            // Update the portion of A00 that is in our current panel with 
            // w01Last and a01Last using two gers. We do not need their bottom
            // entries. We trash the lower triangle of our panel of A since we 
            // are only doing slightly more work and we can replace it
            // afterwards.
            auto a01Last_MC_STAR_Top = a01Last_MC_STAR( ind0, IR(0,1) );
            auto w01Last_MC_STAR_Top = w01Last_MC_STAR( ind0, IR(0,1) );
            auto a01Last_MR_STAR_TopPan = 
                a01Last_MR_STAR( IR(off,kA), IR(0,1) );
            auto w01Last_MR_STAR_TopPan = 
                w01Last_MR_STAR( IR(off,kA), IR(0,1) );
            const F* a01_MC_STAR_Buffer = a01Last_MC_STAR_Top.Buffer();
            const F* w01_MC_STAR_Buffer = w01Last_MC_STAR_Top.Buffer();
            const F* a01_MR_STAR_Buffer = a01Last_MR_STAR_TopPan.Buffer();
            const F* w01_MR_STAR_Buffer = w01Last_MR_STAR_TopPan.Buffer();
            F* A00PanBuffer = A00Pan.Buffer();
            const Int localHeight = A00Pan.LocalHeight();
            const Int localWidth = A00Pan.LocalWidth();
            const Int lDim = A00Pan.LDim();
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const F delta = Conj(a01_MR_STAR_Buffer[jLoc]);
                const F gamma = Conj(w01_MR_STAR_Buffer[jLoc]);
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    A00PanBuffer[iLoc+jLoc*lDim] -=
                        w01_MC_STAR_Buffer[iLoc]*delta +
                        a01_MC_STAR_Buffer[iLoc]*gamma;
            }
        }

        // Form the local portions of (A00 a01) into p01[MC,* ] and q01[MR,* ]:
        //   p01[MC,* ] := triu(A00)[MC,MR] a01[MR,* ]
        //   q01[MR,* ] := triu(A00,+1)'[MR,MC] a01[MC,* ]
        symv::LocalColAccumulate
        ( UPPER, F(1),
          A00, a01_MC_STAR, a01_MR_STAR, p01_MC_STAR, q01_MR_STAR, true,
          ctrl );

        x21_MR_STAR.AlignWith( A02T );
        y21_MR_STAR.AlignWith( A02T );
        Zeros( x21_MR_STAR, A02.Width(), 1 );
        Zeros( y21_MR_STAR, A02.Width(), 1 );
        LocalGemv( ADJOINT, F(1), W02T, a01T_MC_STAR, F(0), x21_MR_STAR );
        LocalGemv( ADJOINT, F(1), A02T, a01T_MC_STAR, F(0), y21_MR_STAR );

        // Combine the AllReduce column summations of x21[MR,* ] and y21[MR,* ]
        {
            const Int x21LocalHeight = x21_MR_STAR.LocalHeight();
            const Int y21LocalHeight = y21_MR_STAR.LocalHeight();
            const Int reduceSize = x21LocalHeight+y21LocalHeight;
            vector<F> colSumSendBuffer(reduceSize),
                      colSumRecvBuffer(reduceSize);
            MemCopy
            ( colSumSendBuffer.data(), x21_MR_STAR.Buffer(), x21LocalHeight );
            MemCopy
            ( &colSumSendBuffer[x21LocalHeight],
              y21_MR_STAR.Buffer(), y21LocalHeight );
            mpi::AllReduce
            ( colSumSendBuffer.data(), colSumRecvBuffer.data(),
              reduceSize, g.ColComm() );
            MemCopy
            ( x21_MR_STAR.Buffer(), colSumRecvBuffer.data(), x21LocalHeight );
            MemCopy
            ( y21_MR_STAR.Buffer(), 
              &colSumRecvBuffer[x21LocalHeight], y21LocalHeight );
        }

        LocalGemv( NORMAL, F(-1), A02T, x21_MR_STAR, F(1), p01T_MC_STAR );
        LocalGemv( NORMAL, F(-1), W02T, y21_MR_STAR, F(1), p01T_MC_STAR );

        // Fast transpose the unsummed q01[MR,* ] -> q01[MC,* ], so that
        // it needs to be summed over process rows instead of process 
        // columns. We immediately add it onto p01[MC,* ], which also needs
        // to be summed within process rows.
        if( onDiagonal )
        {
            const Int a01LocalHeight = a01.LocalHeight();
            F* p01_MC_STAR_Buffer = p01_MC_STAR.Buffer();
            const F* q01_MR_STAR_Buffer = q01_MR_STAR.Buffer();
            for( Int i=0; i<a01LocalHeight; ++i )
                p01_MC_STAR_Buffer[i] += q01_MR_STAR_Buffer[i];
        }
        else
        {
            // Pairwise exchange with the transpose process
            const Int sendSize = A00.LocalWidth();
            const Int recvSize = A00.LocalHeight();
            vector<F> recvBuffer(recvSize);
            mpi::SendRecv
            ( q01_MR_STAR.Buffer(), sendSize, transposeRank,
              recvBuffer.data(),    recvSize, transposeRank, g.VCComm() );

            // Unpack the recv buffer directly onto p01[MC,* ]
            F* p01_MC_STAR_Buffer = p01_MC_STAR.Buffer();
            for( Int i=0; i<recvSize; ++i )
                p01_MC_STAR_Buffer[i] += recvBuffer[i];
        }

        if( k > 0 )
        {
            // This is not the last iteration of the panel factorization, 
            // Reduce to one p01[MC,* ] to the next process column.
            const Int a01LocalHeight = a01.LocalHeight();

            const Int nextProcessRow = (alpha11.ColAlign()+r-1) % r;
            const Int nextProcessCol = (alpha11.RowAlign()+r-1) % r;

            vector<F> reduceToOneRecvBuffer(a01LocalHeight);
            mpi::Reduce
            ( p01_MC_STAR.Buffer(), reduceToOneRecvBuffer.data(),
              a01LocalHeight, nextProcessCol, g.RowComm() );
            if( g.Col() == nextProcessCol )
            {
                // Finish computing w01. During its computation, ensure that 
                // every process has a copy of the last element of the w01.
                // We know a priori that the last element of a01 is one.
                const F* a01_MC_STAR_Buffer = a01_MC_STAR.Buffer();
                F myDotProduct = blas::Dot
                    ( a01LocalHeight, reduceToOneRecvBuffer.data(), 1, 
                                      a01_MC_STAR_Buffer,           1 );
                F sendBuffer[2], recvBuffer[2];
                sendBuffer[0] = myDotProduct;
                sendBuffer[1] = ( g.Row()==nextProcessRow ? 
                                  reduceToOneRecvBuffer[a01LocalHeight-1] : 0 );
                mpi::AllReduce( sendBuffer, recvBuffer, 2, g.ColComm() );
                F dotProduct = recvBuffer[0];

                // Set up for the next iteration by filling in the values for:
                // - w01LastBuffer
                // - w01LastBottomEntry
                F scale = dotProduct*tau / F(2);
                for( Int i=0; i<a01LocalHeight; ++i )
                    w01LastBuffer[i] = Conj(tau)*
                        ( reduceToOneRecvBuffer[i]-
                          scale*a01_MC_STAR_Buffer[i] );
                w01LastBottomEntry = Conj(tau)*( recvBuffer[1]-scale );
            }
        }
        else
        {
            // This is the last iteration, our last task is to finish forming
            // w01[MC,* ] and w01[MR,* ]
            const Int a01LocalHeight = a01.LocalHeight();

            // AllReduce sum p01[MC,* ] over process rows
            vector<F> allReduceRecvBuffer(a01LocalHeight);
            mpi::AllReduce
            ( p01_MC_STAR.Buffer(), allReduceRecvBuffer.data(), 
              a01LocalHeight, g.RowComm() );

            // Finish computing w01. During its computation, ensure that 
            // every process has a copy of the last element of the w01.
            // We know a priori that the last element of a01 is one.
            const F* a01_MC_STAR_Buffer = a01_MC_STAR.Buffer();
            F myDotProduct = blas::Dot
                ( a01LocalHeight, allReduceRecvBuffer.data(), 1, 
                                  a01_MC_STAR_Buffer,         1 );
            const F dotProduct = mpi::AllReduce( myDotProduct, g.ColComm() );

            // Grab views into W[MC,* ] and W[MR,* ]
            auto w01_MC_STAR = W_MC_STAR( ind0, ind1-off );
            auto w01_MR_STAR = W_MR_STAR( ind0, ind1-off );

            // Store w01[MC,* ]
            F scale = dotProduct*tau / F(2);
            F* w01_MC_STAR_Buffer = w01_MC_STAR.Buffer();
            for( Int i=0; i<a01LocalHeight; ++i )
                w01_MC_STAR_Buffer[i] = Conj(tau)*
                    ( allReduceRecvBuffer[i]-scale*a01_MC_STAR_Buffer[i] );

            // Fast transpose w01[MC,* ] -> w01[MR,* ]
            if( onDiagonal )
            {
                MemCopy
                ( w01_MR_STAR.Buffer(),
                  w01_MC_STAR.Buffer(), a01LocalHeight );
            }
            else
            {
                // Pairwise exchange with the transpose process
                const Int sendSize = A00.LocalHeight();
                const Int recvSize = A00.LocalWidth();
                mpi::SendRecv
                ( w01_MC_STAR.Buffer(), sendSize, transposeRank,
                  w01_MR_STAR.Buffer(), recvSize, transposeRank, g.VCComm() );
            }
        }
    }

    SetRealPartOfDiagonal( expandedABR, e, 1 );
}

} // namespace herm_tridiag
} // namespace El

#endif // ifndef EL_HERMITIANTRIDIAG_UPANSQUARE_HPP
