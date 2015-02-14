/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HERMITIANTRIDIAG_LPANSQUARE_HPP
#define EL_HERMITIANTRIDIAG_LPANSQUARE_HPP

namespace El {
namespace herm_tridiag {

// TODO: Fuse Trmvs
// TODO: Custom kernel for W22 update
// TODO: Reuse a small number of preallocated vectors

template<typename F>
void LPanSquare
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
        CallStackEntry cse("herm_tridiag::LPanSquare");
        AssertSameGrids( A, W, t );
        if( n != A.Width() )
            LogicError("A must be square");
        if( n != W.Height() )
            LogicError("A and W must be the same height");
        if( n <= nW )
            LogicError("W must be a column panel");
        if( W.ColAlign() != A.ColAlign() || W.RowAlign() != A.RowAlign() )
            LogicError("W and A must be aligned");
        if( t.Height() != nW || t.Width() != 1 )
            LogicError
            ("t must be a column vector of the same length as W's width");
        if( !A.DiagonalAlignedWith(t,-1) )
            LogicError("t is not aligned with A's subdiagonal");
    )
    typedef Base<F> Real;
    // Find the process holding our transposed data
    const Grid& g = A.Grid();
    const Int r = g.Height();
    const Int transposeRow = Mod( A.ColAlign()+A.RowShift(), r );
    const Int transposeCol = Mod( A.RowAlign()+A.ColShift(), r );
    const Int transposeRank = transposeRow + r*transposeCol;
    const bool onDiagonal = ( transposeRank == g.VCRank() );

    // Create a distributed matrix for storing the subdiagonal
    DistMatrix<Real,MD,STAR> e(g);
    e.SetRoot( A.DiagonalRoot(-1) );
    e.AlignCols( A.DiagonalAlign(-1) );
    e.Resize( nW, 1 );

    vector<F> w21LastBuffer(A.Height()/r+1);
    DistMatrix<F> w21Last(g);
    DistMatrix<F,MC,STAR> a21_MC_STAR(g), p21_MC_STAR(g), 
                          a21Last_MC_STAR(g), w21Last_MC_STAR(g);
    DistMatrix<F,MR,STAR> q21_MR_STAR(g), a21_MR_STAR(g),
                          x01_MR_STAR(g), y01_MR_STAR(g),
                          a21Last_MR_STAR(g), w21Last_MR_STAR(g);

    F tau = 0;
    F w21LastFirstEntry = 0;
    for( Int k=0; k<nW; ++k )
    {
        const Range<Int> ind0( 0,   k   ),
                         ind1( k,   k+1 ),
                         indB( k,   n   ), indR( k, n ),
                         ind2( k+1, n   );
           
        auto A00     = A( ind0, ind0 );
        auto alpha11 = A( ind1, ind1 );
        auto aB1     = A( indB, ind1 );
        auto ABR     = A( indB, indR );
        auto a21     = A( ind2, ind1 );
        auto A22     = A( ind2, ind2 );

        auto alpha21T = A( IR(k+1,k+2), ind1 );
        auto A20B     = A( IR(nW,n),    ind0 );
        auto a21B     = A( IR(k+2,n),   ind1 );

        auto W00 = W( ind0, ind0       );
        auto W22 = W( ind2, IR(k+1,nW) );

        auto W20B = W( IR(nW,n), ind0 );

        auto tau1     = t( ind1, IR(0,1) );
        auto epsilon1 = e( ind1, IR(0,1) );

        a21_MC_STAR.AlignWith( A22 );
        a21_MR_STAR.AlignWith( A22 );
        p21_MC_STAR.AlignWith( A22 );
        q21_MR_STAR.AlignWith( A22 );
        a21_MC_STAR.Resize( n-(k+1), 1 );
        a21_MR_STAR.Resize( n-(k+1), 1 );
        Zeros( p21_MC_STAR, n-(k+1), 1 );
        Zeros( q21_MR_STAR, n-(k+1), 1 );

        // View the portions of a21[MC,* ] and p21[MC,* ] below the current
        // panel's square
        auto a21B_MC_STAR = a21_MC_STAR( IR(nW-(k+1),n-(k+1)), IR(0,1) );
        auto p21B_MC_STAR = p21_MC_STAR( IR(nW-(k+1),n-(k+1)), IR(0,1) );

        if( k > 0 )
        {
            // TODO: Move these and make them auto
            View( a21Last_MC_STAR,  B_MC_STAR, indB, ind1-1 ); 
            View( a21Last_MR_STAR,  B_MR_STAR, indB, ind1-1 );
            View( w21Last,          W,         indB, ind1-1 );
        }
 
        const bool thisIsMyCol = ( g.Col() == alpha11.RowAlign() );
        if( thisIsMyCol )
        {
            if( k > 0 )
            {
                // Finish updating the current column with two axpy's
                const Int aB1LocalHeight = aB1.LocalHeight();
                F* aB1Buffer = aB1.Buffer();
                const F* a21Last_MC_STAR_Buffer = a21Last_MC_STAR.Buffer();
                for( Int i=0; i<aB1LocalHeight; ++i )
                    aB1Buffer[i] -=
                        w21LastBuffer[i] + 
                        a21Last_MC_STAR_Buffer[i]*Conj(w21LastFirstEntry);
            }
            // Compute the Householder reflector
            tau = reflector::Col( alpha21T, a21B );
            tau1.Set( 0, 0, tau );
        }

        // Store the subdiagonal value and turn a21 into a proper scaled 
        // reflector by explicitly placing the implicit one in its first entry.
        GetRealPartOfDiagonal( alpha21T, epsilon1 );
        alpha21T.Set( 0, 0, F(1) );

        // If this is the first iteration, have each member of the owning 
        // process column broadcast tau and a21 within its process row. 
        // Otherwise, also add w21 into the broadcast.
        if( k == 0 )
        {
            const Int a21LocalHeight = a21.LocalHeight();
            vector<F> rowBroadcastBuffer(a21LocalHeight+1);
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
            // Store a21[MC,* ] into B[MC,* ]
            const Int B_MC_STAR_Off = B_MC_STAR.LocalHeight()-a21LocalHeight;
            MemCopy
            ( B_MC_STAR.Buffer(B_MC_STAR_Off,0), 
              rowBroadcastBuffer.data(),
              B_MC_STAR.LocalHeight()-B_MC_STAR_Off );
            // Store tau
            tau = rowBroadcastBuffer[a21LocalHeight];
            
            // Take advantage of the square process grid in order to form 
            // a21[MR,* ] from a21[MC,* ]
            if( onDiagonal )
                MemCopy
                ( a21_MR_STAR.Buffer(), a21_MC_STAR.Buffer(), a21LocalHeight );
            else
                mpi::SendRecv
                ( a21_MC_STAR.Buffer(), A22.LocalHeight(), transposeRank,
                  a21_MR_STAR.Buffer(), A22.LocalWidth(),  transposeRank, 
                  g.VCComm() );

            // Store a21[MR,* ]
            const Int B_MR_STAR_Off = 
                B_MR_STAR.LocalHeight()-a21_MR_STAR.LocalHeight();
            MemCopy
            ( B_MR_STAR.Buffer(B_MR_STAR_Off,A00.Width()),
              a21_MR_STAR.Buffer(),
              B_MR_STAR.LocalHeight()-B_MR_STAR_Off );
        }
        else
        {
            const Int a21LocalHeight = a21.LocalHeight();
            const Int w21LastLocalHeight = aB1.LocalHeight();
            vector<F> 
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
            ( a21_MC_STAR.Buffer(), rowBroadcastBuffer.data(), a21LocalHeight );
            // Store a21[MC,* ] into B[MC,* ]
            const Int B_MC_STAR_Off = B_MC_STAR.LocalHeight()-a21LocalHeight;
            MemCopy
            ( B_MC_STAR.Buffer(B_MC_STAR_Off,A00.Width()), 
              rowBroadcastBuffer.data(),
              B_MC_STAR.LocalHeight()-B_MC_STAR_Off );
            // Store w21Last[MC,* ] into its DistMatrix class
            w21Last_MC_STAR.AlignWith( alpha11 );
            w21Last_MC_STAR.Resize( a21.Height()+1, 1 );
            MemCopy
            ( w21Last_MC_STAR.Buffer(), 
              &rowBroadcastBuffer[a21LocalHeight], w21LastLocalHeight );
            // Store the bottom part of w21Last[MC,* ] into WB[MC,* ] and, 
            // if necessary, w21.
            const Int W_MC_STAR_Off = 
                W_MC_STAR.LocalHeight()-w21LastLocalHeight;
            MemCopy
            ( W_MC_STAR.Buffer(W_MC_STAR_Off,A00.Width()-1),
              &rowBroadcastBuffer[a21LocalHeight],
              W_MC_STAR.LocalHeight()-W_MC_STAR_Off );
            if( g.Col() == w21Last.RowAlign() )
            {
                MemCopy
                ( w21Last.Buffer(),
                  &rowBroadcastBuffer[a21LocalHeight], w21LastLocalHeight );
            }
            // Store tau
            tau = rowBroadcastBuffer[a21LocalHeight+w21LastLocalHeight];

            // Take advantage of the square process grid in order to quickly
            // form a21[MR,* ] and w21Last[MR,* ] from their [MC,* ] 
            // counterparts
            w21Last_MR_STAR.AlignWith( ABR );
            w21Last_MR_STAR.Resize( w21Last.Height(), 1 );
            if( onDiagonal )
            {
                MemCopy
                ( a21_MR_STAR.Buffer(), a21_MC_STAR.Buffer(), a21LocalHeight );
                MemCopy
                ( w21Last_MR_STAR.Buffer(),
                  w21Last_MC_STAR.Buffer(), w21LastLocalHeight );
            }
            else
            {
                const Int sendSize = A22.LocalHeight()+ABR.LocalHeight();
                const Int recvSize = A22.LocalWidth()+ABR.LocalWidth();
                vector<F> sendBuffer(sendSize), recvBuffer(recvSize);

                // Pack the send buffer
                MemCopy
                ( sendBuffer.data(), a21_MC_STAR.Buffer(), A22.LocalHeight() );
                MemCopy
                ( &sendBuffer[A22.LocalHeight()],
                  w21Last_MC_STAR.Buffer(), ABR.LocalHeight() );

                // Pairwise exchange
                mpi::SendRecv
                ( sendBuffer.data(), sendSize, transposeRank,
                  recvBuffer.data(), recvSize, transposeRank, g.VCComm() );

                // Unpack the recv buffer
                MemCopy
                ( a21_MR_STAR.Buffer(), recvBuffer.data(), A22.LocalWidth() );
                MemCopy
                ( w21Last_MR_STAR.Buffer(),
                  &recvBuffer[A22.LocalWidth()], ABR.LocalWidth() );
            }

            // Store w21Last[MR,* ]
            const Int W_MR_STAR_Off = 
                W_MR_STAR.LocalHeight()-w21Last_MR_STAR.LocalHeight();
            MemCopy
            ( W_MR_STAR.Buffer(W_MR_STAR_Off,A00.Width()-1),
              w21Last_MR_STAR.Buffer(),
              W_MR_STAR.LocalHeight()-W_MR_STAR_Off );
            // Store a21[MR,* ]
            const Int B_MR_STAR_Off = 
                B_MR_STAR.LocalHeight()-a21_MR_STAR.LocalHeight();
            MemCopy
            ( B_MR_STAR.Buffer(B_MR_STAR_Off,A00.Width()),
              a21_MR_STAR.Buffer(),
              B_MR_STAR.LocalHeight()-B_MR_STAR_Off );

            // Update the portion of A22 that is in our current panel with 
            // w21Last and a21Last using two gers. We do not need their top 
            // entries. We trash the upper triangle of our panel of A since we 
            // are only doing slightly more work and we can replace it
            // afterwards.
            // TODO: Create a custom kernel for this operation
            auto a21Last_MC_STAR_Bottom = a21Last_MC_STAR( IR(1,n-k), IR(0,1) );
            auto w21Last_MC_STAR_Bottom = w21Last_MC_STAR( IR(1,n-k), IR(0,1) );
            auto a21Last_MR_STAR_Bottom = a21Last_MR_STAR( IR(1,n-k), IR(0,1) );
            auto w21Last_MR_STAR_Bottom = w21Last_MR_STAR( IR(1,n-k), IR(0,1) );
            const F* a21_MC_STAR_Buffer = a21Last_MC_STAR_Bottom.Buffer();
            const F* w21_MC_STAR_Buffer = w21Last_MC_STAR_Bottom.Buffer();
            const F* a21_MR_STAR_Buffer = a21Last_MR_STAR_Bottom.Buffer();
            const F* w21_MR_STAR_Buffer = w21Last_MR_STAR_Bottom.Buffer();
            F* A22Buffer = A22.Buffer();
            const Int localHeight = W22.LocalHeight();
            const Int localWidth = W22.LocalWidth();
            const Int lDim = A22.LDim();
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const F delta = Conj(a21_MR_STAR_Buffer[jLoc]);
                const F gamma = Conj(w21_MR_STAR_Buffer[jLoc]);
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    A22Buffer[iLoc+jLoc*lDim] -=
                        w21_MC_STAR_Buffer[iLoc]*delta +
                        a21_MC_STAR_Buffer[iLoc]*gamma;
            }
        }

        // Form the local portions of (A22 a21) into p21[MC,* ] and q21[MR,* ]:
        //   p21[MC,* ] := tril(A22)[MC,MR] a21[MR,* ]
        //   q21[MR,* ] := tril(A22,-1)'[MR,MC] a21[MC,* ]
        symv::LocalColAccumulate
        ( LOWER, F(1),
          A22, a21_MC_STAR, a21_MR_STAR, p21_MC_STAR, q21_MR_STAR, true,
          ctrl );

        x01_MR_STAR.AlignWith( W20B );
        y01_MR_STAR.AlignWith( W20B );
        Zeros( x01_MR_STAR, W20B.Width(), 1 );
        Zeros( y01_MR_STAR, W20B.Width(), 1 );
        LocalGemv( ADJOINT, F(1), W20B, a21B_MC_STAR, F(0), x01_MR_STAR );
        LocalGemv( ADJOINT, F(1), A20B, a21B_MC_STAR, F(0), y01_MR_STAR );

        // Combine the AllReduce column summations of x01[MR,* ], y01[MR,* ]
        {
            const Int x01LocalHeight = x01_MR_STAR.LocalHeight();
            vector<F> colSumSendBuffer(2*x01LocalHeight),
                      colSumRecvBuffer(2*x01LocalHeight);
            MemCopy
            ( colSumSendBuffer.data(), 
              x01_MR_STAR.Buffer(), x01LocalHeight );
            MemCopy
            ( &colSumSendBuffer[x01LocalHeight],
              y01_MR_STAR.Buffer(), x01LocalHeight );
            mpi::AllReduce
            ( colSumSendBuffer.data(), 
              colSumRecvBuffer.data(), 2*x01LocalHeight, g.ColComm() );
            MemCopy
            ( x01_MR_STAR.Buffer(), 
              colSumRecvBuffer.data(), x01LocalHeight );
            MemCopy
            ( y01_MR_STAR.Buffer(), 
              &colSumRecvBuffer[x01LocalHeight], x01LocalHeight );
        }

        LocalGemv( NORMAL, F(-1), A20B, x01_MR_STAR, F(1), p21B_MC_STAR );
        LocalGemv( NORMAL, F(-1), W20B, y01_MR_STAR, F(1), p21B_MC_STAR );

        // Fast transpose the unsummed q21[MR,* ] -> q21[MC,* ], so that
        // it needs to be summed over process rows instead of process 
        // columns. We immediately add it onto p21[MC,* ], which also needs
        // to be summed within process rows.
        if( onDiagonal )
        {
            blas::Axpy
            ( a21.LocalHeight(), F(1),
              q21_MR_STAR.LockedBuffer(), 1,
              p21_MC_STAR.Buffer(),       1 );
        }
        else
        {
            // Pairwise exchange with the transpose process
            const Int sendSize = A22.LocalWidth();
            const Int recvSize = A22.LocalHeight();
            vector<F> recvBuffer(recvSize);
            mpi::SendRecv
            ( q21_MR_STAR.Buffer(), sendSize, transposeRank, 
              recvBuffer.data(),    recvSize, transposeRank, 
              g.VCComm() );

            // Unpack the recv buffer directly onto p21[MC,* ]
            blas::Axpy
            ( recvSize, F(1), recvBuffer.data(), 1, p21_MC_STAR.Buffer(), 1 );
        }

        if( W22.Width() > 0 )
        {
            // This is not the last iteration of the panel factorization, 
            // Reduce to one p21[MC,* ] to the next process column.
            const Int a21LocalHeight = a21.LocalHeight();

            const Int nextProcessRow = (alpha11.ColAlign()+1) % r;
            const Int nextProcessCol = (alpha11.RowAlign()+1) % r;

            vector<F> reduceToOneRecvBuffer(a21LocalHeight);
            mpi::Reduce
            ( p21_MC_STAR.Buffer(), reduceToOneRecvBuffer.data(),
              a21LocalHeight, nextProcessCol, g.RowComm() );
            if( g.Col() == nextProcessCol )
            {
                // Finish computing w21. During its computation, ensure that 
                // every process has a copy of the first element of the w21.
                // We know a priori that the first element of a21 is one.
                const F* a21_MC_STAR_Buffer = a21_MC_STAR.Buffer();
                F myDotProduct = blas::Dot
                    ( a21LocalHeight, reduceToOneRecvBuffer.data(), 1,
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
                for( Int i=0; i<a21LocalHeight; ++i )
                    w21LastBuffer[i] = Conj(tau)*
                        ( reduceToOneRecvBuffer[i]-
                          scale*a21_MC_STAR_Buffer[i] );
                w21LastFirstEntry = Conj(tau)*( recvBuffer[1]-scale );
            }
        }
        else
        {
            // This is the last iteration, so we just need to form w21[MC,* ]
            // and w21[MR,* ]. 
            const Int a21LocalHeight = a21.LocalHeight();

            // AllReduce sum p21[MC,* ] over process rows
            vector<F> allReduceRecvBuffer(a21LocalHeight);
            mpi::AllReduce
            ( p21_MC_STAR.Buffer(), allReduceRecvBuffer.data(),
              a21LocalHeight, g.RowComm() );

            // Finish computing w21.
            const F* a21_MC_STAR_Buffer = a21_MC_STAR.Buffer();
            F myDotProduct = blas::Dot
                ( a21LocalHeight, allReduceRecvBuffer.data(), 1,
                                  a21_MC_STAR_Buffer,         1 );
            const F dotProduct = mpi::AllReduce( myDotProduct, g.ColComm() );

            // Grab views into W[MC,* ] and W[MR,* ]
            auto w21_MC_STAR = W_MC_STAR( ind2, ind1 );
            auto w21_MR_STAR = W_MR_STAR( ind2, ind1 );

            // Store w21[MC,* ]
            F scale = dotProduct*tau / F(2);
            F* w21_MC_STAR_Buffer = w21_MC_STAR.Buffer();
            for( Int i=0; i<a21LocalHeight; ++i )
                w21_MC_STAR_Buffer[i] = Conj(tau)*
                    ( allReduceRecvBuffer[i]-scale*a21_MC_STAR_Buffer[i] );

            // Fast transpose w21[MC,* ] -> w21[MR,* ]
            if( onDiagonal )
            {
                MemCopy
                ( w21_MR_STAR.Buffer(),
                  w21_MC_STAR.Buffer(), a21LocalHeight );
            }
            else
            {
                // Pairwise exchange with the transpose process
                mpi::SendRecv
                ( w21_MC_STAR.Buffer(), A22.LocalHeight(), transposeRank, 
                  w21_MR_STAR.Buffer(), A22.LocalWidth(),  transposeRank, 
                  g.VCComm() );
            }
        }
    }

    // View the portion of A that e is the subdiagonal of, then place e into it
    auto expandedATL = A( IR(0,nW+1), IR(0,nW+1) );
    SetRealPartOfDiagonal( expandedATL, e, -1 );
}

} // namespace herm_tridiag
} // namespace El

#endif // ifndef EL_HERMITIANTRIDIAG_LPANSQUARE_HPP
