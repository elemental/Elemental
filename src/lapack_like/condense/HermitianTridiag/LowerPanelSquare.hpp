/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HERMITIANTRIDIAG_LOWER_PANEL_SQUARE_HPP
#define EL_HERMITIANTRIDIAG_LOWER_PANEL_SQUARE_HPP

namespace El {
namespace herm_tridiag {

// TODO: Fuse Trmvs
// TODO: Custom kernel for W22 update
// TODO: Reuse a small number of preallocated vectors for buffers

template<typename F>
void LowerPanelSquare
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

    vector<F> w21LastBuf;
    FastResize( w21LastBuf, A.Height()/r+1 );

    DistMatrix<F> w21Last(g);
    DistMatrix<F,MC,STAR> a21_MC(g), p21_MC(g), 
                          a21Last_MC(g), w21Last_MC(g);
    DistMatrix<F,MR,STAR> q21_MR(g), a21_MR(g),
                          x01_MR(g), y01_MR(g),
                          a21Last_MR(g), w21Last_MR(g);

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

        auto alpha21T = A( IR(k+1),     ind1 );
        auto A20B     = A( IR(nW,END),  ind0 );
        auto a21B     = A( IR(k+2,END), ind1 );

        auto W00 = W( ind0, ind0        );
        auto W22 = W( ind2, IR(k+1,END) );

        auto W20B = W( IR(nW,END), ind0 );

        auto tau1     = t( ind1, ALL );
        auto epsilon1 = e( ind1, ALL );

        a21_MC.AlignWith( A22 );
        a21_MR.AlignWith( A22 );
        p21_MC.AlignWith( A22 );
        q21_MR.AlignWith( A22 );
        a21_MC.Resize( n-(k+1), 1 );
        a21_MR.Resize( n-(k+1), 1 );
        Zeros( p21_MC, n-(k+1), 1 );
        Zeros( q21_MR, n-(k+1), 1 );

        // View the portions of a21[MC] and p21[MC] below the current
        // panel's square
        auto a21B_MC = a21_MC( IR(nW-(k+1),END), ALL );
        auto p21B_MC = p21_MC( IR(nW-(k+1),END), ALL );

        if( k > 0 )
        {
            // TODO: Move these and make them auto
            View( a21Last_MC, B_MC_STAR, indB, ind1-1 ); 
            View( a21Last_MR, B_MR_STAR, indB, ind1-1 );
            View( w21Last,    W,         indB, ind1-1 );
        }
 
        const bool thisIsMyCol = ( g.Col() == alpha11.RowAlign() );
        if( thisIsMyCol )
        {
            if( k > 0 )
            {
                // Finish updating the current column with two axpy's
                const Int aB1LocalHeight = aB1.LocalHeight();
                F* aB1Buf = aB1.Buffer();
                const F* a21Last_MC_Buf = a21Last_MC.Buffer();
                for( Int i=0; i<aB1LocalHeight; ++i )
                    aB1Buf[i] -=
                      w21LastBuf[i] + 
                      a21Last_MC_Buf[i]*Conj(w21LastFirstEntry);
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
            vector<F> rowBcastBuf;
            FastResize( rowBcastBuf, a21LocalHeight+1 );

            if( thisIsMyCol )
            {
                // Pack the broadcast buffer with a21 and tau
                MemCopy( rowBcastBuf.data(), a21.Buffer(), a21LocalHeight );
                rowBcastBuf[a21LocalHeight] = tau;
            }
            // Broadcast a21 and tau across the process row
            mpi::Broadcast
            ( rowBcastBuf.data(), 
              a21LocalHeight+1, a21.RowAlign(), g.RowComm() );
            // Store a21[MC] into its DistMatrix class and also store a copy
            // for the next iteration
            MemCopy( a21_MC.Buffer(), rowBcastBuf.data(), a21LocalHeight );
            // Store a21[MC] into B[MC,* ]
            const Int B_MC_STAR_Off = B_MC_STAR.LocalHeight()-a21LocalHeight;
            MemCopy
            ( B_MC_STAR.Buffer(B_MC_STAR_Off,0), 
              rowBcastBuf.data(),
              B_MC_STAR.LocalHeight()-B_MC_STAR_Off );
            // Store tau
            tau = rowBcastBuf[a21LocalHeight];
            
            // Take advantage of the square process grid in order to form 
            // a21[MR] from a21[MC]
            if( onDiagonal )
                MemCopy( a21_MR.Buffer(), a21_MC.Buffer(), a21LocalHeight );
            else
                mpi::SendRecv
                ( a21_MC.Buffer(), A22.LocalHeight(), transposeRank,
                  a21_MR.Buffer(), A22.LocalWidth(),  transposeRank, 
                  g.VCComm() );

            // Store a21[MR]
            const Int B_MR_STAR_Off = 
                B_MR_STAR.LocalHeight()-a21_MR.LocalHeight();
            MemCopy
            ( B_MR_STAR.Buffer(B_MR_STAR_Off,A00.Width()),
              a21_MR.Buffer(),
              B_MR_STAR.LocalHeight()-B_MR_STAR_Off );
        }
        else
        {
            const Int a21LocalHeight = a21.LocalHeight();
            const Int w21LastLocalHeight = aB1.LocalHeight();
            vector<F> rowBcastBuf;
            FastResize( rowBcastBuf, a21LocalHeight+w21LastLocalHeight+1 );
            if( thisIsMyCol ) 
            {
                // Pack the broadcast buffer with a21, w21Last, and tau
                MemCopy
                ( rowBcastBuf.data(), a21.Buffer(), a21LocalHeight );
                MemCopy
                ( &rowBcastBuf[a21LocalHeight], 
                  w21LastBuf.data(), w21LastLocalHeight );
                rowBcastBuf[a21LocalHeight+w21LastLocalHeight] = tau;
            }
            // Broadcast a21, w21Last, and tau across the process row
            mpi::Broadcast
            ( rowBcastBuf.data(), 
              a21LocalHeight+w21LastLocalHeight+1, 
              a21.RowAlign(), g.RowComm() );
            // Store a21[MC] into its DistMatrix class 
            MemCopy( a21_MC.Buffer(), rowBcastBuf.data(), a21LocalHeight );
            // Store a21[MC] into B[MC,* ]
            const Int B_MC_STAR_Off = B_MC_STAR.LocalHeight()-a21LocalHeight;
            MemCopy
            ( B_MC_STAR.Buffer(B_MC_STAR_Off,A00.Width()), 
              rowBcastBuf.data(),
              B_MC_STAR.LocalHeight()-B_MC_STAR_Off );
            // Store w21Last[MC] into its DistMatrix class
            w21Last_MC.AlignWith( alpha11 );
            w21Last_MC.Resize( a21.Height()+1, 1 );
            MemCopy
            ( w21Last_MC.Buffer(), 
              &rowBcastBuf[a21LocalHeight], w21LastLocalHeight );
            // Store the bottom part of w21Last[MC] into WB[MC,* ] and, 
            // if necessary, w21.
            const Int W_MC_STAR_Off = 
                W_MC_STAR.LocalHeight()-w21LastLocalHeight;
            MemCopy
            ( W_MC_STAR.Buffer(W_MC_STAR_Off,A00.Width()-1),
              &rowBcastBuf[a21LocalHeight],
              W_MC_STAR.LocalHeight()-W_MC_STAR_Off );
            if( g.Col() == w21Last.RowAlign() )
            {
                MemCopy
                ( w21Last.Buffer(),
                  &rowBcastBuf[a21LocalHeight], w21LastLocalHeight );
            }
            // Store tau
            tau = rowBcastBuf[a21LocalHeight+w21LastLocalHeight];

            // Take advantage of the square process grid in order to quickly
            // form a21[MR] and w21Last[MR] from their [MC] counterparts
            w21Last_MR.AlignWith( ABR );
            w21Last_MR.Resize( w21Last.Height(), 1 );
            if( onDiagonal )
            {
                MemCopy( a21_MR.Buffer(), a21_MC.Buffer(), a21LocalHeight );
                MemCopy
                ( w21Last_MR.Buffer(),
                  w21Last_MC.Buffer(), w21LastLocalHeight );
            }
            else
            {
                const Int sendSize = A22.LocalHeight()+ABR.LocalHeight();
                const Int recvSize = A22.LocalWidth()+ABR.LocalWidth();
                vector<F> sendBuf, recvBuf;
                FastResize( sendBuf, sendSize );
                FastResize( recvBuf, recvSize );

                // Pack the send buffer
                MemCopy( sendBuf.data(), a21_MC.Buffer(), A22.LocalHeight() );
                MemCopy
                ( &sendBuf[A22.LocalHeight()],
                  w21Last_MC.Buffer(), ABR.LocalHeight() );

                // Pairwise exchange
                mpi::SendRecv
                ( sendBuf.data(), sendSize, transposeRank,
                  recvBuf.data(), recvSize, transposeRank, g.VCComm() );

                // Unpack the recv buffer
                MemCopy
                ( a21_MR.Buffer(), recvBuf.data(), A22.LocalWidth() );
                MemCopy
                ( w21Last_MR.Buffer(),
                  &recvBuf[A22.LocalWidth()], ABR.LocalWidth() );
            }

            // Store w21Last[MR]
            const Int W_MR_STAR_Off = 
                W_MR_STAR.LocalHeight()-w21Last_MR.LocalHeight();
            MemCopy
            ( W_MR_STAR.Buffer(W_MR_STAR_Off,A00.Width()-1),
              w21Last_MR.Buffer(),
              W_MR_STAR.LocalHeight()-W_MR_STAR_Off );
            // Store a21[MR]
            const Int B_MR_STAR_Off = 
                B_MR_STAR.LocalHeight()-a21_MR.LocalHeight();
            MemCopy
            ( B_MR_STAR.Buffer(B_MR_STAR_Off,A00.Width()),
              a21_MR.Buffer(),
              B_MR_STAR.LocalHeight()-B_MR_STAR_Off );

            // Update the portion of A22 that is in our current panel with 
            // w21Last and a21Last using two gers. We do not need their top 
            // entries. We trash the upper triangle of our panel of A since we 
            // are only doing slightly more work and we can replace it
            // afterwards.
            // TODO: Create a custom kernel for this operation
            auto a21Last_MC_Bottom = a21Last_MC( IR(1,n-k), ALL );
            auto w21Last_MC_Bottom = w21Last_MC( IR(1,n-k), ALL );
            auto a21Last_MR_Bottom = a21Last_MR( IR(1,n-k), ALL );
            auto w21Last_MR_Bottom = w21Last_MR( IR(1,n-k), ALL );
            auto A22Sub = A22( IR(0,W22.Height()), IR(0,W22.Width()) );
            Ger2Sub
            ( a21Last_MC_Bottom,
              w21Last_MC_Bottom,
              a21Last_MR_Bottom,
              w21Last_MR_Bottom,
              A22Sub );
        }

        // Form the local portions of (A22 a21) into p21[MC] and q21[MR]:
        //   p21[MC] := tril(A22)[MC,MR] a21[MR]
        //   q21[MR] := tril(A22,-1)'[MR,MC] a21[MC]
        symv::LocalColAccumulate
        ( LOWER, F(1), A22, a21_MC, a21_MR, p21_MC, q21_MR, true, ctrl );

        x01_MR.AlignWith( W20B );
        y01_MR.AlignWith( W20B );
        Zeros( x01_MR, W20B.Width(), 1 );
        Zeros( y01_MR, W20B.Width(), 1 );
        LocalGemv( ADJOINT, F(1), W20B, a21B_MC, F(0), x01_MR );
        LocalGemv( ADJOINT, F(1), A20B, a21B_MC, F(0), y01_MR );

        // Combine the AllReduce column summations of x01[MR], y01[MR]
        {
            const Int x01LocalHeight = x01_MR.LocalHeight();
            vector<F> colSumSendBuf, colSumRecvBuf;
            FastResize( colSumSendBuf, 2*x01LocalHeight );
            FastResize( colSumRecvBuf, 2*x01LocalHeight );
            MemCopy
            ( colSumSendBuf.data(), 
              x01_MR.Buffer(), x01LocalHeight );
            MemCopy
            ( &colSumSendBuf[x01LocalHeight],
              y01_MR.Buffer(), x01LocalHeight );
            mpi::AllReduce
            ( colSumSendBuf.data(), 
              colSumRecvBuf.data(), 2*x01LocalHeight, g.ColComm() );
            MemCopy
            ( x01_MR.Buffer(), 
              colSumRecvBuf.data(), x01LocalHeight );
            MemCopy
            ( y01_MR.Buffer(), 
              &colSumRecvBuf[x01LocalHeight], x01LocalHeight );
        }

        LocalGemv( NORMAL, F(-1), A20B, x01_MR, F(1), p21B_MC );
        LocalGemv( NORMAL, F(-1), W20B, y01_MR, F(1), p21B_MC );

        // Fast transpose the unsummed q21[MR] -> q21[MC], so that
        // it needs to be summed over process rows instead of process 
        // columns. We immediately add it onto p21[MC], which also needs
        // to be summed within process rows.
        if( onDiagonal )
        {
            blas::Axpy
            ( a21.LocalHeight(), F(1),
              q21_MR.LockedBuffer(), 1,
              p21_MC.Buffer(),       1 );
        }
        else
        {
            // Pairwise exchange with the transpose process
            const Int sendSize = A22.LocalWidth();
            const Int recvSize = A22.LocalHeight();
            //vector<F> recvBuf(recvSize);
            vector<F> recvBuf;
            FastResize( recvBuf, recvSize );
            mpi::SendRecv
            ( q21_MR.Buffer(), sendSize, transposeRank, 
              recvBuf.data(),  recvSize, transposeRank, 
              g.VCComm() );

            // Unpack the recv buffer directly onto p21[MC]
            blas::Axpy( recvSize, F(1), recvBuf.data(), 1, p21_MC.Buffer(), 1 );
        }

        if( W22.Width() > 0 )
        {
            // This is not the last iteration of the panel factorization, 
            // Reduce to one p21[MC] to the next process column.
            const Int a21LocalHeight = a21.LocalHeight();

            const Int nextProcessRow = (alpha11.ColAlign()+1) % r;
            const Int nextProcessCol = (alpha11.RowAlign()+1) % r;

            //vector<F> reduceToOneRecvBuf(a21LocalHeight);
            vector<F> reduceToOneRecvBuf;
            FastResize( reduceToOneRecvBuf, a21LocalHeight );
            mpi::Reduce
            ( p21_MC.Buffer(), reduceToOneRecvBuf.data(),
              a21LocalHeight, nextProcessCol, g.RowComm() );
            if( g.Col() == nextProcessCol )
            {
                // Finish computing w21. During its computation, ensure that 
                // every process has a copy of the first element of the w21.
                // We know a priori that the first element of a21 is one.
                const F* a21_MC_Buf = a21_MC.Buffer();
                F myDotProduct = blas::Dot
                    ( a21LocalHeight, reduceToOneRecvBuf.data(), 1,
                                      a21_MC_Buf,                1 );
                F sendBuf[2], recvBuf[2];
                sendBuf[0] = myDotProduct;
                sendBuf[1] = ( g.Row()==nextProcessRow ?
                                  reduceToOneRecvBuf[0] : 0 );
                mpi::AllReduce( sendBuf, recvBuf, 2, g.ColComm() );
                F dotProduct = recvBuf[0];

                // Set up for the next iteration by filling in the values for:
                // - w21LastBuf
                // - w21LastFirstEntry
                F scale = dotProduct*tau / F(2);
                for( Int i=0; i<a21LocalHeight; ++i )
                    w21LastBuf[i] = Conj(tau)*
                      ( reduceToOneRecvBuf[i]-scale*a21_MC_Buf[i] );
                w21LastFirstEntry = Conj(tau)*( recvBuf[1]-scale );
            }
        }
        else
        {
            // This is the last iteration, so we just need to form w21[MC]
            // and w21[MR]. 
            const Int a21LocalHeight = a21.LocalHeight();

            // AllReduce sum p21[MC] over process rows
            vector<F> allReduceRecvBuf;
            FastResize( allReduceRecvBuf, a21LocalHeight );
            mpi::AllReduce
            ( p21_MC.Buffer(), allReduceRecvBuf.data(),
              a21LocalHeight, g.RowComm() );

            // Finish computing w21.
            const F* a21_MC_Buf = a21_MC.Buffer();
            F myDotProduct = blas::Dot
                ( a21LocalHeight, allReduceRecvBuf.data(), 1,
                                  a21_MC_Buf,              1 );
            const F dotProduct = mpi::AllReduce( myDotProduct, g.ColComm() );

            // Grab views into W[MC,* ] and W[MR,* ]
            auto w21_MC = W_MC_STAR( ind2, ind1 );
            auto w21_MR = W_MR_STAR( ind2, ind1 );

            // Store w21[MC]
            F scale = dotProduct*tau / F(2);
            F* w21_MC_Buf = w21_MC.Buffer();
            for( Int i=0; i<a21LocalHeight; ++i )
                w21_MC_Buf[i] = Conj(tau)*
                    ( allReduceRecvBuf[i]-scale*a21_MC_Buf[i] );

            // Fast transpose w21[MC] -> w21[MR]
            if( onDiagonal )
            {
                MemCopy
                ( w21_MR.Buffer(),
                  w21_MC.Buffer(), a21LocalHeight );
            }
            else
            {
                // Pairwise exchange with the transpose process
                mpi::SendRecv
                ( w21_MC.Buffer(), A22.LocalHeight(), transposeRank, 
                  w21_MR.Buffer(), A22.LocalWidth(),  transposeRank, 
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

#endif // ifndef EL_HERMITIANTRIDIAG_LOWER_PANEL_SQUARE_HPP
