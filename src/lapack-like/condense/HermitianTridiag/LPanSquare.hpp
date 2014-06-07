/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HERMITIANTRIDIAG_LPANSQUARE_HPP
#define EL_HERMITIANTRIDIAG_LPANSQUARE_HPP

namespace El {
namespace herm_tridiag {

template<typename F>
void LPanSquare
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
        CallStackEntry cse("herm_tridiag::LPanSquare");
        if( A.Grid() != W.Grid() || W.Grid() != t.Grid() )
            LogicError("A, W, and t must be distributed over the same grid");
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
            LogicError("t is not aligned with A's subdiagonal");
    )
    typedef Base<F> Real;
    // Find the process holding our transposed data
    const Grid& g = A.Grid();
    const Int r = g.Height();
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

    // Create a distributed matrix for storing the subdiagonal
    DistMatrix<Real,MD,STAR> e(g);
    e.SetRoot( A.DiagonalRoot(-1) );
    e.AlignCols( A.DiagonalAlign(-1) );
    e.Resize( nW, 1 );

    std::vector<F> w21LastBuffer(A.Height()/r+1);
    DistMatrix<F> w21Last(g);
    DistMatrix<F,MC,STAR> a21_MC_STAR(g), a21B_MC_STAR(g),
                          p21_MC_STAR(g), p21B_MC_STAR(g),
                          a21Last_MC_STAR(g), w21Last_MC_STAR(g);
    DistMatrix<F,MR,STAR> q21_MR_STAR(g), a21_MR_STAR(g),
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
        auto ABR      = ViewRange( A, k,   k,   n,   n   );
        auto tau1     = View( t, k, 0, 1, 1 );
        auto epsilon1 = View( e, k, 0, 1, 1 );

        a21_MC_STAR.AlignWith( A22 );
        a21_MR_STAR.AlignWith( A22 );
        p21_MC_STAR.AlignWith( A22 );
        q21_MR_STAR.AlignWith( A22 );
        a21_MC_STAR.Resize( n-(k+1), 1 );
        a21_MR_STAR.Resize( n-(k+1), 1 );
        p21_MC_STAR.Resize( n-(k+1), 1 );
        q21_MR_STAR.Resize( n-(k+1), 1 );

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
            
            // Take advantage of the square process grid in order to form 
            // a21[MR,* ] from a21[MC,* ]
            if( onDiagonal )
            {
                MemCopy
                ( a21_MR_STAR.Buffer(), a21_MC_STAR.Buffer(), a21LocalHeight );
            }
            else
            {
                // Pairwise exchange
                const Int sendSize = A22.LocalHeight();
                const Int recvSize = A22.LocalWidth();
                mpi::SendRecv
                ( a21_MC_STAR.Buffer(), sendSize, transposeRank,
                  a21_MR_STAR.Buffer(), recvSize, transposeRank, g.VCComm() );
            }
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
            ( a21_MC_STAR.Buffer(), rowBroadcastBuffer.data(), a21LocalHeight );
            // Store a21[MC,* ] into APan[MC,* ]
            const Int APan_MC_STAR_Offset = 
                APan_MC_STAR.LocalHeight()-a21LocalHeight;
            MemCopy
            ( APan_MC_STAR.Buffer(APan_MC_STAR_Offset,A00.Width()), 
              rowBroadcastBuffer.data(),
              APan_MC_STAR.LocalHeight()-APan_MC_STAR_Offset );
            // Store w21Last[MC,* ] into its DistMatrix class
            w21Last_MC_STAR.AlignWith( alpha11 );
            w21Last_MC_STAR.Resize( a21.Height()+1, 1 );
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
                std::vector<F> sendBuffer(sendSize), recvBuffer(recvSize);

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
        if( A22.ColShift() > A22.RowShift() )
        {
            // We are below the diagonal, so we can multiply without an 
            // offset for tril(A22)[MC,MR] and tril(A22,-1)'[MR,MC]
            if( A22.LocalHeight() != 0 )
            {
                MemCopy
                ( p21_MC_STAR.Buffer(),
                  a21_MR_STAR.Buffer(), A22.LocalHeight() );
                blas::Trmv
                ( 'L', 'N', 'N', A22.LocalHeight(), A22.Buffer(), A22.LDim(), 
                  p21_MC_STAR.Buffer(), 1 );
            }
            if( A22.LocalWidth() != 0 )
            {
                // Our local portion of q21[MR,* ] might be one entry longer 
                // than A22.LocalHeight(), so go ahead and set the last entry 
                // to 0.
                F* q21_MR_STAR_Buffer = q21_MR_STAR.Buffer();
                q21_MR_STAR_Buffer[A22.LocalWidth()-1] = 0;
                MemCopy
                ( q21_MR_STAR_Buffer, a21_MC_STAR.Buffer(), A22.LocalHeight() );
                blas::Trmv
                ( 'L', 'C', 'N', A22.LocalHeight(), A22.Buffer(), A22.LDim(),
                  q21_MR_STAR_Buffer, 1 );
            }
        }
        else if( A22.ColShift() < A22.RowShift() )
        {
            // We are above the diagonal, so we need to use an offset of +1
            // for both tril(A22)[MC,MR] and tril(A22,-1)'[MR,MC]
            const F* a21_MC_STAR_Buffer = a21_MC_STAR.Buffer();
            const F* a21_MR_STAR_Buffer = a21_MR_STAR.Buffer();
            const F* A22Buffer = A22.Buffer();
            if( A22.LocalHeight() != 0 )
            {
                // The first entry of p21[MC,* ] will always be zero due to 
                // the forced offset.
                F* p21_MC_STAR_Buffer = p21_MC_STAR.Buffer();
                p21_MC_STAR_Buffer[0] = 0;
                MemCopy
                ( &p21_MC_STAR_Buffer[1],
                  a21_MR_STAR_Buffer, A22.LocalHeight()-1 );
                blas::Trmv
                ( 'L', 'N', 'N', A22.LocalHeight()-1, &A22Buffer[1], A22.LDim(),
                  &p21_MC_STAR_Buffer[1], 1 );
            }
            if( A22.LocalWidth() != 0 )
            {
                // The last entry of q21[MR,* ] will be zero if the local
                // height and width are equal.
                F* q21_MR_STAR_Buffer = q21_MR_STAR.Buffer();
                q21_MR_STAR_Buffer[A22.LocalWidth()-1] = 0;
                MemCopy
                ( q21_MR_STAR_Buffer,
                  &a21_MC_STAR_Buffer[1], A22.LocalHeight()-1 );
                blas::Trmv
                ( 'L', 'C', 'N', A22.LocalHeight()-1, &A22Buffer[1], A22.LDim(),
                  q21_MR_STAR_Buffer, 1 );
            }
        }
        else
        {
            // We are on the diagonal, so we only need an offset of +1 for
            // tril(A22,-1)'[MR,MC]
            if( A22.LocalHeight() != 0 )
            {
                MemCopy
                ( p21_MC_STAR.Buffer(),
                  a21_MR_STAR.Buffer(), A22.LocalHeight() );
                blas::Trmv
                ( 'L', 'N', 'N', A22.LocalHeight(), A22.Buffer(), A22.LDim(), 
                  p21_MC_STAR.Buffer(), 1 );

                // The last entry of q21[MR,* ] will be zero if the local 
                // height and width are equal.
                const F* a21_MC_STAR_Buffer = a21_MC_STAR.Buffer();
                const F* A22Buffer = A22.Buffer();
                F* q21_MR_STAR_Buffer = q21_MR_STAR.Buffer();
                q21_MR_STAR_Buffer[A22.LocalWidth()-1] = 0;
                MemCopy
                ( q21_MR_STAR_Buffer,
                  &a21_MC_STAR_Buffer[1], A22.LocalHeight()-1 );
                blas::Trmv
                ( 'L', 'C', 'N', A22.LocalHeight()-1, &A22Buffer[1], A22.LDim(),
                  q21_MR_STAR_Buffer, 1 );
            }
        }

        x01_MR_STAR.AlignWith( W20B );
        y01_MR_STAR.AlignWith( W20B );
        Zeros( x01_MR_STAR, W20B.Width(), 1 );
        Zeros( y01_MR_STAR, W20B.Width(), 1 );
        LocalGemv( ADJOINT, F(1), W20B, a21B_MC_STAR, F(0), x01_MR_STAR );
        LocalGemv( ADJOINT, F(1), A20B, a21B_MC_STAR, F(0), y01_MR_STAR );

        // Combine the AllReduce column summations of x01[MR,* ] and 
        // y01[MR,* ]
        {
            const Int x01LocalHeight = x01_MR_STAR.LocalHeight();
            std::vector<F> colSumSendBuffer(2*x01LocalHeight),
                           colSumRecvBuffer(2*x01LocalHeight);
            MemCopy
            ( colSumSendBuffer.data(), 
              x01_MR_STAR.Buffer(), x01LocalHeight );
            MemCopy
            ( &colSumSendBuffer[x01LocalHeight],
              y01_MR_STAR.Buffer(), x01LocalHeight );
            mpi::AllReduce
            ( colSumSendBuffer.data(), colSumRecvBuffer.data(), 
              2*x01LocalHeight, g.ColComm() );
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
            const Int a21LocalHeight = a21.LocalHeight();
            F* p21_MC_STAR_Buffer = p21_MC_STAR.Buffer();
            const F* q21_MR_STAR_Buffer = q21_MR_STAR.Buffer();
            for( Int i=0; i<a21LocalHeight; ++i )
                p21_MC_STAR_Buffer[i] += q21_MR_STAR_Buffer[i];
        }
        else
        {
            // Pairwise exchange with the transpose process
            const Int sendSize = A22.LocalWidth();
            const Int recvSize = A22.LocalHeight();
            std::vector<F> recvBuffer(recvSize);
            mpi::SendRecv
            ( q21_MR_STAR.Buffer(), sendSize, transposeRank, 
              recvBuffer.data(),    recvSize, transposeRank, g.VCComm() );

            // Unpack the recv buffer directly onto p21[MC,* ]
            F* p21_MC_STAR_Buffer = p21_MC_STAR.Buffer();
            for( Int i=0; i<recvSize; ++i )
                p21_MC_STAR_Buffer[i] += recvBuffer[i];
        }

        if( W22.Width() > 0 )
        {
            // This is not the last iteration of the panel factorization, 
            // Reduce to one p21[MC,* ] to the next process column.
            const Int a21LocalHeight = a21.LocalHeight();

            const Int nextProcessRow = (alpha11.ColAlign()+1) % r;
            const Int nextProcessCol = (alpha11.RowAlign()+1) % r;

            std::vector<F> reduceToOneRecvBuffer(a21LocalHeight);
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
            std::vector<F> allReduceRecvBuffer(a21LocalHeight);
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
            auto w21_MC_STAR = ViewRange( W_MC_STAR, k+1, k, n, k+1 );
            auto w21_MR_STAR = ViewRange( W_MR_STAR, k+1, k, n, k+1 );

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
                const Int sendSize = A22.LocalHeight();
                const Int recvSize = A22.LocalWidth();
                mpi::SendRecv
                ( w21_MC_STAR.Buffer(), sendSize, transposeRank, 
                  w21_MR_STAR.Buffer(), recvSize, transposeRank, g.VCComm() );
            }
        }
    }

    // View the portion of A that e is the subdiagonal of, then place e into it
    auto expandedATL = View( A, 0, 0, nW+1, nW+1 );
    expandedATL.SetRealPartOfDiagonal( e, -1 );
}

} // namespace herm_tridiag
} // namespace El

#endif // ifndef EL_HERMITIANTRIDIAG_LPANSQUARE_HPP
