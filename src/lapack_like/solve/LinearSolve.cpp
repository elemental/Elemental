/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

namespace lu {

template<typename F>
void Panel( Matrix<F>& APan, Matrix<Int>& p1 );

template<typename F>
void Panel
( DistMatrix<F,  STAR,STAR>& A11,
  DistMatrix<F,  MC,  STAR>& A21,
  DistMatrix<Int,STAR,STAR>& p1 );

} // namespace lu

// Short-circuited form of LU factorization with partial pivoting
template<typename F> 
inline void
RowEchelon( Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("RowEchelon");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
    )

    Matrix<Int> p1Piv;
    Matrix<Int> p1, p1Inv;

    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int minDimA = Min(mA,nA);
    const Int nB = B.Width();
    const Int bsize = Blocksize();
    for( Int k=0; k<minDimA; k+=bsize )
    {
        const Int nb = Min(bsize,minDimA-k);
        const Range<Int> ind1( k, k+nb ),
                         indB( k, mA   ),
                         ind2Vert( k+nb, mA ), ind2Horz( k+nb, nA );
        auto A11 = A( ind1,     ind1     );
        auto A12 = A( ind1,     ind2Horz );
        auto A21 = A( ind2Vert, ind1     );
        auto A22 = A( ind2Vert, ind2Horz ); 
        auto AB2 = A( indB,     ind2Horz );
        auto B1  = B( ind1,     IR(0,nB)    );
        auto B2  = B( ind2Vert, IR(0,nB)    );
        auto BB  = B( indB,     IR(0,nB)    );

        lu::Panel( AB2, p1Piv );
        PivotsToPartialPermutation( p1Piv, p1, p1Inv ); 
        PermuteRows( BB, p1, p1Inv );

        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11, A12 );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11, B1 );

        Gemm( NORMAL, NORMAL, F(-1), A21, A12, F(1), A22 );
        Gemm( NORMAL, NORMAL, F(-1), A21, B1,  F(1), B2 );
    }
}

// Short-circuited form of LU factorization with partial pivoting
template<typename F> 
inline void
RowEchelon( DistMatrix<F>& A, DistMatrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("RowEchelon");
        AssertSameGrids( A, B );
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
    )
    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int minDimA = Min(mA,nA);
    const Int nB = B.Width();
    const Int bsize = Blocksize();
    const Grid& g = A.Grid();

    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g), B1_STAR_VR(g);
    DistMatrix<F,STAR,MR  > A12_STAR_MR(g), B1_STAR_MR(g);
    DistMatrix<F,MC,  STAR> A21_MC_STAR(g);
    DistMatrix<Int,STAR,STAR> p1Piv_STAR_STAR(g);

    DistMatrix<Int,VC,STAR> p1(g), p1Inv(g);

    // In case B's columns are not aligned with A's
    const bool BAligned = ( B.ColShift() == A.ColShift() );
    DistMatrix<F,MC,STAR> A21_MC_STAR_B(g);

    for( Int k=0; k<minDimA; k+=bsize )
    {
        const Int nb = Min(bsize,minDimA-k);
        const Range<Int> ind1( k, k+nb ),
                         indB( k, mA   ),
                         ind2Vert( k+nb, mA ), ind2Horz( k+nb, nA );
        auto A11 = A( ind1,     ind1     );
        auto A12 = A( ind1,     ind2Horz );
        auto A21 = A( ind2Vert, ind1     );
        auto A22 = A( ind2Vert, ind2Horz ); 
        auto AB2 = A( indB,     ind2Horz );
        auto B1  = B( ind1,     IR(0,nB)    );
        auto B2  = B( ind2Vert, IR(0,nB)    );
        auto BB  = B( indB,     IR(0,nB)    );

        A11_STAR_STAR = A11;
        A21_MC_STAR.AlignWith( A22 );
        A21_MC_STAR = A21;

        lu::Panel( A11_STAR_STAR, A21_MC_STAR, p1Piv_STAR_STAR );
        PivotsToPartialPermutation( p1Piv_STAR_STAR, p1, p1Inv );
        PermuteRows( AB2, p1, p1Inv );
        PermuteRows( BB,  p1, p1Inv );

        A12_STAR_VR.AlignWith( A22 );
        A12_STAR_VR = A12;
        B1_STAR_VR.AlignWith( B1 );
        B1_STAR_VR = B1;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, UNIT, F(1), A11_STAR_STAR, A12_STAR_VR );
        LocalTrsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11_STAR_STAR, B1_STAR_VR );

        A12_STAR_MR.AlignWith( A22 );
        A12_STAR_MR = A12_STAR_VR;
        B1_STAR_MR.AlignWith( B1 );
        B1_STAR_MR = B1_STAR_VR;
        LocalGemm( NORMAL, NORMAL, F(-1), A21_MC_STAR, A12_STAR_MR, F(1), A22 );
        if( BAligned )
        {
            LocalGemm
            ( NORMAL, NORMAL, F(-1), A21_MC_STAR, B1_STAR_MR, F(1), B2 );
        }
        else
        {
            A21_MC_STAR_B.AlignWith( B2 );
            A21_MC_STAR_B = A21_MC_STAR;
            LocalGemm
            ( NORMAL, NORMAL, F(-1), A21_MC_STAR_B, B1_STAR_MR, F(1), B2 );
        }

        A11 = A11_STAR_STAR;
        A12 = A12_STAR_MR;
        B1 = B1_STAR_MR;
    }
}

template<typename F> 
void LinearSolve( Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("LinearSolve"))
    // Perform Gaussian elimination
    RowEchelon( A, B );
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
}

template<typename F> 
void LinearSolve
( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& BPre )
{
    DEBUG_ONLY(CallStackEntry cse("LinearSolve"))
    // Perform Gaussian elimination

    // NOTE: Since only the upper triangle of the factorization is formed,
    //       we could usually get away with A only being a Write proxy.
    auto APtr = ReadWriteProxy<F,MC,MR>( &APre ); auto& A = *APtr;
    auto BPtr = ReadWriteProxy<F,MC,MR>( &BPre ); auto& B = *BPtr;

    RowEchelon( A, B );
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
}

template<typename F>
void LinearSolve( const SparseMatrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("LinearSolve"))
    LogicError("Sequential sparse-direct solvers not yet supported");
}

template<typename F>
void LinearSolve( const DistSparseMatrix<F>& A, DistMultiVec<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("LinearSolve"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = B.Width();
    const Int numLocalEntriesA = A.NumLocalEntries();
    if( m != n )
        LogicError("Cannot solve a linear system with a non-square matrix");
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size(comm);

    // Form the Hermitian quasi-semidefinite system
    // 
    //     | 0 A^H | | X | = | 0 |,
    //     | A 0   | | Y |   | B |
    //
    // where the top-left corner is chosen as positive-semidefinite and the
    // bottom-right corner is chosen as negative semi-definite.
    //
    // NOTE: If A is nonsingular, then Y must be zero up to rounding error.

    // Form J = [0, A^H; A, 0]
    // =======================
    DistSparseMatrix<F> J(comm);
    Zeros( J, 2*n, 2*n );
    {
        // Compute metadata
        // ----------------
        std::vector<int> sendCounts(commSize,0);
        for( Int e=0; e<numLocalEntriesA; ++e )
        {
            const Int i = A.Row(e);
            const Int j = A.Col(e);
            // Sending A
            ++sendCounts[ J.RowOwner(i+n) ];
            // Sending A^H
            ++sendCounts[ J.RowOwner(j) ];
        }
        std::vector<int> recvCounts(commSize);
        mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
        std::vector<int> sendOffsets(commSize), recvOffsets(commSize);
        const int totalSend = Scan( sendCounts, sendOffsets );
        const int totalRecv = Scan( recvCounts, recvOffsets );
        // Pack
        // ----
        std::vector<Int> sSendBuf(totalSend), tSendBuf(totalSend);
        std::vector<F> vSendBuf(totalSend);
        auto offsets = sendOffsets;
        for( Int e=0; e<numLocalEntriesA; ++e )
        {
            const Int i = A.Row(e);
            const Int j = A.Col(e);
            const F value = A.Value(e);

            // Sending A
            int owner = J.RowOwner(i+n);
            sSendBuf[offsets[owner]] = i+n;
            tSendBuf[offsets[owner]] = j;
            vSendBuf[offsets[owner]] = value;
            // Sending A^H
            owner = J.RowOwner(j);
            sSendBuf[offsets[owner]] = j;
            tSendBuf[offsets[owner]] = i+n;
            vSendBuf[offsets[owner]] = Conj(value);

            ++offsets[owner];
        }
        // Exchange
        // --------
        std::vector<Int> sRecvBuf(totalRecv), tRecvBuf(totalRecv);
        std::vector<F> vRecvBuf(totalRecv);
        mpi::AllToAll
        ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        mpi::AllToAll
        ( tSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          tRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        mpi::AllToAll
        ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        // Unpack
        // ------
        J.Reserve( totalRecv );
        for( Int e=0; e<totalRecv; ++e )
            J.QueueLocalUpdate
            ( sRecvBuf[e]-J.FirstLocalRow(), tRecvBuf[e], vRecvBuf[e] );
        J.MakeConsistent();
    }

    // Form D = [0; B]
    // ===============
    DistMultiVec<F> D(comm);
    Zeros( D, 2*n, k );
    {
        // Compute metadata
        // ----------------
        std::vector<int> sendCounts(commSize,0);
        for( Int iLoc=0; iLoc<B.LocalHeight(); ++iLoc )
        {
            const Int i = B.GlobalRow(iLoc);
            sendCounts[ D.RowOwner(i+n) ] += k;
        }
        std::vector<int> recvCounts(commSize);
        mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
        std::vector<int> sendOffsets(commSize), recvOffsets(commSize);
        const int totalSend = Scan( sendCounts, sendOffsets );
        const int totalRecv = Scan( recvCounts, recvOffsets );
        // Pack
        // ----
        std::vector<Int> sSendBuf(totalSend), tSendBuf(totalSend);
        std::vector<F> vSendBuf(totalSend);
        auto offsets = sendOffsets;
        for( Int iLoc=0; iLoc<B.LocalHeight(); ++iLoc )
        {
            const Int i = B.GlobalRow(iLoc);
            int owner = D.RowOwner(i+n);
            
            for( Int j=0; j<k; ++j )
            {
                const F value = B.GetLocal(iLoc,j);
                sSendBuf[offsets[owner]] = i+n;
                tSendBuf[offsets[owner]] = j;
                vSendBuf[offsets[owner]] = value;
                ++offsets[owner];
            }
        }
        // Exchange
        // --------
        std::vector<Int> sRecvBuf(totalRecv), tRecvBuf(totalRecv);
        std::vector<F> vRecvBuf(totalRecv);
        mpi::AllToAll
        ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        mpi::AllToAll
        ( tSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          tRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        mpi::AllToAll
        ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        // Unpack
        // ------
        for( Int e=0; e<totalRecv; ++e )
            D.UpdateLocal
            ( sRecvBuf[e]-D.FirstLocalRow(), tRecvBuf[e], vRecvBuf[e] );
    }

    // Compute the dynamically-regularized quasi-semidefinite fact of J
    // ================================================================
    DistMap map, invMap;
    DistSymmInfo info;
    DistSeparatorTree sepTree;
    DistSymmFrontTree<F> JFrontTree;
    DistMultiVec<Real> regCand(comm), reg(comm);
    DistNodalMultiVec<Real> regCandNodal, regNodal;
    const Real minReductionFactor = 2;
    const Int maxRefineIts = 50; 
    bool print = true;
    const Real epsilon = lapack::MachineEpsilon<Real>();
    const Real regMag = Pow(epsilon,Real(0.5));
    const Real pivTol = MaxNorm(J)*epsilon;
    regCand.Resize( 2*n, 1 );
    for( Int iLoc=0; iLoc<regCand.LocalHeight(); ++iLoc )
    {
        const Int i = regCand.GlobalRow(iLoc);
        if( i < n )
            regCand.SetLocal( iLoc, 0, regMag );
        else
            regCand.SetLocal( iLoc, 0, -regMag );
    }
    NestedDissection( J.LockedDistGraph(), map, sepTree, info );
    map.FormInverse( invMap );
    JFrontTree.Initialize( J, map, sepTree, info );
    regCandNodal.Pull( invMap, info, regCand );
    regNodal.Pull( invMap, info, reg );
    RegularizedQSDLDL
    ( info, JFrontTree, pivTol, regCandNodal, regNodal, LDL_1D );
    regNodal.Push( invMap, info, reg );

    // Successively solve each of the k linear systems
    // ===============================================
    // TODO: Extend the iterative refinement to handle multiple right-hand sides
    DistMultiVec<F> u(comm);
    Zeros( u, 2*n, 1 );
    Matrix<F>& DLoc = D.Matrix();
    Matrix<F>& uLoc = u.Matrix();
    for( Int j=0; j<k; ++j )
    {
        auto dLoc = DLoc( IR(0,2*n), IR(j,j+1) ); 
        Copy( dLoc, uLoc );
        reg_qsd_ldl::SolveAfter
        ( J, reg, invMap, info, JFrontTree, u,
          minReductionFactor, maxRefineIts, print );
        Copy( uLoc, dLoc );
    }

    // Extract X from [X; Y]
    // =====================
    {
        // Compute metadata
        // ----------------
        std::vector<int> sendCounts(commSize,0);
        for( Int iLoc=0; iLoc<D.LocalHeight(); ++iLoc )
        {
            const Int i = D.GlobalRow(iLoc);
            if( i < n )
                sendCounts[ B.RowOwner(i) ] += k;
            else
                break;
        }
        std::vector<int> recvCounts(commSize);
        mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
        std::vector<int> sendOffsets(commSize), recvOffsets(commSize);
        const int totalSend = Scan( sendCounts, sendOffsets );
        const int totalRecv = Scan( recvCounts, recvOffsets );
        // Pack
        // ----
        std::vector<Int> sSendBuf(totalSend), tSendBuf(totalSend);
        std::vector<F> vSendBuf(totalSend);
        auto offsets = sendOffsets;
        for( Int iLoc=0; iLoc<D.LocalHeight(); ++iLoc )
        {
            const Int i = D.GlobalRow(iLoc);
            if( i < n )
            {
                int owner = B.RowOwner(i);
                for( Int j=0; j<k; ++j )
                {
                    const F value = D.GetLocal(iLoc,j);
                    sSendBuf[offsets[owner]] = i;
                    tSendBuf[offsets[owner]] = j;
                    vSendBuf[offsets[owner]] = value;
                    ++offsets[owner];
                }
            }
            else
                break;
        }
        // Exchange
        // --------
        std::vector<Int> sRecvBuf(totalRecv), tRecvBuf(totalRecv);
        std::vector<F> vRecvBuf(totalRecv);
        mpi::AllToAll
        ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        mpi::AllToAll
        ( tSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          tRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        mpi::AllToAll
        ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        // Unpack
        // ------
        for( Int e=0; e<totalRecv; ++e )
            B.UpdateLocal
            ( sRecvBuf[e]-B.FirstLocalRow(), tRecvBuf[e], vRecvBuf[e] );
    }
}

#define PROTO(F) \
  template void LinearSolve( Matrix<F>& A, Matrix<F>& B ); \
  template void LinearSolve \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B ); \
  template void LinearSolve( const SparseMatrix<F>& A, Matrix<F>& B ); \
  template void LinearSolve \
  ( const DistSparseMatrix<F>& A, DistMultiVec<F>& B );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
