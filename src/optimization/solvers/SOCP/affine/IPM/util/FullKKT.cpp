/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "../util.hpp"

namespace El {
namespace socp {
namespace affine {

// Given the Lagrangian
//
//   L(x,s;y,z) = c^T x + y^T (A x - b) + z^T (G x + s - h) + 
//                tau Phi(inv(W)^T s),
//
// where Phi is the standard barrier function of the product cone, tau is a 
// scaling of the barrier, and W is some automorphism of the product cone
// (in our case, it is arises from the Nesterov-Todd scaling), critical points
// must satisfy the first-order optimality conditions:
//
//   D_x L = c + A^T y + G^T z = 0,
//
//   D_y L = A x - b = 0,
//
//   D_z L = G x + s - h = 0, and
//
//   D_s L = z - tau inv(W) inv(inv(W)^T s) = 0,
//
// where the last condition can be simplified to the form
//
//   (inv(W)^T s) o (W z) = tau e,
//
// where "o" is the symmetric product of the Jordan algebra and "e" is the
// corresponding identity element. A Newton step from the point (x,y,z,s)
// therefore corresponds to solving the system
//
//   | 0 A^T     G   | | dx |   | -r_c            |
//   | A 0       0   | | dy |   | -r_b            |,
//   | G 0    -W^T W | | dz | = | -r_h + W^T r_mu |
//
// where W = W^T is a block-diagonal matrix, with each diagonal block 
// equalling -Q_{sqrt(w_i)}, where Q_a is the quadratic representation of 
// a member of the Jordan algebra associated with the SOC:
//   
//   Q_a = 2 a a^T - det(a) R,
//
// where R is the reflection operator diag(1,-I). In our case, we choose the 
// point  w to be the (unique) Nesterov-Todd scaling point satisfying
//
//   Q_w z = s,
//
// or, equivalently,
//
//   Hess(Phi(w)) s = z,
//
// where Phi is the barrier function for the product of second-order cones,
//
//   Phi(w) = -(1/2) sum_i ln det(w_i).
//
// Furthermore, the residuals are defined by:
//
//   r_c  = A^T y + G^T z + c,
//   r_b  = A x - b,
//   r_h  = G x + s - h,
//   r_mu = l - tau inv(l),
//
// where 
//
//   l = W z = inv(W) s
//
// is the Nesterov-Todd scaling of the pair (s,z).
//
// However, a large amount of fill-in is incurred if member cones are large,
// as the the linear operator implied by -W^2 is block-diagonal, with
// each diagonal block corresponding to a member cone. Thankfully, this matrix
// is known [citations!] to be a symmetric rank-two correction to a 
// definite diagonal matrix, and so a sparse embedding can be applied.
// Such a technique will be applied for sufficiently-large subcones in future
// versions of this solver (in the manner described in the ECOS paper of 
// Domahidi et al.).
//

template<typename Real>
void KKT
( const Matrix<Real>& A, 
  const Matrix<Real>& G,
  const Matrix<Real>& s, 
  const Matrix<Real>& z,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  const Matrix<Int>& labels,
        Matrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CSE cse("socp::affine::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    Matrix<Real> w, wDets;
    SOCNesterovTodd( s, z, w, orders, firstInds );
    SOCDets( w, wDets, orders, firstInds );

    Zeros( J, n+m+k, n+m+k );
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd); auto Jxz = J(xInd,zInd); 
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd); auto Jyz = J(yInd,zInd); 
    auto Jzx = J(zInd,xInd); auto Jzy = J(zInd,yInd); auto Jzz = J(zInd,zInd); 

    // Jyx := A
    // ========
    Jyx = A;
    if( !onlyLower )
        Transpose( A, Jxy ); 

    // Jzx := G
    // ========
    Jzx = G;
    if( !onlyLower )
        Transpose( G, Jxz );

    // Jzz := -W^2
    // ===========
    for( Int i=0; i<k; )
    {
        const Int order = orders.Get(i,0);
        const Int firstInd = firstInds.Get(i,0);
        if( i != firstInd )
            LogicError("Inconsistency between firstInds and orders");
        
        auto Jzzi = Jzz(IR(i,i+order),IR(i,i+order));
        auto wi = w(IR(i,i+order),ALL);
        const Real wiDet = wDets.Get(i,0);
        // Jzzi := det(w_i) R - 2 w w^T
        ShiftDiagonal( Jzzi, -wiDet );
        Jzzi.Update( 0, 0, 2*wiDet );
        Syr( LOWER, Real(-2), wi, Jzzi );
        if( !onlyLower )
            MakeSymmetric( LOWER, Jzzi );

        i += order;
    }
}

// TODO: Add tunable cutoffs for both the sparse embedding and parallel
//       scheduling of the Jordan operations
template<typename Real>
void KKT
( const AbstractDistMatrix<Real>& A,    
  const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& s, 
  const AbstractDistMatrix<Real>& z,
  const AbstractDistMatrix<Int>& ordersPre,
  const AbstractDistMatrix<Int>& firstIndsPre,
  const AbstractDistMatrix<Int>& labels,
        AbstractDistMatrix<Real>& JPre, bool onlyLower )
{
    DEBUG_ONLY(CSE cse("socp::affine::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();
    const Grid& g = A.Grid();

    const Int cutoffPar = 1000;

    ProxyCtrl proxCtrl;
    proxCtrl.colConstrain = true;
    proxCtrl.colAlign = 0;

    auto ordersPtr    = ReadProxy<Int,VC,STAR>(&ordersPre,proxCtrl);
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,proxCtrl);
    auto& orders    = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    auto JPtr = WriteProxy<Real,MC,MR>(&JPre); 
    auto& J = *JPtr;

    DistMatrix<Real,VC,STAR> w(g), wDets(g);
    SOCNesterovTodd( s, z, w, orders, firstInds, cutoffPar );
    SOCDets( w, wDets, orders, firstInds, cutoffPar );
    SOCBroadcast( wDets, orders, firstInds, cutoffPar );

    Zeros( J, n+m+k, n+m+k );
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd); auto Jxz = J(xInd,zInd);
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd); auto Jyz = J(yInd,zInd);
    auto Jzx = J(zInd,xInd); auto Jzy = J(zInd,yInd); auto Jzz = J(zInd,zInd);

    // Jyx := A
    // ========
    Jyx = A;
    if( !onlyLower )
        Transpose( A, Jxy );

    // Jzx := G
    // ========
    Jzx = G;
    if( !onlyLower )
        Transpose( G, Jxz );

    // Jzz := -W^2
    // ===========
    // TODO: Create SOCExplicitQuadratic?
    DistMatrix<Int,MC,STAR> orders_MC_STAR(g), firstInds_MC_STAR(g);
    DistMatrix<Real,MC,STAR> w_MC_STAR(g);
    DistMatrix<Real,MR,STAR> w_MR_STAR(g);
    DistMatrix<Real,MC,STAR> wDets_MC_STAR(g);
    orders_MC_STAR.AlignWith( Jzz );
    firstInds_MC_STAR.AlignWith( Jzz );
    w_MC_STAR.AlignWith( Jzz );
    w_MR_STAR.AlignWith( Jzz );
    wDets_MC_STAR.AlignWith( Jzz );
    orders_MC_STAR = orders;
    firstInds_MC_STAR = firstInds;
    w_MC_STAR = w;
    w_MR_STAR = w;
    wDets_MC_STAR = wDets;
    const Int localWidthJzz = Jzz.LocalWidth();
    const Int localHeightJzz = Jzz.LocalHeight();
    // For now, perform this in quadratic time, which does not effect the
    // asymptotic complexity of this routine, as simply zeroing the original
    // matrix has the same asymptotic cost. A more efficient parallel 
    // implementation will require more thought (and likely substantially more
    // code).
    for( Int iLoc=0; iLoc<localHeightJzz; ++iLoc )
    {
        const Int i = Jzz.GlobalRow(iLoc);
        const Int order = orders_MC_STAR.GetLocal(iLoc,0);
        const Int firstInd = firstInds_MC_STAR.GetLocal(iLoc,0);
        const Real omega_i = w_MC_STAR.GetLocal(iLoc,0);
        const Real wDet = wDets_MC_STAR.GetLocal(iLoc,0);
        // TODO: Restrict the range of this loop (perhaps via Length)
        for( Int jLoc=0; jLoc<localWidthJzz; ++jLoc ) 
        {
            const Int j = Jzz.GlobalCol(jLoc);
            const Real omega_j = w_MR_STAR.GetLocal(jLoc,0);

            // Handle the diagonal update, det(w) R
            if( i == firstInd && j == firstInd ) 
                Jzz.UpdateLocal( iLoc, jLoc, wDet );
            else if( i == j )
                Jzz.UpdateLocal( iLoc, jLoc, -wDet );

            // Handle the rank-one update, -2 w w^T
            if( (!onlyLower || i >= j) && 
                j >= firstInd && j < firstInd+order )
                Jzz.UpdateLocal( iLoc, jLoc, -2*omega_i*omega_j ); 
        }
    }
}

template<typename Real>
void KKT
( const SparseMatrix<Real>& A, 
  const SparseMatrix<Real>& G,
  const Matrix<Real>& s, 
  const Matrix<Real>& z,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  const Matrix<Int>& labels,
        SparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CSE cse("socp::affine::KKT"))

    Matrix<Real> w, wDets;
    SOCNesterovTodd( s, z, w, orders, firstInds );
    SOCDets( w, wDets, orders, firstInds );

    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    Zeros( J, n+m+k, n+m+k );
    if( onlyLower )
    {
        // Count the number of entries to queue in the lower triangle
        // ----------------------------------------------------------
        Int numEntries = A.NumEntries() + G.NumEntries();
        for( Int i=0; i<k; )
        {
            if( i != firstInds.Get(i,0) )
                LogicError("Inconsistency between firstInds and orders");

            const Int order = orders.Get(i,0);
            numEntries += (order*(order+1))/2;
            
            i += order;
        }

        // Queue the nonzeros
        // ------------------
        J.Reserve( numEntries );
        for( Int e=0; e<A.NumEntries(); ++e )
            J.QueueUpdate( A.Row(e)+n, A.Col(e), A.Value(e) );
        for( Int e=0; e<G.NumEntries(); ++e )
            J.QueueUpdate( G.Row(e)+n+m, G.Col(e), G.Value(e) );
        for( Int i=0; i<k; ++i )
        {
            const Int firstInd = firstInds.Get(i,0);
            const Real omega_i = w.Get(i,0);
            const Real wDet = wDets.Get(firstInd,0);

            // diag(det(w) R - 2 w w^T)
            if( i == firstInd )
                J.QueueUpdate( i, i, wDet-2*omega_i*omega_i );
            else
                J.QueueUpdate( i, i, -wDet-2*omega_i*omega_i );

            // offdiag(-2 w w^T)
            for( Int j=firstInd; j<i; ++j )
                J.QueueUpdate( i, j, -2*omega_i*w.Get(j,0) );
        }
    }
    else
    {
        // Count the number of entries to queue in the lower triangle
        // ----------------------------------------------------------
        Int numEntries = 2*A.NumEntries() + 2*G.NumEntries();
        for( Int i=0; i<k; )
        {
            const Int order = orders.Get(i,0);
            const Int firstInd = firstInds.Get(i,0);
            if( i != firstInd )
                LogicError("Inconsistency between firstInds and orders");
            numEntries += order*order;
            i += order;
        }

        // Queue the nonzeros
        // ------------------
        J.Reserve( numEntries );
        for( Int e=0; e<A.NumEntries(); ++e )
        {
            J.QueueUpdate( A.Row(e)+n, A.Col(e),   A.Value(e) );
            J.QueueUpdate( A.Col(e),   A.Row(e)+n, A.Value(e) );
        }
        for( Int e=0; e<G.NumEntries(); ++e )
        {
            J.QueueUpdate( G.Row(e)+n+m, G.Col(e),     G.Value(e) );
            J.QueueUpdate( G.Col(e),     G.Row(e)+n+m, G.Value(e) );
        }
        for( Int i=0; i<k; ++i )
        {
            const Int order = orders.Get(i,0);
            const Int firstInd = firstInds.Get(i,0);
            const Real omega_i = w.Get(i,0);
            const Real wDet = wDets.Get(firstInd,0);

            // diag(det(w) R - 2 w w^T)
            if( i == firstInd )
                J.QueueUpdate( i+n+m, i+n+m, wDet-2*omega_i*omega_i );
            else
                J.QueueUpdate( i+n+m, i+n+m, -wDet-2*omega_i*omega_i );

            // offdiag(-2 w w^T)
            for( Int j=firstInd; j<firstInd+order; ++j )
                if( i != j )
                    J.QueueUpdate( i+n+m, j+n+m, -2*omega_i*w.Get(j,0) );
        }

    }
    J.ProcessQueues();
}

template<typename Real>
void KKT
( const DistSparseMatrix<Real>& A, 
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& s,     
  const DistMultiVec<Real>& z,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  const DistMultiVec<Int>& labels,
        DistSparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CSE cse("socp::affine::KKT"))
    mpi::Comm comm = s.Comm();
    const int commSize = mpi::Size(comm);
    const int commRank = mpi::Rank(comm);

    // TODO: Make this a tunable parameter
    const Int cutoffPar = 1000;

    DistMultiVec<Real> w(comm), wDets(comm);
    SOCNesterovTodd( s, z, w, orders, firstInds, cutoffPar );
    SOCDets( w, wDets, orders, firstInds, cutoffPar );
    SOCBroadcast( wDets, orders, firstInds, cutoffPar );

    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    // Gather all of each member cone that we own a piece of
    // -----------------------------------------------------
    // NOTE: We exploit the fact that each process owns a contiguous chunk
    //       of rows so that at most two member cones need to be sent
    //       (and at most two need to be received).
    vector<int> recvOffs;
    vector<Real> recvBuf;
    const Int wLocalHeight = w.LocalHeight();
    {
        // Count the total number of entries of w to send to each process
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        vector<int> sendCounts(commSize,0);
        for( Int iLoc=0; iLoc<wLocalHeight; )
        {
            const Int i = w.GlobalRow(iLoc);
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int numLocalIndsLeft = wLocalHeight-iLoc;
            const Int numLocalIndsCone = 
              Min(numLocalIndsLeft,order-(i-firstInd));

            const int firstOwner = w.RowOwner(firstInd);
            const int lastOwner = w.RowOwner(firstInd+order-1);
            if( firstOwner != commRank )
                sendCounts[firstOwner] += numLocalIndsCone;
            if( lastOwner != commRank )
                sendCounts[lastOwner] += numLocalIndsCone;
            
            iLoc += numLocalIndsCone; 
        }

        // Pack the entries of w to send to each process
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        vector<int> sendOffs;
        const int totalSend = Scan( sendCounts, sendOffs ); 
        vector<Real> sendBuf(totalSend); 
        auto offs = sendOffs;
        for( Int iLoc=0; iLoc<wLocalHeight; ) 
        {
            const Int i = w.GlobalRow(iLoc);
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int numLocalIndsLeft = wLocalHeight-iLoc;
            const Int numLocalIndsCone = 
              Min(numLocalIndsLeft,order-(i-firstInd));

            const int firstOwner = w.RowOwner(firstInd);
            const int lastOwner = w.RowOwner(firstInd+order-1);
            if( firstOwner != commRank )
                for( Int e=0; e<numLocalIndsCone; ++e )
                    sendBuf[offs[firstOwner]++] = w.GetLocal(iLoc+e,0);
            if( lastOwner != commRank )
                for( Int e=0; e<numLocalIndsCone; ++e )
                    sendBuf[offs[lastOwner]++] = w.GetLocal(iLoc+e,0);
            
            iLoc += numLocalIndsCone; 
        }

        // Receive the entries from each process
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        vector<int> recvCounts(commSize);
        mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
        const int totalRecv = Scan( recvCounts, recvOffs );
        recvBuf.resize( totalRecv );
        mpi::AllToAll
        ( sendBuf.data(), sendCounts.data(), sendOffs.data(),
          recvBuf.data(), recvCounts.data(), recvOffs.data(), comm );
    }

    J.SetComm( comm );
    Zeros( J, n+m+k, n+m+k );
    if( onlyLower )
    {
        // Count the number of entries to queue in the lower triangle
        // ----------------------------------------------------------
        Int numRemoteEntries = A.NumLocalEntries() + G.NumLocalEntries();
        for( Int iLoc=0; iLoc<wLocalHeight; ++iLoc )
        {
            const Int i = w.GlobalRow(iLoc);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            numRemoteEntries += (i-firstInd) + 1;
        }

        // Queue the nonzeros
        // ------------------
        J.Reserve( numRemoteEntries, numRemoteEntries );
        for( Int e=0; e<A.NumLocalEntries(); ++e )
            J.QueueUpdate( A.Row(e)+n, A.Col(e), A.Value(e), false );
        for( Int e=0; e<G.NumLocalEntries(); ++e )
            J.QueueUpdate( G.Row(e)+n+m, G.Col(e), G.Value(e), false );
        Int lastFirstInd = -1;
        vector<Real> wBuf;
        auto offs = recvOffs;
        for( Int iLoc=0; iLoc<wLocalHeight; ++iLoc )
        {
            const Int i = w.GlobalRow(iLoc);
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Real omega_i = w.GetLocal(iLoc,0);
            const Real wDet = wDets.GetLocal(iLoc,0);

            if( firstInd != lastFirstInd )
            {
                lastFirstInd = firstInd;
                wBuf.resize( order );
                for( Int j=firstInd; j<firstInd+order; ++j )
                {
                    const int owner = w.RowOwner(j);
                    if( owner == commRank )
                        wBuf[j-firstInd] = w.GetLocal(w.LocalRow(j),0);
                    else
                        wBuf[j-firstInd] = recvBuf[offs[owner]++]; 
                }
            }

            // diag(det(w) R - 2 w w^T)
            if( i == firstInd )
                J.QueueUpdate( i+n+m, i+n+m, wDet-2*omega_i*omega_i, false );
            else
                J.QueueUpdate( i+n+m, i+n+m, -wDet-2*omega_i*omega_i, false );

            // offdiag(-2 w w^T)
            for( Int j=firstInd; j<i; ++j )
                J.QueueUpdate
                ( i+n+m, j+n+m, -2*omega_i*wBuf[j-firstInd], false );
        }
    }
    else
    {
        // Count the number of entries to queue
        // ------------------------------------
        Int numRemoteEntries = 2*A.NumLocalEntries() + 2*G.NumLocalEntries();
        for( Int iLoc=0; iLoc<wLocalHeight; ++iLoc )
        {
            const Int order = orders.GetLocal(iLoc,0);
            numRemoteEntries += order;
        }

        // Queue the nonzeros
        // ------------------
        J.Reserve( numRemoteEntries, numRemoteEntries );
        for( Int e=0; e<A.NumLocalEntries(); ++e )
        {
            J.QueueUpdate( A.Row(e)+n, A.Col(e),   A.Value(e), false );
            J.QueueUpdate( A.Col(e),   A.Row(e)+n, A.Value(e), false );
        }
        for( Int e=0; e<G.NumLocalEntries(); ++e )
        {
            J.QueueUpdate( G.Row(e)+n+m, G.Col(e),     G.Value(e), false );
            J.QueueUpdate( G.Col(e),     G.Row(e)+n+m, G.Value(e), false );
        }
        Int lastFirstInd = -1;
        vector<Real> wBuf;
        auto offs = recvOffs;
        for( Int iLoc=0; iLoc<wLocalHeight; ++iLoc )
        {
            const Int i = w.GlobalRow(iLoc);
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Real omega_i = w.GetLocal(iLoc,0);
            const Real wDet = wDets.GetLocal(iLoc,0);

            if( firstInd != lastFirstInd )
            {
                lastFirstInd = firstInd;
                wBuf.resize( order );
                for( Int j=firstInd; j<firstInd+order; ++j )
                {
                    const int owner = w.RowOwner(j);
                    if( owner == commRank )
                        wBuf[j-firstInd] = w.GetLocal(w.LocalRow(j),0);
                    else
                        wBuf[j-firstInd] = recvBuf[offs[owner]++]; 
                }
            }

            // diag(det(w) R - 2 w w^T)
            if( i == firstInd )
                J.QueueUpdate( i+n+m, i+n+m, wDet-2*omega_i*omega_i, false );
            else
                J.QueueUpdate( i+n+m, i+n+m, -wDet-2*omega_i*omega_i, false );

            // offdiag(-2 w w^T)
            for( Int j=firstInd; j<firstInd+order; ++j )
                if( j != i )
                    J.QueueUpdate
                    ( i+n+m, j+n+m, -2*omega_i*wBuf[j-firstInd], false );
        }
    }
    J.ProcessQueues();
}

#define PROTO(Real) \
  template void KKT \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& G, \
    const Matrix<Real>& s, \
    const Matrix<Real>& z, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    const Matrix<Int>& labels, \
          Matrix<Real>& J, bool onlyLower ); \
  template void KKT \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& G, \
    const AbstractDistMatrix<Real>& s, \
    const AbstractDistMatrix<Real>& z, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, \
    const AbstractDistMatrix<Int>& labels, \
          AbstractDistMatrix<Real>& J, bool onlyLower ); \
  template void KKT \
  ( const SparseMatrix<Real>& A, \
    const SparseMatrix<Real>& G, \
    const Matrix<Real>& s, \
    const Matrix<Real>& z, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    const Matrix<Int>& labels, \
          SparseMatrix<Real>& J, bool onlyLower ); \
  template void KKT \
  ( const DistSparseMatrix<Real>& A, \
    const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& s, \
    const DistMultiVec<Real>& z, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    const DistMultiVec<Int>& labels, \
          DistSparseMatrix<Real>& J, bool onlyLower );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace affine
} // namespace socp
} // namespace El
