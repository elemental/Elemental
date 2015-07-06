/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename Real>
void LongOnlyPortfolio
( const DistSparseMatrix<Real>& Sigma,
  const DistMultiVec<Real>& c,
        Real gamma,
        DistMultiVec<Real>& x )
{
    DEBUG_ONLY(CSE cse("LongOnlyPortfolio"))
    const Int n = c.Height();
    mpi::Comm comm = c.Comm();

    qp::direct::Ctrl<Real> ctrl;
    ctrl.mehrotraCtrl.print = true;
    ctrl.mehrotraCtrl.qsdCtrl.progress = true;

    // Rather than making a copy of Sigma to form gamma*Sigma, scale c
    // ===============================================================
    auto cScaled = c;
    cScaled *= -1/gamma;

    // Enforce 1^T x = 1
    // ================= 
    DistSparseMatrix<Real> A(comm);
    Ones( A, 1, n );
    DistMultiVec<Real> b(comm);
    Ones( b, 1, 1 );

    DistMultiVec<Real> y(comm), z(comm);
    QP( Sigma, A, b, cScaled, x, y, z, ctrl );
}

template<typename Real>
void LongOnlyPortfolio
( const DistMultiVec<Real>& d,
  const DistSparseMatrix<Real>& F,
  const DistMultiVec<Real>& c,
        Real gamma,
        DistMultiVec<Real>& x )
{
    DEBUG_ONLY(CSE cse("LongOnlyPortfolio"))
    const Int n = c.Height();
    mpi::Comm comm = c.Comm();
    const int commRank = mpi::Rank(comm);
 
    // TODO: Expose this as a control parameter
    const bool useSOCP = true;

    if( useSOCP )
    {
        const Int r = F.Width();

        socp::affine::Ctrl<Real> ctrl;
        ctrl.mehrotraCtrl.time = true;
        ctrl.mehrotraCtrl.print = true;
        ctrl.mehrotraCtrl.qsdCtrl.progress = true;

        // Form cHat = [-c^T, gamma, gamma, 0, 0]
        // ======================================
        DistMultiVec<Real> cHat(comm);
        {
            Zeros( cHat, n+4, 1 );
            const Int localHeight = c.LocalHeight();
            if( commRank == 0 )
            {
                cHat.Reserve( localHeight+2 ); 
                cHat.QueueUpdate( n,   0, gamma );
                cHat.QueueUpdate( n+1, 0, gamma );
            }
            else
                cHat.Reserve( localHeight );

            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const Int i = c.GlobalRow(iLoc);
                cHat.QueueUpdate( i, 0, -c.GetLocal(iLoc,0) );
            }
            cHat.ProcessQueues();
        }

        // Form A = [1^T, 0, 0, 0, 0]
        // ==========================
        DistSparseMatrix<Real> A(comm);
        {
            Zeros( A, 1, n+4 );
            const Int localHeight = A.LocalHeight(); 
            if( localHeight > 0 )
            {
                A.Reserve( n );
                for( Int j=0; j<n; ++j )
                    A.QueueLocalUpdate( 0, j, Real(1) );
                A.ProcessLocalQueues();
            }
        }

        // Form b = 1
        // ==========
        DistMultiVec<Real> b(comm);
        Ones( b, 1, 1 );

        // Form G
        // ======
        // G = |   -I      0  0  0  0 |
        //     |    0      0  0 -1  0 |
        //     | -sqrt(D)  0  0  0  0 |
        //     |    0      0  0  0 -1 |
        //     |   -F^T    0  0  0  0 |
        //     |    0     -1  0  0  0 |
        //     |    0      1  0  0  0 |
        //     |    0      0  0 -2  0 |
        //     |    0      0 -1  0  0 |
        //     |    0      0  1  0  0 |
        //     |    0      0  0  0 -2 |
        DistSparseMatrix<Real> G(comm);
        {
            Zeros( G, 2*n+r+8, n+4 );

            // Count the number of purely local updates 
            // ----------------------------------------
            Int numLocalUpdates = 0;
            const Int GLocalHeight = G.LocalHeight();
            for( Int iLoc=0; iLoc<GLocalHeight; ++iLoc ) 
            {
                const Int i = G.GlobalRow(iLoc);
                if( i <= n )
                    ++numLocalUpdates;
                else if( i == 2*n+1 )
                    ++numLocalUpdates;
                else if( i >= 2*n+r+2 )
                    ++numLocalUpdates;
            }

            // Count the number of (potentially) remote updates
            // ------------------------------------------------
            const Int numRemoteUpdates = d.LocalHeight() + F.NumLocalEntries();

            // Queue the updates
            // -----------------
            G.Reserve( numLocalUpdates+numRemoteUpdates, numRemoteUpdates );
            // Local first 
            // ^^^^^^^^^^^
            for( Int iLoc=0; iLoc<GLocalHeight; ++iLoc )
            {
                const Int i = G.GlobalRow(iLoc);
                if( i < n )
                    G.QueueLocalUpdate( iLoc, i, Real(-1) ); 
                else if( i == n )
                    G.QueueLocalUpdate( iLoc, n+2, Real(-1) );
                else if( i == 2*n+1 )
                    G.QueueLocalUpdate( iLoc, n+3, Real(-1) );
                else if( i == 2*n+r+2 )
                    G.QueueLocalUpdate( iLoc, n,   Real(-1) );
                else if( i == 2*n+r+3 )
                    G.QueueLocalUpdate( iLoc, n,   Real(+1) );
                else if( i == 2*n+r+4 )
                    G.QueueLocalUpdate( iLoc, n+2, Real(-2) );
                else if( i == 2*n+r+5 )
                    G.QueueLocalUpdate( iLoc, n+1, Real(-1) );
                else if( i == 2*n+r+6 )
                    G.QueueLocalUpdate( iLoc, n+1, Real(+1) );
                else if( i == 2*n+r+7 )
                    G.QueueLocalUpdate( iLoc, n+3, Real(-2) );
            }
            // Remote second
            // ^^^^^^^^^^^^^
            const Int dLocalHeight = d.LocalHeight();
            for( Int iLoc=0; iLoc<dLocalHeight; ++iLoc )
            {
                const Int i = d.GlobalRow(iLoc);
                G.QueueUpdate( i+n+1, i, -Sqrt(d.GetLocal(iLoc,0)), false );
            }
            const Int numEntriesF = F.NumLocalEntries(); 
            for( Int e=0; e<numEntriesF; ++e )
            {
                const Int i = F.Row(e);
                const Int j = F.Col(e);
                const Real value = F.Value(e);
                G.QueueUpdate( 2*n+2+j, i, -value, false );
            }
            G.ProcessQueues();
        }

        // Form h = [0; 0; 0; 0; 0; 1; 1; 0; 1; 1; 0]
        // ==========================================
        DistMultiVec<Real> h(comm);
        {
            Zeros( h, 2*n+r+8, 1 );
            const Int localHeight = h.LocalHeight();
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const Int i = h.GlobalRow(iLoc);
                if( i == 2*n+r+2 || i == 2*n+r+3 ||
                    i == 2*n+r+5 || i == 2*n+r+6 )
                    h.SetLocal( iLoc, 0, Real(1) );
            }
        }

        // Form orders and firstInds
        // =========================
        DistMultiVec<Int> orders(comm), firstInds(comm); 
        {
            Zeros( orders, 2*n+r+8, 1 );
            Zeros( firstInds, 2*n+r+8, 1 );
            const Int localHeight = orders.LocalHeight();
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const Int i = orders.GlobalRow(iLoc);
                if( i < n )
                {
                    orders.SetLocal( iLoc, 0, 1 );
                    firstInds.SetLocal( iLoc, 0, i );
                }
                else if( i < 2*n+1 )
                {
                    orders.SetLocal( iLoc, 0, n+1 );
                    firstInds.SetLocal( iLoc, 0, n );
                }
                else if( i < 2*n+r+2 )
                {
                    orders.SetLocal( iLoc, 0, r+1 );
                    firstInds.SetLocal( iLoc, 0, 2*n+1 );
                }
                else if( i < 2*n+r+5 )
                {
                    orders.SetLocal( iLoc, 0, 3 );
                    firstInds.SetLocal( iLoc, 0, 2*n+r+2 );
                }
                else
                {
                    orders.SetLocal( iLoc, 0, 3 );
                    firstInds.SetLocal( iLoc, 0, 2*n+r+5 );
                }
            }
        }

        // Solve the Second-Order Cone problem
        // ===================================
        DistMultiVec<Real> xHat(comm), y(comm), z(comm), s(comm);
        SOCP( A, G, b, cHat, h, orders, firstInds, xHat, y, z, s, ctrl );

        // Extract x from [x; t; s; u; v]
        // ==============================
        x = xHat( IR(0,n), ALL );
    }
    else
    {
        DistSparseMatrix<Real> Sigma(comm);
        Diagonal( Sigma, d );
        Syrk( LOWER, NORMAL, Real(1), F, Real(1), Sigma );
        MakeSymmetric( LOWER, Sigma );
        LongOnlyPortfolio( Sigma, c, gamma, x );
    }
}

#define PROTO(Real) \
  template void LongOnlyPortfolio \
  ( const DistSparseMatrix<Real>& Sigma, \
    const DistMultiVec<Real>& c, \
          Real gamma, \
          DistMultiVec<Real>& x ); \
  template void LongOnlyPortfolio \
  ( const DistMultiVec<Real>& d, \
    const DistSparseMatrix<Real>& F, \
    const DistMultiVec<Real>& c, \
          Real gamma, \
          DistMultiVec<Real>& x );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
