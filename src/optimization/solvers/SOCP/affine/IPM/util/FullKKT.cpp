/*
   Copyright (c) 2009-2016, Jack Poulson
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
// Since (Q_a)^2 = Q_(a^2), we have that 
//
//   -W_i^T W_i = -W_i^2 = -(Q_{sqrt(w_i})^2 = -Q_{w_i}.
//              = -2 w_i w_i^T + det(w_i) R,
//
// which can be rewritten as
//
//   -W_i^2 = det(w_i) (-diag(delta_0,I) - u_i u_i^T + v_i v_i^T),
//
// for some delta_0 > 0, and u_i and v_i satisfying
//
//    diag(delta_0,I) - v_i v_i^T > 0.
//
// Let us temporarily drop the "i" subscript and suppose that det(w) = 1.
// Then, if we expose the first entries of u, v, and w as 
//
//    u = [u_0; u_1], v = [v_0; v_1], w = [w_0, w_1],
//
// and require that u_1, v_1, and w_1 are all proportional to some vector
// f_1, with ||f_1||_2 equal to either 0 or 1, so that
//
//    u_1 = ||u_1||_2 f_1, v_1 = ||v_1||_2 f_1, w_1 = ||w_1||_2 f_1.
//
// Furthermore, we require that v_0 = 0 so that
//
//    Q_w = 2 w w^T - det(w) R 
//
//        = 2 w w^T - R
//
//        = | 2 w_0^2 - 1,  2 w_0 w_1^T     |
//          |   2 w_0 w_1,  I + 2 w_1 w_1^T |
//
//        = | 2 ||w_1||_2^2 + 1,         2 w_0 w_1^T          |
//          |      2 w_0 w_1,     I + 2 ||w_1||_2^2 f_1 f_1^T |
//
//        = | delta_0  0 | + | u_0 | | u_0 u_1^T | - | 0   | | 0 v_1^T |
//          |    0     I |   | u_1 |                 | v_1 |
//
//        = | delta_0 + u_0^2,         u_0 u_1^T          |
//          |      u_0 u_1,     I + u_1 u_1^T - v_1 v_1^T |
//
//        = | delta_0 + u_0^2,                   u_0 u_1^T                 |.
//          |      u_0 u_1,     I + (||u_1||_2^2 - 2||v_1||_2^2) f_1 f_1^T |
//
// Then the expansion of a subcone's portion of the  third block row of the 
// KKT system, say
//
//   G dx - W^T W dz = -r_h + W^T r_mu,
//
// can be expanded into the form
//
//  | G | dx + |  -D    v  -u | |  dz | = | -r_h + W^T r_mu |,
//  | 0 |      |  v^T  -1   0 | | eta |   |        0        |
//  | 0 |      | -u^T   0   1 | |  xi |   |        0        |
//
// using the definition D = diag(delta_0,I), xi = u^T dz, and eta = v^T dz.
// It is straight-forward to verify that the requirements that 
//
//     delta_0 > 0 and D - v v^T > 0
//
// imply the defining parameter ||u_1||_2^2 lies within the interval
//
//     (1/(2 psi_w^2 + 1)) (4 psi_w^4 + 4 psi_w^2, 4 psi_w^4 + 4 psi_w^2 + 1)
//
// with psi_w defined to be ||w_1||_2. We therefore choose psi_u^2 = ||u_1||_2^2
// to lie at the central point,
//
//     psi_u^2 = (4 psi_w^4 + 4 psi_w^2 + 1/2) / (2 psi_w^2 + 1),
//
// and set 
//
//     u_0 = 2 w_0 psi_w / psi_u,
//
//     psi_v = sqrt(psi_u^2 - 2 psi_w^2), and
//
//     delta_0 = 2 psi_w^2 + 1 - u_0^2.
//
// In cases where det(w) != 1, we can use the more general formulation
//
//  | G | dx + det(w) |  -D    v  -u | |  dz | = | -r_h + W^T r_mu |,
//  | 0 |             |  v^T  -1   0 | | eta |   |        0        |
//  | 0 |             | -u^T   0   1 | |  xi |   |        0        |
//
// with {D,u,v} being computed from the normalized Nesterov-Todd scaling point
// w / sqrt(det(w)), i.e.,
//
//     psi_u = sqrt((4 psi_w^4 + 4 psi_w^2 + 1/2) / (2 psi_w^2 + 1)),
//
//     u_0 = 2 w_0 psi_w / psi_u,
//
//     psi_v = sqrt(psi_u^2 - 2 psi_w^2), and
//
//     delta_0 = 2 psi_w^2 + 1 - u_0^2,
//
// with psi_w = || w_1 ||_2 / sqrt(det(w)).
//
// For the sake of completeness, the residuals are defined by:
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

template<typename Real>
void KKT
( const Matrix<Real>& A, 
  const Matrix<Real>& G,
  const Matrix<Real>& w,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& J, 
  bool onlyLower )
{
    DEBUG_ONLY(CSE cse("socp::affine::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    Matrix<Real> wDets;
    soc::Dets( w, wDets, orders, firstInds );

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

template<typename Real>
void KKT
( const ElementalMatrix<Real>& A,    
  const ElementalMatrix<Real>& G,
  const ElementalMatrix<Real>& w,
  const ElementalMatrix<Int>& ordersPre,
  const ElementalMatrix<Int>& firstIndsPre,
        ElementalMatrix<Real>& JPre, 
  bool onlyLower, Int cutoffPar )
{
    DEBUG_ONLY(CSE cse("socp::affine::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();
    const Grid& g = A.Grid();

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrixReadProxy<Int,Int,VC,STAR>
      ordersProx( ordersPre, ctrl ),
      firstIndsProx( firstIndsPre, ctrl );
    auto& orders = ordersProx.GetLocked();
    auto& firstInds = firstIndsProx.GetLocked();

    DistMatrixWriteProxy<Real,Real,MC,MR> JProx( JPre );
    auto& J = JProx.Get();

    DistMatrix<Real,VC,STAR> wDets(g);
    soc::Dets( w, wDets, orders, firstInds, cutoffPar );
    cone::Broadcast( wDets, orders, firstInds, cutoffPar );

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
    // TODO: Create soc::ExplicitQuadratic?
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

// TODO: Incorporate {gamma,delta,beta}
template<typename Real>
void KKT
( const SparseMatrix<Real>& A, 
  const SparseMatrix<Real>& G,
  const Matrix<Real>& w,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  const Matrix<Int>& origToSparseOrders,
  const Matrix<Int>& origToSparseFirstInds,
        Int kSparse,
        SparseMatrix<Real>& J, 
  bool onlyLower )
{
    DEBUG_ONLY(CSE cse("socp::affine::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    // NOTE: The following computation is a bit redundant, and the lower norms
    //       are only needed for sufficiently large cones.
    Matrix<Real> wDets, wLowers;
    soc::Dets( w, wDets, orders, firstInds );
    soc::LowerNorms( w, wLowers, orders, firstInds );
    cone::Broadcast( wDets, orders, firstInds );
    cone::Broadcast( wLowers, orders, firstInds );

    Zeros( J, n+m+kSparse, n+m+kSparse );
    if( onlyLower )
    {
        // Count the number of entries to queue in the lower triangle
        // ----------------------------------------------------------
        Int numEntries = A.NumEntries() + G.NumEntries();
        for( Int i=0; i<k; ++i )
        {
            const Int order = orders.Get(i,0);
            const Int firstInd = firstInds.Get(i,0);
            const Int sparseOrder = origToSparseOrders.Get(i,0);

            if( order == sparseOrder )
            {
                numEntries += (i-firstInd) + 1;
            }
            else
            {
                if( i == firstInd )
                {
                    // An entry of D and u (the first entry of v is 0), as well
                    // as the (scaled) -1 and +1 diagonal entries for the 
                    // v and u auxiliary variables
                    numEntries += 4;
                }
                else
                {
                    // An entry of D, u, and v
                    numEntries += 3;
                }
            }
        }

        // Queue the nonzeros
        // ------------------
        J.Reserve( numEntries );
        for( Int e=0; e<A.NumEntries(); ++e )
            J.QueueUpdate( n+A.Row(e), A.Col(e), A.Value(e) );
        for( Int e=0; e<G.NumEntries(); ++e )
        {
            const Int i = G.Row(e);
            const Int firstInd = firstInds.Get(i,0);
            const Int sparseFirstInd = origToSparseFirstInds.Get(i,0);
            const Int iSparse = i + (sparseFirstInd-firstInd);
            J.QueueUpdate( n+m+iSparse, G.Col(e), G.Value(e) );
        }
        for( Int i=0; i<k; ++i )
        {
            const Int order = orders.Get(i,0);
            const Int firstInd = firstInds.Get(i,0);
            const Int sparseOrder = origToSparseOrders.Get(i,0);
            const Int sparseFirstInd = origToSparseFirstInds.Get(i,0);

            const Int sparseOff = sparseFirstInd-firstInd;
            const Int iSparse = i+sparseOff;

            const Real omega_i = w.Get(i,0);
            const Real wDet = wDets.Get(i,0);
            if( order == sparseOrder )
            {
                // diag(det(w) R - 2 w w^T)
                if( i == firstInd )
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+iSparse, +wDet-2*omega_i*omega_i );
                else
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+iSparse, -wDet-2*omega_i*omega_i );

                // offdiag(-2 w w^T)
                for( Int j=firstInd; j<i; ++j )
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+(j+sparseOff), -2*omega_i*w.Get(j,0) );
            }
            else
            {
                const Int coneOff = n+m+sparseFirstInd;
                const Real wLower = wLowers.Get(i,0);
                const Real wPsi = wLower / Sqrt(wDet);
                const Real wPsiSq = wPsi*wPsi;
                const Real uPsi = 
                  Sqrt((4*wPsiSq*wPsiSq+4*wPsiSq+Real(1)/Real(2))/(2*wPsiSq+1));
                // Apply a pseudoinverse of sorts instead of 1/wPsi
                const Real wPsiPinv = 
                  ( wPsi < limits::Epsilon<Real>() ? Real(1) : 1/wPsi );
                // NOTE: This includes the outer wDet factor
                const Real psiMap = wPsiPinv*omega_i*Sqrt(wDet);
                if( i == firstInd )
                {
                    // Queue up an entry of D and u, and then the (scaled) 
                    // -1 and +1
                    const Real u0 = 2*(omega_i/Sqrt(wDet))*wPsi / uPsi;
                    const Real delta0 = 2*wPsiSq + 1 - u0*u0;
                    J.QueueUpdate
                    ( coneOff,         coneOff,         -wDet*delta0 );
                    J.QueueUpdate
                    ( coneOff+order,   coneOff+order,   -wDet        );
                    J.QueueUpdate
                    ( coneOff+order+1, coneOff,         -wDet*u0     );
                    J.QueueUpdate
                    ( coneOff+order+1, coneOff+order+1,  wDet        );
                }
                else
                {
                    // Queue up an entry of D, u, and v
                    const Real vPsi = Sqrt(uPsi*uPsi-2*wPsiSq);
                    J.QueueUpdate( n+m+iSparse,     n+m+iSparse, -wDet        );
                    J.QueueUpdate( coneOff+order,   n+m+iSparse,  vPsi*psiMap );
                    J.QueueUpdate( coneOff+order+1, n+m+iSparse, -uPsi*psiMap );
                }
            }
        }
    }
    else
    {
        // Count the number of entries to queue
        // ------------------------------------
        Int numEntries = 2*A.NumEntries() + 2*G.NumEntries();
        for( Int i=0; i<k; ++i )
        {
            const Int order = orders.Get(i,0);
            const Int firstInd = firstInds.Get(i,0);
            const Int sparseOrder = origToSparseOrders.Get(i,0);

            if( order == sparseOrder )
            {
                numEntries += order;
            }
            else
            {
                if( i == firstInd )
                {
                    // An entry of D and a symmetric update with an entry of u,
                    // as well as the (scaled) -1 and +1 diagonal entries for 
                    // the v and u auxiliary variables
                    numEntries += 5;
                }
                else
                {
                    // An entry of D and symmetric updates with an entry of
                    // u and v
                    numEntries += 5;
                }
            }
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
            const Int i = G.Row(e);
            const Int firstInd = firstInds.Get(i,0);
            const Int sparseFirstInd = origToSparseFirstInds.Get(i,0);
            const Int iSparse = i + (sparseFirstInd-firstInd);
            J.QueueUpdate( n+m+iSparse, G.Col(e),    G.Value(e) );
            J.QueueUpdate( G.Col(e),    n+m+iSparse, G.Value(e) );
        }
        for( Int i=0; i<k; ++i )
        {
            const Int order = orders.Get(i,0);
            const Int firstInd = firstInds.Get(i,0);
            const Int sparseOrder = origToSparseOrders.Get(i,0);
            const Int sparseFirstInd = origToSparseFirstInds.Get(i,0);

            const Int sparseOff = sparseFirstInd-firstInd;
            const Int iSparse = i+sparseOff;

            const Real omega_i = w.Get(i,0);
            const Real wDet = wDets.Get(i,0);
            if( order == sparseOrder )
            {
                // diag(det(w) R - 2 w w^T)
                if( i == firstInd )
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+iSparse, +wDet-2*omega_i*omega_i );
                else
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+iSparse, -wDet-2*omega_i*omega_i );

                // offdiag(-2 w w^T)
                for( Int j=firstInd; j<firstInd+order; ++j )
                    if( j != i )
                        J.QueueUpdate
                        ( n+m+iSparse, n+m+(j+sparseOff), 
                          -2*omega_i*w.Get(j,0) );
            }
            else
            {
                const Int coneOff = n+m+sparseFirstInd;
                const Real wLower = wLowers.Get(i,0);
                const Real wPsi = wLower / Sqrt(wDet);
                const Real wPsiSq = wPsi*wPsi;
                const Real uPsi = 
                  Sqrt((4*wPsiSq*wPsiSq+4*wPsiSq+Real(1)/Real(2))/(2*wPsiSq+1));
                // Apply a pseudoinverse of sorts instead of 1/wPsi
                const Real wPsiPinv =
                  ( wPsi < limits::Epsilon<Real>() ? Real(1) : 1/wPsi );
                // NOTE: This includes the outer wDet factor
                const Real psiMap = wPsiPinv*omega_i*Sqrt(wDet);
                if( i == firstInd )
                {
                    // Queue up an entry of D, a symmetric update with u, and 
                    // then the (scaled) -1 and +1
                    const Real u0 = 2*(omega_i/Sqrt(wDet))*wPsi / uPsi;
                    const Real delta0 = 2*wPsiSq + 1 - u0*u0;
                    J.QueueUpdate
                    ( coneOff,         coneOff,         -wDet*delta0 );
                    J.QueueUpdate
                    ( coneOff+order,   coneOff+order,   -wDet        );
                    J.QueueUpdate
                    ( coneOff+order+1, coneOff,         -wDet*u0     );
                    J.QueueUpdate
                    ( coneOff,         coneOff+order+1, -wDet*u0     );
                    J.QueueUpdate
                    ( coneOff+order+1, coneOff+order+1,  wDet        );
                }
                else
                {
                    // Queue up an entry of D and symmetric updates with 
                    // u and v
                    const Real vPsi = Sqrt(uPsi*uPsi-2*wPsiSq);
                    J.QueueUpdate
                    ( n+m+iSparse,     n+m+iSparse,     -wDet        ); 
                    J.QueueUpdate
                    ( coneOff+order,   n+m+iSparse,      vPsi*psiMap );
                    J.QueueUpdate
                    ( n+m+iSparse,     coneOff+order,    vPsi*psiMap );
                    J.QueueUpdate
                    ( coneOff+order+1, n+m+iSparse,     -uPsi*psiMap );
                    J.QueueUpdate
                    ( n+m+iSparse,     coneOff+order+1, -uPsi*psiMap );
                }
            }
        }
    }
    J.ProcessQueues();
}

template<typename Real>
void StaticKKT
( const SparseMatrix<Real>& A, 
  const SparseMatrix<Real>& G,
        Real gamma,
        Real delta,
        Real beta,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  const Matrix<Int>& origToSparseOrders,
  const Matrix<Int>& origToSparseFirstInds,
        Int kSparse,
        SparseMatrix<Real>& J, 
  bool onlyLower )
{
    DEBUG_ONLY(CSE cse("socp::affine::StaticKKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();
    const Int numEntriesA = A.NumEntries();
    const Int numEntriesG = G.NumEntries();

    Zeros( J, n+m+kSparse, n+m+kSparse );
    if( onlyLower )
    {
        Int numEntries = numEntriesA + numEntriesG + m + n + kSparse;
        // Add on the number of nonzero off-diagonal entries of the Hessian
        // of the barrier will eventually exist
        for( Int i=0; i<k; ++i )
        {
            const Int order = orders.Get(i,0);
            const Int firstInd = firstInds.Get(i,0);
            const Int sparseOrder = origToSparseOrders.Get(i,0);

            if( order == sparseOrder )
            {
                numEntries += i-firstInd;
            }
            else
            {
                if( i == firstInd )
                {
                    // An entry of u (the corresponding entry of v is zero)
                    numEntries += 1;
                }
                else
                {
                    // An entry of u and v
                    numEntries += 2;
                }
            }
        }

        J.Reserve( numEntries );
        for( Int i=0; i<n; ++i )
            J.QueueUpdate( i, i, gamma*gamma );
        for( Int e=0; e<numEntriesA; ++e )
            J.QueueUpdate( n+A.Row(e), A.Col(e), A.Value(e) );
        for( Int i=0; i<m; ++i )
            J.QueueUpdate( i+n, i+n, -delta*delta );
        for( Int e=0; e<numEntriesG; ++e )
        {
            const Int i = G.Row(e);
            const Int firstInd = firstInds.Get(i,0);
            const Int sparseFirstInd = origToSparseFirstInds.Get(i,0);
            const Int iSparse = i + (sparseFirstInd-firstInd);
            J.QueueUpdate( n+m+iSparse, G.Col(e), G.Value(e) );
        }
        for( Int i=0; i<k; ++i )
        {
            const Int order = orders.Get(i,0);
            const Int firstInd = firstInds.Get(i,0);
            const Int sparseOrder = origToSparseOrders.Get(i,0);
            const Int sparseFirstInd = origToSparseFirstInds.Get(i,0);
            const Int iSparse = i + (sparseFirstInd-firstInd);

            if( order == sparseOrder )
            {
                for( Int jSparse=sparseFirstInd; jSparse<iSparse; ++jSparse )
                    J.QueueUpdate( iSparse, jSparse, Real(0) );
                J.QueueUpdate( iSparse, iSparse, -beta*beta );
            }
            else
            {
                if( i == firstInd )
                {
                    // An entry of u (the corresponding entry of v is zero)
                    J.QueueUpdate( sparseFirstInd+order+1, iSparse, Real(0) );
                    // The first diagonal entry and those of u and v
                    J.QueueUpdate
                    ( sparseFirstInd,       sparseFirstInd,       -beta*beta );
                    J.QueueUpdate
                    ( sparseFirstInd+order, sparseFirstInd+order, -beta*beta );
                    J.QueueUpdate
                    ( sparseFirstInd+order+1, sparseFirstInd+order+1, 
                      beta*beta );
                }
                else
                {
                    // An entry of u and v
                    J.QueueUpdate( sparseFirstInd+order,   iSparse, Real(0) );
                    J.QueueUpdate( sparseFirstInd+order+1, iSparse, Real(0) );
                }
            }
        }
    }
    else
    {
        Int numEntries = 2*numEntriesA + 2*numEntriesG + m + n + kSparse;
        for( Int i=0; i<k; ++i )
        {
            const Int order = orders.Get(i,0);
            const Int firstInd = firstInds.Get(i,0);
            const Int sparseOrder = origToSparseOrders.Get(i,0);

            if( order == sparseOrder )
            {
                numEntries += order-1;
            }
            else
            {
                if( i == firstInd )
                {
                    // A symmetric update with an entry of u
                    numEntries += 2;
                }
                else
                {
                    // Symmetric updates with an entry each of u and v
                    numEntries += 4;
                }
            }
        }

        J.Reserve( numEntries );
        for( Int i=0; i<n; ++i )
            J.QueueUpdate( i, i, gamma*gamma );
        for( Int e=0; e<numEntriesA; ++e )
        {
            J.QueueUpdate( A.Row(e)+n, A.Col(e),   A.Value(e) );
            J.QueueUpdate( A.Col(e),   A.Row(e)+n, A.Value(e) );
        }
        for( Int i=0; i<m; ++i )
            J.QueueUpdate( i+n, i+n, -delta*delta );
        for( Int e=0; e<numEntriesG; ++e )
        {
            const Int i = G.Row(e);
            const Int firstInd = firstInds.Get(i,0);
            const Int sparseFirstInd = origToSparseFirstInds.Get(i,0);
            const Int iSparse = i + (sparseFirstInd-firstInd);
            J.QueueUpdate( n+m+iSparse, G.Col(e),    G.Value(e) );
            J.QueueUpdate( G.Col(e),    n+m+iSparse, G.Value(e) );
        }
        for( Int i=0; i<k; ++i )
        {
            const Int order = orders.Get(i,0);
            const Int firstInd = firstInds.Get(i,0);
            const Int sparseOrder = origToSparseOrders.Get(i,0);
            const Int sparseFirstInd = origToSparseFirstInds.Get(i,0);
            const Int iSparse = i + (sparseFirstInd-firstInd);

            if( order == sparseOrder )
            {
                for( Int jSparse=sparseFirstInd; 
                         jSparse<sparseFirstInd+order; ++jSparse )
                {
                    if( jSparse == iSparse )
                        J.QueueUpdate( iSparse, jSparse, -beta*beta );
                    else
                        J.QueueUpdate( iSparse, jSparse, Real(0) );
                }
            }
            else
            {
                if( i == firstInd )
                {
                    // An entry of u (the corresponding entry of v is zero)
                    J.QueueUpdate( sparseFirstInd+order+1, iSparse, Real(0) );
                    J.QueueUpdate( iSparse, sparseFirstInd+order+1, Real(0) );
                    // The first diagonal entry and those of u and v
                    J.QueueUpdate
                    ( sparseFirstInd,       sparseFirstInd,       -beta*beta );
                    J.QueueUpdate
                    ( sparseFirstInd+order, sparseFirstInd+order, -beta*beta );
                    J.QueueUpdate
                    ( sparseFirstInd+order+1, sparseFirstInd+order+1, 
                      beta*beta );
                }
                else
                {
                    // An entry of u and v
                    J.QueueUpdate( sparseFirstInd+order, iSparse, Real(0) );
                    J.QueueUpdate( iSparse, sparseFirstInd+order, Real(0) );
                    J.QueueUpdate( sparseFirstInd+order+1, iSparse, Real(0) );
                    J.QueueUpdate( iSparse, sparseFirstInd+order+1, Real(0) );
                }
            }
        }
    }
    J.ProcessQueues();
    J.FreezeSparsity();
}

template<typename Real>
void FinishKKT
( Int m, Int n, 
  const Matrix<Real>& w,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  const Matrix<Int>& origToSparseOrders,
  const Matrix<Int>& origToSparseFirstInds,
        Int kSparse,
        SparseMatrix<Real>& J, 
  bool onlyLower )
{
    DEBUG_ONLY(CSE cse("socp::affine::FinishKKT"))
    const Int k = w.Height();

    // NOTE: The following computation is a bit redundant, and the lower norms
    //       are only needed for sufficiently large cones.
    Matrix<Real> wDets, wLowers;
    soc::Dets( w, wDets, orders, firstInds );
    soc::LowerNorms( w, wLowers, orders, firstInds );
    cone::Broadcast( wDets, orders, firstInds );
    cone::Broadcast( wLowers, orders, firstInds );

    Zeros( J, n+m+kSparse, n+m+kSparse );
    if( onlyLower )
    {
        // Count the number of entries to queue in the lower triangle
        // ----------------------------------------------------------
        Int numEntries = 0;
        for( Int i=0; i<k; ++i )
        {
            const Int order = orders.Get(i,0);
            const Int firstInd = firstInds.Get(i,0);
            const Int sparseOrder = origToSparseOrders.Get(i,0);

            if( order == sparseOrder )
            {
                numEntries += (i-firstInd) + 1;
            }
            else
            {
                if( i == firstInd )
                {
                    // An entry of D and u (the first entry of v is 0), as well
                    // as the (scaled) -1 and +1 diagonal entries for the 
                    // v and u auxiliary variables
                    numEntries += 4;
                }
                else
                {
                    // An entry of D, u, and v
                    numEntries += 3;
                }
            }
        }

        // Queue the nonzeros
        // ------------------
        if( !J.FrozenSparsity() )
            J.Reserve( J.NumEntries()+numEntries );
        for( Int i=0; i<k; ++i )
        {
            const Int order = orders.Get(i,0);
            const Int firstInd = firstInds.Get(i,0);
            const Int sparseOrder = origToSparseOrders.Get(i,0);
            const Int sparseFirstInd = origToSparseFirstInds.Get(i,0);

            const Int sparseOff = sparseFirstInd-firstInd;
            const Int iSparse = i+sparseOff;

            const Real omega_i = w.Get(i,0);
            const Real wDet = wDets.Get(i,0);
            if( order == sparseOrder )
            {
                // diag(det(w) R - 2 w w^T)
                if( i == firstInd )
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+iSparse, +wDet-2*omega_i*omega_i );
                else
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+iSparse, -wDet-2*omega_i*omega_i );

                // offdiag(-2 w w^T)
                for( Int j=firstInd; j<i; ++j )
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+(j+sparseOff), -2*omega_i*w.Get(j,0) );
            }
            else
            {
                const Int coneOff = n+m+sparseFirstInd;
                const Real wLower = wLowers.Get(i,0);
                const Real wPsi = wLower / Sqrt(wDet);
                const Real wPsiSq = wPsi*wPsi;
                const Real uPsi = 
                  Sqrt((4*wPsiSq*wPsiSq+4*wPsiSq+Real(1)/Real(2))/(2*wPsiSq+1));
                // Apply a pseudoinverse of sorts instead of 1/wPsi
                const Real wPsiPinv = 
                  ( wPsi < limits::Epsilon<Real>() ? Real(1) : 1/wPsi );
                // NOTE: This includes the outer wDet factor
                const Real psiMap = wPsiPinv*omega_i*Sqrt(wDet);
                if( i == firstInd )
                {
                    // Queue up an entry of D and u, and then the (scaled) 
                    // -1 and +1
                    const Real u0 = 2*(omega_i/Sqrt(wDet))*wPsi / uPsi;
                    const Real delta0 = 2*wPsiSq + 1 - u0*u0;
                    J.QueueUpdate
                    ( coneOff,         coneOff,         -wDet*delta0 );
                    J.QueueUpdate
                    ( coneOff+order,   coneOff+order,   -wDet        );
                    J.QueueUpdate
                    ( coneOff+order+1, coneOff,         -wDet*u0     );
                    J.QueueUpdate
                    ( coneOff+order+1, coneOff+order+1,  wDet        );
                }
                else
                {
                    // Queue up an entry of D, u, and v
                    const Real vPsi = Sqrt(uPsi*uPsi-2*wPsiSq);
                    J.QueueUpdate( n+m+iSparse,     n+m+iSparse, -wDet        );
                    J.QueueUpdate( coneOff+order,   n+m+iSparse,  vPsi*psiMap );
                    J.QueueUpdate( coneOff+order+1, n+m+iSparse, -uPsi*psiMap );
                }
            }
        }
    }
    else
    {
        // Count the number of entries to queue
        // ------------------------------------
        Int numEntries = 0;
        for( Int i=0; i<k; ++i )
        {
            const Int order = orders.Get(i,0);
            const Int firstInd = firstInds.Get(i,0);
            const Int sparseOrder = origToSparseOrders.Get(i,0);

            if( order == sparseOrder )
            {
                numEntries += order;
            }
            else
            {
                if( i == firstInd )
                {
                    // An entry of D and a symmetric update with an entry of u,
                    // as well as the (scaled) -1 and +1 diagonal entries for 
                    // the v and u auxiliary variables
                    numEntries += 5;
                }
                else
                {
                    // An entry of D and symmetric updates with an entry of
                    // u and v
                    numEntries += 5;
                }
            }
        }

        // Queue the nonzeros
        // ------------------
        if( !J.FrozenSparsity() )
            J.Reserve( J.NumEntries()+numEntries );
        for( Int i=0; i<k; ++i )
        {
            const Int order = orders.Get(i,0);
            const Int firstInd = firstInds.Get(i,0);
            const Int sparseOrder = origToSparseOrders.Get(i,0);
            const Int sparseFirstInd = origToSparseFirstInds.Get(i,0);

            const Int sparseOff = sparseFirstInd-firstInd;
            const Int iSparse = i+sparseOff;

            const Real omega_i = w.Get(i,0);
            const Real wDet = wDets.Get(i,0);
            if( order == sparseOrder )
            {
                // diag(det(w) R - 2 w w^T)
                if( i == firstInd )
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+iSparse, +wDet-2*omega_i*omega_i );
                else
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+iSparse, -wDet-2*omega_i*omega_i );

                // offdiag(-2 w w^T)
                for( Int j=firstInd; j<firstInd+order; ++j )
                    if( j != i )
                        J.QueueUpdate
                        ( n+m+iSparse, n+m+(j+sparseOff), 
                          -2*omega_i*w.Get(j,0) );
            }
            else
            {
                const Int coneOff = n+m+sparseFirstInd;
                const Real wLower = wLowers.Get(i,0);
                const Real wPsi = wLower / Sqrt(wDet);
                const Real wPsiSq = wPsi*wPsi;
                const Real uPsi = 
                  Sqrt((4*wPsiSq*wPsiSq+4*wPsiSq+Real(1)/Real(2))/(2*wPsiSq+1));
                // Apply a pseudoinverse of sorts instead of 1/wPsi
                const Real wPsiPinv =
                  ( wPsi < limits::Epsilon<Real>() ? Real(1) : 1/wPsi );
                // NOTE: This includes the outer wDet factor
                const Real psiMap = wPsiPinv*omega_i*Sqrt(wDet);
                if( i == firstInd )
                {
                    // Queue up an entry of D, a symmetric update with u, and 
                    // then the (scaled) -1 and +1
                    const Real u0 = 2*(omega_i/Sqrt(wDet))*wPsi / uPsi;
                    const Real delta0 = 2*wPsiSq + 1 - u0*u0;
                    J.QueueUpdate
                    ( coneOff,         coneOff,         -wDet*delta0 );
                    J.QueueUpdate
                    ( coneOff+order,   coneOff+order,   -wDet        );
                    J.QueueUpdate
                    ( coneOff+order+1, coneOff,         -wDet*u0     );
                    J.QueueUpdate
                    ( coneOff,         coneOff+order+1, -wDet*u0     );
                    J.QueueUpdate
                    ( coneOff+order+1, coneOff+order+1,  wDet        );
                }
                else
                {
                    // Queue up an entry of D and symmetric updates with 
                    // u, and v
                    const Real vPsi = Sqrt(uPsi*uPsi-2*wPsiSq);
                    J.QueueUpdate
                    ( n+m+iSparse,     n+m+iSparse,     -wDet        ); 
                    J.QueueUpdate
                    ( coneOff+order,   n+m+iSparse,      vPsi*psiMap );
                    J.QueueUpdate
                    ( n+m+iSparse,     coneOff+order,    vPsi*psiMap );
                    J.QueueUpdate
                    ( coneOff+order+1, n+m+iSparse,     -uPsi*psiMap );
                    J.QueueUpdate
                    ( n+m+iSparse,     coneOff+order+1, -uPsi*psiMap );
                }
            }
        }
    }
    J.ProcessQueues();
}

// TODO: Incorporate {gamma,delta,beta}
template<typename Real>
void KKT
( const DistSparseMatrix<Real>& A, 
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  const DistMultiVec<Int>& origToSparseOrders,
  const DistMultiVec<Int>& origToSparseFirstInds,
        Int kSparse,
        DistSparseMatrix<Real>& J, 
  bool onlyLower, Int cutoffPar )
{
    DEBUG_ONLY(CSE cse("socp::affine::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = w.Comm();
    const int commSize = mpi::Size(comm);
    const int commRank = mpi::Rank(comm);

    // NOTE: The following computation is a bit redundant, and the lower norms
    //       are only needed for sufficiently large cones.
    DistMultiVec<Real> wDets(comm), wLowers(comm);
    soc::Dets( w, wDets, orders, firstInds, cutoffPar );
    soc::LowerNorms( w, wLowers, orders, firstInds, cutoffPar );
    cone::Broadcast( wDets, orders, firstInds, cutoffPar );
    cone::Broadcast( wLowers, orders, firstInds, cutoffPar );

    // Gather all of each non-sparsified member cone that we own a piece of
    // --------------------------------------------------------------------
    // NOTE: We exploit the fact that each process owns a contiguous chunk
    //       of rows 
    vector<int> recvOffs;
    vector<Real> recvBuf;
    const Int wLocalHeight = w.LocalHeight();
    {
        // TODO: Rewrite this using QueuePull!
        // ===================================

        // Count the total number of entries of w to send to each process
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        vector<int> sendCounts(commSize,0);
        for( Int iLoc=0; iLoc<wLocalHeight; )
        {
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseOrder = origToSparseOrders.GetLocal(iLoc,0);

            const Int i = w.GlobalRow(iLoc);
            const Int numLocalIndsLeft = wLocalHeight-iLoc;
            const Int numLocalIndsCone = 
              Min(numLocalIndsLeft,order-(i-firstInd));

            if( order == sparseOrder )
            {
                const int firstOwner = w.RowOwner(firstInd);
                const int lastOwner = w.RowOwner(firstInd+order-1);
                for( Int owner=firstOwner; owner<=lastOwner; ++owner )
                    if( owner != commRank )
                        sendCounts[owner] += numLocalIndsCone;
            }
            
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
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseOrder = origToSparseOrders.GetLocal(iLoc,0);

            const Int i = w.GlobalRow(iLoc);
            const Int numLocalIndsLeft = wLocalHeight-iLoc;
            const Int numLocalIndsCone = 
              Min(numLocalIndsLeft,order-(i-firstInd));

            if( order == sparseOrder )
            {
                const int firstOwner = w.RowOwner(firstInd);
                const int lastOwner = w.RowOwner(firstInd+order-1);
                for( Int owner=firstOwner; owner<=lastOwner; ++owner )
                    if( owner != commRank )
                        for( Int e=0; e<numLocalIndsCone; ++e )
                            sendBuf[offs[owner]++] = w.GetLocal(iLoc+e,0);
            }
            
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
    Zeros( J, n+m+kSparse, n+m+kSparse );
    if( onlyLower )
    {
        // Count the number of entries to queue in the lower triangle
        // ----------------------------------------------------------
        Int numRemoteEntries = A.NumLocalEntries() + G.NumLocalEntries();
        for( Int iLoc=0; iLoc<wLocalHeight; ++iLoc )
        {
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseOrder = origToSparseOrders.GetLocal(iLoc,0);

            const Int i = w.GlobalRow(iLoc);
            if( order == sparseOrder )
            {
                numRemoteEntries += (i-firstInd) + 1;
            }
            else
            {
                if( i == firstInd )
                {
                    // An entry of D and u (the first entry of v is 0), as well
                    // as the (scaled) -1 and +1 diagonal entries for the 
                    // v and u auxiliary variables
                    numRemoteEntries += 4;
                }
                else
                {
                    // An entry of D, u, and v
                    numRemoteEntries += 3;
                }
            }
        }

        // Queue the nonzeros
        // ------------------
        J.Reserve( numRemoteEntries, numRemoteEntries );
        for( Int e=0; e<A.NumLocalEntries(); ++e )
            J.QueueUpdate( n+A.Row(e), A.Col(e), A.Value(e) );
        for( Int e=0; e<G.NumLocalEntries(); ++e )
        {
            const Int i = G.Row(e);
            const Int iLoc = G.LocalRow(i);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseFirstInd = origToSparseFirstInds.GetLocal(iLoc,0);
            const Int iSparse = i + (sparseFirstInd-firstInd);
            J.QueueUpdate( n+m+iSparse, G.Col(e), G.Value(e) );
        }
        Int lastFirstInd = -1;
        vector<Real> wBuf;
        auto offs = recvOffs;
        for( Int iLoc=0; iLoc<wLocalHeight; ++iLoc )
        {
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseOrder = origToSparseOrders.GetLocal(iLoc,0);
            const Int sparseFirstInd = origToSparseFirstInds.GetLocal(iLoc,0);

            const Int i = w.GlobalRow(iLoc);
            const Int sparseOff = sparseFirstInd-firstInd;
            const Int iSparse = i+sparseOff;

            const Real omega_i = w.GetLocal(iLoc,0);
            const Real wDet = wDets.GetLocal(iLoc,0);
            if( order == sparseOrder )
            {
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
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+iSparse, +wDet-2*omega_i*omega_i );
                else
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+iSparse, -wDet-2*omega_i*omega_i );

                // offdiag(-2 w w^T)
                for( Int j=firstInd; j<i; ++j )
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+(j+sparseOff), 
                      -2*omega_i*wBuf[j-firstInd] );
            }
            else
            {
                const Int coneOff = n+m+sparseFirstInd;
                const Real wLower = wLowers.GetLocal(iLoc,0);
                const Real wPsi = wLower / Sqrt(wDet);
                const Real wPsiSq = wPsi*wPsi;
                const Real uPsi = 
                  Sqrt((4*wPsiSq*wPsiSq+4*wPsiSq+Real(1)/Real(2))/(2*wPsiSq+1));
                // Apply a pseudoinverse of sorts instead of 1/wPsi
                const Real wPsiPinv = 
                  ( wPsi < limits::Epsilon<Real>() ? Real(1) : 1/wPsi );
                // NOTE: This includes the outer wDet factor
                const Real psiMap = wPsiPinv*omega_i*Sqrt(wDet);
                if( i == firstInd )
                {
                    // Queue up an entry of D and u, and then the (scaled) 
                    // -1 and +1
                    const Real u0 = 2*(omega_i/Sqrt(wDet))*wPsi / uPsi;
                    const Real delta0 = 2*wPsiSq + 1 - u0*u0;
                    J.QueueUpdate( coneOff,         coneOff,     -wDet*delta0 );
                    J.QueueUpdate( coneOff+order,   coneOff+order,   -wDet );
                    J.QueueUpdate( coneOff+order+1, coneOff,         -wDet*u0 );
                    J.QueueUpdate( coneOff+order+1, coneOff+order+1,  wDet );
                }
                else
                {
                    // Queue up an entry of D, u, and v
                    const Real vPsi = Sqrt(uPsi*uPsi-2*wPsiSq);
                    J.QueueUpdate( n+m+iSparse,     n+m+iSparse, -wDet );
                    J.QueueUpdate( coneOff+order,   n+m+iSparse,  vPsi*psiMap );
                    J.QueueUpdate( coneOff+order+1, n+m+iSparse, -uPsi*psiMap );
                }
            }
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
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseOrder = origToSparseOrders.GetLocal(iLoc,0);

            const Int i = w.GlobalRow(iLoc);
            if( order == sparseOrder )
            {
                numRemoteEntries += order;
            }
            else
            {
                if( i == firstInd )
                {
                    // An entry of D and a symmetric update with an entry of u,
                    // as well as the (scaled) -1 and +1 diagonal entries for 
                    // the v and u auxiliary variables
                    numRemoteEntries += 5;
                }
                else
                {
                    // An entry of D and symmetric updates with an entry of
                    // u and v
                    numRemoteEntries += 5;
                }
            }
        }

        // Queue the nonzeros
        // ------------------
        J.Reserve( numRemoteEntries, numRemoteEntries );
        for( Int e=0; e<A.NumLocalEntries(); ++e )
        {
            J.QueueUpdate( A.Row(e)+n, A.Col(e),   A.Value(e) );
            J.QueueUpdate( A.Col(e),   A.Row(e)+n, A.Value(e) );
        }
        for( Int e=0; e<G.NumLocalEntries(); ++e )
        {
            const Int i = G.Row(e);
            const Int iLoc = G.LocalRow(i);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseFirstInd = origToSparseFirstInds.GetLocal(iLoc,0);
            const Int iSparse = i + (sparseFirstInd-firstInd);
            J.QueueUpdate( n+m+iSparse, G.Col(e),    G.Value(e) );
            J.QueueUpdate( G.Col(e),    n+m+iSparse, G.Value(e) );
        }
        Int lastFirstInd = -1;
        vector<Real> wBuf;
        auto offs = recvOffs;
        for( Int iLoc=0; iLoc<wLocalHeight; ++iLoc )
        {
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseOrder = origToSparseOrders.GetLocal(iLoc,0);
            const Int sparseFirstInd = origToSparseFirstInds.GetLocal(iLoc,0);

            const Int i = w.GlobalRow(iLoc);
            const Int sparseOff = sparseFirstInd-firstInd;
            const Int iSparse = i+sparseOff;

            const Real omega_i = w.GetLocal(iLoc,0);
            const Real wDet = wDets.GetLocal(iLoc,0);
            if( order == sparseOrder )
            {
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
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+iSparse, +wDet-2*omega_i*omega_i );
                else
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+iSparse, -wDet-2*omega_i*omega_i );

                // offdiag(-2 w w^T)
                for( Int j=firstInd; j<firstInd+order; ++j )
                    if( j != i )
                        J.QueueUpdate
                        ( n+m+iSparse, n+m+(j+sparseOff), 
                          -2*omega_i*wBuf[j-firstInd] );
            }
            else
            {
                const Int coneOff = n+m+sparseFirstInd;
                const Real wLower = wLowers.GetLocal(iLoc,0);
                const Real wPsi = wLower / Sqrt(wDet);
                const Real wPsiSq = wPsi*wPsi;
                const Real uPsi = 
                  Sqrt((4*wPsiSq*wPsiSq+4*wPsiSq+Real(1)/Real(2))/(2*wPsiSq+1));
                // Apply a pseudoinverse of sorts instead of 1/wPsi
                const Real wPsiPinv =
                  ( wPsi < limits::Epsilon<Real>() ? Real(1) : 1/wPsi );
                // NOTE: This includes the outer wDet factor
                const Real psiMap = wPsiPinv*omega_i*Sqrt(wDet);
                if( i == firstInd )
                {
                    // Queue up an entry of D, a symmetric update with u, and 
                    // then the (scaled) -1 and +1
                    const Real u0 = 2*(omega_i/Sqrt(wDet))*wPsi / uPsi;
                    const Real delta0 = 2*wPsiSq + 1 - u0*u0;
                    J.QueueUpdate
                    ( coneOff,         coneOff,         -wDet*delta0 );
                    J.QueueUpdate
                    ( coneOff+order,   coneOff+order,   -wDet );
                    J.QueueUpdate
                    ( coneOff+order+1, coneOff,         -wDet*u0 );
                    J.QueueUpdate
                    ( coneOff,         coneOff+order+1, -wDet*u0 );
                    J.QueueUpdate
                    ( coneOff+order+1, coneOff+order+1,  wDet );
                }
                else
                {
                    // Queue up an entry of D and symmetric updates with 
                    // u, and v
                    const Real vPsi = Sqrt(uPsi*uPsi-2*wPsiSq);
                    J.QueueUpdate
                    ( n+m+iSparse,     n+m+iSparse,     -wDet ); 
                    J.QueueUpdate
                    ( coneOff+order,   n+m+iSparse,      vPsi*psiMap );
                    J.QueueUpdate
                    ( n+m+iSparse,     coneOff+order,    vPsi*psiMap );
                    J.QueueUpdate
                    ( coneOff+order+1, n+m+iSparse,     -uPsi*psiMap );
                    J.QueueUpdate
                    ( n+m+iSparse,     coneOff+order+1, -uPsi*psiMap );
                }
            }
        }
    }
    J.ProcessQueues();
}

template<typename Real>
void StaticKKT
( const DistSparseMatrix<Real>& A, 
  const DistSparseMatrix<Real>& G,
        Real gamma,
        Real delta,
        Real beta,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  const DistMultiVec<Int>& origToSparseOrders,
  const DistMultiVec<Int>& origToSparseFirstInds,
        Int kSparse,
        DistSparseMatrix<Real>& J, 
  bool onlyLower )
{
    DEBUG_ONLY(CSE cse("socp::affine::StaticKKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int coneLocalHeight = orders.LocalHeight();
    const Int numEntriesA = A.NumLocalEntries();
    const Int numEntriesG = G.NumLocalEntries();

    mpi::Comm comm = A.Comm();
    J.SetComm( comm );
    Zeros( J, n+m+kSparse, n+m+kSparse );
    const Int JLocalHeight = J.LocalHeight();

    Int numLocalEntries = 0;
    for( Int iLoc=0; iLoc<JLocalHeight; ++iLoc )
    {
        const Int i = J.GlobalRow(iLoc);
        if( i < n+m )    
            ++numLocalEntries;
        else
            break;
    }

    if( onlyLower )
    {
        Int numRemoteEntries = numEntriesA + numEntriesG;
        for( Int iLoc=0; iLoc<coneLocalHeight; ++iLoc )
        {
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseOrder = origToSparseOrders.GetLocal(iLoc,0);

            const Int i = orders.GlobalRow(iLoc);
            if( order == sparseOrder )
            {
                numRemoteEntries += (i-firstInd) + 1;
            }
            else
            {
                if( i == firstInd )
                {
                    // An entry of D and u (the first entry of v is 0), as well
                    // as the (scaled) -1 and +1 diagonal entries for the 
                    // v and u auxiliary variables
                    numRemoteEntries += 4;
                }
                else
                {
                    // An entry of D, u, and v
                    numRemoteEntries += 3;
                }
            }
        }

        J.Reserve( numLocalEntries+numRemoteEntries, numRemoteEntries );

        for( Int iLoc=0; iLoc<JLocalHeight; ++iLoc )
        {
            const Int i = J.GlobalRow(iLoc);
            if( i < n )        J.QueueLocalUpdate( iLoc, i,  gamma*gamma );
            else if( i < n+m ) J.QueueLocalUpdate( iLoc, i, -delta*delta );
            else break;
        }

        for( Int e=0; e<numEntriesA; ++e )
            J.QueueUpdate( n+A.Row(e), A.Col(e), A.Value(e) );

        for( Int e=0; e<numEntriesG; ++e )
        {
            const Int i = G.Row(e);
            const Int iLoc = G.LocalRow(i);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseFirstInd = origToSparseFirstInds.GetLocal(iLoc,0);
            const Int iSparse = i + (sparseFirstInd-firstInd);
            J.QueueUpdate( n+m+iSparse, G.Col(e), G.Value(e) );
        }

        for( Int iLoc=0; iLoc<coneLocalHeight; ++iLoc )
        {
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseOrder = origToSparseOrders.GetLocal(iLoc,0);
            const Int sparseFirstInd = origToSparseFirstInds.GetLocal(iLoc,0);

            const Int i = orders.GlobalRow(iLoc);
            const Int sparseOff = sparseFirstInd-firstInd;
            const Int iSparse = i+sparseOff;

            if( order == sparseOrder )
            {
                // offdiag(-2 w w^T)
                for( Int j=firstInd; j<i; ++j )
                    J.QueueUpdate( n+m+iSparse, n+m+(j+sparseOff), Real(0) );

                // diag(det(w) R - 2 w w^T)
                J.QueueUpdate( n+m+iSparse, n+m+iSparse, -beta*beta );
            }
            else
            {
                const Int coneOff = n+m+sparseFirstInd;
                if( i == firstInd )
                {
                    // Queue up an entry of D and u, and then the (scaled) 
                    // -1 and +1
                    J.QueueUpdate
                    ( coneOff,         coneOff,       -beta*beta );
                    J.QueueUpdate
                    ( coneOff+order,   coneOff+order, -beta*beta );
                    J.QueueUpdate
                    ( coneOff+order+1, coneOff,         Real(0) );
                    J.QueueUpdate
                    ( coneOff+order+1, coneOff+order+1, beta*beta );
                }
                else
                {
                    // Queue up an entry of D, u, and v
                    J.QueueUpdate( n+m+iSparse,     n+m+iSparse, -beta*beta );
                    J.QueueUpdate( coneOff+order,   n+m+iSparse,  Real(0) );
                    J.QueueUpdate( coneOff+order+1, n+m+iSparse,  Real(0) );
                }
            }
        }
    }
    else
    {
        Int numRemoteEntries = 2*numEntriesA + 2*numEntriesG;
        for( Int iLoc=0; iLoc<coneLocalHeight; ++iLoc )
        {
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseOrder = origToSparseOrders.GetLocal(iLoc,0);

            const Int i = orders.GlobalRow(iLoc);
            if( order == sparseOrder )
            {
                numRemoteEntries += order;
            }
            else
            {
                if( i == firstInd )
                {
                    // An entry of D and u (the first entry of v is 0), as well
                    // as the (scaled) -1 and +1 diagonal entries for the 
                    // v and u auxiliary variables
                    numRemoteEntries += 5;
                }
                else
                {
                    // An entry of D, u, and v (and the transposes for u and v)
                    numRemoteEntries += 5;
                }
            }
        }

        J.Reserve( numLocalEntries+numRemoteEntries, numRemoteEntries );

        for( Int iLoc=0; iLoc<JLocalHeight; ++iLoc )
        {
            const Int i = J.GlobalRow(iLoc);
            if( i < n )        J.QueueLocalUpdate( iLoc, i,  gamma*gamma );
            else if( i < n+m ) J.QueueLocalUpdate( iLoc, i, -delta*delta );
            else break;
        }

        for( Int e=0; e<numEntriesA; ++e )
        {
            J.QueueUpdate( A.Row(e)+n, A.Col(e),   A.Value(e) );
            J.QueueUpdate( A.Col(e),   A.Row(e)+n, A.Value(e) );
        }
        for( Int e=0; e<numEntriesG; ++e )
        {
            const Int i = G.Row(e);
            const Int iLoc = G.LocalRow(i);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseFirstInd = origToSparseFirstInds.GetLocal(iLoc,0);
            const Int iSparse = i + (sparseFirstInd-firstInd);
            J.QueueUpdate( n+m+iSparse, G.Col(e),    G.Value(e) );
            J.QueueUpdate( G.Col(e),    n+m+iSparse, G.Value(e) );
        }

        for( Int iLoc=0; iLoc<coneLocalHeight; ++iLoc )
        {
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseOrder = origToSparseOrders.GetLocal(iLoc,0);
            const Int sparseFirstInd = origToSparseFirstInds.GetLocal(iLoc,0);

            const Int i = orders.GlobalRow(iLoc);
            const Int sparseOff = sparseFirstInd-firstInd;
            const Int iSparse = i+sparseOff;

            if( order == sparseOrder )
            {
                for( Int j=firstInd; j<firstInd+order; ++j )
                {
                    if( j == i )
                        // diag(det(w) R - 2 w w^T)
                        J.QueueUpdate( n+m+iSparse, n+m+iSparse, -beta*beta );
                    else
                        // offdiag(-2 w w^T)
                        J.QueueUpdate
                        ( n+m+iSparse, n+m+(j+sparseOff), Real(0) );
                }
            }
            else
            {
                const Int coneOff = n+m+sparseFirstInd;
                if( i == firstInd )
                {
                    // Queue up an entry of D and u, and then the (scaled) 
                    // -1 and +1
                    J.QueueUpdate
                    ( coneOff,         coneOff,       -beta*beta );
                    J.QueueUpdate
                    ( coneOff+order,   coneOff+order, -beta*beta );
                    J.QueueUpdate
                    ( coneOff+order+1, coneOff,         Real(0) );
                    J.QueueUpdate
                    ( coneOff, coneOff+order+1,         Real(0) );
                    J.QueueUpdate
                    ( coneOff+order+1, coneOff+order+1, beta*beta );
                }
                else
                {
                    // Queue up an entry of D, u, and v
                    J.QueueUpdate( n+m+iSparse,     n+m+iSparse,  -beta*beta );
                    J.QueueUpdate( coneOff+order,   n+m+iSparse,     Real(0) );
                    J.QueueUpdate( n+m+iSparse,     coneOff+order,   Real(0) );
                    J.QueueUpdate( coneOff+order+1, n+m+iSparse,     Real(0) );
                    J.QueueUpdate( n+m+iSparse,     coneOff+order+1, Real(0) );
                }
            }
        }
    }
    J.ProcessQueues();
    J.FreezeSparsity();
}

template<typename Real>
void FinishKKT
( Int m, Int n,
  const DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  const DistMultiVec<Int>& origToSparseOrders,
  const DistMultiVec<Int>& origToSparseFirstInds,
        Int kSparse,
        DistSparseMatrix<Real>& J, 
  bool onlyLower, Int cutoffPar )
{
    DEBUG_ONLY(CSE cse("socp::affine::FinishKKT"))
    mpi::Comm comm = w.Comm();
    const int commSize = mpi::Size(comm);
    const int commRank = mpi::Rank(comm);

    // NOTE: The following computation is a bit redundant, and the lower norms
    //       are only needed for sufficiently large cones.
    DistMultiVec<Real> wDets(comm), wLowers(comm);
    soc::Dets( w, wDets, orders, firstInds, cutoffPar );
    soc::LowerNorms( w, wLowers, orders, firstInds, cutoffPar );
    cone::Broadcast( wDets, orders, firstInds, cutoffPar );
    cone::Broadcast( wLowers, orders, firstInds, cutoffPar );

    // Gather all of each non-sparsified member cone that we own a piece of
    // --------------------------------------------------------------------
    // NOTE: We exploit the fact that each process owns a contiguous chunk
    //       of rows 
    vector<int> recvOffs;
    vector<Real> recvBuf;
    const Int wLocalHeight = w.LocalHeight();
    {
        // TODO: Reimplement this using the new QueuePull interface!
        // =========================================================

        // Count the total number of entries of w to send to each process
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        vector<int> sendCounts(commSize,0);
        for( Int iLoc=0; iLoc<wLocalHeight; )
        {
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseOrder = origToSparseOrders.GetLocal(iLoc,0);

            const Int i = w.GlobalRow(iLoc);
            const Int numLocalIndsLeft = wLocalHeight-iLoc;
            const Int numLocalIndsCone = 
              Min(numLocalIndsLeft,order-(i-firstInd));

            if( order == sparseOrder )
            {
                const int firstOwner = w.RowOwner(firstInd);
                const int lastOwner = w.RowOwner(firstInd+order-1);
                for( Int owner=firstOwner; owner<=lastOwner; ++owner )
                    if( owner != commRank )
                        sendCounts[owner] += numLocalIndsCone;
            }
            
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
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseOrder = origToSparseOrders.GetLocal(iLoc,0);

            const Int i = w.GlobalRow(iLoc);
            const Int numLocalIndsLeft = wLocalHeight-iLoc;
            const Int numLocalIndsCone = 
              Min(numLocalIndsLeft,order-(i-firstInd));

            if( order == sparseOrder )
            {
                const int firstOwner = w.RowOwner(firstInd);
                const int lastOwner = w.RowOwner(firstInd+order-1);
                for( Int owner=firstOwner; owner<=lastOwner; ++owner )
                    if( owner != commRank )
                        for( Int e=0; e<numLocalIndsCone; ++e )
                            sendBuf[offs[owner]++] = w.GetLocal(iLoc+e,0);
            }
            
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

    if( onlyLower )
    {
        // Count the number of entries to queue in the lower triangle
        // ----------------------------------------------------------
        Int numRemoteEntries = 0;
        for( Int iLoc=0; iLoc<wLocalHeight; ++iLoc )
        {
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseOrder = origToSparseOrders.GetLocal(iLoc,0);

            const Int i = w.GlobalRow(iLoc);
            if( order == sparseOrder )
            {
                numRemoteEntries += (i-firstInd) + 1;
            }
            else
            {
                if( i == firstInd )
                {
                    // An entry of D and u (the first entry of v is 0), as well
                    // as the (scaled) -1 and +1 diagonal entries for the 
                    // v and u auxiliary variables
                    numRemoteEntries += 4;
                }
                else
                {
                    // An entry of D, u, and v
                    numRemoteEntries += 3;
                }
            }
        }

        // Queue the nonzeros
        // ------------------
        if( !J.FrozenSparsity() )
            J.Reserve( J.NumLocalEntries()+numRemoteEntries, numRemoteEntries );
        Int lastFirstInd = -1;
        vector<Real> wBuf;
        auto offs = recvOffs;
        for( Int iLoc=0; iLoc<wLocalHeight; ++iLoc )
        {
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseOrder = origToSparseOrders.GetLocal(iLoc,0);
            const Int sparseFirstInd = origToSparseFirstInds.GetLocal(iLoc,0);

            const Int i = w.GlobalRow(iLoc);
            const Int sparseOff = sparseFirstInd-firstInd;
            const Int iSparse = i+sparseOff;

            const Real omega_i = w.GetLocal(iLoc,0);
            const Real wDet = wDets.GetLocal(iLoc,0);
            if( order == sparseOrder )
            {
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
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+iSparse, +wDet-2*omega_i*omega_i );
                else
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+iSparse, -wDet-2*omega_i*omega_i );

                // offdiag(-2 w w^T)
                for( Int j=firstInd; j<i; ++j )
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+(j+sparseOff), 
                      -2*omega_i*wBuf[j-firstInd] );
            }
            else
            {
                const Int coneOff = n+m+sparseFirstInd;
                const Real wLower = wLowers.GetLocal(iLoc,0);
                const Real wPsi = wLower / Sqrt(wDet);
                const Real wPsiSq = wPsi*wPsi;
                const Real uPsi = 
                  Sqrt((4*wPsiSq*wPsiSq+4*wPsiSq+Real(1)/Real(2))/(2*wPsiSq+1));
                // Apply a pseudoinverse of sorts instead of 1/wPsi
                const Real wPsiPinv = 
                  ( wPsi < limits::Epsilon<Real>() ? Real(1) : 1/wPsi );
                // NOTE: This includes the outer wDet factor
                const Real psiMap = wPsiPinv*omega_i*Sqrt(wDet);
                if( i == firstInd )
                {
                    // Queue up an entry of D and u, and then the (scaled) 
                    // -1 and +1
                    const Real u0 = 2*(omega_i/Sqrt(wDet))*wPsi / uPsi;
                    const Real delta0 = 2*wPsiSq + 1 - u0*u0;
                    J.QueueUpdate( coneOff,         coneOff,     -wDet*delta0 );
                    J.QueueUpdate( coneOff+order,   coneOff+order,   -wDet );
                    J.QueueUpdate( coneOff+order+1, coneOff,         -wDet*u0 );
                    J.QueueUpdate( coneOff+order+1, coneOff+order+1,  wDet );
                }
                else
                {
                    // Queue up an entry of D, u, and v
                    const Real vPsi = Sqrt(uPsi*uPsi-2*wPsiSq);
                    J.QueueUpdate( n+m+iSparse,     n+m+iSparse, -wDet );
                    J.QueueUpdate( coneOff+order,   n+m+iSparse,  vPsi*psiMap );
                    J.QueueUpdate( coneOff+order+1, n+m+iSparse, -uPsi*psiMap );
                }
            }
        }
    }
    else
    {
        // Count the number of entries to queue
        // ------------------------------------
        Int numRemoteEntries = 0;
        for( Int iLoc=0; iLoc<wLocalHeight; ++iLoc )
        {
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseOrder = origToSparseOrders.GetLocal(iLoc,0);

            const Int i = w.GlobalRow(iLoc);
            if( order == sparseOrder )
            {
                numRemoteEntries += order;
            }
            else
            {
                if( i == firstInd )
                {
                    // An entry of D and a symmetric update with an entry of u,
                    // as well as the (scaled) -1 and +1 diagonal entries for 
                    // the v and u auxiliary variables
                    numRemoteEntries += 5;
                }
                else
                {
                    // An entry of D and symmetric updates with an entry of
                    // u and v
                    numRemoteEntries += 5;
                }
            }
        }

        // Queue the nonzeros
        // ------------------
        if( !J.FrozenSparsity() )
            J.Reserve( J.NumLocalEntries()+numRemoteEntries, numRemoteEntries );
        Int lastFirstInd = -1;
        vector<Real> wBuf;
        auto offs = recvOffs;
        for( Int iLoc=0; iLoc<wLocalHeight; ++iLoc )
        {
            const Int order = orders.GetLocal(iLoc,0);
            const Int firstInd = firstInds.GetLocal(iLoc,0);
            const Int sparseOrder = origToSparseOrders.GetLocal(iLoc,0);
            const Int sparseFirstInd = origToSparseFirstInds.GetLocal(iLoc,0);

            const Int i = w.GlobalRow(iLoc);
            const Int sparseOff = sparseFirstInd-firstInd;
            const Int iSparse = i+sparseOff;

            const Real omega_i = w.GetLocal(iLoc,0);
            const Real wDet = wDets.GetLocal(iLoc,0);
            if( order == sparseOrder )
            {
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
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+iSparse, +wDet-Real(2)*omega_i*omega_i );
                else
                    J.QueueUpdate
                    ( n+m+iSparse, n+m+iSparse, -wDet-Real(2)*omega_i*omega_i );

                // offdiag(-2 w w^T)
                for( Int j=firstInd; j<firstInd+order; ++j )
                    if( j != i )
                        J.QueueUpdate
                        ( n+m+iSparse, n+m+(j+sparseOff), 
                          -2*omega_i*wBuf[j-firstInd] );
            }
            else
            {
                const Int coneOff = n+m+sparseFirstInd;
                const Real wLower = wLowers.GetLocal(iLoc,0);
                const Real wPsi = wLower / Sqrt(wDet);
                const Real wPsiSq = wPsi*wPsi;
                const Real uPsi = 
                  Sqrt((4*wPsiSq*wPsiSq+4*wPsiSq+Real(1)/Real(2))/(2*wPsiSq+1));
                // Apply a pseudoinverse of sorts instead of 1/wPsi
                const Real wPsiPinv =
                  ( wPsi < limits::Epsilon<Real>() ? Real(1) : 1/wPsi );
                // NOTE: This includes the outer wDet factor
                const Real psiMap = wPsiPinv*omega_i*Sqrt(wDet);
                if( i == firstInd )
                {
                    // Queue up an entry of D, a symmetric update with u, and 
                    // then the (scaled) -1 and +1
                    const Real u0 = 2*(omega_i/Sqrt(wDet))*wPsi / uPsi;
                    const Real delta0 = 2*wPsiSq + 1 - u0*u0;
                    J.QueueUpdate( coneOff,         coneOff,     -wDet*delta0 );
                    J.QueueUpdate( coneOff+order,   coneOff+order,   -wDet );
                    J.QueueUpdate( coneOff+order+1, coneOff,         -wDet*u0 );
                    J.QueueUpdate( coneOff,         coneOff+order+1, -wDet*u0 );
                    J.QueueUpdate( coneOff+order+1, coneOff+order+1,  wDet );
                }
                else
                {
                    // Queue up an entry of D and symmetric updates with 
                    // u, and v
                    const Real vPsi = Sqrt(uPsi*uPsi-2*wPsiSq);
                    J.QueueUpdate
                    ( n+m+iSparse,     n+m+iSparse,     -wDet ); 
                    J.QueueUpdate
                    ( coneOff+order,   n+m+iSparse,      vPsi*psiMap );
                    J.QueueUpdate
                    ( n+m+iSparse,     coneOff+order,    vPsi*psiMap );
                    J.QueueUpdate
                    ( coneOff+order+1, n+m+iSparse,     -uPsi*psiMap );
                    J.QueueUpdate
                    ( n+m+iSparse,     coneOff+order+1, -uPsi*psiMap );
                }
            }
        }
    }
    J.ProcessQueues();
}

template<typename Real>
void KKTRHS
( const Matrix<Real>& rc,
  const Matrix<Real>& rb,
  const Matrix<Real>& rh,
  const Matrix<Real>& rmu,
  const Matrix<Real>& wRoot,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& d )
{
    DEBUG_ONLY(CallStackEntry cse("socp::affine::KKTRHS"))
    const Int n = rc.Height();
    const Int m = rb.Height();
    const Int k = rh.Height();
    Zeros( d, n+m+k, 1 );

    auto dx = d(IR(0,n),ALL);
    dx = rc;
    dx *= -1;

    auto dy = d(IR(n,n+m),ALL);
    dy = rb;
    dy *= -1;

    auto dz = d(IR(n+m,n+m+k),ALL);
    dz = rmu;
    soc::ApplyQuadratic( wRoot, dz, orders, firstInds );
    dz -= rh;
}

template<typename Real>
void KKTRHS
( const ElementalMatrix<Real>& rc,  
  const ElementalMatrix<Real>& rb,
  const ElementalMatrix<Real>& rh,  
  const ElementalMatrix<Real>& rmu,
  const ElementalMatrix<Real>& wRoot,
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds,
        ElementalMatrix<Real>& dPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("qp::affine::KKTRHS"))

    DistMatrixWriteProxy<Real,Real,MC,MR> dProx( dPre );
    auto& d = dProx.Get();

    const Int n = rc.Height();
    const Int m = rb.Height();
    const Int k = rh.Height();
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    Zeros( d, n+m+k, 1 );

    auto dx = d(xInd,ALL);
    dx = rc;
    dx *= -1;

    auto dy = d(yInd,ALL);
    dy = rb;
    dy *= -1;

    auto dz = d(zInd,ALL);
    dz = rmu;
    soc::ApplyQuadratic( wRoot, dz, orders, firstInds, cutoff );
    dz -= rh;
}

template<typename Real>
void KKTRHS
( const Matrix<Real>& rc, 
  const Matrix<Real>& rb, 
  const Matrix<Real>& rh, 
  const Matrix<Real>& rmu,
  const Matrix<Real>& wRoot, 
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  const Matrix<Int>& origToSparseFirstInds,
        Int kSparse,
        Matrix<Real>& d )
{
    DEBUG_ONLY(CSE cse("qp::affine::KKTRHS"))
    const Int n = rc.Height();
    const Int m = rb.Height();
    Zeros( d, n+m+kSparse, 1 );

    Matrix<Real> W_rmu;
    soc::ApplyQuadratic( wRoot, rmu, W_rmu, orders, firstInds );

    auto dx = d( IR(0,n),             ALL );
    auto dy = d( IR(n,n+m),           ALL );
    auto dz = d( IR(n+m,n+m+kSparse), ALL );
    dx = rc;
    dx *= -1;
    dy = rb;
    dy *= -1;
    for( Int i=0; i<kSparse; ++i )
    {
        const Int firstInd = firstInds.Get(i,0);
        const Int sparseFirstInd = origToSparseFirstInds.Get(i,0);
        const Int iSparse = i + (sparseFirstInd-firstInd);
        Real value = W_rmu.Get(i,0) - rh.Get(i,0);
        d.Set( n+m+iSparse, 0, value );
    }
}

template<typename Real>
void KKTRHS
( const DistMultiVec<Real>& rc, 
  const DistMultiVec<Real>& rb, 
  const DistMultiVec<Real>& rh, 
  const DistMultiVec<Real>& rmu,
  const DistMultiVec<Real>& wRoot, 
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  const DistMultiVec<Int>& origToSparseFirstInds,
        Int kSparse,
        DistMultiVec<Real>& d,
  Int cutoffPar )
{
    DEBUG_ONLY(CSE cse("qp::affine::KKTRHS"))
    const Int n = rc.Height();
    const Int m = rb.Height();
    d.SetComm( rc.Comm() );
    Zeros( d, n+m+kSparse, 1 );

    DistMultiVec<Real> W_rmu( rc.Comm() );
    soc::ApplyQuadratic( wRoot, rmu, W_rmu, orders, firstInds, cutoffPar );

    Int numEntries = rc.LocalHeight() + rb.LocalHeight() + rmu.LocalHeight();
    d.Reserve( numEntries );
    for( Int iLoc=0; iLoc<rc.LocalHeight(); ++iLoc )
    {
        Int i = rc.GlobalRow(iLoc);
        Real value = -rc.GetLocal(iLoc,0);
        d.QueueUpdate( i, 0, value );
    }
    for( Int iLoc=0; iLoc<rb.LocalHeight(); ++iLoc )
    {
        Int i = n + rb.GlobalRow(iLoc);
        Real value = -rb.GetLocal(iLoc,0);
        d.QueueUpdate( i, 0, value );
    }
    for( Int iLoc=0; iLoc<rmu.LocalHeight(); ++iLoc )
    {
        const Int firstInd = firstInds.GetLocal(iLoc,0);
        const Int sparseFirstInd = origToSparseFirstInds.GetLocal(iLoc,0);
        const Int i = rmu.GlobalRow(iLoc);
        const Int iSparse = i + (sparseFirstInd-firstInd);
        Real value = W_rmu.GetLocal(iLoc,0) - rh.GetLocal(iLoc,0);
        d.QueueUpdate( n+m+iSparse, 0, value );
    }
    d.ProcessQueues();
}

template<typename Real>
void ExpandSolution
( Int m, Int n,
  const Matrix<Real>& d, 
  const Matrix<Real>& rmu,
  const Matrix<Real>& wRoot,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& dx, 
        Matrix<Real>& dy,
        Matrix<Real>& dz, 
        Matrix<Real>& ds )
{
    DEBUG_ONLY(CSE cse("qp::affine::ExpandSolution"))
    const Int k = wRoot.Height();
    qp::affine::ExpandCoreSolution( m, n, k, d, dx, dy, dz );
    // ds := - W^T ( rmu + W dz )
    // ==========================
    ds = dz;
    soc::ApplyQuadratic( wRoot, ds, orders, firstInds );
    ds += rmu;
    soc::ApplyQuadratic( wRoot, ds, orders, firstInds );
    ds *= -1;
}

template<typename Real>
void ExpandSolution
( Int m, Int n,
  const ElementalMatrix<Real>& d, 
  const ElementalMatrix<Real>& rmu,
  const ElementalMatrix<Real>& wRoot, 
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds,
        ElementalMatrix<Real>& dx, 
        ElementalMatrix<Real>& dy,
        ElementalMatrix<Real>& dz, 
        ElementalMatrix<Real>& ds,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("qp::affine::ExpandSolution"))
    const Int k = wRoot.Height();
    qp::affine::ExpandCoreSolution( m, n, k, d, dx, dy, dz );
    // ds := - W^T ( rmu + W dz )
    // ==========================
    ds = dz;
    soc::ApplyQuadratic( wRoot, ds, orders, firstInds, cutoff );
    ds += rmu;
    soc::ApplyQuadratic( wRoot, ds, orders, firstInds, cutoff );
    ds *= -1;
}

template<typename Real>
void ExpandSolution
( Int m, Int n,
  const Matrix<Real>& d, 
  const Matrix<Real>& rmu,
  const Matrix<Real>& wRoot,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  const Matrix<Int>& sparseOrders,
  const Matrix<Int>& sparseFirstInds,
  const Matrix<Int>& sparseToOrigOrders,
  const Matrix<Int>& sparseToOrigFirstInds,
        Matrix<Real>& dx, 
        Matrix<Real>& dy,
        Matrix<Real>& dz, 
        Matrix<Real>& ds )
{
    DEBUG_ONLY(CSE cse("qp::affine::ExpandSolution"))
    const Int k = wRoot.Height();
    const Int kSparse = d.Height() - (n+m);

    // Extract dx
    // ==========
    dx = d( IR(0,n), ALL );

    // Extract dy
    // ==========
    dy = d( IR(n,n+m), ALL );

    // Extract dz
    // ==========
    auto dzSparse = d( IR(n+m,END), ALL );
    Zeros( dz, k, 1 );
    for( Int iSparse=0; iSparse<kSparse; ++iSparse )
    {
        const Int order = sparseToOrigOrders.Get(iSparse,0);
        const Int firstInd = sparseToOrigFirstInds.Get(iSparse,0);
        const Int sparseFirstInd = sparseFirstInds.Get(iSparse,0);
        const Int i = iSparse - (sparseFirstInd-firstInd);
        if( i < firstInd+order )
            dz.Set( i, 0, dzSparse.Get(iSparse,0) );
    }

    // ds := - W^T ( rmu + W dz )
    // ==========================
    ds = dz;
    soc::ApplyQuadratic( wRoot, ds, orders, firstInds );
    ds += rmu;
    soc::ApplyQuadratic( wRoot, ds, orders, firstInds );
    ds *= -1;
}

template<typename Real>
void ExpandSolution
( Int m, Int n,
  const DistMultiVec<Real>& d, 
  const DistMultiVec<Real>& rmu,
  const DistMultiVec<Real>& wRoot,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  const DistMultiVec<Int>& sparseOrders,
  const DistMultiVec<Int>& sparseFirstInds,
  const DistMultiVec<Int>& sparseToOrigOrders,
  const DistMultiVec<Int>& sparseToOrigFirstInds,
        DistMultiVec<Real>& dx, 
        DistMultiVec<Real>& dy,
        DistMultiVec<Real>& dz, 
        DistMultiVec<Real>& ds,
  Int cutoffPar )
{
    DEBUG_ONLY(CSE cse("qp::affine::ExpandSolution"))
    const Int k = wRoot.Height();

    // Extract dx
    // ==========
    dx = d( IR(0,n), ALL );

    // Extract dy
    // ==========
    dy = d( IR(n,n+m), ALL );

    // Extract dz
    // ==========
    auto dzSparse = d( IR(n+m,END), ALL );
    Zeros( dz, k, 1 );
    Int numQueues = 0;
    const Int localHeight = dzSparse.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int iSparse = dzSparse.GlobalRow(iLoc);
        const Int order = sparseToOrigOrders.GetLocal(iLoc,0);
        const Int firstInd = sparseToOrigFirstInds.GetLocal(iLoc,0);
        const Int sparseFirstInd = sparseFirstInds.GetLocal(iLoc,0);
        const Int i = iSparse - (sparseFirstInd-firstInd);
        if( i < firstInd+order )
            ++numQueues;
    }
    dz.Reserve( numQueues );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int iSparse = dzSparse.GlobalRow(iLoc);
        const Int order = sparseToOrigOrders.GetLocal(iLoc,0);
        const Int firstInd = sparseToOrigFirstInds.GetLocal(iLoc,0);
        const Int sparseFirstInd = sparseFirstInds.GetLocal(iLoc,0);
        const Int i = iSparse - (sparseFirstInd-firstInd);
        if( i < firstInd+order )
            dz.QueueUpdate( i, 0, dzSparse.GetLocal(iLoc,0) );
    }
    dz.ProcessQueues();

    // ds := - W^T ( rmu + W dz )
    // ==========================
    ds = dz;
    soc::ApplyQuadratic( wRoot, ds, orders, firstInds, cutoffPar );
    ds += rmu;
    soc::ApplyQuadratic( wRoot, ds, orders, firstInds, cutoffPar );
    ds *= -1;
}

#define PROTO(Real) \
  template void KKT \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& G, \
    const Matrix<Real>& w, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
          Matrix<Real>& J, \
    bool onlyLower ); \
  template void KKT \
  ( const ElementalMatrix<Real>& A, \
    const ElementalMatrix<Real>& G, \
    const ElementalMatrix<Real>& w, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, \
          ElementalMatrix<Real>& J, \
    bool onlyLower, Int cutoff ); \
  template void KKT \
  ( const SparseMatrix<Real>& A, \
    const SparseMatrix<Real>& G, \
    const Matrix<Real>& w, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    const Matrix<Int>& origToSparseOrders, \
    const Matrix<Int>& origToSparseFirstInds, \
          Int kSparse, \
          SparseMatrix<Real>& J, \
    bool onlyLower ); \
  template void StaticKKT \
  ( const SparseMatrix<Real>& A, \
    const SparseMatrix<Real>& G, \
          Real gamma, \
          Real delta, \
          Real beta, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    const Matrix<Int>& origToSparseOrders, \
    const Matrix<Int>& origToSparseFirstInds, \
          Int kSparse, \
          SparseMatrix<Real>& J, \
    bool onlyLower ); \
  template void FinishKKT \
  ( Int m, Int n, \
    const Matrix<Real>& w, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    const Matrix<Int>& origToSparseOrders, \
    const Matrix<Int>& origToSparseFirstInds, \
          Int kSparse, \
          SparseMatrix<Real>& J, \
    bool onlyLower ); \
  template void KKT \
  ( const DistSparseMatrix<Real>& A, \
    const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& w, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    const DistMultiVec<Int>& origToSparseOrders, \
    const DistMultiVec<Int>& origToSparseFirstInds, \
          Int kSparse, \
          DistSparseMatrix<Real>& J, \
    bool onlyLower, Int cutoffPar ); \
  template void StaticKKT \
  ( const DistSparseMatrix<Real>& A, \
    const DistSparseMatrix<Real>& G, \
          Real gamma, \
          Real delta, \
          Real beta, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    const DistMultiVec<Int>& origToSparseOrders, \
    const DistMultiVec<Int>& origToSparseFirstInds, \
          Int kSparse, \
          DistSparseMatrix<Real>& J, \
    bool onlyLower ); \
  template void FinishKKT \
  ( Int m, Int n, \
    const DistMultiVec<Real>& w, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    const DistMultiVec<Int>& origToSparseOrders, \
    const DistMultiVec<Int>& origToSparseFirstInds, \
          Int kSparse, \
          DistSparseMatrix<Real>& J, \
    bool onlyLower, Int cutoffPar ); \
  template void KKTRHS \
  ( const Matrix<Real>& rc, \
    const Matrix<Real>& rb, \
    const Matrix<Real>& rh, \
    const Matrix<Real>& rmu, \
    const Matrix<Real>& wRoot, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
          Matrix<Real>& d ); \
  template void KKTRHS \
  ( const ElementalMatrix<Real>& rc, \
    const ElementalMatrix<Real>& rb, \
    const ElementalMatrix<Real>& rh, \
    const ElementalMatrix<Real>& rmu, \
    const ElementalMatrix<Real>& wRoot, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, \
          ElementalMatrix<Real>& d, \
    Int cutoff ); \
  template void KKTRHS \
  ( const Matrix<Real>& rc, \
    const Matrix<Real>& rb, \
    const Matrix<Real>& rh, \
    const Matrix<Real>& rmu, \
    const Matrix<Real>& wRoot, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    const Matrix<Int>& origToSparseFirstInds, \
          Int kSparse, \
          Matrix<Real>& d ); \
  template void KKTRHS \
  ( const DistMultiVec<Real>& rc, \
    const DistMultiVec<Real>& rb, \
    const DistMultiVec<Real>& rh, \
    const DistMultiVec<Real>& rmu, \
    const DistMultiVec<Real>& wRoot, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    const DistMultiVec<Int>& origToSparseFirstInds, \
          Int kSparse, \
          DistMultiVec<Real>& d, \
    Int cutoffPar ); \
  template void ExpandSolution \
  ( Int m, Int n, \
    const Matrix<Real>& d, \
    const Matrix<Real>& rmu, \
    const Matrix<Real>& wRoot, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
          Matrix<Real>& dx, \
          Matrix<Real>& dy, \
          Matrix<Real>& dz, \
          Matrix<Real>& ds ); \
  template void ExpandSolution \
  ( Int m, Int n, \
    const ElementalMatrix<Real>& d, \
    const ElementalMatrix<Real>& rmu, \
    const ElementalMatrix<Real>& wRoot, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, \
          ElementalMatrix<Real>& dx, \
          ElementalMatrix<Real>& dy, \
          ElementalMatrix<Real>& dz, \
          ElementalMatrix<Real>& ds, \
          Int cutoff ); \
  template void ExpandSolution \
  ( Int m, Int n, \
    const Matrix<Real>& d, \
    const Matrix<Real>& rmu, \
    const Matrix<Real>& wRoot, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    const Matrix<Int>& sparseOrders, \
    const Matrix<Int>& sparseFirstInds, \
    const Matrix<Int>& sparseToOrigOrders, \
    const Matrix<Int>& sparseToOrigFirstInds, \
          Matrix<Real>& dx, \
          Matrix<Real>& dy, \
          Matrix<Real>& dz, \
          Matrix<Real>& ds ); \
  template void ExpandSolution \
  ( Int m, Int n, \
    const DistMultiVec<Real>& d, \
    const DistMultiVec<Real>& rmu, \
    const DistMultiVec<Real>& wRoot, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    const DistMultiVec<Int>& sparseOrders, \
    const DistMultiVec<Int>& sparseFirstInds, \
    const DistMultiVec<Int>& sparseToOrigOrders, \
    const DistMultiVec<Int>& sparseToOrigFirstInds, \
          DistMultiVec<Real>& dx, \
          DistMultiVec<Real>& dy, \
          DistMultiVec<Real>& dz, \
          DistMultiVec<Real>& ds, \
    Int cutoffPar );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace affine
} // namespace socp
} // namespace El
