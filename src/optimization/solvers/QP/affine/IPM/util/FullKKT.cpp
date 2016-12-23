/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace qp {
namespace affine {

// The full KKT system is of the form
//
//   | Q A^T      G     | | x |   |        -c             |
//   | A 0        0     | | y |   |         b             |,
//   | G 0    -(z <> s) | | z | = | -z <> (s o z + tau e) |
//
// and the particular system solved is of the form
//
//   | Q A^T      G     | | dx |   |     -rc        |
//   | A 0        0     | | dy |   |     -rb        |,
//   | G 0    -(z <> s) | | dz | = | -rh + z <> rmu |
//
// where
//
//   rc  = Q x + A^T y + G^T z + c,
//   rb  = A x - b,
//   rh  = G x + s - h,
//   rmu = s o z - tau e

template<typename Real>
void KKT
( const Matrix<Real>& Q,
  const Matrix<Real>& A,
  const Matrix<Real>& G,
  const Matrix<Real>& s,
  const Matrix<Real>& z,
        Matrix<Real>& J,
  bool onlyLower )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    Zeros( J, n+m+k, n+m+k );
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd); auto Jxz = J(xInd,zInd);
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd); auto Jyz = J(yInd,zInd);
    auto Jzx = J(zInd,xInd); auto Jzy = J(zInd,yInd); auto Jzz = J(zInd,zInd);

    // Jxx := Q
    // ========
    Jxx = Q;

    // Jyx := A
    // ========
    Jyx = A;

    // Jzx := G
    // ========
    Jzx = G;

    // Jzz := - z <> s
    // ===============
    Matrix<Real> t;
    t = s;
    DiagonalSolve( LEFT, NORMAL, z, t );
    t *= -1;
    Diagonal( Jzz, t );

    if( !onlyLower )
    {
        // Jxy := A^T
        // ==========
        Transpose( A, Jxy );

        // Jxz := G^T
        // ==========
        Transpose( G, Jxz );
    }
}

template<typename Real>
void KKT
( const ElementalMatrix<Real>& Q,
  const ElementalMatrix<Real>& A,
  const ElementalMatrix<Real>& G,
  const ElementalMatrix<Real>& s,
  const ElementalMatrix<Real>& z,
        ElementalMatrix<Real>& JPre,
  bool onlyLower )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    DistMatrixWriteProxy<Real,Real,MC,MR> JProx( JPre );
    auto& J = JProx.Get();

    Zeros( J, n+m+k, n+m+k );
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd); auto Jxz = J(xInd,zInd);
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd); auto Jyz = J(yInd,zInd);
    auto Jzx = J(zInd,xInd); auto Jzy = J(zInd,yInd); auto Jzz = J(zInd,zInd);

    // Jxx := Q
    // ========
    Jxx = Q;

    // Jyx := A
    // ========
    Jyx = A;

    // Jzx := G
    // ========
    Jzx = G;

    // Jzz := - z <> s
    // ===============
    DistMatrix<Real,MC,STAR> t(s);
    DiagonalSolve( LEFT, NORMAL, z, t );
    t *= -1;
    Diagonal( Jzz, t );

    if( !onlyLower )
    {
        // Jxy := A^T
        // ==========
        Transpose( A, Jxy );

        // Jxz := G
        // ========
        Transpose( G, Jxz );
    }
}

template<typename Real>
void KKT
( const SparseMatrix<Real>& Q,
  const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
  const Matrix<Real>& s,
  const Matrix<Real>& z,
        SparseMatrix<Real>& J,
  bool onlyLower )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    Zeros( J, n+m+k, n+m+k );
    const Int numEntriesQ = Q.NumEntries();
    const Int numEntriesA = A.NumEntries();
    const Int numEntriesG = G.NumEntries();
    // Count the number of entries of Q that we'll use
    Int numUsedEntriesQ;
    if( onlyLower )
    {
        numUsedEntriesQ = 0;
        for( Int e=0; e<numEntriesQ; ++e )
            if( Q.Row(e) >= Q.Col(e) )
                ++numUsedEntriesQ;
    }
    else
        numUsedEntriesQ = numEntriesQ;

    if( onlyLower )
        J.Reserve( numUsedEntriesQ + numEntriesA + numEntriesG + k );
    else
        J.Reserve( numUsedEntriesQ + 2*numEntriesA + 2*numEntriesG + k );

    // Jxx = Q
    // =======
    for( Int e=0; e<numEntriesQ; ++e )
    {
        const Int i = Q.Row(e);
        const Int j = Q.Col(e);
        if( i >= j || !onlyLower )
            J.QueueUpdate( i, j, Q.Value(e) );
    }

    // Jyx = A (and Jxy = A^T)
    // =======================
    for( Int e=0; e<numEntriesA; ++e )
    {
        J.QueueUpdate( n+A.Row(e), A.Col(e), A.Value(e) );
        if( !onlyLower )
            J.QueueUpdate( A.Col(e), n+A.Row(e), A.Value(e) );
    }

    // Jzx = G (and Jxz = G^T)
    // =======================
    for( Int e=0; e<numEntriesG; ++e )
    {
        J.QueueUpdate( n+m+G.Row(e), G.Col(e), G.Value(e) );
        if( !onlyLower )
            J.QueueUpdate( G.Col(e), n+m+G.Row(e), G.Value(e) );
    }

    // Jzz = -z <> s
    // =============
    for( Int e=0; e<k; ++e )
        J.QueueUpdate( n+m+e, n+m+e, -s(e)/z(e) );

    J.ProcessQueues();
}

template<typename Real>
void StaticKKT
( const SparseMatrix<Real>& Q,
  const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
        Real gamma,
        Real delta,
        Real beta,
        SparseMatrix<Real>& J,
  bool onlyLower )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    Zeros( J, n+m+k, n+m+k );
    const Int numEntriesQ = Q.NumEntries();
    const Int numEntriesA = A.NumEntries();
    const Int numEntriesG = G.NumEntries();
    // Count the number of entries of Q that we'll use
    Int numUsedEntriesQ;
    if( onlyLower )
    {
        numUsedEntriesQ = 0;
        for( Int e=0; e<numEntriesQ; ++e )
            if( Q.Row(e) >= Q.Col(e) )
                ++numUsedEntriesQ;
    }
    else
        numUsedEntriesQ = numEntriesQ;

    if( onlyLower )
        J.Reserve( numUsedEntriesQ + numEntriesA + numEntriesG + n+m+k );
    else
        J.Reserve( numUsedEntriesQ + 2*numEntriesA + 2*numEntriesG + n+m+k );

    // Jxx = Q + gamma^2*I
    // ===================
    for( Int e=0; e<numEntriesQ; ++e )
    {
        const Int i = Q.Row(e);
        const Int j = Q.Col(e);
        if( i >= j || !onlyLower )
            J.QueueUpdate( i, j, Q.Value(e) );
    }
    for( Int i=0; i<n; ++i )
        J.QueueUpdate( i, i, gamma*gamma );

    // Jyy = -delta^2*I
    // ================
    for( Int i=0; i<m; ++i )
        J.QueueUpdate( i+n, i+n, -delta*delta );

    // Jyx = A (and Jxy = A^T)
    // =======================
    for( Int e=0; e<numEntriesA; ++e )
    {
        J.QueueUpdate( n+A.Row(e), A.Col(e), A.Value(e) );
        if( !onlyLower )
            J.QueueUpdate( A.Col(e), n+A.Row(e), A.Value(e) );
    }

    // Jzx = G (and Jxz = G^T)
    // =======================
    for( Int e=0; e<numEntriesG; ++e )
    {
        J.QueueUpdate( n+m+G.Row(e), G.Col(e), G.Value(e) );
        if( !onlyLower )
            J.QueueUpdate( G.Col(e), n+m+G.Row(e), G.Value(e) );
    }

    // Jzz = -beta^2*I
    // ===============
    for( Int i=0; i<k; ++i )
        J.QueueUpdate( i+m+n, i+m+n, -beta*beta );

    J.ProcessQueues();
    J.FreezeSparsity();
}

template<typename Real>
void FinishKKT
( Int m, Int n,
  const Matrix<Real>& s,
  const Matrix<Real>& z,
        SparseMatrix<Real>& J )
{
    EL_DEBUG_CSE
    const Int k = s.Height();

    // Jzz = -z <> s
    // =============
    if( !J.FrozenSparsity() )
        J.Reserve( J.NumEntries()+k );
    for( Int e=0; e<k; ++e )
        J.QueueUpdate( n+m+e, n+m+e, -s(e)/z(e) );
    J.ProcessQueues();
}

template<typename Real>
void KKT
( const DistSparseMatrix<Real>& Q,
  const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& s,
  const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J,
  bool onlyLower )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();
    const Int numEntriesQ = Q.NumLocalEntries();
    const Int numEntriesA = A.NumLocalEntries();
    const Int numEntriesG = G.NumLocalEntries();
    auto& sLoc = s.LockedMatrix();
    auto& zLoc = z.LockedMatrix();

    J.SetGrid( A.Grid() );
    Zeros( J, n+m+k, n+m+k );

    // Compute the number of entries to send
    // =====================================
    Int numEntries = numEntriesA + numEntriesG + s.LocalHeight();
    if( !onlyLower )
        numEntries += numEntriesA + numEntriesG;
    for( Int e=0; e<numEntriesQ; ++e )
    {
        const Int i = Q.Row(e);
        const Int j = Q.Col(e);
        if( i >= j || !onlyLower )
            ++numEntries;
    }

    // Queue and process the entries
    // =============================
    J.Reserve( numEntries, numEntries );
    // Pack Q
    // ------
    for( Int e=0; e<numEntriesQ; ++e )
    {
        const Int i = Q.Row(e);
        const Int j = Q.Col(e);
        if( i >= j || !onlyLower )
            J.QueueUpdate( i, j, Q.Value(e) );
    }
    // Pack A
    // ------
    for( Int e=0; e<numEntriesA; ++e )
    {
        const Int i = A.Row(e) + n;
        const Int j = A.Col(e);
        J.QueueUpdate( i, j, A.Value(e) );
        if( !onlyLower )
            J.QueueUpdate( j, i, A.Value(e) );
    }
    // Pack G
    // ------
    for( Int e=0; e<numEntriesG; ++e )
    {
        const Int i = G.Row(e) + n + m;
        const Int j = G.Col(e);
        J.QueueUpdate( i, j, G.Value(e) );
        if( !onlyLower )
            J.QueueUpdate( j, i, G.Value(e) );
    }
    // Pack -z <> s
    // ------------
    for( Int iLoc=0; iLoc<s.LocalHeight(); ++iLoc )
    {
        const Int i = m+n + s.GlobalRow(iLoc);
        const Real value = -sLoc(iLoc)/zLoc(iLoc);
        J.QueueUpdate( i, i, value );
    }
    J.ProcessQueues();
}

template<typename Real>
void StaticKKT
( const DistSparseMatrix<Real>& Q,
  const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
        Real gamma,
        Real delta,
        Real beta,
        DistSparseMatrix<Real>& J,
  bool onlyLower )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();
    const Int numEntriesQ = Q.NumLocalEntries();
    const Int numEntriesA = A.NumLocalEntries();
    const Int numEntriesG = G.NumLocalEntries();

    J.SetGrid( A.Grid() );
    Zeros( J, n+m+k, n+m+k );
    const Int localHeightJ = J.LocalHeight();

    // Compute the number of entries to send
    // =====================================
    Int numEntries = numEntriesA + numEntriesG;
    if( !onlyLower )
        numEntries += numEntriesA + numEntriesG;
    for( Int e=0; e<numEntriesQ; ++e )
    {
        const Int i = Q.Row(e);
        const Int j = Q.Col(e);
        if( i >= j || !onlyLower )
            ++numEntries;
    }

    // Queue and process the entries
    // =============================
    J.Reserve( numEntries+localHeightJ, numEntries );
    // Pack Q
    // ------
    for( Int e=0; e<numEntriesQ; ++e )
    {
        const Int i = Q.Row(e);
        const Int j = Q.Col(e);
        if( i >= j || !onlyLower )
            J.QueueUpdate( i, j, Q.Value(e) );
    }
    // Pack A
    // ------
    for( Int e=0; e<numEntriesA; ++e )
    {
        const Int i = A.Row(e) + n;
        const Int j = A.Col(e);
        J.QueueUpdate( i, j, A.Value(e) );
        if( !onlyLower )
            J.QueueUpdate( j, i, A.Value(e) );
    }
    // Pack G
    // ------
    for( Int e=0; e<numEntriesG; ++e )
    {
        const Int i = G.Row(e) + n + m;
        const Int j = G.Col(e);
        J.QueueUpdate( i, j, G.Value(e) );
        if( !onlyLower )
            J.QueueUpdate( j, i, G.Value(e) );
    }
    // Add the regularization
    // ----------------------
    for( Int iLoc=0; iLoc<localHeightJ; ++iLoc )
    {
        const Int i = J.GlobalRow(iLoc);
        if( i < n )        J.QueueLocalUpdate( iLoc, i,  gamma*gamma );
        else if( i < n+m ) J.QueueLocalUpdate( iLoc, i, -delta*delta );
        else               J.QueueLocalUpdate( iLoc, i, -beta*beta );
    }
    J.ProcessQueues();
    J.FreezeSparsity();
}

template<typename Real>
void FinishKKT
( Int m, Int n,
  const DistMultiVec<Real>& s,
  const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J )
{
    EL_DEBUG_CSE
    auto& sLoc = s.LockedMatrix();
    auto& zLoc = z.LockedMatrix();

    // Pack -z <> s
    // ------------
    const Int numEntries = s.LocalHeight();
    if( !J.FrozenSparsity() )
        J.Reserve( numEntries, numEntries );
    for( Int iLoc=0; iLoc<s.LocalHeight(); ++iLoc )
    {
        const Int i = m+n + s.GlobalRow(iLoc);
        const Real value = -sLoc(iLoc)/zLoc(iLoc);
        J.QueueUpdate( i, i, value );
    }
    J.ProcessQueues();
}

template<typename Real>
void KKTRHS
( const Matrix<Real>& rc,
  const Matrix<Real>& rb,
  const Matrix<Real>& rh,
  const Matrix<Real>& rmu,
  const Matrix<Real>& z,
        Matrix<Real>& d )
{
    EL_DEBUG_CSE
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
    DiagonalSolve( LEFT, NORMAL, z, dz );
    dz -= rh;
}

template<typename Real>
void KKTRHS
( const ElementalMatrix<Real>& rc,
  const ElementalMatrix<Real>& rb,
  const ElementalMatrix<Real>& rh,
  const ElementalMatrix<Real>& rmu,
  const ElementalMatrix<Real>& z,
        ElementalMatrix<Real>& dPre )
{
    EL_DEBUG_CSE

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
    DiagonalSolve( LEFT, NORMAL, z, dz );
    dz -= rh;
}

template<typename Real>
void KKTRHS
( const DistMultiVec<Real>& rc,
  const DistMultiVec<Real>& rb,
  const DistMultiVec<Real>& rh,
  const DistMultiVec<Real>& rmu,
  const DistMultiVec<Real>& z,
        DistMultiVec<Real>& d )
{
    EL_DEBUG_CSE
    const Int n = rc.Height();
    const Int m = rb.Height();
    const int k = rh.Height();
    auto& rcLoc = rc.LockedMatrix();
    auto& rbLoc = rb.LockedMatrix();
    auto& rhLoc = rh.LockedMatrix();
    auto& rmuLoc = rmu.LockedMatrix();
    auto& zLoc = z.LockedMatrix();

    d.SetGrid( rc.Grid() );
    Zeros( d, n+m+k, 1 );

    Int numEntries = rc.LocalHeight() + rb.LocalHeight() + rmu.LocalHeight();
    d.Reserve( numEntries );
    for( Int iLoc=0; iLoc<rc.LocalHeight(); ++iLoc )
    {
        Int i = rc.GlobalRow(iLoc);
        Real value = -rcLoc(iLoc);
        d.QueueUpdate( i, 0, value );
    }
    for( Int iLoc=0; iLoc<rb.LocalHeight(); ++iLoc )
    {
        Int i = n + rb.GlobalRow(iLoc);
        Real value = -rbLoc(iLoc);
        d.QueueUpdate( i, 0, value );
    }
    for( Int iLoc=0; iLoc<rmu.LocalHeight(); ++iLoc )
    {
        Int i = n+m + rmu.GlobalRow(iLoc);
        Real value = rmuLoc(iLoc)/zLoc(iLoc) - rhLoc(iLoc);
        d.QueueUpdate( i, 0, value );
    }
    d.ProcessQueues();
}

template<typename Real>
void ExpandCoreSolution
( Int m, Int n, Int k,
  const Matrix<Real>& d,
        Matrix<Real>& dx,
        Matrix<Real>& dy,
        Matrix<Real>& dz )
{
    EL_DEBUG_CSE
    if( d.Height() != n+m+k || d.Width() != 1 )
        LogicError("Right-hand side was the wrong size");
    dx = d(IR(0,  n    ),ALL);
    dy = d(IR(n,  n+m  ),ALL);
    dz = d(IR(n+m,n+m+k),ALL);
}

template<typename Real>
void ExpandCoreSolution
( Int m, Int n, Int k,
  const ElementalMatrix<Real>& dPre,
        ElementalMatrix<Real>& dx,
        ElementalMatrix<Real>& dy,
        ElementalMatrix<Real>& dz )
{
    EL_DEBUG_CSE

    DistMatrixReadProxy<Real,Real,MC,MR> dProx( dPre );
    auto& d = dProx.GetLocked();

    if( d.Height() != n+m+k || d.Width() != 1 )
        LogicError("Right-hand side was the wrong size");

    dx = d(IR(0,  n    ),ALL);
    dy = d(IR(n,  n+m  ),ALL);
    dz = d(IR(n+m,n+m+k),ALL);
}

template<typename Real>
void ExpandCoreSolution
( Int m, Int n, Int k,
  const DistMultiVec<Real>& d,
        DistMultiVec<Real>& dx,
        DistMultiVec<Real>& dy,
        DistMultiVec<Real>& dz )
{
    EL_DEBUG_CSE
    if( d.Height() != n+m+k || d.Width() != 1 )
        LogicError("Right-hand side was the wrong size");

    const Grid& grid = d.Grid();
    dx.SetGrid( grid );
    dy.SetGrid( grid );
    dz.SetGrid( grid );
    Zeros( dx, n, 1 );
    Zeros( dy, m, 1 );
    Zeros( dz, k, 1 );

    // Count the number of entries to send to each piece
    // =================================================
    Int dxNumEntries=0, dyNumEntries=0, dzNumEntries=0;
    for( Int iLoc=0; iLoc<d.LocalHeight(); ++iLoc )
    {
        const Int i = d.GlobalRow(iLoc);
        if( i < n )
            ++dxNumEntries;
        else if( i < n+m )
            ++dyNumEntries;
        else
            ++dzNumEntries;
    }

    // Pack and process the entries
    // ============================
    dx.Reserve( dxNumEntries );
    dy.Reserve( dyNumEntries );
    dz.Reserve( dzNumEntries );
    auto& dLoc = d.LockedMatrix();
    for( Int iLoc=0; iLoc<d.LocalHeight(); ++iLoc )
    {
        const Int i = d.GlobalRow(iLoc);
        const Real value = dLoc(iLoc);
        if( i < n )
            dx.QueueUpdate( i, 0, value );
        else if( i < n+m )
            dy.QueueUpdate( i-n, 0, value );
        else
            dz.QueueUpdate( i-(n+m), 0, value );
    }
    dx.ProcessQueues();
    dy.ProcessQueues();
    dz.ProcessQueues();
}

template<typename Real>
void ExpandSolution
( Int m, Int n,
  const Matrix<Real>& d,
  const Matrix<Real>& rmu,
  const Matrix<Real>& s,
  const Matrix<Real>& z,
        Matrix<Real>& dx,
        Matrix<Real>& dy,
        Matrix<Real>& dz,
        Matrix<Real>& ds )
{
    EL_DEBUG_CSE
    const Int k = s.Height();
    ExpandCoreSolution( m, n, k, d, dx, dy, dz );
    // ds := - z <> ( rmu + s o dz )
    // =============================
    ds = dz;
    DiagonalScale( LEFT, NORMAL, s, ds );
    ds += rmu;
    DiagonalSolve( LEFT, NORMAL, z, ds );
    ds *= -1;
}

template<typename Real>
void ExpandSolution
( Int m, Int n,
  const ElementalMatrix<Real>& d,
  const ElementalMatrix<Real>& rmu,
  const ElementalMatrix<Real>& s,
  const ElementalMatrix<Real>& z,
        ElementalMatrix<Real>& dx,
        ElementalMatrix<Real>& dy,
        ElementalMatrix<Real>& dz,
        ElementalMatrix<Real>& ds )
{
    EL_DEBUG_CSE
    const int k = s.Height();
    ExpandCoreSolution( m, n, k, d, dx, dy, dz );
    // ds := - z <> ( rmu + s o dz )
    // =============================
    ds = dz;
    DiagonalScale( LEFT, NORMAL, s, ds );
    ds += rmu;
    DiagonalSolve( LEFT, NORMAL, z, ds );
    ds *= -1;
}

template<typename Real>
void ExpandSolution
( Int m, Int n,
  const DistMultiVec<Real>& d,
  const DistMultiVec<Real>& rmu,
  const DistMultiVec<Real>& s,
  const DistMultiVec<Real>& z,
        DistMultiVec<Real>& dx,
        DistMultiVec<Real>& dy,
        DistMultiVec<Real>& dz,
        DistMultiVec<Real>& ds )
{
    EL_DEBUG_CSE
    const Int k = s.Height();
    ExpandCoreSolution( m, n, k, d, dx, dy, dz );
    // ds := - z <> ( rmu + s o dz )
    // =============================
    ds = dz;
    DiagonalScale( LEFT, NORMAL, s, ds );
    ds += rmu;
    DiagonalSolve( LEFT, NORMAL, z, ds );
    ds *= -1;
}

#define PROTO(Real) \
  template void KKT \
  ( const Matrix<Real>& Q, \
    const Matrix<Real>& A, \
    const Matrix<Real>& G, \
    const Matrix<Real>& s, \
    const Matrix<Real>& z, \
          Matrix<Real>& J, \
    bool onlyLower ); \
  template void KKT \
  ( const ElementalMatrix<Real>& Q, \
    const ElementalMatrix<Real>& A, \
    const ElementalMatrix<Real>& G, \
    const ElementalMatrix<Real>& s, \
    const ElementalMatrix<Real>& z, \
          ElementalMatrix<Real>& J, \
    bool onlyLower ); \
  template void KKT \
  ( const SparseMatrix<Real>& Q, \
    const SparseMatrix<Real>& A, \
    const SparseMatrix<Real>& G, \
    const Matrix<Real>& s, \
    const Matrix<Real>& z, \
          SparseMatrix<Real>& J, \
    bool onlyLower ); \
  template void StaticKKT \
  ( const SparseMatrix<Real>& Q, \
    const SparseMatrix<Real>& A, \
    const SparseMatrix<Real>& G, \
          Real gamma, \
          Real delta, \
          Real beta, \
          SparseMatrix<Real>& J, \
    bool onlyLower ); \
  template void FinishKKT \
  ( Int m, Int n, \
    const Matrix<Real>& s, \
    const Matrix<Real>& z, \
          SparseMatrix<Real>& J ); \
  template void KKT \
  ( const DistSparseMatrix<Real>& Q, \
    const DistSparseMatrix<Real>& A, \
    const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& s, \
    const DistMultiVec<Real>& z, \
          DistSparseMatrix<Real>& J, \
    bool onlyLower ); \
  template void StaticKKT \
  ( const DistSparseMatrix<Real>& Q, \
    const DistSparseMatrix<Real>& A, \
    const DistSparseMatrix<Real>& G, \
          Real gamma, \
          Real delta, \
          Real beta, \
          DistSparseMatrix<Real>& J, \
    bool onlyLower ); \
  template void FinishKKT \
  ( Int m, Int n, \
    const DistMultiVec<Real>& s, \
    const DistMultiVec<Real>& z, \
          DistSparseMatrix<Real>& J ); \
  template void KKTRHS \
  ( const Matrix<Real>& rc, \
    const Matrix<Real>& rb, \
    const Matrix<Real>& rh, \
    const Matrix<Real>& rmu, \
    const Matrix<Real>& z, \
          Matrix<Real>& d ); \
  template void KKTRHS \
  ( const ElementalMatrix<Real>& rc, \
    const ElementalMatrix<Real>& rb, \
    const ElementalMatrix<Real>& rh, \
    const ElementalMatrix<Real>& rmu, \
    const ElementalMatrix<Real>& z, \
          ElementalMatrix<Real>& d ); \
  template void KKTRHS \
  ( const DistMultiVec<Real>& rc, \
    const DistMultiVec<Real>& rb, \
    const DistMultiVec<Real>& rh, \
    const DistMultiVec<Real>& rmu, \
    const DistMultiVec<Real>& z, \
          DistMultiVec<Real>& d ); \
  template void ExpandCoreSolution \
  ( Int m, Int n, Int k, \
    const Matrix<Real>& d, \
          Matrix<Real>& dx, \
          Matrix<Real>& dy, \
          Matrix<Real>& dz ); \
  template void ExpandCoreSolution \
  ( Int m, Int n, Int k, \
    const ElementalMatrix<Real>& d, \
          ElementalMatrix<Real>& dx, \
          ElementalMatrix<Real>& dy, \
          ElementalMatrix<Real>& dz ); \
  template void ExpandCoreSolution \
  ( Int m, Int n, Int k, \
    const DistMultiVec<Real>& d, \
          DistMultiVec<Real>& dx, \
          DistMultiVec<Real>& dy, \
          DistMultiVec<Real>& dz ); \
  template void ExpandSolution \
  ( Int m, Int n, \
    const Matrix<Real>& d, \
    const Matrix<Real>& rmu, \
    const Matrix<Real>& s, \
    const Matrix<Real>& z, \
          Matrix<Real>& dx, \
          Matrix<Real>& dy, \
          Matrix<Real>& dz, \
          Matrix<Real>& ds ); \
  template void ExpandSolution \
  ( Int m, Int n, \
    const ElementalMatrix<Real>& d, \
    const ElementalMatrix<Real>& rmu, \
    const ElementalMatrix<Real>& s, \
    const ElementalMatrix<Real>& z, \
          ElementalMatrix<Real>& dx, \
          ElementalMatrix<Real>& dy, \
          ElementalMatrix<Real>& dz, \
          ElementalMatrix<Real>& ds ); \
  template void ExpandSolution \
  ( Int m, Int n, \
    const DistMultiVec<Real>& d, \
    const DistMultiVec<Real>& rmu, \
    const DistMultiVec<Real>& s, \
    const DistMultiVec<Real>& z, \
          DistMultiVec<Real>& dx, \
          DistMultiVec<Real>& dy, \
          DistMultiVec<Real>& dz, \
          DistMultiVec<Real>& ds );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace affine
} // namespace qp
} // namespace El
