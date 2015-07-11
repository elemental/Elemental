/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Find the Nesterov-Todd scaling point w such that 
//
//   Q_w z = s,
//
// or, equivalently, for each subcone,
//
//   Hess(-(1/2) ln det(w)) s = z,
//
// where -(1/2) ln det(w) is the barrier function for the Second-Order Cone.

namespace {

template<typename Real>
void ClassicalNT
( const Matrix<Real>& s, 
  const Matrix<Real>& z,
        Matrix<Real>& w,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("ClassicalNT"))
    typedef Promote<Real> PReal;

    Matrix<PReal> sProm, zProm;
    Copy( s, sProm );
    Copy( z, zProm );

    Matrix<PReal> sRoot;
    SOCSquareRoot( sProm, sRoot, orders, firstInds );

    // a := Q_{sqrt(s)}(z)
    // -------------------
    Matrix<PReal> a;
    SOCApplyQuadratic( sRoot, zProm, a, orders, firstInds );

    // a := inv(sqrt(a)) = inv(sqrt((Q_{sqrt(s)}(z))))
    // -----------------------------------------------
    Matrix<PReal> b; 
    SOCSquareRoot( a, b, orders, firstInds );
    SOCInverse( b, a, orders, firstInds );

    // w := Q_{sqrt(s)}(a)
    // -------------------
    Matrix<PReal> wProm;
    SOCApplyQuadratic( sRoot, a, wProm, orders, firstInds );
    Copy( wProm, w );
}

template<typename Real>
void ClassicalNT
( const AbstractDistMatrix<Real>& sPre, 
  const AbstractDistMatrix<Real>& zPre,
        AbstractDistMatrix<Real>& wPre,
  const AbstractDistMatrix<Int>& ordersPre, 
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("ClassicalNT"))
    typedef Promote<Real> PReal;
    AssertSameGrids( sPre, zPre, wPre, ordersPre, firstIndsPre );

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto sPtr = ReadProxy<PReal,VC,STAR>(&sPre,ctrl); 
    auto zPtr = ReadProxy<PReal,VC,STAR>(&zPre,ctrl);
    auto wPtr = WriteProxy<PReal,VC,STAR>(&wPre,ctrl);
    auto ordersPtr = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& s = *sPtr;
    auto& z = *zPtr;
    auto& w = *wPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    DistMatrix<PReal,VC,STAR> sRoot(s.Grid());
    SOCSquareRoot( s, sRoot, orders, firstInds, cutoff );

    // a := Q_{sqrt(s)}(z)
    // -------------------
    DistMatrix<PReal,VC,STAR> a(z.Grid());
    SOCApplyQuadratic( sRoot, z, a, orders, firstInds, cutoff );

    // a := inv(sqrt(a)) = inv(sqrt((Q_{sqrt(s)}(z))))
    // -----------------------------------------------
    DistMatrix<PReal,VC,STAR> b(z.Grid()); 
    SOCSquareRoot( a, b, orders, firstInds, cutoff );
    SOCInverse( b, a, orders, firstInds, cutoff );

    // w := Q_{sqrt(s)}(a)
    // -------------------
    SOCApplyQuadratic( sRoot, a, w, orders, firstInds, cutoff );
}

template<typename Real>
void ClassicalNT
( const DistMultiVec<Real>& s, 
  const DistMultiVec<Real>& z,
        DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("ClassicalNT"))
    typedef Promote<Real> PReal;
    mpi::Comm comm = s.Comm();

    DistMultiVec<PReal> sProm(comm), zProm(comm);
    Copy( s, sProm );
    Copy( z, zProm );

    DistMultiVec<PReal> sRoot(comm);
    SOCSquareRoot( sProm, sRoot, orders, firstInds, cutoff );

    // a := Q_{sqrt(s)}(z)
    // -------------------
    DistMultiVec<PReal> a(comm);
    SOCApplyQuadratic( sRoot, zProm, a, orders, firstInds, cutoff );

    // a := inv(sqrt(a)) = inv(sqrt((Q_{sqrt(s)}(z))))
    // -----------------------------------------------
    DistMultiVec<PReal> b(comm); 
    SOCSquareRoot( a, b, orders, firstInds, cutoff );
    SOCInverse( b, a, orders, firstInds, cutoff );

    // w := Q_{sqrt(s)}(a)
    // -------------------
    DistMultiVec<PReal> wProm(comm);
    SOCApplyQuadratic( sRoot, a, wProm, orders, firstInds, cutoff );
    Copy( wProm, w );
}

// See Section 4.2 of 
// http://www.seas.ucla.edu/~vandenbe/publications/coneprog.pdf

template<typename Real>
void VandenbergheNT
( const Matrix<Real>& s, 
  const Matrix<Real>& z,
        Matrix<Real>& w,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("VandenbergheNT"))
    typedef Promote<Real> PReal;
    const Int n = s.Height();

    // Promote the floating-point precision
    // ====================================
    Matrix<PReal> sProm, zProm;
    Copy( s, sProm );
    Copy( z, zProm );

    // Normalize with respect to the Jordan determinant
    // ================================================
    Matrix<PReal> sDets, zDets;
    SOCDets( sProm, sDets, orders, firstInds );
    SOCDets( zProm, zDets, orders, firstInds );
    ConeBroadcast( sDets, orders, firstInds );
    ConeBroadcast( zDets, orders, firstInds );
    for( Int i=0; i<n; ++i )
    {
        sProm.Set( i, 0, sProm.Get(i,0)/Sqrt(sDets.Get(i,0)) );
        zProm.Set( i, 0, zProm.Get(i,0)/Sqrt(zDets.Get(i,0)) );
    }

    // Compute the 'gamma' coefficients
    // ================================
    Matrix<PReal> gammas;
    SOCDots( zProm, sProm, gammas, orders, firstInds );
    ConeBroadcast( gammas, orders, firstInds );
    for( Int i=0; i<n; ++i )
        gammas.Set( i, 0, Sqrt((PReal(1)+gammas.Get(i,0))/PReal(2)) );

    // Compute the normalized scaling point
    // ====================================
    auto wProm = zProm;
    SOCReflect( wProm, orders, firstInds );
    wProm += sProm;
    DiagonalSolve( LEFT, NORMAL, gammas, wProm );
    wProm *= PReal(1)/PReal(2);

    // Rescale the scaling point
    // =========================
    for( Int i=0; i<n; ++i )
    {
        const PReal sDet = sDets.Get(i,0);
        const PReal zDet = zDets.Get(i,0);
        const PReal scale = Pow(sDet,PReal(0.25))/Pow(zDet,PReal(0.25));
        wProm.Set( i, 0, wProm.Get(i,0)*scale );
    }
    Copy( wProm, w );
}

template<typename Real>
void VandenbergheNT
( const AbstractDistMatrix<Real>& sPre, 
  const AbstractDistMatrix<Real>& zPre,
        AbstractDistMatrix<Real>& wPre,
  const AbstractDistMatrix<Int>& ordersPre, 
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("VandenbergheNT"))
    typedef Promote<Real> PReal;
    const Grid& grid = sPre.Grid();
    AssertSameGrids( sPre, zPre, wPre, ordersPre, firstIndsPre );

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrix<PReal,VC,STAR> s(grid), z(grid);
    Copy( sPre, s );
    Copy( zPre, z );

    auto wPtr = WriteProxy<PReal,VC,STAR>(&wPre,ctrl);
    auto ordersPtr = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& w = *wPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    // Normalize with respect to the Jordan determinant
    // ================================================
    const Int nLocal = s.LocalHeight();
    DistMatrix<PReal,VC,STAR> sDets(grid), zDets(grid);
    SOCDets( s, sDets, orders, firstInds, cutoff );
    SOCDets( z, zDets, orders, firstInds, cutoff );
    ConeBroadcast( sDets, orders, firstInds, cutoff );
    ConeBroadcast( zDets, orders, firstInds, cutoff );
    for( Int iLoc=0; iLoc<nLocal; ++iLoc )
    {
        s.SetLocal
        ( iLoc, 0, s.GetLocal(iLoc,0)/Sqrt(sDets.GetLocal(iLoc,0)) );
        z.SetLocal
        ( iLoc, 0, z.GetLocal(iLoc,0)/Sqrt(zDets.GetLocal(iLoc,0)) );
    }

    // Compute the 'gamma' coefficients
    // ================================
    DistMatrix<PReal,VC,STAR> gammas(grid);
    SOCDots( z, s, gammas, orders, firstInds, cutoff );
    ConeBroadcast( gammas, orders, firstInds, cutoff );
    for( Int iLoc=0; iLoc<nLocal; ++iLoc )
        gammas.SetLocal
        ( iLoc, 0, Sqrt((PReal(1)+gammas.GetLocal(iLoc,0))/PReal(2)) );

    // Compute the normalized scaling point
    // ====================================
    w = z;
    SOCReflect( w, orders, firstInds );
    w += s;
    DiagonalSolve( LEFT, NORMAL, gammas, w );
    w *= PReal(1)/PReal(2);

    // Rescale the scaling point
    // =========================
    for( Int iLoc=0; iLoc<nLocal; ++iLoc )
    {
        const PReal sDet = sDets.GetLocal(iLoc,0);
        const PReal zDet = zDets.GetLocal(iLoc,0);
        const PReal scale = Pow(sDet,PReal(0.25))/Pow(zDet,PReal(0.25));
        w.SetLocal( iLoc, 0, w.GetLocal(iLoc,0)*scale );
    }
}

template<typename Real>
void VandenbergheNT
( const DistMultiVec<Real>& s, 
  const DistMultiVec<Real>& z,
        DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("VandenbergheNT"))
    typedef Promote<Real> PReal;
    mpi::Comm comm = s.Comm();

    DistMultiVec<PReal> sProm(comm), zProm(comm);
    Copy( s, sProm );
    Copy( z, zProm );

    // Normalize with respect to the Jordan determinant
    // ================================================
    const Int nLocal = sProm.LocalHeight();
    DistMultiVec<PReal> sDets(comm), zDets(comm);
    SOCDets( sProm, sDets, orders, firstInds, cutoff );
    SOCDets( zProm, zDets, orders, firstInds, cutoff );
    ConeBroadcast( sDets, orders, firstInds, cutoff );
    ConeBroadcast( zDets, orders, firstInds, cutoff );
    for( Int iLoc=0; iLoc<nLocal; ++iLoc )
    {
        sProm.SetLocal
        ( iLoc, 0, sProm.GetLocal(iLoc,0)/Sqrt(sDets.GetLocal(iLoc,0)) );
        zProm.SetLocal
        ( iLoc, 0, zProm.GetLocal(iLoc,0)/Sqrt(zDets.GetLocal(iLoc,0)) );
    }

    // Compute the 'gamma' coefficients
    // ================================
    DistMultiVec<PReal> gammas(comm);
    SOCDots( zProm, sProm, gammas, orders, firstInds, cutoff );
    ConeBroadcast( gammas, orders, firstInds, cutoff );
    for( Int iLoc=0; iLoc<nLocal; ++iLoc )
        gammas.SetLocal
        ( iLoc, 0, Sqrt((PReal(1)+gammas.GetLocal(iLoc,0))/PReal(2)) );

    // Compute the normalized scaling point
    // ====================================
    auto wProm = zProm;
    SOCReflect( wProm, orders, firstInds );
    wProm += sProm;
    DiagonalSolve( LEFT, NORMAL, gammas, wProm );
    wProm *= PReal(1)/PReal(2);

    // Rescale the scaling point
    // =========================
    for( Int iLoc=0; iLoc<nLocal; ++iLoc )
    {
        const PReal sDet = sDets.GetLocal(iLoc,0);
        const PReal zDet = zDets.GetLocal(iLoc,0);
        const PReal scale = Pow(sDet,PReal(0.25))/Pow(zDet,PReal(0.25));
        wProm.SetLocal( iLoc, 0, wProm.GetLocal(iLoc,0)*scale );
    }
    Copy( wProm, w );
}

} // anonymous namespace

template<typename Real>
void SOCNesterovTodd
( const Matrix<Real>& s, 
  const Matrix<Real>& z,
        Matrix<Real>& w,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCNesterovTodd"))
    const bool useClassical = false;
    if( useClassical )
        ClassicalNT( s, z, w, orders, firstInds );
    else
        VandenbergheNT( s, z, w, orders, firstInds );
}

template<typename Real>
void SOCNesterovTodd
( const AbstractDistMatrix<Real>& s, 
  const AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& w,
  const AbstractDistMatrix<Int>& orders, 
  const AbstractDistMatrix<Int>& firstInds,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCNesterovTodd"))
    const bool useClassical = false;
    if( useClassical )
        ClassicalNT( s, z, w, orders, firstInds, cutoff );
    else
        VandenbergheNT( s, z, w, orders, firstInds, cutoff );
}

template<typename Real>
void SOCNesterovTodd
( const DistMultiVec<Real>& s, 
  const DistMultiVec<Real>& z,
        DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCNesterovTodd"))
    const bool useClassical = false;
    if( useClassical )
        ClassicalNT( s, z, w, orders, firstInds, cutoff );
    else
        VandenbergheNT( s, z, w, orders, firstInds, cutoff );
}

#define PROTO(Real) \
  template void SOCNesterovTodd \
  ( const Matrix<Real>& s, \
    const Matrix<Real>& z, \
          Matrix<Real>& w, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void SOCNesterovTodd \
  ( const AbstractDistMatrix<Real>& s, \
    const AbstractDistMatrix<Real>& z, \
          AbstractDistMatrix<Real>& w, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, Int cutoff ); \
  template void SOCNesterovTodd \
  ( const DistMultiVec<Real>& s, \
    const DistMultiVec<Real>& z, \
          DistMultiVec<Real>& w, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
