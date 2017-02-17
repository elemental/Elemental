/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace soc {

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

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
void ClassicalNT
( const Matrix<Real>& s,
  const Matrix<Real>& z,
        Matrix<Real>& w,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds )
{
    EL_DEBUG_CSE
    typedef Promote<Real> PReal;

    Matrix<PReal> sProm, zProm;
    Copy( s, sProm );
    Copy( z, zProm );

    Matrix<PReal> sRoot;
    soc::SquareRoot( sProm, sRoot, orders, firstInds );

    // a := Q_{sqrt(s)}(z)
    // -------------------
    Matrix<PReal> a;
    soc::ApplyQuadratic( sRoot, zProm, a, orders, firstInds );

    // a := inv(sqrt(a)) = inv(sqrt((Q_{sqrt(s)}(z))))
    // -----------------------------------------------
    Matrix<PReal> b;
    soc::SquareRoot( a, b, orders, firstInds );
    soc::Inverse( b, a, orders, firstInds );

    // w := Q_{sqrt(s)}(a)
    // -------------------
    Matrix<PReal> wProm;
    soc::ApplyQuadratic( sRoot, a, wProm, orders, firstInds );
    Copy( wProm, w );
}

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
void ClassicalNT
( const AbstractDistMatrix<Real>& sPre,
  const AbstractDistMatrix<Real>& zPre,
        AbstractDistMatrix<Real>& wPre,
  const AbstractDistMatrix<Int>& ordersPre,
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    EL_DEBUG_CSE
    typedef Promote<Real> PReal;
    AssertSameGrids( sPre, zPre, wPre, ordersPre, firstIndsPre );

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrixReadProxy<Real,PReal,VC,STAR>
      sProx( sPre, ctrl ),
      zProx( zPre, ctrl );
    DistMatrixWriteProxy<Real,PReal,VC,STAR>
      wProx( wPre, ctrl );
    DistMatrixReadProxy<Int,Int,VC,STAR>
      ordersProx( ordersPre, ctrl ),
      firstIndsProx( firstIndsPre, ctrl );
    auto& s = sProx.GetLocked();
    auto& z = zProx.GetLocked();
    auto& w = wProx.Get();
    auto& orders = ordersProx.GetLocked();
    auto& firstInds = firstIndsProx.GetLocked();

    DistMatrix<PReal,VC,STAR> sRoot(s.Grid());
    soc::SquareRoot( s, sRoot, orders, firstInds, cutoff );

    // a := Q_{sqrt(s)}(z)
    // -------------------
    DistMatrix<PReal,VC,STAR> a(z.Grid());
    soc::ApplyQuadratic( sRoot, z, a, orders, firstInds, cutoff );

    // a := inv(sqrt(a)) = inv(sqrt((Q_{sqrt(s)}(z))))
    // -----------------------------------------------
    DistMatrix<PReal,VC,STAR> b(z.Grid());
    soc::SquareRoot( a, b, orders, firstInds, cutoff );
    soc::Inverse( b, a, orders, firstInds, cutoff );

    // w := Q_{sqrt(s)}(a)
    // -------------------
    soc::ApplyQuadratic( sRoot, a, w, orders, firstInds, cutoff );
}

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
void ClassicalNT
( const DistMultiVec<Real>& s,
  const DistMultiVec<Real>& z,
        DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff )
{
    EL_DEBUG_CSE
    typedef Promote<Real> PReal;
    const Grid& grid = s.Grid();

    DistMultiVec<PReal> sProm(grid), zProm(grid);
    Copy( s, sProm );
    Copy( z, zProm );

    DistMultiVec<PReal> sRoot(grid);
    soc::SquareRoot( sProm, sRoot, orders, firstInds, cutoff );

    // a := Q_{sqrt(s)}(z)
    // -------------------
    DistMultiVec<PReal> a(grid);
    soc::ApplyQuadratic( sRoot, zProm, a, orders, firstInds, cutoff );

    // a := inv(sqrt(a)) = inv(sqrt((Q_{sqrt(s)}(z))))
    // -----------------------------------------------
    DistMultiVec<PReal> b(grid);
    soc::SquareRoot( a, b, orders, firstInds, cutoff );
    soc::Inverse( b, a, orders, firstInds, cutoff );

    // w := Q_{sqrt(s)}(a)
    // -------------------
    DistMultiVec<PReal> wProm(grid);
    soc::ApplyQuadratic( sRoot, a, wProm, orders, firstInds, cutoff );
    Copy( wProm, w );
}

// See Section 4.2 of
// http://www.seas.ucla.edu/~vandenbe/publications/coneprog.pdf

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
void VandenbergheNT
( const Matrix<Real>& s,
  const Matrix<Real>& z,
        Matrix<Real>& w,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds )
{
    EL_DEBUG_CSE
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
    soc::Dets( sProm, sDets, orders, firstInds );
    soc::Dets( zProm, zDets, orders, firstInds );
    cone::Broadcast( sDets, orders, firstInds );
    cone::Broadcast( zDets, orders, firstInds );
    for( Int i=0; i<n; ++i )
    {
        sProm(i) /= Sqrt(sDets(i));
        zProm(i) /= Sqrt(zDets(i));
    }

    // Compute the 'gamma' coefficients
    // ================================
    Matrix<PReal> gammas;
    soc::Dots( zProm, sProm, gammas, orders, firstInds );
    cone::Broadcast( gammas, orders, firstInds );
    for( Int i=0; i<n; ++i )
        gammas(i) = Sqrt((PReal(1)+gammas(i))/PReal(2));

    // Compute the normalized scaling point
    // ====================================
    auto wProm = zProm;
    soc::Reflect( wProm, orders, firstInds );
    wProm += sProm;
    DiagonalSolve( LEFT, NORMAL, gammas, wProm );
    wProm *= PReal(1)/PReal(2);

    // Rescale the scaling point
    // =========================
    for( Int i=0; i<n; ++i )
    {
        const PReal sDet = sDets(i);
        const PReal zDet = zDets(i);
        const PReal scale = Pow(sDet,PReal(0.25))/Pow(zDet,PReal(0.25));
        wProm(i) *= scale;
    }
    Copy( wProm, w );
}

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
void VandenbergheNT
( const AbstractDistMatrix<Real>& sPre,
  const AbstractDistMatrix<Real>& zPre,
        AbstractDistMatrix<Real>& wPre,
  const AbstractDistMatrix<Int>& ordersPre,
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    EL_DEBUG_CSE
    typedef Promote<Real> PReal;
    const Grid& grid = sPre.Grid();
    AssertSameGrids( sPre, zPre, wPre, ordersPre, firstIndsPre );

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrix<PReal,VC,STAR> s(grid), z(grid);
    Copy( sPre, s );
    Copy( zPre, z );

    DistMatrixWriteProxy<Real,PReal,VC,STAR>
      wProx( wPre, ctrl );
    DistMatrixReadProxy<Int,Int,VC,STAR>
      ordersProx( ordersPre, ctrl ),
      firstIndsProx( firstIndsPre, ctrl );
    auto& w = wProx.Get();
    auto& orders = ordersProx.GetLocked();
    auto& firstInds = firstIndsProx.GetLocked();

    // Normalize with respect to the Jordan determinant
    // ================================================
    const Int nLocal = s.LocalHeight();
    DistMatrix<PReal,VC,STAR> sDets(grid), zDets(grid);
    soc::Dets( s, sDets, orders, firstInds, cutoff );
    soc::Dets( z, zDets, orders, firstInds, cutoff );
    cone::Broadcast( sDets, orders, firstInds, cutoff );
    cone::Broadcast( zDets, orders, firstInds, cutoff );
    auto& sLoc = s.Matrix();
    auto& zLoc = z.Matrix();
    auto& sDetsLoc = sDets.LockedMatrix();
    auto& zDetsLoc = zDets.LockedMatrix();
    for( Int iLoc=0; iLoc<nLocal; ++iLoc )
    {
        sLoc(iLoc) /= Sqrt(sDetsLoc(iLoc));
        zLoc(iLoc) /= Sqrt(zDetsLoc(iLoc));
    }

    // Compute the 'gamma' coefficients
    // ================================
    DistMatrix<PReal,VC,STAR> gammas(grid);
    soc::Dots( z, s, gammas, orders, firstInds, cutoff );
    cone::Broadcast( gammas, orders, firstInds, cutoff );
    auto& gammasLoc = gammas.Matrix();
    for( Int iLoc=0; iLoc<nLocal; ++iLoc )
        gammasLoc(iLoc) = Sqrt((PReal(1)+gammasLoc(iLoc))/PReal(2));

    // Compute the normalized scaling point
    // ====================================
    w = z;
    soc::Reflect( w, orders, firstInds );
    w += s;
    DiagonalSolve( LEFT, NORMAL, gammas, w );
    w *= PReal(1)/PReal(2);

    // Rescale the scaling point
    // =========================
    auto& wLoc = w.Matrix();
    for( Int iLoc=0; iLoc<nLocal; ++iLoc )
    {
        const PReal sDet = sDetsLoc(iLoc);
        const PReal zDet = zDetsLoc(iLoc);
        const PReal scale = Pow(sDet,PReal(0.25))/Pow(zDet,PReal(0.25));
        wLoc(iLoc) *= scale;
    }
}

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
void VandenbergheNT
( const DistMultiVec<Real>& s,
  const DistMultiVec<Real>& z,
        DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    EL_DEBUG_CSE
    typedef Promote<Real> PReal;
    const Grid& grid = s.Grid();

    DistMultiVec<PReal> sProm(grid), zProm(grid);
    Copy( s, sProm );
    Copy( z, zProm );

    // Normalize with respect to the Jordan determinant
    // ================================================
    const Int nLocal = sProm.LocalHeight();
    DistMultiVec<PReal> sDets(grid), zDets(grid);
    soc::Dets( sProm, sDets, orders, firstInds, cutoff );
    soc::Dets( zProm, zDets, orders, firstInds, cutoff );
    cone::Broadcast( sDets, orders, firstInds, cutoff );
    cone::Broadcast( zDets, orders, firstInds, cutoff );
    auto& sPromLoc = sProm.Matrix();
    auto& zPromLoc = zProm.Matrix();
    auto& sDetsLoc = sDets.LockedMatrix();
    auto& zDetsLoc = zDets.LockedMatrix();
    for( Int iLoc=0; iLoc<nLocal; ++iLoc )
    {
        sPromLoc(iLoc) /= Sqrt(sDetsLoc(iLoc));
        zPromLoc(iLoc) /= Sqrt(zDetsLoc(iLoc));
    }

    // Compute the 'gamma' coefficients
    // ================================
    DistMultiVec<PReal> gammas(grid);
    soc::Dots( zProm, sProm, gammas, orders, firstInds, cutoff );
    cone::Broadcast( gammas, orders, firstInds, cutoff );
    auto& gammasLoc = gammas.Matrix();
    for( Int iLoc=0; iLoc<nLocal; ++iLoc )
        gammasLoc(iLoc) = Sqrt((PReal(1)+gammasLoc(iLoc))/PReal(2));

    // Compute the normalized scaling point
    // ====================================
    auto wProm = zProm;
    soc::Reflect( wProm, orders, firstInds );
    wProm += sProm;
    DiagonalSolve( LEFT, NORMAL, gammas, wProm );
    wProm *= PReal(1)/PReal(2);

    // Rescale the scaling point
    // =========================
    auto& wPromLoc = wProm.Matrix();
    for( Int iLoc=0; iLoc<nLocal; ++iLoc )
    {
        const PReal sDet = sDetsLoc(iLoc);
        const PReal zDet = zDetsLoc(iLoc);
        const PReal scale = Pow(sDet,PReal(0.25))/Pow(zDet,PReal(0.25));
        wPromLoc(iLoc) *= scale;
    }
    Copy( wProm, w );
}

} // anonymous namespace

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void NesterovTodd
( const Matrix<Real>& s,
  const Matrix<Real>& z,
        Matrix<Real>& w,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds )
{
    EL_DEBUG_CSE
    const bool useClassical = false;
    if( useClassical )
        ClassicalNT( s, z, w, orders, firstInds );
    else
        VandenbergheNT( s, z, w, orders, firstInds );
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void NesterovTodd
( const AbstractDistMatrix<Real>& s,
  const AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& w,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
  Int cutoff )
{
    EL_DEBUG_CSE
    const bool useClassical = false;
    if( useClassical )
        ClassicalNT( s, z, w, orders, firstInds, cutoff );
    else
        VandenbergheNT( s, z, w, orders, firstInds, cutoff );
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void NesterovTodd
( const DistMultiVec<Real>& s,
  const DistMultiVec<Real>& z,
        DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    EL_DEBUG_CSE
    const bool useClassical = false;
    if( useClassical )
        ClassicalNT( s, z, w, orders, firstInds, cutoff );
    else
        VandenbergheNT( s, z, w, orders, firstInds, cutoff );
}

#define PROTO(Real) \
  template void NesterovTodd \
  ( const Matrix<Real>& s, \
    const Matrix<Real>& z, \
          Matrix<Real>& w, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void NesterovTodd \
  ( const AbstractDistMatrix<Real>& s, \
    const AbstractDistMatrix<Real>& z, \
          AbstractDistMatrix<Real>& w, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, \
    Int cutoff ); \
  template void NesterovTodd \
  ( const DistMultiVec<Real>& s, \
    const DistMultiVec<Real>& z, \
          DistMultiVec<Real>& w, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace soc
} // namespace El
