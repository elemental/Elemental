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
    cScaled *= 1/gamma;

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
 
    // TODO: Expose this as a control parameter
    const bool useSOCP = false;

    if( useSOCP )
    {
        LogicError("This option is not yet supported");
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
