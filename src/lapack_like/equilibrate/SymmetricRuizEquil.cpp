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
inline Real DampScaling( Real alpha )
{
    static const Real tol = Pow(Epsilon<Real>(),Real(0.33));
    if( alpha == Real(0) )
        return 1;
    else 
        return Max(alpha,tol);
}

template<typename Real>
inline Real SquareRootScaling( Real alpha )
{
    return Sqrt(alpha);
}

template<typename F>
void SymmetricRuizEquil
( Matrix<F>& A, 
  Matrix<Base<F>>& d, 
  Int maxIter, bool progress )
{
    DEBUG_ONLY(CSE cse("SymmetricRuizEquil"))
    LogicError("This routine is not yet written");
}

template<typename F>
void SymmetricRuizEquil
( ElementalMatrix<F>& APre, 
  ElementalMatrix<Base<F>>& dPre,
  Int maxIter, bool progress )
{
    DEBUG_ONLY(CSE cse("SymmetricRuizEquil"))
    typedef Base<F> Real;

    ElementalProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;
    auto APtr = ReadWriteProxy<F,MC,MR>(&APre,control);
    auto dPtr = WriteProxy<Real,MC,STAR>(&dPre,control); 
    auto& A = *APtr;
    auto& d = *dPtr;

    const Int n = A.Height();
    Ones( d, n, 1 );

    DistMatrix<Real,MR,STAR> scales(A.Grid());
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Rescale the columns (and rows)
        // ------------------------------
        ColumnMaxNorms( A, scales );
        EntrywiseMap( scales, function<Real(Real)>(DampScaling<Real>) );
        EntrywiseMap( scales, function<Real(Real)>(SquareRootScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, d );
        // TODO: Replace with SymmetricDiagonalSolve
        DiagonalSolve( RIGHT, NORMAL, scales, A );
        DiagonalSolve( LEFT, NORMAL, scales, A );
    }
    SetIndent( indent );
}

template<typename F>
void SymmetricRuizEquil
( SparseMatrix<F>& A, 
  Matrix<Base<F>>& d,
  Int maxIter, bool progress )
{
    DEBUG_ONLY(CSE cse("SymmetricRuizEquil"))
    typedef Base<F> Real;
    const Int n = A.Height();
    Ones( d, n, 1 );

    Matrix<Real> scales;
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Rescale the columns (and rows)
        // ------------------------------
        ColumnMaxNorms( A, scales );
        EntrywiseMap( scales, function<Real(Real)>(DampScaling<Real>) );
        EntrywiseMap( scales, function<Real(Real)>(SquareRootScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, d );
        SymmetricDiagonalSolve( scales, A );
    }
    SetIndent( indent );
}

template<typename F>
void SymmetricRuizEquil
( DistSparseMatrix<F>& A, 
  DistMultiVec<Base<F>>& d, 
  Int maxIter, bool progress )
{
    DEBUG_ONLY(CSE cse("SymmetricRuizEquil"))
    typedef Base<F> Real;
    const Int n = A.Height();
    mpi::Comm comm = A.Comm();
    d.SetComm( comm );
    Ones( d, n, 1 );

    DistMultiVec<Real> scales(comm);
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Rescale the columns (and rows)
        // ------------------------------
        ColumnMaxNorms( A, scales );
        EntrywiseMap( scales, function<Real(Real)>(DampScaling<Real>) );
        EntrywiseMap( scales, function<Real(Real)>(SquareRootScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, d );
        SymmetricDiagonalSolve( scales, A );
    }
    SetIndent( indent );
}

#define PROTO(F) \
  template void SymmetricRuizEquil \
  ( Matrix<F>& A, \
    Matrix<Base<F>>& d, \
    Int maxIter, bool progress ); \
  template void SymmetricRuizEquil \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<Base<F>>& d, \
    Int maxIter, bool progress ); \
  template void SymmetricRuizEquil \
  ( SparseMatrix<F>& A, \
    Matrix<Base<F>>& d, \
    Int maxIter, bool progress ); \
  template void SymmetricRuizEquil \
  ( DistSparseMatrix<F>& A, \
    DistMultiVec<Base<F>>& d, \
    Int maxIter, bool progress );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
