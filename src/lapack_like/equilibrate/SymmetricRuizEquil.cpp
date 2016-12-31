/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename Real>
Real DampScaling( const Real& alpha )
{
    static const Real tol = Pow(limits::Epsilon<Real>(),Real(0.33));
    if( alpha == Real(0) )
        return 1;
    else
        return Max(alpha,tol);
}

template<typename Real>
Real SquareRootScaling( const Real& alpha )
{
    return Sqrt(alpha);
}

template<typename Field>
void SymmetricRuizEquil
( Matrix<Field>& A,
  Matrix<Base<Field>>& d,
  Int maxIter, bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = A.Height();
    Ones( d, n, 1 );

    Matrix<Real> scales;
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Rescale the columns (and rows)
        // ------------------------------
        ColumnMaxNorms( A, scales );
        EntrywiseMap( scales, MakeFunction(DampScaling<Real>) );
        EntrywiseMap( scales, MakeFunction(SquareRootScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, d );
        // TODO(poulson): Replace with SymmetricDiagonalSolve
        DiagonalSolve( RIGHT, NORMAL, scales, A );
        DiagonalSolve( LEFT, NORMAL, scales, A );
    }
    SetIndent( indent );
}

template<typename Field>
void SymmetricRuizEquil
( AbstractDistMatrix<Field>& APre,
  AbstractDistMatrix<Base<Field>>& dPre,
  Int maxIter, bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    ElementalProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;

    DistMatrixReadWriteProxy<Field,Field,MC,MR> AProx( APre, control );
    DistMatrixWriteProxy<Real,Real,MC,STAR> dProx( dPre, control );
    auto& A = AProx.Get();
    auto& d = dProx.Get();

    const Int n = A.Height();
    Ones( d, n, 1 );

    DistMatrix<Real,MR,STAR> scales(A.Grid());
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Rescale the columns (and rows)
        // ------------------------------
        ColumnMaxNorms( A, scales );
        EntrywiseMap( scales, MakeFunction(DampScaling<Real>) );
        EntrywiseMap( scales, MakeFunction(SquareRootScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, d );
        // TODO(poulson): Replace with SymmetricDiagonalSolve
        DiagonalSolve( RIGHT, NORMAL, scales, A );
        DiagonalSolve( LEFT, NORMAL, scales, A );
    }
    SetIndent( indent );
}

template<typename Field>
void SymmetricRuizEquil
( SparseMatrix<Field>& A,
  Matrix<Base<Field>>& d,
  Int maxIter, bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = A.Height();
    Ones( d, n, 1 );

    Matrix<Real> scales;
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Rescale the columns (and rows)
        // ------------------------------
        ColumnMaxNorms( A, scales );
        EntrywiseMap( scales, MakeFunction(DampScaling<Real>) );
        EntrywiseMap( scales, MakeFunction(SquareRootScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, d );
        SymmetricDiagonalSolve( scales, A );
    }
    SetIndent( indent );
}

template<typename Field>
void SymmetricRuizEquil
( DistSparseMatrix<Field>& A,
  DistMultiVec<Base<Field>>& d,
  Int maxIter, bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = A.Height();
    const Grid& grid = A.Grid();
    d.SetGrid( grid );
    Ones( d, n, 1 );

    DistMultiVec<Real> scales(grid);
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Rescale the columns (and rows)
        // ------------------------------
        ColumnMaxNorms( A, scales );
        EntrywiseMap( scales, MakeFunction(DampScaling<Real>) );
        EntrywiseMap( scales, MakeFunction(SquareRootScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, d );
        SymmetricDiagonalSolve( scales, A );
    }
    SetIndent( indent );
}

#define PROTO(Field) \
  template void SymmetricRuizEquil \
  ( Matrix<Field>& A, \
    Matrix<Base<Field>>& d, \
    Int maxIter, bool progress ); \
  template void SymmetricRuizEquil \
  ( AbstractDistMatrix<Field>& A, \
    AbstractDistMatrix<Base<Field>>& d, \
    Int maxIter, bool progress ); \
  template void SymmetricRuizEquil \
  ( SparseMatrix<Field>& A, \
    Matrix<Base<Field>>& d, \
    Int maxIter, bool progress ); \
  template void SymmetricRuizEquil \
  ( DistSparseMatrix<Field>& A, \
    DistMultiVec<Base<Field>>& d, \
    Int maxIter, bool progress );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
