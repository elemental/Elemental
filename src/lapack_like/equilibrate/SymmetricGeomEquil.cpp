/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include "./Util.hpp"

// The following routines are adaptations of the approach uses by
// Saunders et al. (originally recommended by Joseph Fourer) for iteratively
// rescaling the columns and rows by their approximate geometric means in order
// to better scale the original problem. After this iteration is finished,
// the columns or rows are rescaled so that their maximum entry has magnitude
// one (unless the row/column is identically zero).
//
// The implementation of Saunders et al. is commonly referred to as either
// gmscale.m or gmscal.m.

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
void SymmetricGeomEquil
( Matrix<Field>& A, Matrix<Base<Field>>& d, bool progress )
{
    EL_DEBUG_CSE
    // TODO(poulson): Ensure A is symmetric
    typedef Base<Field> Real;
    const Int n = A.Height();
    Ones( d, n, 1 );

    // TODO(poulson): Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 10;
    const Real damp = Real(1)/Real(1000);
    const Real relTol = Real(9)/Real(10);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbs = MaxAbsLoc( A );
    const Real maxAbsVal = maxAbs.value;
    if( maxAbsVal == Real(0) )
        return;
    const Real minAbsVal = MinAbsNonzero( A, maxAbsVal );
    Real ratio = maxAbsVal / minAbsVal;
    if( progress )
        cout << "    Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
             << ratio << endl;

    Matrix<Real> scales;
    Zeros( scales, n, 1 );
    const Real sqrtDamp = Sqrt(damp);
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically equilibrate the columns
        // TODO(poulson): Reimplement this as a standalone routine?
        for( Int j=0; j<n; ++j )
        {
            auto aCol = A( ALL, IR(j) );
            auto maxColAbs = VectorMaxAbsLoc( aCol );
            const Real maxColAbsVal = maxColAbs.value;
            const Real minColAbsVal = MinAbsNonzero( aCol, maxColAbsVal );
            const Real propScale = Sqrt(minColAbsVal*maxColAbsVal);
            const Real scale = Max(propScale,sqrtDamp*maxColAbsVal);
            scales(j) = scale;
        }
        EntrywiseMap( scales, MakeFunction(DampScaling<Real>) );
        EntrywiseMap( scales, MakeFunction(SquareRootScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, d );
        SymmetricDiagonalSolve( scales, A );

        auto newMaxAbs = MaxAbsLoc( A );
        const Real newMaxAbsVal = newMaxAbs.value;
        const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress )
            cout << "    New ratio is " << newMaxAbsVal << "/"
                 << newMinAbsVal << "=" << newRatio << endl;
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }

    // Normalize the maximum absolute values towards one
    for( Int normIter=0; normIter<3; ++normIter )
    {
        for( Int j=0; j<n; ++j )
        {
            Real maxAbs = 1;
            for( Int i=0; i<n; ++i )
                maxAbs = Max( Abs(A(i,j)), maxAbs );
            scales(j) = Sqrt(maxAbs);
        }
        DiagonalScale( LEFT, NORMAL, scales, d );
        SymmetricDiagonalSolve( scales, A );
    }
    auto newMaxAbs = MaxAbsLoc( A );
    const Real newMaxAbsVal = newMaxAbs.value;
    const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
    const Real newRatio = newMaxAbsVal / newMinAbsVal;
    if( progress )
        cout << "    Final ratio is " << newMaxAbsVal << "/"
             << newMinAbsVal << "=" << newRatio << endl;
}

template<typename Field>
void SymmetricGeomEquil
( AbstractDistMatrix<Field>& APre,
  AbstractDistMatrix<Base<Field>>& dPre,
  bool progress )
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
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    Ones( d, n, 1 );

    // TODO(poulson): Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 10;
    const Real relTol = Real(9)/Real(10);

    // TODO(poulson): Incorporate damping
    //const Real damp = Real(1)/Real(1000);
    //const Real sqrtDamp = Sqrt(damp);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbs = MaxAbsLoc( A );
    const Real maxAbsVal = maxAbs.value;
    if( maxAbsVal == Real(0) )
        return;
    const Real minAbsVal = MinAbsNonzero( A, maxAbsVal );
    Real ratio = maxAbsVal / minAbsVal;
    if( progress && A.Grid().Rank() == 0 )
        cout << "    Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
             << ratio << endl;

    DistMatrix<Real,MR,STAR> scales(A.Grid());
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically equilibrate the columns
        // -------------------------------------
        // TODO(poulson): Remove GeometricColumnScaling
        GeometricColumnScaling( A, scales );
        EntrywiseMap( scales, MakeFunction(DampScaling<Real>) );
        EntrywiseMap( scales, MakeFunction(SquareRootScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, d );
        // TODO(poulson): SymmetricDiagonalSolve
        DiagonalSolve( RIGHT, NORMAL, scales, A );
        DiagonalSolve( LEFT, NORMAL, scales, A );

        auto newMaxAbs = MaxAbsLoc( A );
        const Real newMaxAbsVal = newMaxAbs.value;
        const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress && A.Grid().Rank() == 0 )
            cout << "    New ratio is " << newMaxAbsVal << "/"
                 << newMinAbsVal << "=" << newRatio << endl;
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }

    // Normalize the maximum absolute values towards one
    scales.AlignWith( A );
    Zeros( scales, n, 1 );
    auto& ALoc = A.Matrix();
    auto& scalesLoc = scales.Matrix();
    for( Int normIter=0; normIter<3; ++normIter )
    {
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            Real maxAbs = 1;
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                maxAbs = Max( Abs(ALoc(iLoc,jLoc)), maxAbs );
            const Real scale = Sqrt(maxAbs);
            scalesLoc(jLoc) = scale;
        }
        mpi::AllReduce( scales.Buffer(), nLocal, mpi::MAX, A.ColComm() );
        DiagonalScale( LEFT, NORMAL, scales, d );
        // TODO(poulson): SymmetricDiagonalSolve
        DiagonalSolve( RIGHT, NORMAL, scales, A );
        DiagonalSolve( LEFT, NORMAL, scales, A );
    }
    auto newMaxAbs = MaxAbsLoc( A );
    const Real newMaxAbsVal = newMaxAbs.value;
    const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
    const Real newRatio = newMaxAbsVal / newMinAbsVal;
    if( progress && A.Grid().Rank() == 0 )
        cout << "    Final ratio is " << newMaxAbsVal << "/"
             << newMinAbsVal << "=" << newRatio << endl;
}

template<typename Field>
void SymmetricGeomEquil
( SparseMatrix<Field>& A,
  Matrix<Base<Field>>& d,
  bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = A.Height();
    Ones( d, n, 1 );

    // TODO(poulson): Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 10;
    const Real damp = Real(1)/Real(1000);
    const Real relTol = Real(9)/Real(10);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbs = MaxAbsLoc( A );
    const Real maxAbsVal = maxAbs.value;
    if( maxAbsVal == Real(0) )
        return;
    const Real minAbsVal = MinAbsNonzero( A, maxAbsVal );
    Real ratio = maxAbsVal / minAbsVal;
    if( progress )
        cout << "    Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
             << ratio << endl;

    const Real sqrtDamp = Sqrt(damp);
    Matrix<Real> scales(n,1), rowMaxAbs, rowMinAbs;
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Symmetrically geometrically rescale
        // -----------------------------------
        RowMaxNorms( A, rowMaxAbs );
        RowMinAbsNonzero( A, rowMaxAbs, rowMinAbs );
        for( Int j=0; j<n; ++j )
        {
            const Real maxAbs = rowMaxAbs(j);
            const Real minAbs = rowMinAbs(j);
            const Real propScale = Sqrt(minAbs*maxAbs);
            const Real scale = Max(propScale,sqrtDamp*maxAbs);
            scales(j) = scale;
        }
        EntrywiseMap( scales, MakeFunction(DampScaling<Real>) );
        EntrywiseMap( scales, MakeFunction(SquareRootScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, d );
        SymmetricDiagonalSolve( scales, A );

        // Determine whether we are done or not
        // ------------------------------------
        auto newMaxAbs = MaxAbsLoc( A );
        const Real newMaxAbsVal = newMaxAbs.value;
        const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress )
            cout << "    New ratio is " << newMaxAbsVal << "/"
                 << newMinAbsVal << "=" << newRatio << endl;
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }

    // Normalize the maximum absolute values towards one
    for( Int normIter=0; normIter<3; ++normIter )
    {
        for( Int i=0; i<n; ++i )
        {
            Real maxAbs = 1;
            const Int numConnections = A.NumConnections(i);
            const Int entryOff = A.RowOffset(i);
            for( Int k=0; k<numConnections; ++k )
                maxAbs = Max( Abs(A.Value(entryOff+k)), maxAbs );
            scales(i) = Sqrt(maxAbs);
        }
        DiagonalScale( LEFT, NORMAL, scales, d );
        SymmetricDiagonalSolve( scales, A );
    }
    auto newMaxAbs = MaxAbsLoc( A );
    const Real newMaxAbsVal = newMaxAbs.value;
    const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
    const Real newRatio = newMaxAbsVal / newMinAbsVal;
    if( progress )
        cout << "    Final ratio is " << newMaxAbsVal << "/"
             << newMinAbsVal << "=" << newRatio << endl;
}

template<typename Field>
void SymmetricGeomEquil
( DistSparseMatrix<Field>& A,
  DistMultiVec<Base<Field>>& d,
  bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = A.Height();
    const Grid& grid = A.Grid();
    const int commRank = grid.Rank();
    d.SetGrid( grid );
    Ones( d, n, 1 );

    // TODO(poulson): Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 10;
    const Real damp = Real(1)/Real(1000);
    const Real relTol = Real(9)/Real(10);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbs = MaxAbsLoc( A );
    const Real maxAbsVal = maxAbs.value;
    if( maxAbsVal == Real(0) )
        return;
    const Real minAbsVal = MinAbsNonzero( A, maxAbsVal );
    Real ratio = maxAbsVal / minAbsVal;
    if( progress && commRank == 0 )
        cout << "    Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
             << ratio << endl;

    DistMultiVec<Real> scales(grid), rowMaxAbs(grid), rowMinAbs(grid);
    Zeros( scales, n, 1 );
    auto& scalesLoc = scales.Matrix();
    const Real sqrtDamp = Sqrt(damp);
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Symmetrically geometrically rescale
        // -----------------------------------
        RowMaxNorms( A, rowMaxAbs );
        RowMinAbsNonzero( A, rowMaxAbs, rowMinAbs );
        auto& rowMaxAbsLoc = rowMaxAbs.Matrix();
        auto& rowMinAbsLoc = rowMinAbs.Matrix();
        const Int localHeight = scales.LocalHeight();
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Real maxAbs = rowMaxAbsLoc(iLoc);
            const Real minAbs = rowMinAbsLoc(iLoc);
            const Real propScale = Sqrt(minAbs*maxAbs);
            const Real scale = Max(propScale,sqrtDamp*maxAbs);
            scalesLoc(iLoc) = scale;
        }
        EntrywiseMap( scales, MakeFunction(DampScaling<Real>) );
        EntrywiseMap( scales, MakeFunction(SquareRootScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, d );
        SymmetricDiagonalSolve( scales, A );

        // Determine whether we are done or not
        // ------------------------------------
        auto newMaxAbs = MaxAbsLoc( A );
        const Real newMaxAbsVal = newMaxAbs.value;
        const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress && commRank == 0 )
            cout << "    New ratio is " << newMaxAbsVal << "/"
                 << newMinAbsVal << "=" << newRatio << endl;
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }

    // Normalize the maximum absolute values towards one
    const Int localHeight = A.LocalHeight();
    for( Int normIter=0; normIter<3; ++normIter )
    {
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            Real maxAbs = 1;
            const Int numConnections = A.NumConnections(iLoc);
            const Int entryOff = A.RowOffset(iLoc);
            for( Int k=0; k<numConnections; ++k )
                maxAbs = Max( Abs(A.Value(entryOff+k)), maxAbs );
            scalesLoc(iLoc) = Sqrt(maxAbs);
        }
        DiagonalScale( LEFT, NORMAL, scales, d );
        SymmetricDiagonalSolve( scales, A );
    }
    auto newMaxAbs = MaxAbsLoc( A );
    const Real newMaxAbsVal = newMaxAbs.value;
    const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
    const Real newRatio = newMaxAbsVal / newMinAbsVal;
    if( progress && commRank == 0 )
        cout << "    Final ratio is " << newMaxAbsVal << "/"
             << newMinAbsVal << "=" << newRatio << endl;
}

#define PROTO(Field) \
  template void SymmetricGeomEquil \
  ( Matrix<Field>& A, Matrix<Base<Field>>& d, bool progress ); \
  template void SymmetricGeomEquil \
  ( AbstractDistMatrix<Field>& A, \
    AbstractDistMatrix<Base<Field>>& d, bool progress ); \
  template void SymmetricGeomEquil \
  ( SparseMatrix<Field>& A, Matrix<Base<Field>>& d, bool progress ); \
  template void SymmetricGeomEquil \
  ( DistSparseMatrix<Field>& A, DistMultiVec<Base<Field>>& d, bool progress );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
