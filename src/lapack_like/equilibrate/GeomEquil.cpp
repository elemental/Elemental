/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F>
inline Base<F> MinAbsNonzero( Matrix<F>& A, Base<F> upperBound )
{
    DEBUG_ONLY(CallStackEntry cse("MinAbsNonzero"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Real minAbs = upperBound;
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            const Real alphaAbs = Abs(A.Get(i,j));
            if( alphaAbs > Real(0) && alphaAbs < minAbs )
                minAbs = alphaAbs;
        }
    }
    return minAbs;
}

template<typename F>
void GeomEquil( Matrix<F>& A, Matrix<F>& dRow, Matrix<F>& dCol )
{
    DEBUG_ONLY(CallStackEntry cse("GeomEquil"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Ones( dRow, m, 1 );
    Ones( dCol, n, 1 );

    // TODO: Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 10; 
    const Real damp = Real(1)/Real(1000);
    const Real relTol = Real(9)/Real(10);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbs = MaxAbs( A );
    Real maxAbsVal = maxAbs.value;
    if( maxAbsVal == Real(0) )
        return;
    Real minAbsVal = MinAbsNonzero( A, maxAbsVal );
    Real ratio = maxAbsVal / minAbsVal;

    const Real sqrtDamp = Sqrt(damp);
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically equilibrate the columns
        for( Int j=0; j<n; ++j )
        {
            auto aCol = A( IR(0,m), IR(j,j+1) );
            auto maxColAbs = VectorMaxAbs( aCol );
            const Real maxColAbsVal = maxColAbs.value;
            const Real minColAbsVal = MinAbsNonzero( aCol, maxColAbsVal );
            if( maxColAbsVal > Real(0) )
            {
                const Real propScale = Sqrt(minColAbsVal*maxColAbsVal);
                const Real scale = Max(propScale,sqrtDamp*maxColAbsVal);
                Scale( Real(1)/scale, aCol );         
                dCol.Set( j, 0, scale*dCol.Get(j,0) );
            }
        }

        // Geometrically equilibrate the rows
        for( Int i=0; i<m; ++i )
        {
            auto aRow = A( IR(i,i+1), IR(0,n) );
            auto maxRowAbs = VectorMaxAbs( aRow );
            const Real maxRowAbsVal = maxRowAbs.value;
            const Real minRowAbsVal = MinAbsNonzero( aRow, maxRowAbsVal );
            if( maxRowAbsVal > Real(0) )
            {
                const Real propScale = Sqrt(minRowAbsVal*maxRowAbsVal);
                const Real scale = Max(propScale,sqrtDamp*maxRowAbsVal);
                Scale( Real(1)/scale, aRow );         
                dRow.Set( i, 0, scale*dRow.Get(i,0) );
            }
        }

        auto newMaxAbs = MaxAbs( A );
        Real newMaxAbsVal = newMaxAbs.value;
        Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }

    // Scale each column so that its maximum entry is 1 or 0
    for( Int j=0; j<n; ++j )
    {
        auto aCol = A( IR(0,m), IR(j,j+1) );
        auto maxColAbs = VectorMaxAbs( aCol );
        const Real maxColAbsVal = maxColAbs.value;
        const Real minColAbsVal = MinAbsNonzero( aCol, maxColAbsVal );
        if( maxColAbsVal > Real(0) )
        {
            Scale( Real(1)/maxColAbsVal, aCol );         
            dCol.Set( j, 0, maxColAbsVal*dCol.Get(j,0) );
        }
    }
}

template<typename F>
void GeomEquil
( AbstractDistMatrix<F>& A, 
  AbstractDistMatrix<F>& dRow, AbstractDistMatrix<F>& dCol )
{
    DEBUG_ONLY(CallStackEntry cse("GeomEquil"))
    LogicError("This routine is not yet written");
}

template<typename F>
void GeomEquil
( SparseMatrix<F>& A, Matrix<F>& dRow, Matrix<F>& dCol )
{
    DEBUG_ONLY(CallStackEntry cse("GeomEquil"))
    LogicError("This routine is not yet written");
}

template<typename F>
void GeomEquil
( DistSparseMatrix<F>& A, DistMultiVec<F>& dRow, DistMultiVec<F>& dCol )
{
    DEBUG_ONLY(CallStackEntry cse("GeomEquil"))
    LogicError("This routine is not yet written");
}

#define PROTO(F) \
  template void GeomEquil \
  ( Matrix<F>& A, Matrix<F>& dRow, Matrix<F>& dCol ); \
  template void GeomEquil \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& dRow, AbstractDistMatrix<F>& dCol ); \
  template void GeomEquil \
  ( SparseMatrix<F>& A, Matrix<F>& dRow, Matrix<F>& dCol ); \
  template void GeomEquil \
  ( DistSparseMatrix<F>& A, DistMultiVec<F>& dRow, DistMultiVec<F>& dCol );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
