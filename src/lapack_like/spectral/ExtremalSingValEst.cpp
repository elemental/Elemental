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
pair<Base<F>,Base<F>>
ExtremalSingValEst( const SparseMatrix<F>& A, Int basisSize )
{
    DEBUG_ONLY(CallStackEntry cse("ExtremalSingValEst"))
    typedef Base<F> Real;
    const Int minDim = Min(A.Height(),A.Width());

    Matrix<Real> T;
    ProductLanczos( A, T, basisSize );
    const Int k = T.Height();
    if( k == 0 )
    {
        return pair<Real,Real>(0,0);
    }
    else if( k == 1 )
    {
        const Real onlySing = Sqrt( Max(T.Get(0,0),0) );
        return pair<Real,Real>(onlySing,onlySing);
    }

    Matrix<Real> d, dSub;
    d = GetDiagonal( T );
    dSub = GetDiagonal( T, -1 );
    
    Matrix<Real> w;
    HermitianTridiagEig( d, dSub, w, ASCENDING );
    
    pair<Real,Real> extremal;
    extremal.first = Sqrt( Max(w.Get(0,0),0) );
    extremal.second = Sqrt( Max(w.Get(minDim-1,0),0) );
    return extremal;
}

template<typename F>
pair<Base<F>,Base<F>>
ExtremalSingValEst( const DistSparseMatrix<F>& A, Int basisSize )
{
    DEBUG_ONLY(CallStackEntry cse("ExtremalSingValEst"))
    typedef Base<F> Real;
    const Int minDim = Min(A.Height(),A.Width());

    Matrix<Real> T;
    ProductLanczos( A, T, basisSize );
    const Int k = T.Height();
    if( k == 0 )
    {
        return pair<Real,Real>(0,0);
    }
    else if( k == 1 )
    {
        const Real onlySing = Sqrt( Max(T.Get(0,0),0) );
        return pair<Real,Real>(onlySing,onlySing);
    }

    Matrix<Real> d, dSub;
    d = GetDiagonal( T );
    dSub = GetDiagonal( T, -1 );
    
    Matrix<Real> w;
    HermitianTridiagEig( d, dSub, w, ASCENDING );
    
    pair<Real,Real> extremal;
    extremal.first = Sqrt( Max(w.Get(0,0),0) );
    extremal.second = Sqrt( Max(w.Get(minDim-1,0),0) );
    return extremal;
}

#define PROTO(F) \
  template pair<Base<F>,Base<F>> ExtremalSingValEst \
  ( const SparseMatrix<F>& A, Int basisSize ); \
  template pair<Base<F>,Base<F>> ExtremalSingValEst \
  ( const DistSparseMatrix<F>& A, Int basisSize );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
