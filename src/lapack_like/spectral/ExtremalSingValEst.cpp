/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename F>
pair<Base<F>,Base<F>>
ExtremalSingValEst( const SparseMatrix<F>& A, Int basisSize )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    Matrix<Real> T;
    ProductLanczos( A, T, basisSize );
    const Int k = T.Height();
    if( k == 0 )
        return pair<Real,Real>(0,0);

    Matrix<Real> d, dSub;
    d = GetDiagonal( T );
    dSub = GetDiagonal( T, -1 );

    Matrix<Real> w;
    HermitianTridiagEig( d, dSub, w );

    pair<Real,Real> extremal;
    extremal.first = Sqrt( Max(w(0),Real(0)) );
    extremal.second = Sqrt( Max(w(k-1),Real(0)) );
    return extremal;
}

template<typename F>
pair<Base<F>,Base<F>>
ExtremalSingValEst( const DistSparseMatrix<F>& A, Int basisSize )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Grid& grid = A.Grid();
    DistMatrix<Real,STAR,STAR> T(grid);
    ProductLanczos( A, T, basisSize );
    const Int k = T.Height();
    if( k == 0 )
        return pair<Real,Real>(0,0);

    auto d = GetDiagonal( T.Matrix() );
    auto dSub = GetDiagonal( T.Matrix(), -1 );

    Matrix<Real> w;
    HermitianTridiagEig( d, dSub, w );

    pair<Real,Real> extremal;
    extremal.first = Sqrt( Max(w(0),Real(0)) );
    extremal.second = Sqrt( Max(w(k-1),Real(0)) );
    return extremal;
}

template<typename F>
pair<Base<F>,Base<F>>
HermitianExtremalSingValEst( const SparseMatrix<F>& A, Int basisSize )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    Matrix<Real> T;
    Lanczos( A, T, basisSize );
    const Int k = T.Height();
    if( k == 0 )
        return pair<Real,Real>(0,0);

    Matrix<Real> d, dSub;
    d = GetDiagonal( T );
    dSub = GetDiagonal( T, -1 );

    Matrix<Real> w;
    HermitianTridiagEig( d, dSub, w );

    pair<Real,Real> extremal;
    extremal.second = MaxNorm(w);
    extremal.first = extremal.second;
    for( Int i=0; i<k; ++i )
        extremal.first = Min(extremal.first,Abs(w(i)));
    return extremal;
}

template<typename F>
pair<Base<F>,Base<F>>
HermitianExtremalSingValEst( const DistSparseMatrix<F>& A, Int basisSize )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Grid& grid = A.Grid();

    DistMatrix<Real,STAR,STAR> T(grid);
    Lanczos( A, T, basisSize );
    const Int k = T.Height();
    if( k == 0 )
        return pair<Real,Real>(0,0);

    auto d = GetDiagonal( T.Matrix() );
    auto dSub = GetDiagonal( T.Matrix(), -1 );

    Matrix<Real> w;
    HermitianTridiagEig( d, dSub, w );

    pair<Real,Real> extremal;
    extremal.second = MaxNorm(w);
    extremal.first = extremal.second;
    for( Int i=0; i<k; ++i )
        extremal.first = Min(extremal.first,Abs(w(i)));

    return extremal;
}

#define PROTO(F) \
  template pair<Base<F>,Base<F>> ExtremalSingValEst \
  ( const SparseMatrix<F>& A, Int basisSize ); \
  template pair<Base<F>,Base<F>> ExtremalSingValEst \
  ( const DistSparseMatrix<F>& A, Int basisSize ); \
  template pair<Base<F>,Base<F>> HermitianExtremalSingValEst \
  ( const SparseMatrix<F>& A, Int basisSize ); \
  template pair<Base<F>,Base<F>> HermitianExtremalSingValEst \
  ( const DistSparseMatrix<F>& A, Int basisSize );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
