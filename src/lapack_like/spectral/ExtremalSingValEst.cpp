/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

namespace extremal_sing_val {

template<typename F,typename=EnableIf<IsBlasScalar<F>>>
pair<Base<F>,Base<F>>
Helper( const SparseMatrix<F>& A, Int basisSize )
{
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
    HermitianTridiagEig( d, dSub, w, ASCENDING );
    
    pair<Real,Real> extremal;
    extremal.first = Sqrt( Max(w.Get(0,0),Real(0)) );
    extremal.second = Sqrt( Max(w.Get(k-1,0),Real(0)) );
    return extremal;
}

template<typename F,typename=EnableIf<IsBlasScalar<F>>>
pair<Base<F>,Base<F>>
Helper( const DistSparseMatrix<F>& A, Int basisSize )
{
    typedef Base<F> Real;
    Grid grid( A.Comm() );
    DistMatrix<Real,STAR,STAR> T(grid);
    ProductLanczos( A, T, basisSize );
    const Int k = T.Height();
    if( k == 0 )
        return pair<Real,Real>(0,0);

    auto d = GetDiagonal( T.Matrix() );
    auto dSub = GetDiagonal( T.Matrix(), -1 );
    
    Matrix<Real> w;
    HermitianTridiagEig( d, dSub, w, ASCENDING );
    
    pair<Real,Real> extremal;
    extremal.first = Sqrt( Max(w.Get(0,0),Real(0)) );
    extremal.second = Sqrt( Max(w.Get(k-1,0),Real(0)) );
    return extremal;
}

template<typename F,typename=EnableIf<IsBlasScalar<F>>>
pair<Base<F>,Base<F>>
HermitianHelper( const SparseMatrix<F>& A, Int basisSize )
{
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
        extremal.first = Min(extremal.first,Abs(w.Get(i,0)));
    return extremal;
}

template<typename F,typename=EnableIf<IsBlasScalar<F>>>
pair<Base<F>,Base<F>>
HermitianHelper( const DistSparseMatrix<F>& A, Int basisSize )
{
    typedef Base<F> Real;
    Grid grid( A.Comm() );

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
        extremal.first = Min(extremal.first,Abs(w.Get(i,0)));

    return extremal;
}

// Since we don't currently have non-standard eigensolvers, and this is
// just an estimate, convert to and from double-precision
template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
pair<Real,Real>
Helper( const SparseMatrix<Real>& A, Int basisSize )
{
    SparseMatrix<double> ADbl;
    Copy( A, ADbl );
    auto pairDbl = Helper( ADbl, basisSize );
    // NOTE: It seems the conversion could be implicitly performed
    return std::make_pair( Real(pairDbl.first), Real(pairDbl.second) );
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
pair<Real,Real>
Helper( const SparseMatrix<Complex<Real>>& A, Int basisSize )
{
    SparseMatrix<Complex<double>> ADbl;
    Copy( A, ADbl );
    auto pairDbl = Helper( ADbl, basisSize );
    return std::make_pair( Real(pairDbl.first), Real(pairDbl.second) );
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
pair<Real,Real>
Helper( const DistSparseMatrix<Real>& A, Int basisSize )
{
    DistSparseMatrix<double> ADbl(A.Comm());
    Copy( A, ADbl );
    auto pairDbl = Helper( ADbl, basisSize );
    return std::make_pair( Real(pairDbl.first), Real(pairDbl.second) );
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
pair<Real,Real>
Helper( const DistSparseMatrix<Complex<Real>>& A, Int basisSize )
{
    DistSparseMatrix<Complex<double>> ADbl(A.Comm());
    Copy( A, ADbl );
    auto pairDbl = Helper( ADbl, basisSize );
    return std::make_pair( Real(pairDbl.first), Real(pairDbl.second) );
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
pair<Real,Real>
HermitianHelper( const SparseMatrix<Real>& A, Int basisSize )
{
    SparseMatrix<double> ADbl;
    Copy( A, ADbl );
    auto pairDbl = HermitianHelper( ADbl, basisSize );
    return std::make_pair( Real(pairDbl.first), Real(pairDbl.second) );
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
pair<Real,Real>
HermitianHelper( const SparseMatrix<Complex<Real>>& A, Int basisSize )
{
    SparseMatrix<Complex<double>> ADbl;
    Copy( A, ADbl );
    auto pairDbl = HermitianHelper( ADbl, basisSize );
    return std::make_pair( Real(pairDbl.first), Real(pairDbl.second) );
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
pair<Real,Real>
HermitianHelper( const DistSparseMatrix<Real>& A, Int basisSize )
{
    DistSparseMatrix<double> ADbl(A.Comm());
    Copy( A, ADbl );
    auto pairDbl = HermitianHelper( ADbl, basisSize );
    return std::make_pair( Real(pairDbl.first), Real(pairDbl.second) );
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
pair<Real,Real>
HermitianHelper( const DistSparseMatrix<Complex<Real>>& A, Int basisSize )
{
    DistSparseMatrix<Complex<double>> ADbl(A.Comm());
    Copy( A, ADbl );
    auto pairDbl = HermitianHelper( ADbl, basisSize );
    return std::make_pair( Real(pairDbl.first), Real(pairDbl.second) );
}

} // namespace extrmal_sing_val

template<typename F>
pair<Base<F>,Base<F>>
ExtremalSingValEst( const SparseMatrix<F>& A, Int basisSize )
{
    DEBUG_ONLY(CSE cse("ExtremalSingValEst"))
    return extremal_sing_val::Helper( A, basisSize );
}

template<typename F>
pair<Base<F>,Base<F>>
ExtremalSingValEst( const DistSparseMatrix<F>& A, Int basisSize )
{
    DEBUG_ONLY(CSE cse("ExtremalSingValEst"))
    return extremal_sing_val::Helper( A, basisSize );
}

template<typename F>
pair<Base<F>,Base<F>>
HermitianExtremalSingValEst( const SparseMatrix<F>& A, Int basisSize )
{
    DEBUG_ONLY(CSE cse("HermitianExtremalSingValEst"))
    return extremal_sing_val::HermitianHelper( A, basisSize );
}

template<typename F>
pair<Base<F>,Base<F>>
HermitianExtremalSingValEst( const DistSparseMatrix<F>& A, Int basisSize )
{
    DEBUG_ONLY(CSE cse("HermitianExtremalSingValEst"))
    return extremal_sing_val::HermitianHelper( A, basisSize );
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
#include "El/macros/Instantiate.h"

} // namespace El
