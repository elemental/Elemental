/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

namespace eig {

template<typename Real>
void Helper
( Matrix<Real>& A,
  Matrix<Complex<Real>>& w,
  Matrix<Complex<Real>>& X )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;

    Matrix<Real> Q;
    Schur( A, w, Q );

    Matrix<C> ACpx;
    schur::RealToComplex( A, Q, ACpx, X );

    Matrix<C> R;
    TriangEig( ACpx, R );

    Trmm( RIGHT, UPPER, NORMAL, NON_UNIT, C(1), R, X );
}

template<typename Real>
void Helper
( Matrix<Complex<Real>>& A,
  Matrix<Complex<Real>>& w,
  Matrix<Complex<Real>>& X )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;

    Schur( A, w, X );

    Matrix<C> R;
    TriangEig( A, R );

    Trmm( RIGHT, UPPER, NORMAL, NON_UNIT, C(1), R, X );
}

template<typename Real>
void Helper
( AbstractDistMatrix<Real>& APre,
  AbstractDistMatrix<Complex<Real>>& w,
  AbstractDistMatrix<Complex<Real>>& XPre )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;

    DistMatrixReadProxy<Real,Real,MC,MR> AProxy( APre );
    DistMatrixWriteProxy<C,C,MC,MR> XProxy( XPre );
    auto& A = AProxy.Get();
    auto& X = XProxy.Get();

    const Grid& g = A.Grid();
    DistMatrix<Real> Q(g);
    Schur( A, w, Q );

    DistMatrix<C> ACpx(g);
    schur::RealToComplex( A, Q, ACpx, X );

    DistMatrix<C> R(g);
    TriangEig( ACpx, R );

    Trmm( RIGHT, UPPER, NORMAL, NON_UNIT, C(1), R, X );
}

template<typename Real>
void Helper
( AbstractDistMatrix<Complex<Real>>& APre,
  AbstractDistMatrix<Complex<Real>>& w,
  AbstractDistMatrix<Complex<Real>>& XPre )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;

    DistMatrixReadProxy<C,C,MC,MR> AProxy( APre );
    DistMatrixWriteProxy<C,C,MC,MR> XProxy( XPre );
    auto& A = AProxy.Get();
    auto& X = XProxy.Get();

    Schur( A, w, X );

    DistMatrix<C> R( A.Grid() );
    TriangEig( A, R );

    Trmm( RIGHT, UPPER, NORMAL, NON_UNIT, C(1), R, X );
}

} // namespace eig

template<typename Field>
void Eig
( Matrix<Field>& A,
  Matrix<Complex<Base<Field>>>& w,
  Matrix<Complex<Base<Field>>>& X )
{
    EL_DEBUG_CSE
    eig::Helper( A, w, X );
}

template<typename Field>
void Eig
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Complex<Base<Field>>>& w,
  AbstractDistMatrix<Complex<Base<Field>>>& X )
{
    EL_DEBUG_CSE
    eig::Helper( A, w, X );
}

#define PROTO(Field) \
  template void Eig \
  ( Matrix<Field>& A, \
    Matrix<Complex<Base<Field>>>& w, \
    Matrix<Complex<Base<Field>>>& X ); \
  template void Eig \
  ( AbstractDistMatrix<Field>& A, \
    AbstractDistMatrix<Complex<Base<Field>>>& w, \
    AbstractDistMatrix<Complex<Base<Field>>>& X );

#define EL_NO_INT_PROTO
#include <El/macros/Instantiate.h>

} // namespace El
