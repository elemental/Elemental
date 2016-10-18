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
void ImageAndKernel
( const Matrix<F>& B,
        Matrix<F>& M,
        Matrix<F>& K )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    const Real eps = limits::Epsilon<Real>();

    SVDCtrl<Real> ctrl;
    ctrl.bidiagSVDCtrl.approach = FULL_SVD;
    Matrix<F> U, V;
    Matrix<Real> s;
    SVD( B, U, s, V, ctrl );

    const Int numSingVals = s.Height();
    const Real twoNorm = ( numSingVals==0 ? Real(0) : s(0) );

    // TODO: Incorporate a user-defined (relative) threshold
    const Real relTol = Max(m,n)*eps;
    const Real tol = twoNorm*relTol;

    Int rank = numSingVals;
    for( Int j=0; j<numSingVals; ++j )
    {
        if( s(j) <= tol )
        {
            rank = j;
            break;
        }
    }

    auto UL = U( ALL, IR(0,rank) );
    auto VR = V( ALL, IR(rank,END) );
    M = UL;
    K = VR;
}

template<typename F>
void ImageAndKernel
( const ElementalMatrix<F>& B,
        ElementalMatrix<F>& M,
        ElementalMatrix<F>& K )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    const Real eps = limits::Epsilon<Real>();
    const Grid& g = B.Grid();

    SVDCtrl<Real> ctrl;
    ctrl.bidiagSVDCtrl.approach = FULL_SVD;
    DistMatrix<F> U(g), V(g);
    DistMatrix<Real,STAR,STAR> s(g);
    SVD( B, U, s, V, ctrl );

    const Int numSingVals = s.Height();
    const Real twoNorm = ( numSingVals==0 ? Real(0) : s.Get(0,0) );

    // TODO: Incorporate a user-defined (relative) threshold
    const Real relTol = Max(m,n)*eps;
    const Real tol = twoNorm*relTol;

    Int rank = numSingVals;
    auto& sLoc = s.Matrix();
    for( Int j=0; j<numSingVals; ++j )
    {
        if( sLoc(j) <= tol )
        {
            rank = j;
            break;
        }
    }

    auto UL = U( ALL, IR(0,rank) );
    auto VR = V( ALL, IR(rank,END) );
    Copy( UL, M );
    Copy( VR, K );
}

template<typename F>
void Image
( const Matrix<F>& B,
        Matrix<F>& M )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    const Real eps = limits::Epsilon<Real>();

    SVDCtrl<Real> ctrl;
    ctrl.bidiagSVDCtrl.approach = COMPACT_SVD;
    ctrl.bidiagSVDCtrl.wantV = false;
    // TODO(poulson): Let the user change these defaults
    ctrl.bidiagSVDCtrl.tolType = RELATIVE_TO_MAX_SING_VAL_TOL;
    ctrl.bidiagSVDCtrl.tol = Max(m,n)*eps;
    Matrix<F> V;
    Matrix<Real> s;
    SVD( B, M, s, V, ctrl );
}

template<typename F>
void Image
( const ElementalMatrix<F>& B,
        ElementalMatrix<F>& M )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    const Real eps = limits::Epsilon<Real>();
    const Grid& g = B.Grid();

    SVDCtrl<Real> ctrl;
    ctrl.bidiagSVDCtrl.approach = COMPACT_SVD;
    ctrl.bidiagSVDCtrl.wantV = false;
    // TODO(poulson): Let the user change these defaults
    ctrl.bidiagSVDCtrl.tolType = RELATIVE_TO_MAX_SING_VAL_TOL;
    ctrl.bidiagSVDCtrl.tol = Max(m,n)*eps;
    DistMatrix<F> V(g);
    DistMatrix<Real,STAR,STAR> s(g);
    SVD( B, M, s, V, ctrl );
}

template<typename F>
void Kernel
( const Matrix<F>& B,
        Matrix<F>& K )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    const Real eps = limits::Epsilon<Real>();

    SVDCtrl<Real> ctrl;
    ctrl.bidiagSVDCtrl.approach = FULL_SVD;
    ctrl.bidiagSVDCtrl.wantU = false;
    Matrix<F> U, V;
    Matrix<Real> s;
    SVD( B, U, s, V, ctrl );

    const Int numSingVals = s.Height();
    const Real twoNorm = ( numSingVals==0 ? Real(0) : s(0) );

    // TODO: Incorporate a user-defined (relative) threshold
    const Real relTol = Max(m,n)*eps;
    const Real tol = twoNorm*relTol;

    Int rank = numSingVals;
    for( Int j=0; j<numSingVals; ++j )
    {
        if( s(j) <= tol )
        {
            rank = j;
            break;
        }
    }

    auto VR = V( ALL, IR(rank,END) );
    K = VR;
}

template<typename F>
void Kernel
( const ElementalMatrix<F>& B,
        ElementalMatrix<F>& K )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    const Real eps = limits::Epsilon<Real>();
    const Grid& g = B.Grid();

    SVDCtrl<Real> ctrl;
    ctrl.bidiagSVDCtrl.approach = FULL_SVD;
    ctrl.bidiagSVDCtrl.wantU = false;
    DistMatrix<F> U(g), V(g);
    DistMatrix<Real,STAR,STAR> s(g);
    SVD( B, U, s, V, ctrl );

    const Int numSingVals = s.Height();
    const Real twoNorm = ( numSingVals==0 ? Real(0) : s.Get(0,0) );

    // TODO: Incorporate a user-defined (relative) threshold
    const Real relTol = Max(m,n)*eps;
    const Real tol = twoNorm*relTol;

    Int rank = numSingVals;
    auto& sLoc = s.Matrix();
    for( Int j=0; j<numSingVals; ++j )
    {
        if( sLoc(j) <= tol )
        {
            rank = j;
            break;
        }
    }

    auto VR = V( ALL, IR(rank,END) );
    Copy( VR, K );
}

#define PROTO(F) \
  template void ImageAndKernel \
  ( const Matrix<F>& B, \
          Matrix<F>& M, \
          Matrix<F>& K ); \
  template void ImageAndKernel \
  ( const ElementalMatrix<F>& B, \
          ElementalMatrix<F>& M, \
          ElementalMatrix<F>& K ); \
  template void Image \
  ( const Matrix<F>& B, \
          Matrix<F>& M ); \
  template void Image \
  ( const ElementalMatrix<F>& B, \
          ElementalMatrix<F>& M ); \
  template void Kernel \
  ( const Matrix<F>& B, \
          Matrix<F>& K ); \
  template void Kernel \
  ( const ElementalMatrix<F>& B, \
          ElementalMatrix<F>& K );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
