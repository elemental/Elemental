/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename Field>
void ImageAndKernel
( const Matrix<Field>& B,
        Matrix<Field>& M,
        Matrix<Field>& K )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    const Real eps = limits::Epsilon<Real>();

    SVDCtrl<Real> ctrl;
    ctrl.bidiagSVDCtrl.approach = FULL_SVD;
    Matrix<Field> U, V;
    Matrix<Real> s;
    SVD( B, U, s, V, ctrl );

    const Int numSingVals = s.Height();
    const Real twoNorm = ( numSingVals==0 ? Real(0) : s(0) );

    // TODO(poulson): Incorporate a user-defined (relative) threshold
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

template<typename Field>
void ImageAndKernel
( const AbstractDistMatrix<Field>& B,
        AbstractDistMatrix<Field>& M,
        AbstractDistMatrix<Field>& K )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    const Real eps = limits::Epsilon<Real>();
    const Grid& g = B.Grid();

    SVDCtrl<Real> ctrl;
    ctrl.bidiagSVDCtrl.approach = FULL_SVD;
    DistMatrix<Field> U(g), V(g);
    DistMatrix<Real,STAR,STAR> s(g);
    SVD( B, U, s, V, ctrl );

    const Int numSingVals = s.Height();
    const Real twoNorm = ( numSingVals==0 ? Real(0) : s.Get(0,0) );

    // TODO(poulson): Incorporate a user-defined (relative) threshold
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

template<typename Field>
void Image
( const Matrix<Field>& B,
        Matrix<Field>& M )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    const Real eps = limits::Epsilon<Real>();

    SVDCtrl<Real> ctrl;
    ctrl.bidiagSVDCtrl.approach = COMPACT_SVD;
    ctrl.bidiagSVDCtrl.wantV = false;
    // TODO(poulson): Let the user change these defaults
    ctrl.bidiagSVDCtrl.tolType = RELATIVE_TO_MAX_SING_VAL_TOL;
    ctrl.bidiagSVDCtrl.tol = Max(m,n)*eps;
    Matrix<Field> V;
    Matrix<Real> s;
    SVD( B, M, s, V, ctrl );
}

template<typename Field>
void Image
( const AbstractDistMatrix<Field>& B,
        AbstractDistMatrix<Field>& M )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
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
    DistMatrix<Field> V(g);
    DistMatrix<Real,STAR,STAR> s(g);
    SVD( B, M, s, V, ctrl );
}

template<typename Field>
void Kernel
( const Matrix<Field>& B,
        Matrix<Field>& K )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    const Real eps = limits::Epsilon<Real>();

    SVDCtrl<Real> ctrl;
    ctrl.bidiagSVDCtrl.approach = FULL_SVD;
    ctrl.bidiagSVDCtrl.wantU = false;
    Matrix<Field> U, V;
    Matrix<Real> s;
    SVD( B, U, s, V, ctrl );

    const Int numSingVals = s.Height();
    const Real twoNorm = ( numSingVals==0 ? Real(0) : s(0) );

    // TODO(poulson): Incorporate a user-defined (relative) threshold
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

template<typename Field>
void Kernel
( const AbstractDistMatrix<Field>& B,
        AbstractDistMatrix<Field>& K )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    const Real eps = limits::Epsilon<Real>();
    const Grid& g = B.Grid();

    SVDCtrl<Real> ctrl;
    ctrl.bidiagSVDCtrl.approach = FULL_SVD;
    ctrl.bidiagSVDCtrl.wantU = false;
    DistMatrix<Field> U(g), V(g);
    DistMatrix<Real,STAR,STAR> s(g);
    SVD( B, U, s, V, ctrl );

    const Int numSingVals = s.Height();
    const Real twoNorm = ( numSingVals==0 ? Real(0) : s.Get(0,0) );

    // TODO(poulson): Incorporate a user-defined (relative) threshold
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

#define PROTO(Field) \
  template void ImageAndKernel \
  ( const Matrix<Field>& B, \
          Matrix<Field>& M, \
          Matrix<Field>& K ); \
  template void ImageAndKernel \
  ( const AbstractDistMatrix<Field>& B, \
          AbstractDistMatrix<Field>& M, \
          AbstractDistMatrix<Field>& K ); \
  template void Image \
  ( const Matrix<Field>& B, \
          Matrix<Field>& M ); \
  template void Image \
  ( const AbstractDistMatrix<Field>& B, \
          AbstractDistMatrix<Field>& M ); \
  template void Kernel \
  ( const Matrix<Field>& B, \
          Matrix<Field>& K ); \
  template void Kernel \
  ( const AbstractDistMatrix<Field>& B, \
          AbstractDistMatrix<Field>& K );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
