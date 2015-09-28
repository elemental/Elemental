/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// TODO: Better convergence criterions. E.g., accept a relative tolerance
//       in addition to the maximum number of iterations.

template<typename Real>
void NMF
( const Matrix<Real>& A, 
        Matrix<Real>& X, Matrix<Real>& Y,
  const NMFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("NMF"))

    Matrix<Real> AAdj, XAdj, YAdj;
    Adjoint( A, AAdj );

    for( Int iter=0; iter<ctrl.maxIter; ++iter )
    {
        NNLS( X, A, YAdj, ctrl.nnlsCtrl );
        Adjoint( YAdj, Y );
        NNLS( Y, AAdj, XAdj, ctrl.nnlsCtrl );
        Adjoint( XAdj, X );
    }
}

template<typename Real>
void NMF
( const ElementalMatrix<Real>& APre, 
        ElementalMatrix<Real>& XPre, 
        ElementalMatrix<Real>& YPre,
  const NMFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("NMF"))

    auto APtr = ReadProxy<Real,MC,MR>( &APre );      auto& A = *APtr;
    auto XPtr = ReadWriteProxy<Real,MC,MR>( &XPre ); auto& X = *XPtr;
    auto YPtr = WriteProxy<Real,MC,MR>( &YPre );     auto& Y = *YPtr;

    DistMatrix<Real> AAdj(A.Grid()), XAdj(A.Grid()), YAdj(A.Grid());
    Adjoint( A, AAdj );

    for( Int iter=0; iter<ctrl.maxIter; ++iter )
    {
        NNLS( X, A, YAdj, ctrl.nnlsCtrl );
        Adjoint( YAdj, Y );
        NNLS( Y, AAdj, XAdj, ctrl.nnlsCtrl );
        Adjoint( XAdj, X );
    }
}

#define PROTO(Real) \
  template void NMF \
  ( const Matrix<Real>& A, \
          Matrix<Real>& X, \
          Matrix<Real>& Y, \
    const NMFCtrl<Real>& ctrl ); \
  template void NMF \
  ( const ElementalMatrix<Real>& A, \
          ElementalMatrix<Real>& X, \
          ElementalMatrix<Real>& Y, \
    const NMFCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
