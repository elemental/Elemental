/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename Real>
void NMF( const Matrix<Real>& A, Matrix<Real>& X, Matrix<Real>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("NMF"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");

    Matrix<Real> AAdj, XAdj, YAdj;
    Adjoint( A, AAdj );

    const Int maxIter = 20;
    for( Int iter=0; iter<maxIter; ++iter )
    {
        NonNegativeLeastSquares( X, A, YAdj );
        Adjoint( YAdj, Y );
        NonNegativeLeastSquares( Y, AAdj, XAdj );
        Adjoint( XAdj, X );
    }
}

template<typename Real>
void NMF
( const AbstractDistMatrix<Real>& APre, AbstractDistMatrix<Real>& XPre, 
        AbstractDistMatrix<Real>& YPre )
{
    DEBUG_ONLY(CallStackEntry cse("NonNegativeLeastSquares"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");

    auto APtr = ReadProxy<Real,MC,MR>( &APre );      auto& A = *APtr;
    auto XPtr = ReadWriteProxy<Real,MC,MR>( &XPre ); auto& X = *XPtr;
    auto YPtr = WriteProxy<Real,MC,MR>( &YPre );     auto& Y = *YPtr;

    DistMatrix<Real> AAdj(A.Grid()), XAdj(A.Grid()), YAdj(A.Grid());
    Adjoint( A, AAdj );

    const Int maxIter = 20;
    for( Int iter=0; iter<maxIter; ++iter )
    {
        NonNegativeLeastSquares( X, A, YAdj );
        Adjoint( YAdj, Y );
        NonNegativeLeastSquares( Y, AAdj, XAdj );
        Adjoint( XAdj, X );
    }
}

#define PROTO(Real) \
  template void NMF \
  ( const Matrix<Real>& A, Matrix<Real>& X, Matrix<Real>& Y ); \
  template void NMF \
  ( const AbstractDistMatrix<Real>& A, AbstractDistMatrix<Real>& X, \
          AbstractDistMatrix<Real>& Y );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
