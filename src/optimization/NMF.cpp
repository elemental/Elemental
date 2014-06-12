/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include EL_ZEROS_INC

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
void NMF( const DistMatrix<Real>& A, DistMatrix<Real>& X, DistMatrix<Real>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("NonNegativeLeastSquares"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");

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
  ( const DistMatrix<Real>& A, DistMatrix<Real>& X, DistMatrix<Real>& Y );

PROTO(float)
PROTO(double)

} // namespace El
