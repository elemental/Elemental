/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./SVD/Chan.hpp"
#include "./SVD/Thresholded.hpp"

namespace El {

// Grab the full SVD of the general matrix A, A = U diag(s) V^H
// ============================================================
// NOTE: On exit, A is overwritten with U

template<typename F>
void SVD
( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V, const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("SVD"))
    if( ctrl.thresholded )
    {
        svd::Thresholded( A, s, V, ctrl.tol, ctrl.relative );
    }
    else
    {
        if( ctrl.seqQR )
            svd::QRSVD( A, s, V );
        else
            svd::DivideAndConquerSVD( A, s, V );
    }
}

template<typename F>
void SVD
( AbstractDistMatrix<F>& A, AbstractDistMatrix<Base<F>>& s, 
  AbstractDistMatrix<F>& V, const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("SVD"))
    if( ctrl.thresholded )
    {
        if( A.ColDist() == VC && A.RowDist() == STAR )
        {
            auto& ACast = dynamic_cast<DistMatrix<F,VC,STAR>&>( A );
            svd::Thresholded( ACast, s, V, ctrl.tol, ctrl.relative );
        }
        else
            svd::Thresholded( A, s, V, ctrl.tol, ctrl.relative );
    }
    else
        svd::Chan( A, s, V, ctrl.fullChanRatio );
}

// Return the singular values
// ==========================

template<typename F>
void SVD( Matrix<F>& A, Matrix<Base<F>>& s )
{
    DEBUG_ONLY(CallStackEntry cse("SVD"))
    const Int m = A.Height();
    const Int n = A.Width();
    s.Resize( Min(m,n), 1 );
    lapack::SVD( m, n, A.Buffer(), A.LDim(), s.Buffer() );
}

template<typename F>
void SVD
( AbstractDistMatrix<F>& A, AbstractDistMatrix<Base<F>>& s, 
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("SVD"))
    // TODO: Add more options
    svd::Chan( A, s, ctrl.valChanRatio );
}

#define PROTO(F) \
  template void SVD( Matrix<F>& A, Matrix<Base<F>>& s ); \
  template void SVD \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<Base<F>>& s, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template void SVD \
  ( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template void SVD \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<Base<F>>& s, \
    AbstractDistMatrix<F>& V, const SVDCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
