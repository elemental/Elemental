/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

namespace El {

template<typename F>
void LeastSquares
( Orientation orientation,
  const DistSparseMatrix<F>& A, const DistMultiVec<F>& Y, DistMultiVec<F>& X,
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(
        CallStackEntry cse("LeastSquares");
        if( orientation == NORMAL && A.Height() != Y.Height() )
            LogicError("Heights of A and Y must match");
        if( orientation != NORMAL && A.Width() != Y.Height() )
            LogicError("Width of A and height of Y must match");
    )
    DistSparseMatrix<F> C(A.Comm());
    if( orientation == NORMAL )
    {
        const Int n = A.Width();
        Herk( LOWER, ADJOINT, Base<F>(1), A, C );
        MakeHermitian( LOWER, C );
        X.SetComm( Y.Comm() );
        Zeros( X, n, Y.Width() );
        Multiply( ADJOINT, F(1), A, Y, F(0), X ); 
    }
    else if( orientation == ADJOINT || !IsComplex<F>::val )
    {
        const Int n = A.Height();
        Herk( LOWER, NORMAL, Base<F>(1), A, C );
        MakeHermitian( LOWER, C );
        X.SetComm( Y.Comm() );
        Zeros( X, n, Y.Width() );
        Multiply( NORMAL, F(1), A, Y, F(0), X );
    }
    else
    {
        LogicError("Complex transposed option not yet supported");
    }
    HermitianSolve( C, X, ctrl );
}

#define PROTO(F) \
  template void LeastSquares \
  ( Orientation orientation, \
    const DistSparseMatrix<F>& A, const DistMultiVec<F>& Y, \
    DistMultiVec<F>& X, const BisectCtrl& ctrl );
 
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
