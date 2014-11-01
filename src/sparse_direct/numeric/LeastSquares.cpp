/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

namespace El {

// TODO: Switch to a control structure for SymmetricSolve
template<typename F>
void LeastSquares
( const DistSparseMatrix<F>& A, const DistMultiVec<F>& Y, DistMultiVec<F>& X,
  bool sequential, int numDistSeps, int numSeqSteps, int cutoff )
{
    DEBUG_ONLY(CallStackEntry cse("LeastSquares"))
    const Int n = A.Width();
    DistSparseMatrix<F> C(A.Comm());
    C.Resize( n, n );
    Syrk( ADJOINT, F(1), A, F(0), C );
    Zeros( X, n, Y.Width() );
    Multiply( ADJOINT, F(1), A, Y, F(0), X ); 
    HermitianSolve( C, X, sequential, numDistSeps, numSeqSteps, cutoff );
}

#define PROTO(F) \
  template void LeastSquares \
  ( const DistSparseMatrix<F>& A, const DistMultiVec<F>& Y, \
    DistMultiVec<F>& X, bool sequential, int numDistSeps, int numSeqSeps, \
    int cutoff );
 
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
