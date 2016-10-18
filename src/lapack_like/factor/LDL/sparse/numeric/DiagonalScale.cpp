/*
   Copyright (c) 2009-2012, Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin.
   All rights reserved.

   Copyright (c) 2013, Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright (c) 2013-2014, Jack Poulson and 
   The Georgia Institute of Technology.
   All rights reserved.

   Copyright (c) 2014-2015, Jack Poulson and Stanford University.
   All rights reserved.
   
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace ldl {

template<typename F>
void DiagonalScale
( const NodeInfo& info, const Front<F>& front, MatrixNode<F>& X )
{
    DEBUG_CSE

    const Int numChildren = info.children.size();
    for( Int c=0; c<numChildren; ++c )
        DiagonalScale( *info.children[c], *front.children[c], *X.children[c] );

    if( PivotedFactorization(front.type) )
        QuasiDiagonalScale
        ( LEFT, LOWER, front.diag, front.subdiag, 
          X.matrix, front.isHermitian );
    else
        DiagonalScale( LEFT, NORMAL, front.diag, X.matrix );
}

template<typename F>
void DiagonalScale
( const DistNodeInfo& info, const DistFront<F>& front, DistMultiVecNode<F>& X )
{
    DEBUG_CSE

    if( front.child == nullptr )
    {
        DiagonalScale( *info.duplicate, *front.duplicate, *X.duplicate );
        return;
    }
    DiagonalScale( *info.child, *front.child, *X.child );

    if( PivotedFactorization(front.type) )
        QuasiDiagonalScale
        ( LEFT, LOWER, front.diag, front.subdiag, X.matrix, front.isHermitian );
    else
        DiagonalScale( LEFT, NORMAL, front.diag, X.matrix );
}

template<typename F>
void DiagonalScale
( const DistNodeInfo& info, const DistFront<F>& front, DistMatrixNode<F>& X )
{
    DEBUG_CSE

    if( front.child == nullptr )
    {
        DiagonalScale( *info.duplicate, *front.duplicate, *X.duplicate );
        return;
    }
    DiagonalScale( *info.child, *front.child, *X.child );

    if( PivotedFactorization(front.type) )
        QuasiDiagonalScale
        ( LEFT, LOWER, front.diag, front.subdiag, X.matrix, front.isHermitian );
    else
        DiagonalScale( LEFT, NORMAL, front.diag, X.matrix );
}

#define PROTO(F) \
  template void DiagonalScale \
  ( const NodeInfo& info, const Front<F>& front, \
    MatrixNode<F>& X ); \
  template void DiagonalScale \
  ( const DistNodeInfo& info, const DistFront<F>& front, \
    DistMultiVecNode<F>& X ); \
  template void DiagonalScale \
  ( const DistNodeInfo& info, const DistFront<F>& front, \
    DistMatrixNode<F>& X );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace ldl
} // namespace El
