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

namespace El {

template<typename F>
inline void LocalDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMultiVec<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LocalDiagonalSolve"))
    const Int numLocalNodes = info.localNodes.size();
    if( PivotedFactorization(L.frontType) )
    {
        for( Int s=0; s<numLocalNodes; ++s )
            QuasiDiagonalSolve
            ( LEFT, LOWER, L.localFronts[s].diag, L.localFronts[s].subdiag, 
              X.localNodes[s], L.isHermitian );
    }
    else
    {
        for( Int s=0; s<numLocalNodes; ++s )
            DiagonalSolve
            ( LEFT, NORMAL, L.localFronts[s].diag, X.localNodes[s], true );
    }
}

template<typename F>
inline void LocalDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LocalDiagonalSolve"))
    const Int numLocalNodes = info.localNodes.size();
    if( PivotedFactorization(L.frontType) )
    {
        for( Int s=0; s<numLocalNodes; ++s ) 
            QuasiDiagonalSolve
            ( LEFT, LOWER, L.localFronts[s].diag, L.localFronts[s].subdiag, 
              X.localNodes[s], L.isHermitian );
    }
    else
    {
        for( Int s=0; s<numLocalNodes; ++s )
            DiagonalSolve
            ( LEFT, NORMAL, L.localFronts[s].diag, X.localNodes[s], true );
    }
}

template<typename F> 
inline void DistDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMultiVec<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DistDiagonalSolve"))
    const Int numDistNodes = info.distNodes.size();

    if( PivotedFactorization(L.frontType) )
    {
        for( Int s=1; s<numDistNodes; ++s )
        {
            const DistSymmFront<F>& front = L.distFronts[s];
            QuasiDiagonalSolve
            ( LEFT, LOWER, front.diag1d, front.subdiag1d, X.distNodes[s-1],
              L.isHermitian );
        }
    }
    else
    {
        for( Int s=1; s<numDistNodes; ++s )
        {
            const DistSymmFront<F>& front = L.distFronts[s];
            DiagonalSolve( LEFT, NORMAL, front.diag1d, X.distNodes[s-1], true );
        }
    }
}

template<typename F> 
inline void DistDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DistDiagonalSolve"))
    const Int numDistNodes = info.distNodes.size();

    if( PivotedFactorization(L.frontType) )
    {
        for( Int s=1; s<numDistNodes; ++s )
        {
            const DistSymmFront<F>& front = L.distFronts[s];
            QuasiDiagonalSolve
            ( LEFT, LOWER, front.diag1d, front.subdiag1d, X.distNodes[s-1],
              L.isHermitian );
        }
    }
    else
    {
        for( Int s=1; s<numDistNodes; ++s )
        {
            const DistSymmFront<F>& front = L.distFronts[s];
            DiagonalSolve( LEFT, NORMAL, front.diag1d, X.distNodes[s-1], true );
        }
    }
}

template<typename F>
void DiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMultiVec<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalSolve"))
    LocalDiagonalSolve( info, L, X );
    DistDiagonalSolve( info, L, X );
}

template<typename F>
void DiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalSolve"))
    LocalDiagonalSolve( info, L, X );
    DistDiagonalSolve( info, L, X );
}

#define PROTO(F) \
  template void DiagonalSolve \
  ( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, \
    DistNodalMultiVec<F>& X ); \
  template void DiagonalSolve \
  ( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, \
    DistNodalMatrix<F>& X );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
