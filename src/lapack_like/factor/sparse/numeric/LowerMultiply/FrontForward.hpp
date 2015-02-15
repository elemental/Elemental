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
#pragma once
#ifndef EL_FACTOR_NUMERIC_LOWERMULTIPLY_FRONTFORWARD_HPP
#define EL_FACTOR_NUMERIC_LOWERMULTIPLY_FRONTFORWARD_HPP

#include "../LowerSolve/FrontUtil.hpp"

namespace El {

template<typename F>
inline void FrontVanillaLowerForwardMultiply( const Matrix<F>& L, Matrix<F>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontVanillaLowerForwardMultiply");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal multiply:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
    )
    Matrix<F> LT, LB, XT, XB;
    LockedPartitionDown( L, LT, LB, L.Width() );
    PartitionDown( X, XT, XB, L.Width() );

    Gemm( NORMAL, NORMAL, F(1), LB, XT, F(1), XB );
    Trmm( LEFT, LOWER, NORMAL, UNIT, F(1), LT, XT );
}

template<typename F>
inline void 
FrontLowerForwardMultiply( const SymmFront<F>& front, Matrix<F>& W )
{
    DEBUG_ONLY(CallStackEntry cse("FrontLowerForwardMultiply"))
    const SymmFrontType type = front.type;
    if( Unfactored(type) )
        LogicError("Cannot multiply against an unfactored front");
    if( BlockFactorization(type) || PivotedFactorization(type) )
        LogicError("Blocked and pivoted factorizations not supported");

    FrontVanillaLowerForwardMultiply( front.L, W );
}

template<typename F>
inline void FrontVanillaLowerForwardMultiply
( const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
    DEBUG_ONLY(
      CallStackEntry cse("FrontVanillaLowerForwardMultiply");
      if( L.Grid() != X.Grid() )
          LogicError("L and X must be distributed over the same grid");
      if( L.Height() < L.Width() || L.Height() != X.Height() )
          LogicError
          ("Nonconformal multiply:\n",
           DimsString(L,"L"),"\n",DimsString(X,"X"));
      if( L.ColAlign() != X.ColAlign() )
          LogicError("L and X are assumed to be aligned");
    )
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        FrontVanillaLowerForwardMultiply( L.LockedMatrix(), X.Matrix() );
        return;
    }

    // Separate the top and bottom portions of X and L
    const Int snSize = L.Width();
    DistMatrix<F,VC,STAR> LT(g), LB(g), XT(g), XB(g);
    LockedPartitionDown( L, LT, LB, snSize );
    PartitionDown( X, XT, XB, snSize );

    // XB := XB + LB XT
    Gemm( NORMAL, NORMAL, F(1), LB, XT, F(1), XB );

    // XT := LT XT
    Trmm( LEFT, LOWER, NORMAL, UNIT, F(1), LT, XT );
}

template<typename F>
inline void FrontVanillaLowerForwardMultiply
( const DistMatrix<F>& L, DistMatrix<F>& X )
{
    DEBUG_ONLY(
      CallStackEntry cse("FrontVanillaLowerForwardMultiply");
      if( L.Grid() != X.Grid() )
          LogicError("L and X must be distributed over the same grid");
      if( L.Height() < L.Width() || L.Height() != X.Height() )
          LogicError
          ("Nonconformal multiply:\n",
           DimsString(L,"L"),"\n",DimsString(X,"X"));
    )
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        FrontVanillaLowerForwardMultiply( L.LockedMatrix(), X.Matrix() );
        return;
    }

    // Separate the top and bottom portions of X and L
    const Int snSize = L.Width();
    DistMatrix<F> LT(g), LB(g), XT(g), XB(g);
    LockedPartitionDown( L, LT, LB, snSize );
    PartitionDown( X, XT, XB, snSize );

    // XB := XB + LB XT
    Gemm( NORMAL, NORMAL, F(1), LB, XT, F(1), XB );

    // XT := LT XT
    Trmm( LEFT, LOWER, NORMAL, UNIT, F(1), LT, XT );
}

template<typename F>
inline void 
FrontLowerForwardMultiply
( const DistSymmFront<F>& front, DistMatrix<F,VC,STAR>& W )
{
    DEBUG_ONLY(CallStackEntry cse("FrontLowerForwardMultiply"))
    const SymmFrontType type = front.type;

    if( type == LDL_1D )
        FrontVanillaLowerForwardMultiply( front.L1D, W );
    else
        LogicError("Unsupported front type");
}

template<typename F>
inline void 
FrontLowerForwardMultiply( const DistSymmFront<F>& front, DistMatrix<F>& W )
{
    DEBUG_ONLY(CallStackEntry cse("FrontLowerForwardMultiply"))
    const SymmFrontType type = front.type;

    if( type == LDL_2D )
        FrontVanillaLowerForwardMultiply( front.L2D, W );
    else
        LogicError("Unsupported front type");
}

} // namespace El

#endif // ifndef EL_FACTOR_NUMERIC_LOWERMULTIPLY_FRONTFORWARD_HPP
