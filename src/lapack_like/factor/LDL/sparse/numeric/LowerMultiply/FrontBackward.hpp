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
#ifndef EL_FACTOR_LDL_NUMERIC_LOWERMULTIPLY_FRONTBACKWARD_HPP
#define EL_FACTOR_LDL_NUMERIC_LOWERMULTIPLY_FRONTBACKWARD_HPP

namespace El {
namespace ldl {

template<typename F>
void FrontVanillaLowerBackwardMultiply
( const Matrix<F>& L, Matrix<F>& X, bool conjugate )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( L.Height() < L.Width() || L.Height() != X.Height() )
          LogicError
          ("Nonconformal multiply:\n",
           DimsString(L,"L"),"\n",DimsString(X,"X"));
    )
    Matrix<F> LT, LB, XT, XB;
    LockedPartitionDown( L, LT, LB, L.Width() );
    PartitionDown( X, XT, XB, L.Width() );

    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    Trmm( LEFT, LOWER, orientation, UNIT, F(1), LT, XT );
    Gemm( orientation, NORMAL, F(1), LB, XB, F(1), XT );
}

template<typename F>
void FrontLowerBackwardMultiply
( const Front<F>& front, Matrix<F>& W, bool conjugate )
{
    EL_DEBUG_CSE
    auto type = front.type;
    if( Unfactored(type) )
        LogicError("Cannot multiply against an unfactored matrix");

    if( front.sparseLeaf )
    {
        LogicError
        ("Sparse leaves not yet supported in FrontLowerBackwardMultiply");
    }
    else
    {
        if( type == LDL_2D )
            FrontVanillaLowerBackwardMultiply( front.LDense, W, conjugate );
        else
            LogicError("Unsupported front type");
    }
}

template<typename F>
void FrontVanillaLowerBackwardMultiply
( const DistMatrix<F,VC,STAR>& L,
        DistMatrix<F,VC,STAR>& X,
  bool conjugate )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
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
        FrontVanillaLowerBackwardMultiply
        ( L.LockedMatrix(), X.Matrix(), conjugate );
        return;
    }

    DistMatrix<F,VC,STAR> LT(g), LB(g), XT(g), XB(g);
    LockedPartitionDown( L, LT, LB, L.Width() );
    PartitionDown( X, XT, XB, L.Width() );

    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    Trmm( LEFT, LOWER, orientation, UNIT, F(1), LT, XT );

    if( XB.Height() != 0 )
    {
        // Subtract off the parent updates
        DistMatrix<F,STAR,STAR> Z(g);
        LocalGemm( orientation, NORMAL, F(1), LB, XB, Z );
        AxpyContract( F(1), Z, XT );
    }
}

template<typename F>
void FrontVanillaLowerBackwardMultiply
( const DistMatrix<F>& L,
        DistMatrix<F>& X,
  bool conjugate )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
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
        FrontVanillaLowerBackwardMultiply
        ( L.LockedMatrix(), X.Matrix(), conjugate );
        return;
    }

    DistMatrix<F> LT(g), LB(g), XT(g), XB(g);
    LockedPartitionDown( L, LT, LB, L.Width() );
    PartitionDown( X, XT, XB, L.Width() );

    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    Trmm( LEFT, LOWER, orientation, UNIT, F(1), LT, XT );
    Gemm( orientation, NORMAL, F(1), LB, XB, F(1), XT );
}

template<typename F>
void FrontLowerBackwardMultiply
( const DistFront<F>& front, DistMatrix<F>& W, bool conjugate )
{
    EL_DEBUG_CSE
    if( Unfactored(front.type) )
        LogicError("Cannot multiply against an unfactored matrix");

    if( front.type == LDL_2D )
        FrontVanillaLowerBackwardMultiply( front.L2D, W, conjugate );
    else
        LogicError("Unsupported front type");
}

template<typename F>
void FrontLowerBackwardMultiply
( const DistFront<F>& front, DistMatrix<F,VC,STAR>& W, bool conjugate )
{
    EL_DEBUG_CSE
    if( Unfactored(front.type) )
        LogicError("Cannot multiply against an unfactored matrix");

    if( front.type == LDL_1D )
        FrontVanillaLowerBackwardMultiply( front.L1D, W, conjugate );
    else
        LogicError("Unsupported front type");
}

} // namespace ldl
} // namespace El

#endif // ifndef EL_FACTOR_LDL_NUMERIC_LOWERMULTIPLY_FRONTBACKWARD_HPP
