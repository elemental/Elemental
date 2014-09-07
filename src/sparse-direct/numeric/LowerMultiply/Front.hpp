/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SPARSEDIRECT_NUMERIC_LOWERMULTIPLY_FRONT_HPP
#define EL_SPARSEDIRECT_NUMERIC_LOWERMULTIPLY_FRONT_HPP

namespace El {

namespace internal {

template<typename T>
void ModifyForTrmm
( Matrix<T>& D, int diagOff, std::vector<Matrix<T>>& diagonals )
{
    DEBUG_ONLY(CallStackEntry cse("ModifyForTrmm"))
    diagonals.resize( -diagOff );
    for( int i=0; i<-diagOff; ++i )
        diagonals[i].Resize( D.DiagonalLength(-i), 1 );

    const int height = D.Height();
    for( int j=0; j<height; ++j )
    {
        const int length = std::min(-diagOff,height-j);
        for( int i=0; i<length; ++i )
        {
            diagonals[i].Set( j, 0, D.Get(j+i,j) );
            D.Set( j+i, j, T(0) );
        }
    }
}

template<typename T>
void ReplaceAfterTrmm
( Matrix<T>& D, int diagOff,
  const std::vector<Matrix<T>>& diagonals )
{
    DEBUG_ONLY(CallStackEntry cse("ReplaceAfterTrmm"))
    const int height = D.Height();
    for( int j=0; j<height; ++j )
    {
        const int length = std::min(-diagOff,height-j);
        for( int i=0; i<length; ++i )
            D.Set( j+i, j, diagonals[i].Get(j,0) );
    }
}

template<typename T>
void ModifyForTrmm( DistMatrix<T,STAR,STAR>& D, int diagOff )
{
    DEBUG_ONLY(CallStackEntry cse("ModifyForTrmm"))
    const int height = D.Height();
    for( int j=0; j<height; ++j )
    {
        const int length = std::min(-diagOff,height-j);
        MemZero( D.Buffer(j,j), length );
    }
}

} // namespace internal

template<typename T>
inline void FrontLowerMultiply
( Orientation orientation, int diagOff, const Matrix<T>& L, Matrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("FrontLowerMultiply"))
    if( orientation == NORMAL )
        FrontLowerMultiplyNormal( diagOff, L, X );
    else
    {
        const bool conjugate = ( orientation==ADJOINT );
        FrontLowerMultiplyTranspose( diagOff, L, X, conjugate );
    }
}

template<typename T>
inline void FrontLowerMultiply
( Orientation orientation, int diagOff,
  const DistMatrix<T,VC,STAR>& L, DistMatrix<T,VC,STAR>& X )
{
    DEBUG_ONLY(CallStackEntry cse("FrontLowerMultiply"))
    if( orientation == NORMAL )
        FrontLowerMultiplyNormal( diagOff, L, X );
    else
    {
        const bool conjugate = ( orientation==ADJOINT );
        FrontLowerMultiplyTranspose( diagOff, L, X, conjugate );
    }
}

template<typename T>
inline void FrontLowerMultiplyNormal
( int diagOff, const Matrix<T>& L, Matrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontLowerMultiplyNormal");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal multiply:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
        if( diagOff > 0 )
            LogicError("Diagonal offsets cannot be positive");
    )
    // Danger, Will Robinson!
    Matrix<T>* LMod = const_cast<Matrix<T>*>(&L);
    Matrix<T> LT, LB;
    PartitionDown( *LMod, LT, LB, L.Width() );
    Matrix<T> XT, XB;
    PartitionDown( X, XT, XB, L.Width() );

    Gemm( NORMAL, NORMAL, T(1), LB, XT, T(1), XB );

    if( diagOff == 0 )
    {
        Trmm( LEFT, LOWER, NORMAL, NON_UNIT, T(1), LT, XT );
    }
    else
    {
        std::vector<Matrix<T>> diagonals;
        internal::ModifyForTrmm( LT, diagOff, diagonals );
        Trmm( LEFT, LOWER, NORMAL, NON_UNIT, T(1), LT, XT );
        internal::ReplaceAfterTrmm( LT, diagOff, diagonals );
    }
}

// TODO: Simplify this implementation
template<typename T>
inline void FrontLowerMultiplyNormal
( int diagOff, const DistMatrix<T,VC,STAR>& L, DistMatrix<T,VC,STAR>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontLowerMultiplyNormal");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal multiply:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
        if( L.ColAlign() != X.ColAlign() )
            LogicError("L and X are assumed to be aligned");
        if( diagOff > 0 )
            LogicError("Diagonal offsets cannot be positive");
    )
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T,VC,STAR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T,VC,STAR> XT(g),  X0(g),
                          XB(g),  X1(g),
                                  X2(g);

    // Temporary distributions
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,STAR,STAR> X1_STAR_STAR(g);

    // Start the algorithm
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, L.Width() );
    PartitionDown
    ( X, XT,
         XB, L.Width() );
    while( XT.Height() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        RepartitionUp
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );

        //--------------------------------------------------------------------//
        X1_STAR_STAR = X1;
        LocalGemm( NORMAL, NORMAL, T(1), L21, X1_STAR_STAR, T(1), X2 );

        if( diagOff == 0 )
        {
            L11_STAR_STAR = L11;
            LocalTrmm
            ( LEFT, LOWER, NORMAL, NON_UNIT, 
              T(1), L11_STAR_STAR, X1_STAR_STAR );
        }
        else
        {
            L11_STAR_STAR = L11;
            internal::ModifyForTrmm( L11_STAR_STAR, diagOff );
            LocalTrmm
            ( LEFT, LOWER, NORMAL, NON_UNIT, 
              T(1), L11_STAR_STAR, X1_STAR_STAR );
        }
        X1 = X1_STAR_STAR;
        //--------------------------------------------------------------------//

        SlideLockedPartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        SlidePartitionUp
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );
    }
}

template<typename T>
inline void FrontLowerMultiplyTranspose
( int diagOff, const Matrix<T>& L, Matrix<T>& X, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontLowerMultiplyTranspose");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal solve:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
        if( diagOff > 0 )
            LogicError("Diagonal offsets cannot be positive");
    )
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    // Danger, Will Robinson!
    Matrix<T>* LMod = const_cast<Matrix<T>*>(&L);
    Matrix<T> LT, LB;
    PartitionDown( *LMod, LT, LB, L.Width() );
    Matrix<T> XT, XB;
    PartitionDown( X, XT, XB, L.Width() );
    if( diagOff == 0 )
    {
        Trmm( LEFT, LOWER, orientation, NON_UNIT, T(1), LT, XT );
    }
    else
    {
        std::vector<Matrix<T>> diagonals;
        internal::ModifyForTrmm( LT, diagOff, diagonals );
        Trmm( LEFT, LOWER, orientation, NON_UNIT, T(1), LT, XT );
        internal::ReplaceAfterTrmm( LT, diagOff, diagonals );
    }

    Gemm( orientation, NORMAL, T(1), LB, XB, T(1), XT );
}

// TODO: Simplify this implementation
template<typename T>
inline void FrontLowerMultiplyTranspose
( int diagOff, const DistMatrix<T,VC,STAR>& L, DistMatrix<T,VC,STAR>& X,
  bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontLowerMultiplyTranspose");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal multiply:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
        if( L.ColAlign() != X.ColAlign() )
            LogicError("L and X are assumed to be aligned");
        if( diagOff > 0 )
            LogicError("Diagonal offsets cannot be positive");
    )
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T,VC,STAR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T,VC,STAR> XT(g),  X0(g),
                          XB(g),  X1(g),
                                  X2(g);

    // Temporary distributions
    DistMatrix<T,STAR,STAR> X1_STAR_STAR(g);
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,STAR,STAR> Z1_STAR_STAR(g);

    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( LTL.Width() < L.Width() )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,   L00, /**/ L01, L02,
         /*************/  /******************/
               /**/        L10, /**/ L11, L12,
          LBL, /**/ LBR,   L20, /**/ L21, L22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2, L11.Height() );

        //--------------------------------------------------------------------//
        X1_STAR_STAR = X1; // Can this be avoided?
        L11_STAR_STAR = L11;
        if( diagOff == 0 )
        {
            LocalTrmm
            ( LEFT, LOWER, orientation, NON_UNIT, 
              T(1), L11_STAR_STAR, X1_STAR_STAR );
        }
        else
        {
            internal::ModifyForTrmm( L11_STAR_STAR, diagOff );
            LocalTrmm
            ( LEFT, LOWER, orientation, NON_UNIT, 
              T(1), L11_STAR_STAR, X1_STAR_STAR );
        }
        X1 = X1_STAR_STAR;

        Z1_STAR_STAR.Resize( X1.Height(), X1.Width() );
        LocalGemm( orientation, NORMAL, T(1), L21, X2, T(0), Z1_STAR_STAR );
        X1.SumScatterUpdate( T(1), Z1_STAR_STAR );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
}

} // namespace El

#endif // ifndef EL_SPARSEDIRECT_NUMERIC_LOWERMULTIPLY_FRONT_HPP
