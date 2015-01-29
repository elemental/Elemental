/*
   Copyright (c) 2009-2015, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SPARSEDIRECT_NUMERIC_LOWERSOLVE_FRONT_HPP
#define EL_SPARSEDIRECT_NUMERIC_LOWERSOLVE_FRONT_HPP

namespace El {

template<typename F>
inline void FrontLowerForwardSolve( const Matrix<F>& L, Matrix<F>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontLowerForwardSolve");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal solve:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
    )
    Matrix<F> LT, LB, XT, XB;
    LockedPartitionDown( L, LT, LB, L.Width() );
    PartitionDown( X, XT, XB, L.Width() );

    Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), LT, XT, true );
    Gemm( NORMAL, NORMAL, F(-1), LB, XT, F(1), XB );
}

template<typename F>
inline void FrontIntraPivLowerForwardSolve
( const Matrix<F>& L, const Matrix<Int>& p, Matrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("FrontIntraPivLowerForwardSolve"))
    Matrix<F> XT, XB;
    PartitionDown( X, XT, XB, L.Width() );
    PermuteRows( XT, p );
    FrontLowerForwardSolve( L, X );
}

namespace internal {

// TODO: Compress implementation
template<typename F>
inline void ForwardMany
( const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        FrontLowerForwardSolve( L.LockedMatrix(), X.Matrix() );
        return;
    }

    // Matrix views
    DistMatrix<F,VC,STAR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<F,VC,STAR> XT(g),  X0(g),
                          XB(g),  X1(g),
                                  X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> X1_STAR_STAR(g);

    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( LTL.Width() < L.Width() )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2, L11.Height() );

        //--------------------------------------------------------------------//
        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[VC,* ]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VC,* ]

        // X1[* ,* ] := (L11[* ,* ])^-1 X1[* ,* ]
        LocalTrsm
        ( LEFT, LOWER, NORMAL, NON_UNIT, 
          F(1), L11_STAR_STAR, X1_STAR_STAR, true );
        X1 = X1_STAR_STAR;

        // X2[VC,* ] -= L21[VC,* ] X1[* ,* ]
        LocalGemm( NORMAL, NORMAL, F(-1), L21, X1_STAR_STAR, F(1), X2 );
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

template<typename F>
void FormDiagonalBlocks
( const DistMatrix<F,VC,STAR>& L, DistMatrix<F,STAR,STAR>& D, bool conjugate )
{
    const Grid& g = L.Grid();

    const Int height = L.Width();
    const Int blocksize = Blocksize();

    const int commRank = g.VCRank();
    const int commSize = g.Size();

    const Int localHeight = Length(height,commRank,commSize);
    const Int maxLocalHeight = MaxLength(height,commSize);
    const Int portionSize = maxLocalHeight*blocksize;

    std::vector<F> sendBuffer( portionSize );
    const Int colShift = L.ColShift();
    const Int LLDim = L.LDim();
    const F* LBuffer = L.LockedBuffer();
    if( conjugate )
    {
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*commSize;
            const Int block = i / blocksize;
            const Int jStart = block*blocksize;
            const Int b = std::min(height-jStart,blocksize);
            for( Int jOff=0; jOff<b; ++jOff )
                sendBuffer[iLoc*blocksize+jOff] = 
                    Conj(LBuffer[iLoc+(jStart+jOff)*LLDim]);
        }
    }
    else
    {
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*commSize;
            const Int block = i / blocksize;
            const Int jStart = block*blocksize;
            const Int b = std::min(height-jStart,blocksize);
            for( Int jOff=0; jOff<b; ++jOff )
                sendBuffer[iLoc*blocksize+jOff] = 
                    LBuffer[iLoc+(jStart+jOff)*LLDim];
        }
    }

    std::vector<F> recvBuffer( portionSize*commSize );
    mpi::AllGather
    ( &sendBuffer[0], portionSize, &recvBuffer[0], portionSize, g.VCComm() );
    SwapClear( sendBuffer );
    
    D.Resize( blocksize, height );
    F* DBuffer = D.Buffer();
    const Int DLDim = D.LDim();
    for( Int proc=0; proc<commSize; ++proc )
    {
        const F* procRecv = &recvBuffer[proc*portionSize];
        const Int procLocalHeight = Length(height,proc,commSize);
        for( Int iLoc=0; iLoc<procLocalHeight; ++iLoc )
        {
            const Int i = proc + iLoc*commSize;
            for( Int jOff=0; jOff<blocksize; ++jOff )
                DBuffer[jOff+i*DLDim] = procRecv[jOff+iLoc*blocksize];
        }
    }
}

template<typename F>
void AccumulateRHS( const DistMatrix<F,VC,STAR>& X, DistMatrix<F,STAR,STAR>& Z )
{
    const Int height = X.Height();
    const Int width = X.Width();
    Z.Empty();
    Zeros( Z, height, width );

    const Int localHeight = X.LocalHeight();
    const Int colShift = X.ColShift();
    const int commSize = X.Grid().Size();
    const F* XBuffer = X.LockedBuffer();
    F* ZBuffer = Z.Buffer();
    const Int XLDim = X.LDim();
    const Int ZLDim = Z.LDim();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = colShift + iLoc*commSize;
        for( Int j=0; j<width; ++j )
            ZBuffer[i+j*ZLDim] = XBuffer[iLoc+j*XLDim];
    }
    mpi::AllReduce( ZBuffer, ZLDim*width, mpi::SUM, X.Grid().VCComm() );
}

// TODO: Compress implementation
template<typename F>
void ForwardSingle( const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        FrontLowerForwardSolve( L.LockedMatrix(), X.Matrix() );
        return;
    }

    // Matrix views
    DistMatrix<F,VC,STAR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<F,VC,STAR> XT(g),  X0(g),
                          XB(g),  X1(g),
                                  X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11Trans_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> X1_STAR_STAR(g);

    DistMatrix<F,STAR,STAR> D(g);
    FormDiagonalBlocks( L, D, false );

    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( LTL.Width() < L.Width() )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2, L11.Height() );

        LockedView
        ( L11Trans_STAR_STAR, D, 0, L00.Height(), L11.Height(), L11.Height() );

        //--------------------------------------------------------------------//
        AccumulateRHS( X1, X1_STAR_STAR ); // X1[* ,* ] <- X1[VC,* ]

        // X1[* ,* ] := (L11[* ,* ])^-1 X1[* ,* ]
        LocalTrsm
        ( LEFT, UPPER, TRANSPOSE, NON_UNIT, 
          F(1), L11Trans_STAR_STAR, X1_STAR_STAR, true );
        X1 = X1_STAR_STAR;

        // X2[VC,* ] -= L21[VC,* ] X1[* ,* ]
        LocalGemm( NORMAL, NORMAL, F(-1), L21, X1_STAR_STAR, F(1), X2 );
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

template<typename F>
void BackwardMany
( const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X,
  bool conjugate=false )
{
    // TODO: Replace this with modified inline code?
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    trsm::LLTSmall( orientation, NON_UNIT,  L, X, true );
}

// TODO: Compress implementation
template<typename F>
void BackwardSingle
( const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X,
  bool conjugate=false )
{
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<F,VC,STAR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<F,VC,STAR> XT(g),  X0(g),
                          XB(g),  X1(g),
                                  X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> L11Trans_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> Z1_STAR_STAR(g);

    DistMatrix<F,STAR,STAR> D(g);
    FormDiagonalBlocks( L, D, conjugate );
    const Int blocksize = Blocksize();
    const Int firstBlocksize = 
        ( L.Height()%blocksize==0 ?
          blocksize :
          L.Height()%blocksize );

    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    // Start the algorithm
    Int b = firstBlocksize;
    LockedPartitionUpDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionUp
    ( X, XT,
         XB, 0 );
    while( XT.Height() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22, b );

        RepartitionUp
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2, b );

        LockedView( L11Trans_STAR_STAR, D, 0, L00.Height(), b, b );

        //--------------------------------------------------------------------//
        // X1 -= L21' X2
        LocalGemm( orientation, NORMAL, F(-1), L21, X2, Z1_STAR_STAR );
        AddInLocalData( X1, Z1_STAR_STAR );
        El::AllReduce( Z1_STAR_STAR, X1.DistComm() );

        // X1 := L11^-1 X1
        LocalTrsm
        ( LEFT, UPPER, NORMAL, UNIT, F(1), L11Trans_STAR_STAR, Z1_STAR_STAR );
        X1 = Z1_STAR_STAR;
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

        b = blocksize;
    }
}

} // namespace internal

template<typename F>
inline void FrontLowerForwardSolve
( const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X,
  bool singleL11AllGather=true )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontLowerForwardSolve");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal solve:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
        if( L.ColAlign() != X.ColAlign() )
            LogicError("L and X are assumed to be aligned");
    )
    if( singleL11AllGather )
        internal::ForwardSingle( L, X );
    else
        internal::ForwardMany( L, X );
}

template<typename F>
inline void FrontIntraPivLowerForwardSolve
( const DistMatrix<F,VC,STAR>& L, const DistMatrix<Int,VC,STAR>& p, 
  DistMatrix<F,VC,STAR>& X, bool singleL11AllGather=true )
{
    DEBUG_ONLY(CallStackEntry cse("FrontIntraPivLowerForwardSolve"))

    // TODO: Cache the send and recv data for the pivots to avoid p[*,*]
    const Grid& g = L.Grid();
    DistMatrix<F,VC,STAR> XT(g), XB(g);
    PartitionDown( X, XT, XB, L.Width() );
    PermuteRows( XT, p );

    FrontLowerForwardSolve( L, X, singleL11AllGather );
}

template<typename F>
inline void FrontLowerForwardSolve( const DistMatrix<F>& L, DistMatrix<F>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontLowerForwardSolve");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal solve:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
    )
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        FrontLowerForwardSolve( L.LockedMatrix(), X.Matrix() );
        return;
    }

    // Separate the top and bottom portions of X and L
    const Int snSize = L.Width();
    DistMatrix<F> LT(g), LB(g), XT(g), XB(g);
    LockedPartitionDown( L, LT, LB, snSize );
    PartitionDown( X, XT, XB, snSize );

    // XT := LT XT
    // TODO: Replace with TrsmLLNMedium?
    Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), LT, XT );

    // XB := XB - LB XT
    Gemm( NORMAL, NORMAL, F(-1), LB, XT, F(1), XB );
}

template<typename F>
inline void FrontIntraPivLowerForwardSolve
( const DistMatrix<F>& L, const DistMatrix<Int,VC,STAR>& p, DistMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("FrontIntraPivLowerForwardSolve"))

    // TODO: Cache the send and recv data for the pivots to avoid p[*,*]
    const Grid& g = L.Grid();
    DistMatrix<F> XT(g), XB(g);
    PartitionDown( X, XT, XB, L.Width() );
    PermuteRows( XT, p );

    FrontLowerForwardSolve( L, X );
}

template<typename F>
inline void FrontFastLowerForwardSolve
( const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontFastLowerForwardSolve");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal solve:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
        if( L.ColAlign() != X.ColAlign() )
            LogicError("L and X are assumed to be aligned");
    )
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        FrontLowerForwardSolve( L.LockedMatrix(), X.Matrix() );
        return;
    }

    // Separate the top and bottom portions of X and L
    const int snSize = L.Width();
    DistMatrix<F,VC,STAR> LT(g), LB(g), XT(g), XB(g);
    LockedPartitionDown( L, LT, LB, snSize );
    PartitionDown( X, XT, XB, snSize );

    // XT := LT XT
    DistMatrix<F,STAR,STAR> XT_STAR_STAR( XT );
    LocalGemm( NORMAL, NORMAL, F(1), LT, XT_STAR_STAR, F(0), XT );

    // XB := XB - LB XT
    if( LB.Height() != 0 )
    {
        XT_STAR_STAR = XT;
        LocalGemm( NORMAL, NORMAL, F(-1), LB, XT_STAR_STAR, F(1), XB );
    }
}

template<typename F>
inline void FrontFastIntraPivLowerForwardSolve
( const DistMatrix<F,VC,STAR>& L, const DistMatrix<Int,VC,STAR>& p,
  DistMatrix<F,VC,STAR>& X )
{
    DEBUG_ONLY(CallStackEntry cse("FrontFastIntraPivLowerForwardSolve"))

    // TODO: Cache the send and recv data for the pivots to avoid p[*,*]
    const Grid& g = L.Grid();
    DistMatrix<F,VC,STAR> XT(g), XB(g);
    PartitionDown( X, XT, XB, L.Width() );
    PermuteRows( XT, p );

    FrontFastLowerForwardSolve( L, X );
}

template<typename F>
inline void FrontFastLowerForwardSolve
( const DistMatrix<F>& L, DistMatrix<F,VC,STAR>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontFastLowerForwardSolve");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal solve:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
    )
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        FrontLowerForwardSolve( L.LockedMatrix(), X.Matrix() );
        return;
    }

    // Separate the top and bottom portions of X and L
    const int snSize = L.Width();
    DistMatrix<F> LT(g), LB(g);
    LockedPartitionDown( L, LT, LB, snSize );
    DistMatrix<F,VC,STAR> XT(g), XB(g);
    PartitionDown( X, XT, XB, snSize );

    // Get ready for the local multiply
    DistMatrix<F,MR,STAR> XT_MR_STAR(g);
    XT_MR_STAR.AlignWith( LT );

    {
        // ZT[MC,* ] := LT[MC,MR] XT[MR,* ], 
        DistMatrix<F,MC,STAR> ZT_MC_STAR(g);
        ZT_MC_STAR.AlignWith( LT );
        XT_MR_STAR = XT;
        LocalGemm( NORMAL, NORMAL, F(1), LT, XT_MR_STAR, ZT_MC_STAR );

        Contract( ZT_MC_STAR, XT );
    }

    if( LB.Height() != 0 )
    {
        // Set up for the local multiply
        XT_MR_STAR = XT;

        // ZB[MC,* ] := LB[MC,MR] XT[MR,* ]
        DistMatrix<F,MC,STAR> ZB_MC_STAR(g);
        ZB_MC_STAR.AlignWith( LB );
        LocalGemm( NORMAL, NORMAL, F(-1), LB, XT_MR_STAR, ZB_MC_STAR );

        // XB[VC,* ] += ZB[MC,* ]
        AxpyContract( F(1), ZB_MC_STAR, XB );
    }
}

template<typename F>
inline void FrontFastIntraPivLowerForwardSolve
( const DistMatrix<F>& L, const DistMatrix<Int,VC,STAR>& p,
  DistMatrix<F,VC,STAR>& X )
{
    DEBUG_ONLY(CallStackEntry cse("FrontFastIntraPivLowerForwardSolve"))

    // TODO: Cache the send and recv data for the pivots to avoid p[*,*]
    const Grid& g = L.Grid();
    DistMatrix<F,VC,STAR> XT(g), XB(g);
    PartitionDown( X, XT, XB, L.Width() );
    PermuteRows( XT, p );

    FrontFastLowerForwardSolve( L, X );
}

template<typename F>
inline void FrontFastLowerForwardSolve
( const DistMatrix<F>& L, DistMatrix<F>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontFastLowerForwardSolve");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal solve:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
    )
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        FrontLowerForwardSolve( L.LockedMatrix(), X.Matrix() );
        return;
    }

    // Separate the top and bottom portions of X and L
    const int snSize = L.Width();
    DistMatrix<F> LT(g), LB(g), XT(g), XB(g);
    LockedPartitionDown( L, LT, LB, snSize );
    PartitionDown( X, XT, XB, snSize );

    // XT := LT XT
    DistMatrix<F> YT( XT );
    Gemm( NORMAL, NORMAL, F(1), LT, YT, F(0), XT );

    // XB := XB - LB XT
    Gemm( NORMAL, NORMAL, F(-1), LB, XT, F(1), XB );
}

template<typename F>
inline void FrontFastIntraPivLowerForwardSolve
( const DistMatrix<F>& L, const DistMatrix<Int,VC,STAR>& p, DistMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("FrontFastIntraPivLowerForwardSolve"))

    // TODO: Cache the send and recv data for the pivots to avoid p[*,*]
    const Grid& g = L.Grid();
    DistMatrix<F> XT(g), XB(g);
    PartitionDown( X, XT, XB, L.Width() );
    PermuteRows( XT, p );

    FrontFastLowerForwardSolve( L, X );
}

template<typename F>
inline void FrontLowerBackwardSolve
( const Matrix<F>& L, Matrix<F>& X, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontLowerBackwardSolve");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal solve:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
    )
    Matrix<F> LT, LB, XT, XB;
    LockedPartitionDown( L, LT, LB, L.Width() );
    PartitionDown( X, XT, XB, L.Width() );

    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    Gemm( orientation, NORMAL, F(-1), LB, XB, F(1), XT );
    Trsm( LEFT, LOWER, orientation, NON_UNIT, F(1), LT, XT, true );
}

template<typename F>
inline void FrontIntraPivLowerBackwardSolve
( const Matrix<F>& L, const Matrix<Int>& p, Matrix<F>& X, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("FrontIntraPivLowerBackwardSolve"))
    FrontLowerBackwardSolve( L, X, conjugate );
    Matrix<F> XT, XB;
    PartitionDown( X, XT, XB, L.Width() );
    InversePermuteRows( XT, p );
}

template<typename F>
inline void FrontLowerBackwardSolve
( const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X,
  bool conjugate, bool singleL11AllGather=true )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontLowerBackwardSolve");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal solve:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
        if( L.ColAlign() != X.ColAlign() )
            LogicError("L and X are assumed to be aligned");
    )
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        FrontLowerBackwardSolve( L.LockedMatrix(), X.Matrix(), conjugate );
        return;
    }

    DistMatrix<F,VC,STAR> LT(g), LB(g), XT(g), XB(g);
    LockedPartitionDown( L, LT, LB, L.Width() );
    PartitionDown( X, XT, XB, L.Width() );

    if( XB.Height() != 0 )
    {
        // Subtract off the parent updates
        DistMatrix<F,STAR,STAR> Z(g);
        const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
        LocalGemm( orientation, NORMAL, F(-1), LB, XB, Z );
        AxpyContract( F(1), Z, XT );
    }

    if( singleL11AllGather )
        internal::BackwardSingle( LT, XT, conjugate );
    else
        internal::BackwardMany( LT, XT, conjugate );
}

template<typename F>
inline void FrontIntraPivLowerBackwardSolve
( const DistMatrix<F,VC,STAR>& L, const DistMatrix<Int,VC,STAR>& p,
  DistMatrix<F,VC,STAR>& X, bool conjugate, bool singleL11AllGather=true )
{
    DEBUG_ONLY(CallStackEntry cse("FrontIntraPivLowerBackwardSolve"))

    FrontLowerBackwardSolve( L, X, conjugate, singleL11AllGather );

    // TODO: Cache the send and recv data for the pivots to avoid p[*,*]
    const Grid& g = L.Grid();
    DistMatrix<F,VC,STAR> XT(g), XB(g);
    PartitionDown( X, XT, XB, L.Width() );
    InversePermuteRows( XT, p );
}

template<typename F>
inline void FrontLowerBackwardSolve
( const DistMatrix<F>& L, DistMatrix<F>& X, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontLowerBackwardSolve");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal solve:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
    )
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        FrontLowerBackwardSolve( L.LockedMatrix(), X.Matrix(), conjugate );
        return;
    }

    DistMatrix<F> LT(g), LB(g), XT(g), XB(g);
    LockedPartitionDown( L, LT, LB, L.Width() );
    PartitionDown( X, XT, XB, L.Width() );

    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    Gemm( orientation, NORMAL, F(-1), LB, XB, F(1), XT );
    Trsm( LEFT, LOWER, orientation, NON_UNIT, F(1), LT, XT );
}

template<typename F>
inline void FrontIntraPivLowerBackwardSolve
( const DistMatrix<F>& L, const DistMatrix<Int,VC,STAR>& p,
  DistMatrix<F>& X, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("FrontIntraPivLowerBackwardSolve"))

    FrontLowerBackwardSolve( L, X, conjugate );

    // TODO: Cache the send and recv data for the pivots to avoid p[*,*]
    const Grid& g = L.Grid();
    DistMatrix<F> XT(g), XB(g);
    PartitionDown( X, XT, XB, L.Width() );
    InversePermuteRows( XT, p );
}

template<typename F>
inline void FrontFastLowerBackwardSolve
( const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X,
  bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontFastLowerBackwardSolve");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal solve:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
        if( L.ColAlign() != X.ColAlign() )
            LogicError("L and X are assumed to be aligned");
    )
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        FrontLowerBackwardSolve( L.LockedMatrix(), X.Matrix(), conjugate );
        return;
    }

    const int snSize = L.Width();
    DistMatrix<F,VC,STAR> LT(g), LB(g), XT(g), XB(g);
    LockedPartitionDown( L, LT, LB, snSize );
    PartitionDown( X, XT, XB, snSize );

    // XT := XT - LB^{T/H} XB
    DistMatrix<F,STAR,STAR> Z(g);
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    if( XB.Height() != 0 )
    {
        LocalGemm( orientation, NORMAL, F(-1), LB, XB, Z );
        AxpyContract( F(1), Z, XT );
    }

    // XT := LT^{T/H} XT
    LocalGemm( orientation, NORMAL, F(1), LT, XT, Z );
    Contract( Z, XT );
}

template<typename F>
inline void FrontFastIntraPivLowerBackwardSolve
( const DistMatrix<F,VC,STAR>& L, const DistMatrix<Int,VC,STAR>& p,
  DistMatrix<F,VC,STAR>& X, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("FrontFastIntraPivLowerBackwardSolve"))

    FrontFastLowerBackwardSolve( L, X, conjugate );

    // TODO: Cache the send and recv data for the pivots to avoid p[*,*]
    const Grid& g = L.Grid();
    DistMatrix<F,VC,STAR> XT(g), XB(g);
    PartitionDown( X, XT, XB, L.Width() );
    InversePermuteRows( XT, p );
}

template<typename F>
inline void FrontFastLowerBackwardSolve
( const DistMatrix<F>& L, DistMatrix<F,VC,STAR>& X, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontFastLowerBackwardSolve");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal solve:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
    )
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        FrontLowerBackwardSolve( L.LockedMatrix(), X.Matrix(), conjugate );
        return;
    }

    const int snSize = L.Width();
    DistMatrix<F> LT(g), LB(g);
    LockedPartitionDown( L, LT, LB, snSize );
    DistMatrix<F,VC,STAR> XT(g), XB(g);
    PartitionDown( X, XT, XB, snSize );

    DistMatrix<F,MR,STAR> ZT_MR_STAR( g );
    DistMatrix<F,VR,STAR> ZT_VR_STAR( g );
    ZT_MR_STAR.AlignWith( LB );
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    if( XB.Height() != 0 )
    {
        // ZT[MR,* ] := -(LB[MC,MR])^{T/H} XB[MC,* ]
        DistMatrix<F,MC,STAR> XB_MC_STAR( g );
        XB_MC_STAR.AlignWith( LB );
        XB_MC_STAR = XB;
        LocalGemm( orientation, NORMAL, F(-1), LB, XB_MC_STAR, ZT_MR_STAR );

        Contract( ZT_MR_STAR, ZT_VR_STAR );

        // ZT[VC,* ] := ZT[VR,* ]
        DistMatrix<F,VC,STAR> ZT_VC_STAR( g );
        ZT_VC_STAR.AlignWith( XT );
        ZT_VC_STAR = ZT_VR_STAR;

        // XT[VC,* ] += ZT[VC,* ]
        Axpy( F(1), ZT_VC_STAR, XT );
    }

    {
        // ZT[MR,* ] := (LT[MC,MR])^{T/H} XT[MC,* ]
        DistMatrix<F,MC,STAR> XT_MC_STAR( g );
        XT_MC_STAR.AlignWith( LT );
        XT_MC_STAR = XT;
        LocalGemm( orientation, NORMAL, F(1), LT, XT_MC_STAR, ZT_MR_STAR );

        Contract( ZT_MR_STAR, ZT_VR_STAR );

        // XT[VC,* ] := ZT[VR,* ]
        XT = ZT_VR_STAR;
    }
}

template<typename F>
inline void FrontFastIntraPivLowerBackwardSolve
( const DistMatrix<F>& L, const DistMatrix<Int,VC,STAR>& p,
  DistMatrix<F,VC,STAR>& X, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("FrontFastIntraPivLowerBackwardSolve"))

    FrontFastLowerBackwardSolve( L, X, conjugate );

    // TODO: Cache the send and recv data for the pivots to avoid p[*,*]
    const Grid& g = L.Grid();
    DistMatrix<F,VC,STAR> XT(g), XB(g);
    PartitionDown( X, XT, XB, L.Width() );
    InversePermuteRows( XT, p );
}

template<typename F>
inline void FrontFastLowerBackwardSolve
( const DistMatrix<F>& L, DistMatrix<F>& X, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontFastLowerBackwardSolve");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal solve:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
    )
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        FrontLowerBackwardSolve( L.LockedMatrix(), X.Matrix(), conjugate );
        return;
    }

    const int snSize = L.Width();
    DistMatrix<F> LT(g), LB(g), XT(g), XB(g);
    LockedPartitionDown( L, LT, LB, snSize );
    PartitionDown( X, XT, XB, snSize );

    // XT := XT - LB^{T/H} XB
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    Gemm( orientation, NORMAL, F(-1), LB, XB, F(1), XT );

    // XT := LT^{T/H} XT
    DistMatrix<F> Z(XT.Grid());
    Gemm( orientation, NORMAL, F(1), LT, XT, Z );
    XT = Z;
}

template<typename F>
inline void FrontFastIntraPivLowerBackwardSolve
( const DistMatrix<F>& L, const DistMatrix<Int,VC,STAR>& p,
  DistMatrix<F>& X, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("FrontFastIntraPivLowerBackwardSolve"))

    FrontFastLowerBackwardSolve( L, X, conjugate );

    // TODO: Cache the send and recv data for the pivots to avoid p[*,*]
    const Grid& g = L.Grid();
    DistMatrix<F> XT(g), XB(g);
    PartitionDown( X, XT, XB, L.Width() );
    InversePermuteRows( XT, p );
}

} // namespace El

#endif // ifndef EL_SPARSEDIRECT_NUMERIC_LOWERSOLVE_FRONT_HPP
