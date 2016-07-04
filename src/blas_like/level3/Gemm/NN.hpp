/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace gemm {

// Cannon's algorithm
template<typename T>
void Cannon_NN
( T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& BPre,
        AbstractDistMatrix<T>& CPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( APre, BPre, CPre );
      if( APre.Height() != CPre.Height() || BPre.Width() != CPre.Width() ||
          APre.Width() != BPre.Height() )
          LogicError
          ("Nonconformal matrices:\n",
           DimsString(APre,"A"),"\n",DimsString(BPre,"B"),"\n",
           DimsString(CPre,"C"));
    )
    const Grid& g = APre.Grid();
    if( g.Height() != g.Width() )
        LogicError("Process grid must be square for Cannon's");

    // Force A, B, and C to be in [MC,MR] distributions aligned with C
    DistMatrixReadWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& C = CProx.Get();

    ElementalProxyCtrl ctrlA, ctrlB;
    ctrlA.colConstrain = true; ctrlA.colAlign = C.ColAlign();
    ctrlB.rowConstrain = true; ctrlB.rowAlign = C.RowAlign();
    
    DistMatrixReadProxy<T,T,MC,MR> AProx( APre, ctrlA );
    DistMatrixReadProxy<T,T,MC,MR> BProx( BPre, ctrlB );
    auto& A = AProx.GetLocked();
    auto& B = BProx.GetLocked();

    const Int row = g.Row();
    const Int col = g.Col();
    const Int pSqrt = g.Height();
    mpi::Comm rowComm = g.RowComm();
    mpi::Comm colComm = g.ColComm(); 
    if( A.Width() % pSqrt != 0 )
        LogicError("For now, width(A) must be integer multiple of sqrt(p)");

    // Load the initial A and B packages (may want to transpose B...)
    const Int localHeightA = A.LocalHeight();
    const Int localHeightB = B.LocalHeight();
    const Int localWidthA = A.LocalWidth();
    const Int localWidthB = B.LocalWidth();
    Matrix<T> pkgA(localHeightA,localWidthA,localHeightA), 
              pkgB(localHeightB,localWidthB,localHeightB);
    for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
        MemCopy
        ( pkgA.Buffer(0,jLoc), A.LockedBuffer(0,jLoc), localHeightA );
    for( Int jLoc=0; jLoc<localWidthB; ++jLoc )
        MemCopy
        ( pkgB.Buffer(0,jLoc), B.LockedBuffer(0,jLoc), localHeightB );

    // Perform the initial circular shifts so that our A and B packages align
    const Int rowShiftA = A.RowShift();
    const Int colShiftB = B.ColShift();
    const Int leftInitA  = Mod(col-colShiftB,pSqrt);
    const Int rightInitA = Mod(col+colShiftB,pSqrt);
    const Int aboveInitB = Mod(row-rowShiftA,pSqrt);
    const Int belowInitB = Mod(row+rowShiftA,pSqrt);
    const Int pkgSizeA = localHeightA*localWidthA;
    const Int pkgSizeB = localHeightB*localWidthB;
    mpi::SendRecv( pkgA.Buffer(), pkgSizeA, leftInitA, rightInitA, rowComm );
    mpi::SendRecv( pkgB.Buffer(), pkgSizeB, aboveInitB, belowInitB, colComm );

    // Now begin the data flow
    const Int aboveRow = Mod(row-1,pSqrt);
    const Int belowRow = Mod(row+1,pSqrt);
    const Int leftCol  = Mod(col-1,pSqrt);
    const Int rightCol = Mod(col+1,pSqrt);
    for( Int q=0; q<pSqrt; ++q )
    {
        Gemm( NORMAL, NORMAL, alpha, pkgA, pkgB, T(1), C.Matrix() );
        if( q != pSqrt-1 )
        {
            mpi::SendRecv
            ( pkgA.Buffer(), pkgSizeA, leftCol, rightCol, rowComm );
            mpi::SendRecv
            ( pkgB.Buffer(), pkgSizeB, aboveRow, belowRow, colComm );
        }
    }
}

// Normal Normal Gemm that avoids communicating the matrix A
template<typename T>
void SUMMA_NNA
( T alpha,
  const AbstractDistMatrix<T>& APre, 
  const AbstractDistMatrix<T>& BPre,
        AbstractDistMatrix<T>& CPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( APre, BPre, CPre );
      if( APre.Height() != CPre.Height() || BPre.Width() != CPre.Width() ||
          APre.Width() != BPre.Height() )
          LogicError
          ("Nonconformal matrices:\n",
           DimsString(APre,"A"),"\n",DimsString(BPre,"B"),"\n",
           DimsString(CPre,"C"));
    )
    const Int n = CPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> AProx( APre );
    DistMatrixReadProxy<T,T,MC,MR> BProx( BPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.GetLocked();
    auto& C = CProx.Get();

    // Temporary distributions
    DistMatrix<T,VR,STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR> B1Trans_STAR_MR(g);
    DistMatrix<T,MC,STAR> D1_MC_STAR(g);

    B1_VR_STAR.AlignWith( A );
    B1Trans_STAR_MR.AlignWith( A );
    D1_MC_STAR.AlignWith( A );

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        auto B1 = B( ALL, IR(k,k+nb) );
        auto C1 = C( ALL, IR(k,k+nb) );

        // D1[MC,*] := alpha A[MC,MR] B1[MR,*]
        B1_VR_STAR = B1;
        Transpose( B1_VR_STAR, B1Trans_STAR_MR );
        LocalGemm( NORMAL, TRANSPOSE, alpha, A, B1Trans_STAR_MR, D1_MC_STAR );

        // C1[MC,MR] += scattered result of D1[MC,*] summed over grid rows
        AxpyContract( T(1), D1_MC_STAR, C1 );
    }
}

// Normal Normal Gemm that avoids communicating the matrix B
template<typename T>
void SUMMA_NNB
( T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& BPre,
        AbstractDistMatrix<T>& CPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( APre, BPre, CPre );
      if( APre.Height() != CPre.Height() || BPre.Width() != CPre.Width() ||
          APre.Width() != BPre.Height() )
          LogicError
          ("Nonconformal matrices:\n",
           DimsString(APre,"A"),"\n",DimsString(BPre,"B"),"\n",
           DimsString(CPre,"C"));
    )
    const Int m = CPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> AProx( APre );
    DistMatrixReadProxy<T,T,MC,MR> BProx( BPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.GetLocked();
    auto& C = CProx.Get();

    // Temporary distributions
    DistMatrix<T,STAR,MC> A1_STAR_MC(g);
    DistMatrix<T,MR,STAR> D1Trans_MR_STAR(g);

    A1_STAR_MC.AlignWith( B );
    D1Trans_MR_STAR.AlignWith( B );

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);
        auto A1 = A( IR(k,k+nb), ALL );
        auto C1 = C( IR(k,k+nb), ALL );

        // D1^T[MR,* ] := alpha B^T[MR,MC] A1^T[MC,* ]
        A1_STAR_MC = A1;
        LocalGemm
        ( TRANSPOSE, TRANSPOSE, alpha, B, A1_STAR_MC, D1Trans_MR_STAR );

        TransposeAxpyContract( T(1), D1Trans_MR_STAR, C1 );
    }
}

// Normal Normal Gemm that avoids communicating the matrix C
template<typename T>
void SUMMA_NNC
( T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& BPre,
        AbstractDistMatrix<T>& CPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( APre, BPre, CPre );
      if( APre.Height() != CPre.Height() || BPre.Width() != CPre.Width() ||
          APre.Width() != BPre.Height() )
          LogicError
          ("Nonconformal matrices:\n",
           DimsString(APre,"A"),"\n",DimsString(BPre,"B"),"\n",
           DimsString(CPre,"C"));
    )
    const Int sumDim = APre.Width();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> AProx( APre );
    DistMatrixReadProxy<T,T,MC,MR> BProx( BPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.GetLocked();
    auto& C = CProx.Get();

    // Temporary distributions
    DistMatrix<T,MC,STAR> A1_MC_STAR(g);
    DistMatrix<T,MR,STAR> B1Trans_MR_STAR(g); 

    A1_MC_STAR.AlignWith( C );
    B1Trans_MR_STAR.AlignWith( C );

    for( Int k=0; k<sumDim; k+=bsize )
    {
        const Int nb = Min(bsize,sumDim-k);
        auto A1 = A( ALL,        IR(k,k+nb) );
        auto B1 = B( IR(k,k+nb), ALL        );

        // C[MC,MR] += alpha A1[MC,*] (B1^T[MR,*])^T
        //           = alpha A1[MC,*] B1[*,MR]
        A1_MC_STAR = A1; 
        Transpose( B1, B1Trans_MR_STAR );
        LocalGemm
        ( NORMAL, TRANSPOSE, alpha, A1_MC_STAR, B1Trans_MR_STAR, T(1), C );
    }
}

// Normal Normal Gemm for panel-panel dot products
template<typename T>
void SUMMA_NNDot
( T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& BPre,
        AbstractDistMatrix<T>& CPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( APre, BPre, CPre );
      if( APre.Height() != CPre.Height() || BPre.Width() != CPre.Width() ||
          APre.Width() != BPre.Height() )
          LogicError
          ("Nonconformal matrices:\n",
           DimsString(APre,"A"),"\n",DimsString(BPre,"B"),"\n",
           DimsString(CPre,"C"));
    )
    const Int m = CPre.Height();
    const Int n = CPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> AProx( APre );
    DistMatrixReadProxy<T,T,MC,MR> BProx( BPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.GetLocked();
    auto& C = CProx.Get();

    if( A.Height() > B.Width() )
    {
        // Temporary distributions
        DistMatrix<T,STAR,VC> A1_STAR_VC(g);
        DistMatrix<T,VC,STAR> B1_VC_STAR(g);
        DistMatrix<T,STAR,STAR> C11_STAR_STAR(g);

        for( Int kOuter=0; kOuter<m; kOuter+=bsize )
        {
            const Int nbOuter = Min(bsize,m-kOuter);
            const Range<Int> indOuter( kOuter, kOuter+nbOuter );

            auto A1 = A( indOuter, ALL );

            A1_STAR_VC = A1; 
            B1_VC_STAR.AlignWith( A1_STAR_VC );

            for( Int kInner=0; kInner<n; kInner+=bsize )
            {
                const Int nbInner = Min(bsize,n-kInner);
                const Range<Int> indInner( kInner, kInner+nbInner );

                auto B1  = B( ALL,      indInner );
                auto C11 = C( indOuter, indInner );

                B1_VC_STAR = B1;
                LocalGemm
                ( NORMAL, NORMAL, 
                  alpha, A1_STAR_VC, B1_VC_STAR, C11_STAR_STAR );

                AxpyContract( T(1), C11_STAR_STAR, C11 );
            }
        }
    }
    else
    {
        // Temporary distributions
        DistMatrix<T,STAR,VR> A1_STAR_VR(g);
        DistMatrix<T,VR,STAR> B1_VR_STAR(g);
        DistMatrix<T,STAR,STAR> C11_STAR_STAR(g);

        for( Int kOuter=0; kOuter<n; kOuter+=bsize )
        {
            const Int nbOuter = Min(bsize,n-kOuter);
            const Range<Int> indOuter( kOuter, kOuter+nbOuter );

            auto B1 = B( ALL, indOuter );

            B1_VR_STAR = B1;
            A1_STAR_VR.AlignWith( B1_VR_STAR );

            for( Int kInner=0; kInner<m; kInner+=bsize )
            {
                const Int nbInner = Min(bsize,m-kInner);
                const Range<Int> indInner( kInner, kInner+nbInner );

                auto A1  = A( indInner, ALL      );
                auto C11 = C( indInner, indOuter );

                A1_STAR_VR = A1;
                LocalGemm
                ( NORMAL, NORMAL, 
                  alpha, A1_STAR_VR, B1_VR_STAR, C11_STAR_STAR );
                AxpyContract( T(1), C11_STAR_STAR, C11 );
            }
        }
    }
}

template<typename T>
void SUMMA_NN
( T alpha,
  const AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& B,
        AbstractDistMatrix<T>& C,
  GemmAlgorithm alg=GEMM_DEFAULT )
{
    DEBUG_CSE
    const Int m = C.Height();
    const Int n = C.Width();
    const Int sumDim = A.Width();
    const double weightTowardsC = 2.;
    const double weightAwayFromDot = 10.;

    switch( alg )
    {
    case GEMM_DEFAULT:
        if( weightAwayFromDot*m <= sumDim && weightAwayFromDot*n <= sumDim )
            SUMMA_NNDot( alpha, A, B, C );
        else if( m <= n && weightTowardsC*m <= sumDim )
            SUMMA_NNB( alpha, A, B, C );    
        else if( n <= m && weightTowardsC*n <= sumDim )
            SUMMA_NNA( alpha, A, B, C );
        else
            SUMMA_NNC( alpha, A, B, C );
        break;
    case GEMM_SUMMA_A:   SUMMA_NNA( alpha, A, B, C ); break;
    case GEMM_SUMMA_B:   SUMMA_NNB( alpha, A, B, C ); break;
    case GEMM_SUMMA_C:   SUMMA_NNC( alpha, A, B, C ); break;
    case GEMM_SUMMA_DOT: SUMMA_NNDot( alpha, A, B, C ); break;
    default: LogicError("Unsupported Gemm option");
    }
}

} // namespace gemm
} // namespace El
