/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRRK_LOCAL_HPP
#define EL_TRRK_LOCAL_HPP

#include <El/io.hpp>

namespace El {

namespace trrk {

#ifndef EL_RELEASE

template<typename T>
void EnsureConformal
( const DistMatrix<T,MC,STAR>& A, const DistMatrix<T>& C, string name )
{
    if( A.Height() != C.Height() || A.ColAlign() != C.ColAlign() )
        LogicError(name," not conformal with C");
}

template<typename T>
void EnsureConformal
( const DistMatrix<T,STAR,MC>& A, const DistMatrix<T>& C, string name )
{
    if( A.Width() != C.Height() || A.RowAlign() != C.ColAlign() )
        LogicError(name," not conformal with C");
}

template<typename T>
void EnsureConformal
( const DistMatrix<T,MR,STAR>& A, const DistMatrix<T>& C, string name )
{
    if( A.Height() != C.Width() || A.ColAlign() != C.RowAlign() )
        LogicError(name," not conformal with C");
}

template<typename T>
void EnsureConformal
( const DistMatrix<T,STAR,MR>& A, const DistMatrix<T>& C, string name )
{
    if( A.Width() != C.Width() || A.RowAlign() != C.RowAlign() )
        LogicError(name," not conformal with C");
}

template<typename T,Dist UA,Dist VA,Dist UB,Dist VB>
void CheckInput
( const DistMatrix<T,UA,VA>& A, const DistMatrix<T,UB,VB>& B,
  const DistMatrix<T>& C )
{
    AssertSameGrids( A, B, C );
    EnsureConformal( A, C, "A" );
    EnsureConformal( B, C, "B" );
}

// Local C := A B + C
template<typename T>
void CheckInputNN( const Matrix<T>& A, const Matrix<T>& B, const Matrix<T>& C )
{
    if( A.Height() != C.Height() || B.Width()  != C.Width() ||
        A.Width()  != B.Height() || A.Height() != B.Width() )
        LogicError
        ("Nonconformal LocalTrrk:\n",
         DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
}

// Local C := A B^{T/H} + C
template<typename T>
void CheckInputNT
( Orientation orientationOfB,
  const Matrix<T>& A, const Matrix<T>& B, const Matrix<T>& C )
{
    if( orientationOfB == NORMAL )
        LogicError("B must be (Conjugate)Transpose'd");
    if( A.Height() != C.Height() || B.Height() != C.Width() ||
        A.Width()  != B.Width()  || A.Height() != B.Height() )
        LogicError
        ("Nonconformal LocalTrrk:\n",
         DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
}

// Local C := A^{T/H} B + C
template<typename T>
void CheckInputTN
( Orientation orientationOfA,
  const Matrix<T>& A, const Matrix<T>& B, const Matrix<T>& C )
{
    if( orientationOfA == NORMAL )
        LogicError("A must be (Conjugate)Transpose'd");
    if( A.Width() != C.Height() || B.Width() != C.Width() ||
        A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError
        ("Nonconformal LocalTrrk:\n",
         DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
}

// Local C := A^{T/H} B^{T/H} + C
template<typename T>
void CheckInputTT
( Orientation orientationOfA,
  Orientation orientationOfB,
  const Matrix<T>& A, const Matrix<T>& B, const Matrix<T>& C )
{
    if( orientationOfA == NORMAL )
        LogicError("A must be (Conjugate)Transpose'd");
    if( orientationOfB == NORMAL )
        LogicError("B must be (Conjugate)Transpose'd");
    if( A.Width() != C.Height() || B.Height() != C.Width() ||
        A.Height() != B.Width() || A.Width() != B.Height() )
        LogicError
        ("Nonconformal LocalTrrk:\n",
         DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
}

#endif // ifndef EL_RELEASE

// Local C := alpha A B + C
template<typename T>
void TrrkNNKernel
( UpperOrLower uplo, 
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
                 Matrix<T>& C )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(CheckInputNN( A, B, C ))

    const Int half = C.Height()/2;
    const auto indTL = IR(0,half);
    const auto indBR = IR(half,END);

    auto AT = A(indTL,ALL);
    auto AB = A(indBR,ALL);
    auto BL = B(ALL,indTL);
    auto BR = B(ALL,indBR);
    auto CTL = C(indTL,indTL);
    auto CBR = C(indBR,indBR);

    if( uplo == LOWER )
    {
        auto CBL = C(indBR,indTL);
        Gemm( NORMAL, NORMAL, alpha, AB, BL, T(1), CBL );
    }
    else
    {
        auto CTR = C(indTL,indBR);
        Gemm( NORMAL, NORMAL, alpha, AT, BR, T(1), CTR );
    }

    // TODO(poulson): Avoid the temporary copy
    Matrix<T> DTL;
    Gemm( NORMAL, NORMAL, alpha, AT, BL, DTL );
    AxpyTrapezoid( uplo, T(1), DTL, CTL );

    // TODO(poulson): Avoid the temporary copy
    Matrix<T> DBR;
    Gemm( NORMAL, NORMAL, alpha, AB, BR, DBR );
    AxpyTrapezoid( uplo, T(1), DBR, CBR );
}

// Distributed C := alpha A B + C
template<typename T>
void LocalTrrkKernel
( UpperOrLower uplo, 
  T alpha, const DistMatrix<T,MC,  STAR>& A,
           const DistMatrix<T,STAR,MR  >& B,
                 DistMatrix<T>& C )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(CheckInput( A, B, C ))
    const Grid& g = C.Grid();

    const Int half = C.Height()/2;
    const auto indTL = IR(0,half);
    const auto indBR = IR(half,END);

    auto AT = A(indTL,ALL);
    auto AB = A(indBR,ALL);
    auto BL = B(ALL,indTL);
    auto BR = B(ALL,indBR);
    auto CTL = C(indTL,indTL);
    auto CBR = C(indBR,indBR);

    if( uplo == LOWER )
    {
        auto CBL = C(indBR,indTL);
        LocalGemm( NORMAL, NORMAL, alpha, AB, BL, T(1), CBL );
    }
    else
    {
        auto CTR = C(indTL,indBR);
        LocalGemm( NORMAL, NORMAL, alpha, AT, BR, T(1), CTR );
    }

    // TODO(poulson): Avoid the temporary copy
    DistMatrix<T> DTL(g);
    DTL.AlignWith( CTL );
    LocalGemm( NORMAL, NORMAL, alpha, AT, BL, DTL );
    LocalAxpyTrapezoid( uplo, T(1), DTL, CTL );

    // TODO(poulson): Avoid the temporary copy
    DistMatrix<T> DBR(g);
    DBR.AlignWith( CBR );
    LocalGemm( NORMAL, NORMAL, alpha, AB, BR, DBR );
    LocalAxpyTrapezoid( uplo, T(1), DBR, CBR );
}

// Local C := alpha A B^{T/H} + C
template<typename T>
void TrrkNTKernel
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
                 Matrix<T>& C )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(CheckInputNT( orientationOfB, A, B, C ))

    const Int half = C.Height()/2;
    const auto indTL = IR(0,half);
    const auto indBR = IR(half,END);

    auto AT = A(indTL,ALL);
    auto AB = A(indBR,ALL);
    auto BT = B(indTL,ALL);
    auto BB = B(indBR,ALL);
    auto CTL = C(indTL,indTL);
    auto CBR = C(indBR,indBR);

    if( uplo == LOWER )
    {
        auto CBL = C(indBR,indTL);
        Gemm( NORMAL, orientationOfB, alpha, AB, BT, T(1), CBL );
    }
    else
    {
        auto CTR = C(indTL,indBR);
        Gemm( NORMAL, orientationOfB, alpha, AT, BB, T(1), CTR );
    }

    // TODO(poulson): Avoid the temporary copy
    Matrix<T> DTL;
    Gemm( NORMAL, orientationOfB, alpha, AT, BT, DTL );
    AxpyTrapezoid( uplo, T(1), DTL, CTL );

    // TODO(poulson): Avoid the temporary copy
    Matrix<T> DBR;
    Gemm( NORMAL, orientationOfB, alpha, AB, BB, DBR );
    AxpyTrapezoid( uplo, T(1), DBR, CBR );
}

// Distributed C := alpha A B^{T/H} + C
template<typename T>
void LocalTrrkKernel
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,STAR>& A,
           const DistMatrix<T,MR,STAR>& B,
                 DistMatrix<T>& C )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(CheckInput( A, B, C ))
    const Grid& g = C.Grid();

    const Int half = C.Height()/2;
    const auto indTL = IR(0,half);
    const auto indBR = IR(half,END);

    auto AT = A(indTL,ALL);
    auto AB = A(indBR,ALL);
    auto BT = B(indTL,ALL);
    auto BB = B(indBR,ALL);
    auto CTL = C(indTL,indTL);
    auto CBR = C(indBR,indBR);

    if( uplo == LOWER )
    {
        auto CBL = C(indBR,indTL);
        LocalGemm( NORMAL, orientationOfB, alpha, AB, BT, T(1), CBL );
    }
    else
    {
        auto CTR = C(indTL,indBR);
        LocalGemm( NORMAL, orientationOfB, alpha, AT, BB, T(1), CTR );
    }

    // TODO(poulson): Avoid the temporary copy
    DistMatrix<T> DTL(g);
    DTL.AlignWith( CTL );
    LocalGemm( NORMAL, orientationOfB, alpha, AT, BT, DTL );
    LocalAxpyTrapezoid( uplo, T(1), DTL, CTL );

    // TODO(poulson): Avoid the temporary copy
    DistMatrix<T> DBR(g);
    DBR.AlignWith( CBR );
    LocalGemm( NORMAL, orientationOfB, alpha, AB, BB, DBR );
    LocalAxpyTrapezoid( uplo, T(1), DBR, CBR );
}

// Local C := alpha A^{T/H} B + C
template<typename T>
void TrrkTNKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
                 Matrix<T>& C )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(CheckInputTN( orientationOfA, A, B, C ))

    const Int half = C.Height()/2;
    const auto indTL = IR(0,half);
    const auto indBR = IR(half,END);

    auto AL = A(ALL,indTL);
    auto AR = A(ALL,indBR);
    auto BL = B(ALL,indTL);
    auto BR = B(ALL,indBR);
    auto CTL = C(indTL,indTL);
    auto CBR = C(indBR,indBR);

    if( uplo == LOWER )
    {
        auto CBL = C(indBR,indTL);
        Gemm( orientationOfA, NORMAL, alpha, AR, BL, T(1), CBL );
    }
    else
    {
        auto CTR = C(indTL,indBR);
        Gemm( orientationOfA, NORMAL, alpha, AL, BR, T(1), CTR );
    }

    // TODO(poulson): Avoid the temporary copy
    Matrix<T> DTL;
    Gemm( orientationOfA, NORMAL, alpha, AL, BL, DTL );
    AxpyTrapezoid( uplo, T(1), DTL, CTL );

    // TODO(poulson): Avoid the temporary copy
    Matrix<T> DBR;
    Gemm( orientationOfA, NORMAL, alpha, AR, BR, DBR );
    AxpyTrapezoid( uplo, T(1), DBR, CBR );
}

// Distributed C := alpha A^{T/H} B + C
template<typename T>
void LocalTrrkKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T,STAR,MC>& A,
           const DistMatrix<T,STAR,MR>& B,
                 DistMatrix<T>& C )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(CheckInput( A, B, C ))
    const Grid& g = C.Grid();

    const Int half = C.Height()/2;
    const auto indTL = IR(0,half);
    const auto indBR = IR(half,END);

    auto AL = A(ALL,indTL);
    auto AR = A(ALL,indBR);
    auto BL = B(ALL,indTL);
    auto BR = B(ALL,indBR);
    auto CTL = C(indTL,indTL);
    auto CBR = C(indBR,indBR);

    if( uplo == LOWER )
    {
        auto CBL = C(indBR,indTL);
        LocalGemm( orientationOfA, NORMAL, alpha, AR, BL, T(1), CBL );
    }
    else
    {
        auto CTR = C(indTL,indBR);
        LocalGemm( orientationOfA, NORMAL, alpha, AL, BR, T(1), CTR );
    }

    // TODO(poulson): Avoid the temporary copy
    DistMatrix<T> DTL(g);
    DTL.AlignWith( CTL );
    LocalGemm( orientationOfA, NORMAL, alpha, AL, BL, DTL );
    LocalAxpyTrapezoid( uplo, T(1), DTL, CTL );

    // TODO(poulson): Avoid the temporary copy
    DistMatrix<T> DBR(g);
    DBR.AlignWith( CBR );
    LocalGemm( orientationOfA, NORMAL, alpha, AR, BR, DBR );
    LocalAxpyTrapezoid( uplo, T(1), DBR, CBR );
}

// Local C := alpha A^{T/H} B^{T/H} + C
template<typename T>
void TrrkTTKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
                 Matrix<T>& C )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(CheckInputTT( orientationOfA, orientationOfB, A, B, C ))

    const Int half = C.Height()/2;
    const auto indTL = IR(0,half);
    const auto indBR = IR(half,END);

    auto AL = A(ALL,indTL);
    auto AR = A(ALL,indBR);
    auto BT = B(indTL,ALL);
    auto BB = B(indBR,ALL);
    auto CTL = C(indTL,indTL);
    auto CBR = C(indBR,indBR);

    if( uplo == LOWER )
    {
        auto CBL = C(indBR,indTL);
        Gemm( orientationOfA, orientationOfB, alpha, AR, BT, T(1), CBL );
    }
    else
    {
        auto CTR = C(indTL,indBR);
        Gemm( orientationOfA, orientationOfB, alpha, AL, BB, T(1), CTR );
    }

    // TODO(poulson): Avoid the temporary copy
    Matrix<T> DTL;
    Gemm( orientationOfA, orientationOfB, alpha, AL, BT, DTL );
    AxpyTrapezoid( uplo, T(1), DTL, CTL );

    // TODO(poulson): Avoid the temporary copy
    Matrix<T> DBR;
    Gemm( orientationOfA, orientationOfB, alpha, AR, BB, DBR );
    AxpyTrapezoid( uplo, T(1), DBR, CBR );
}

// Distributed C := alpha A^{T/H} B^{T/H} + C
template<typename T>
void LocalTrrkKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,STAR,MC  >& A,
           const DistMatrix<T,MR,  STAR>& B,
                 DistMatrix<T>& C )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(CheckInput( A, B, C ))
    const Grid& g = C.Grid();

    const Int half = C.Height()/2;
    const auto indTL = IR(0,half);
    const auto indBR = IR(half,END);

    auto AL = A(ALL,indTL);
    auto AR = A(ALL,indBR);
    auto BT = B(indTL,ALL);
    auto BB = B(indBR,ALL);
    auto CTL = C(indTL,indTL);
    auto CBR = C(indBR,indBR);

    if( uplo == LOWER )
    {
        auto CBL = C(indBR,indTL);
        LocalGemm( orientationOfA, orientationOfB, alpha, AR, BT, T(1), CBL );
    }
    else
    {
        auto CTR = C(indTL,indBR);
        LocalGemm( orientationOfA, orientationOfB, alpha, AL, BB, T(1), CTR );
    }

    // TODO(poulson): Avoid the temporary copy
    DistMatrix<T> DTL(g);
    DTL.AlignWith( CTL );
    LocalGemm( orientationOfA, orientationOfB, alpha, AL, BT, DTL );
    LocalAxpyTrapezoid( uplo, T(1), DTL, CTL );

    // TODO(poulson): Avoid the temporary copy
    DistMatrix<T> DBR(g);
    DBR.AlignWith( CBR );
    LocalGemm( orientationOfA, orientationOfB, alpha, AR, BB, DBR );
    LocalAxpyTrapezoid( uplo, T(1), DBR, CBR );
}

// Local C := alpha A B + C
template<typename T>
void TrrkNN
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
                 Matrix<T>& C )
{
    EL_DEBUG_CSE
    using namespace trrk;
    EL_DEBUG_ONLY(CheckInputNN( A, B, C ))
    if( C.Height() < LocalTrrkBlocksize<T>() )
    {
        TrrkNNKernel( uplo, alpha, A, B, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        const Int half = C.Height()/2;
        const auto indTL = IR(0,half);
        const auto indBR = IR(half,END);

        auto AT = A(indTL,ALL);
        auto AB = A(indBR,ALL);
        auto BL = B(ALL,indTL);
        auto BR = B(ALL,indBR);
        auto CTL = C(indTL,indTL);
        auto CBR = C(indBR,indBR);

        if( uplo == LOWER )
        {
            auto CBL = C(indBR,indTL);
            Gemm( NORMAL, NORMAL, alpha, AB, BL, T(1), CBL );
        }
        else
        {
            auto CTR = C(indTL,indBR);
            Gemm( NORMAL, NORMAL, alpha, AT, BR, T(1), CTR );
        }

        // Recurse
        TrrkNN( uplo, alpha, AT, BL, CTL );
        TrrkNN( uplo, alpha, AB, BR, CBR );
    }
}

// Local C := alpha A B^{T/H} + C
template<typename T>
void TrrkNT
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
                 Matrix<T>& C )
{
    EL_DEBUG_CSE
    using namespace trrk;
    EL_DEBUG_ONLY(CheckInputNT( orientationOfB, A, B, C ))
    if( C.Height() < LocalTrrkBlocksize<T>() )
    {
        TrrkNTKernel( uplo, orientationOfB, alpha, A, B, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        const Int half = C.Height()/2;
        const auto indTL = IR(0,half);
        const auto indBR = IR(half,END);

        auto AT = A(indTL,ALL);
        auto AB = A(indBR,ALL);
        auto BT = B(indTL,ALL);
        auto BB = B(indBR,ALL);
        auto CTL = C(indTL,indTL);
        auto CBR = C(indBR,indBR);

        if( uplo == LOWER )
        {
            auto CBL = C(indBR,indTL);
            Gemm( NORMAL, orientationOfB, alpha, AB, BT, T(1), CBL );
        }
        else
        {
            auto CTR = C(indTL,indBR);
            Gemm( NORMAL, orientationOfB, alpha, AT, BB, T(1), CTR );
        }

        // Recurse
        TrrkNT( uplo, orientationOfB, alpha, AT, BT, CTL );
        TrrkNT( uplo, orientationOfB, alpha, AB, BB, CBR );
    }
}

// Local C := alpha A^{T/H} B + C
template<typename T>
void TrrkTN
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
                 Matrix<T>& C )
{
    EL_DEBUG_CSE
    using namespace trrk;
    EL_DEBUG_ONLY(CheckInputTN( orientationOfA, A, B, C ))
    if( C.Height() < LocalTrrkBlocksize<T>() )
    {
        TrrkTNKernel( uplo, orientationOfA, alpha, A, B, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        const Int half = C.Height()/2;
        const auto indTL = IR(0,half);
        const auto indBR = IR(half,END);

        auto AL = A(ALL,indTL);
        auto AR = A(ALL,indBR);
        auto BL = B(ALL,indTL);
        auto BR = B(ALL,indBR);
        auto CTL = C(indTL,indTL);
        auto CBR = C(indBR,indBR);

        if( uplo == LOWER )
        {
            auto CBL = C(indBR,indTL);
            Gemm( orientationOfA, NORMAL, alpha, AR, BL, T(1), CBL );
        }
        else
        {
            auto CTR = C(indTL,indBR);
            Gemm( orientationOfA, NORMAL, alpha, AL, BR, T(1), CTR );
        }

        // Recurse
        TrrkTN( uplo, orientationOfA, alpha, AL, BL, CTL );
        TrrkTN( uplo, orientationOfA, alpha, AR, BR, CBR );
    }
}

// Local C := alpha A^{T/H} B^{T/H} + C
template<typename T>
void TrrkTT
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
                 Matrix<T>& C )
{
    EL_DEBUG_CSE
    using namespace trrk;
    EL_DEBUG_ONLY(CheckInputTT( orientationOfA, orientationOfB, A, B, C ))
    if( C.Height() < LocalTrrkBlocksize<T>() )
    {
        TrrkTTKernel( uplo, orientationOfA, orientationOfB, alpha, A, B, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        const Int half = C.Height()/2;
        const auto indTL = IR(0,half);
        const auto indBR = IR(half,END);

        auto AL = A(ALL,indTL);
        auto AR = A(ALL,indBR);
        auto BT = B(indTL,ALL);
        auto BB = B(indBR,ALL);
        auto CTL = C(indTL,indTL);
        auto CBR = C(indBR,indBR);

        if( uplo == LOWER )
        {
            auto CBL = C(indBR,indTL);
            Gemm( orientationOfA, orientationOfB, alpha, AR, BT, T(1), CBL );
        }
        else
        {
            auto CTR = C(indTL,indBR);
            Gemm( orientationOfA, orientationOfB, alpha, AL, BB, T(1), CTR );
        }

        // Recurse
        TrrkTT( uplo, orientationOfA, orientationOfB, alpha, AL, BT, CTL );
        TrrkTT( uplo, orientationOfA, orientationOfB, alpha, AR, BB, CBR );
    }
}

} // namespace trrk

// Distributed C := alpha A B + beta C
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,  STAR>& A,
           const DistMatrix<T,STAR,MR  >& B,
  T beta,        DistMatrix<T>& C )
{
    EL_DEBUG_CSE
    using namespace trrk;
    EL_DEBUG_ONLY(CheckInput( A, B, C ))
    const Grid& g = C.Grid();
    ScaleTrapezoid( beta, uplo, C );

    if( C.Height() < g.Width()*LocalTrrkBlocksize<T>() )
    {
        LocalTrrkKernel( uplo, alpha, A, B, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        const Int half = C.Height()/2;
        const auto indTL = IR(0,half);
        const auto indBR = IR(half,END);

        auto AT = A(indTL,ALL);
        auto AB = A(indBR,ALL);
        auto BL = B(ALL,indTL);
        auto BR = B(ALL,indBR);
        auto CTL = C(indTL,indTL);
        auto CBR = C(indBR,indBR);

        if( uplo == LOWER )
        {
            auto CBL = C(indBR,indTL);
            LocalGemm( NORMAL, NORMAL, alpha, AB, BL, T(1), CBL );
        }
        else
        {
            auto CTR = C(indTL,indBR);
            LocalGemm( NORMAL, NORMAL, alpha, AT, BR, T(1), CTR );
        }

        // Recurse
        LocalTrrk( uplo, alpha, AT, BL, T(1), CTL );
        LocalTrrk( uplo, alpha, AB, BR, T(1), CBR );
    }
}

// Distributed C := alpha A B^{T/H} + beta C
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,STAR>& A,
           const DistMatrix<T,MR,STAR>& B,
  T beta,        DistMatrix<T>& C )
{
    EL_DEBUG_CSE
    using namespace trrk;
    EL_DEBUG_ONLY(CheckInput( A, B, C ))
    const Grid& g = C.Grid();
    ScaleTrapezoid( beta, uplo, C );

    if( C.Height() < g.Width()*LocalTrrkBlocksize<T>() )
    {
        LocalTrrkKernel( uplo, orientationOfB, alpha, A, B, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        const Int half = C.Height()/2;
        const auto indTL = IR(0,half);
        const auto indBR = IR(half,END);

        auto AT = A(indTL,ALL);
        auto AB = A(indBR,ALL);
        auto BT = B(indTL,ALL);
        auto BB = B(indBR,ALL);
        auto CTL = C(indTL,indTL);
        auto CBR = C(indBR,indBR);

        if( uplo == LOWER )
        {
            auto CBL = C(indBR,indTL);
            LocalGemm( NORMAL, orientationOfB, alpha, AB, BT, T(1), CBL );
        }
        else
        {
            auto CTR = C(indTL,indBR);
            LocalGemm( NORMAL, orientationOfB, alpha, AT, BB, T(1), CTR );
        }

        // Recurse
        LocalTrrk( uplo, orientationOfB, alpha, AT, BT, T(1), CTL );
        LocalTrrk( uplo, orientationOfB, alpha, AB, BB, T(1), CBR );
    }
}

// Distributed C := alpha A^{T/H} B + beta C
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T,STAR,MC>& A,
           const DistMatrix<T,STAR,MR>& B,
  T beta,        DistMatrix<T>& C )
{
    EL_DEBUG_CSE
    using namespace trrk;
    EL_DEBUG_ONLY(CheckInput( A, B, C ))
    const Grid& g = C.Grid();
    ScaleTrapezoid( beta, uplo, C );

    if( C.Height() < g.Width()*LocalTrrkBlocksize<T>() )
    {
        LocalTrrkKernel( uplo, orientationOfA, alpha, A, B, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        const Int half = C.Height()/2;
        const auto indTL = IR(0,half);
        const auto indBR = IR(half,END);

        auto AL = A(ALL,indTL);
        auto AR = A(ALL,indBR);
        auto BL = B(ALL,indTL);
        auto BR = B(ALL,indBR);
        auto CTL = C(indTL,indTL);
        auto CBR = C(indBR,indBR);

        if( uplo == LOWER )
        {
            auto CBL = C(indBR,indTL);
            LocalGemm( orientationOfA, NORMAL, alpha, AR, BL, T(1), CBL );
        }
        else
        {
            auto CTR = C(indTL,indBR);
            LocalGemm( orientationOfA, NORMAL, alpha, AL, BR, T(1), CTR );
        }

        // Recurse
        LocalTrrk( uplo, orientationOfA, alpha, AL, BL, T(1), CTL );
        LocalTrrk( uplo, orientationOfA, alpha, AR, BR, T(1), CBR );
    }
}

// Distributed C := alpha A^{T/H} B^{T/H} + beta C
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,STAR,MC  >& A,
           const DistMatrix<T,MR,  STAR>& B,
  T beta,        DistMatrix<T>& C )
{
    EL_DEBUG_CSE
    using namespace trrk;
    EL_DEBUG_ONLY(CheckInput( A, B, C ))
    const Grid& g = C.Grid();
    ScaleTrapezoid( beta, uplo, C );

    if( C.Height() < g.Width()*LocalTrrkBlocksize<T>() )
    {
        LocalTrrkKernel( uplo, orientationOfA, orientationOfB, alpha, A, B, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        const Int half = C.Height()/2;
        const auto indTL = IR(0,half);
        const auto indBR = IR(half,END);

        auto AL = A(ALL,indTL);
        auto AR = A(ALL,indBR);
        auto BT = B(indTL,ALL);
        auto BB = B(indBR,ALL);
        auto CTL = C(indTL,indTL);
        auto CBR = C(indBR,indBR);

        if( uplo == LOWER )
        {
            auto CBL = C(indBR,indTL);
            LocalGemm
            ( orientationOfA, orientationOfB, alpha, AR, BT, T(1), CBL );
        }
        else
        {
            auto CTR = C(indTL,indBR);
            LocalGemm
            ( orientationOfA, orientationOfB, alpha, AL, BB, T(1), CTR );
        }

        // Recurse
        LocalTrrk
        ( uplo, orientationOfA, orientationOfB, alpha, AL, BT, T(1), CTL );
        LocalTrrk
        ( uplo, orientationOfA, orientationOfB, alpha, AR, BB, T(1), CBR );
    }
}

} // namespace El

#endif // ifndef EL_TRRK_LOCAL_HPP
