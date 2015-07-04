/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRRK_LOCAL_HPP
#define EL_TRRK_LOCAL_HPP

namespace El {

namespace trrk {

#ifndef EL_RELEASE

template<typename Ring>
void EnsureConformal
( const DistMatrix<Ring,MC,STAR>& A, const DistMatrix<Ring>& C, string name )
{
    if( A.Height() != C.Height() || A.ColAlign() != C.ColAlign() )
        LogicError(name," not conformal with C");
}

template<typename Ring>
void EnsureConformal
( const DistMatrix<Ring,STAR,MC>& A, const DistMatrix<Ring>& C, string name )
{
    if( A.Width() != C.Height() || A.RowAlign() != C.ColAlign() )
        LogicError(name," not conformal with C");
}

template<typename Ring>
void EnsureConformal
( const DistMatrix<Ring,MR,STAR>& A, const DistMatrix<Ring>& C, string name )
{
    if( A.Height() != C.Width() || A.ColAlign() != C.RowAlign() )
        LogicError(name," not conformal with C");
}

template<typename Ring>
void EnsureConformal
( const DistMatrix<Ring,STAR,MR>& A, const DistMatrix<Ring>& C, string name )
{
    if( A.Width() != C.Width() || A.RowAlign() != C.RowAlign() )
        LogicError(name," not conformal with C");
}

template<typename Ring,Dist UA,Dist VA,Dist UB,Dist VB>
void CheckInput
( const DistMatrix<Ring,UA,VA>& A, const DistMatrix<Ring,UB,VB>& B,
  const DistMatrix<Ring>& C )
{
    AssertSameGrids( A, B, C );
    EnsureConformal( A, C, "A" );
    EnsureConformal( B, C, "B" );
}

// Local C := A B + C
template<typename Ring>
void CheckInputNN
( const Matrix<Ring>& A, const Matrix<Ring>& B, const Matrix<Ring>& C )
{
    if( A.Height() != C.Height() || B.Width()  != C.Width() ||
        A.Width()  != B.Height() || A.Height() != B.Width() )
        LogicError
        ("Nonconformal LocalTrrk:\n",
         DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
}

// Local C := A B^{T/H} + C
template<typename Ring>
void CheckInputNT
( Orientation orientB,
  const Matrix<Ring>& A, const Matrix<Ring>& B, const Matrix<Ring>& C )
{
    if( orientB == NORMAL )
        LogicError("B must be (Conjugate)Transpose'd");
    if( A.Height() != C.Height() || B.Height() != C.Width() ||
        A.Width()  != B.Width()  || A.Height() != B.Height() )
        LogicError
        ("Nonconformal LocalTrrk:\n",
         DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
}

// Local C := A^{T/H} B + C
template<typename Ring>
void CheckInputTN
( Orientation orientA,
  const Matrix<Ring>& A, const Matrix<Ring>& B, const Matrix<Ring>& C )
{
    if( orientA == NORMAL )
        LogicError("A must be (Conjugate)Transpose'd");
    if( A.Width() != C.Height() || B.Width() != C.Width() ||
        A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError
        ("Nonconformal LocalTrrk:\n",
         DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
}

// Local C := A^{T/H} B^{T/H} + C
template<typename Ring>
void CheckInputTT
( Orientation orientA,
  Orientation orientB,
  const Matrix<Ring>& A, const Matrix<Ring>& B, const Matrix<Ring>& C )
{
    if( orientA == NORMAL )
        LogicError("A must be (Conjugate)Transpose'd");
    if( orientB == NORMAL )
        LogicError("B must be (Conjugate)Transpose'd");
    if( A.Width() != C.Height() || B.Height() != C.Width() ||
        A.Height() != B.Width() || A.Width() != B.Height() )
        LogicError
        ("Nonconformal LocalTrrk:\n",
         DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
}

#endif // ifndef EL_RELEASE

// Local C := alpha A B + C
template<typename Ring>
inline void
TrrkNNKernel
( UpperOrLower uplo, 
  Ring alpha, const Matrix<Ring>& A, 
              const Matrix<Ring>& B,
                    Matrix<Ring>& C )
{
    DEBUG_ONLY(
      CSE cse("TrrkNNKernel");
      CheckInputNN( A, B, C );
    )
    Matrix<Ring> AT, AB;
    Matrix<Ring> BL, BR;
    Matrix<Ring> CTL, CTR,
                 CBL, CBR;
    Matrix<Ring> DTL, DBR;

    const Int half = C.Height()/2;
    LockedPartitionDown( A, AT, AB, half );
    LockedPartitionRight( B, BL, BR, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    if( uplo == LOWER )
        Gemm( alpha, AB.N(), BL.N(), Ring(1), CBL );
    else
        Gemm( alpha, AT.N(), BR.N(), Ring(1), CTR );

    Gemm( alpha, AT.N(), BL.N(), DTL );
    AxpyTrapezoid( uplo, Ring(1), DTL, CTL );

    Gemm( alpha, AB.N(), BR.N(), DBR );
    AxpyTrapezoid( uplo, Ring(1), DBR, CBR );
}

// Distributed C := alpha A B + C
template<typename Ring>
inline void
LocalTrrkKernel
( UpperOrLower uplo, 
  Ring alpha, const DistMatrix<Ring,MC,  STAR>& A,
              const DistMatrix<Ring,STAR,MR  >& B,
                    DistMatrix<Ring>& C )
{
    DEBUG_ONLY(
      CSE cse("LocalTrrkKernel");
      CheckInput( A, B, C );
    )
    const Grid& g = C.Grid();

    DistMatrix<Ring,MC,STAR> AT(g), AB(g);
    DistMatrix<Ring,STAR,MR> BL(g), BR(g);
    DistMatrix<Ring> CTL(g), CTR(g),
                     CBL(g), CBR(g);
    DistMatrix<Ring> DTL(g), DBR(g);

    const Int half = C.Height()/2;
    LockedPartitionDown( A, AT, AB, half );
    LockedPartitionRight( B, BL, BR, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    if( uplo == LOWER )
        LocalGemm( alpha, AB.N(), BL.N(), Ring(1), CBL );
    else
        LocalGemm( alpha, AT.N(), BR.N(), Ring(1), CTR );

    DTL.AlignWith( CTL );
    LocalGemm( alpha, AT.N(), BL.N(), DTL );
    AxpyTrapezoid( uplo, Ring(1), DTL, CTL );

    DBR.AlignWith( CBR );
    LocalGemm( alpha, AB.N(), BR.N(), DBR );
    AxpyTrapezoid( uplo, Ring(1), DBR, CBR );
}

// Local C := alpha A B^{T/H} + C
template<typename Ring>
inline void
TrrkNTKernel
( UpperOrLower uplo,
  Orientation orientB,
  Ring alpha, const Matrix<Ring>& A, 
              const Matrix<Ring>& B,
                    Matrix<Ring>& C )
{
    DEBUG_ONLY(
      CSE cse("TrrkNTKernel");
      CheckInputNT( orientB, A, B, C );
    )
    Matrix<Ring> AT, AB;
    Matrix<Ring> BT, BB;
    Matrix<Ring> CTL, CTR,
                 CBL, CBR;
    Matrix<Ring> DTL, DBR;

    const Int half = C.Height()/2;
    LockedPartitionDown( A, AT, AB, half );
    LockedPartitionDown( B, BT, BB, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    if( uplo == LOWER )
        Gemm( alpha, AB.N(), BT.Orient(orientB), Ring(1), CBL );
    else
        Gemm( alpha, AT.N(), BB.Orient(orientB), Ring(1), CTR );

    Gemm( alpha, AT.N(), BT.Orient(orientB), DTL );
    AxpyTrapezoid( uplo, Ring(1), DTL, CTL );

    Gemm( alpha, AB.N(), BB.Orient(orientB), DBR );
    AxpyTrapezoid( uplo, Ring(1), DBR, CBR );
}

// Distributed C := alpha A B^{T/H} + C
template<typename Ring>
inline void
LocalTrrkKernel
( UpperOrLower uplo,
  Orientation orientB,
  Ring alpha, const DistMatrix<Ring,MC,STAR>& A,
              const DistMatrix<Ring,MR,STAR>& B,
                    DistMatrix<Ring>& C )
{
    DEBUG_ONLY(
      CSE cse("LocalTrrkKernel");
      CheckInput( A, B, C );
    )
    const Grid& g = C.Grid();

    DistMatrix<Ring,MC,STAR> AT(g), AB(g);
    DistMatrix<Ring,MR,STAR> BT(g), BB(g);
    DistMatrix<Ring> CTL(g), CTR(g),
                     CBL(g), CBR(g);
    DistMatrix<Ring> DTL(g), DBR(g);

    const Int half = C.Height()/2;
    LockedPartitionDown( A, AT, AB, half );
    LockedPartitionDown( B, BT, BB, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    if( uplo == LOWER )
        LocalGemm( alpha, AB.N(), BT.Orient(orientB), Ring(1), CBL );
    else
        LocalGemm( alpha, AT.N(), BB.Orient(orientB), Ring(1), CTR );

    DTL.AlignWith( CTL );
    LocalGemm( alpha, AT.N(), BT.Orient(orientB), DTL );
    AxpyTrapezoid( uplo, Ring(1), DTL, CTL );

    DBR.AlignWith( CBR );
    LocalGemm( alpha, AB.N(), BB.Orient(orientB), DBR );
    AxpyTrapezoid( uplo, Ring(1), DBR, CBR );
}

// Local C := alpha A^{T/H} B + C
template<typename Ring>
inline void
TrrkTNKernel
( UpperOrLower uplo,
  Orientation orientA,
  Ring alpha, const Matrix<Ring>& A, 
              const Matrix<Ring>& B,
                    Matrix<Ring>& C )
{
    DEBUG_ONLY(
      CSE cse("TrrkTNKernel");
      CheckInputTN( orientA, A, B, C );
    )
    Matrix<Ring> AL, AR;
    Matrix<Ring> BL, BR;
    Matrix<Ring> CTL, CTR,
                 CBL, CBR;
    Matrix<Ring> DTL, DBR;

    const Int half = C.Height()/2;
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionRight( B, BL, BR, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    if( uplo == LOWER )
        Gemm( alpha, AR.Orient(orientA), BL.N(), Ring(1), CBL );
    else
        Gemm( alpha, AL.Orient(orientA), BR.N(), Ring(1), CTR );

    Gemm( alpha, AL.Orient(orientA), BL.N(), DTL );
    AxpyTrapezoid( uplo, Ring(1), DTL, CTL );

    Gemm( alpha, AR.Orient(orientA), BR.N(), DBR );
    AxpyTrapezoid( uplo, Ring(1), DBR, CBR );
}

// Distributed C := alpha A^{T/H} B + C
template<typename Ring>
inline void
LocalTrrkKernel
( UpperOrLower uplo,
  Orientation orientA,
  Ring alpha, const DistMatrix<Ring,STAR,MC>& A,
              const DistMatrix<Ring,STAR,MR>& B,
                    DistMatrix<Ring>& C )
{
    DEBUG_ONLY(
      CSE cse("LocalTrrkKernel");
      CheckInput( A, B, C );
    )
    const Grid& g = C.Grid();

    DistMatrix<Ring,STAR,MC> AL(g), AR(g);
    DistMatrix<Ring,STAR,MR> BL(g), BR(g);
    DistMatrix<Ring> CTL(g), CTR(g),
                     CBL(g), CBR(g);
    DistMatrix<Ring> DTL(g), DBR(g);

    const Int half = C.Height()/2;
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionRight( B, BL, BR, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    if( uplo == LOWER )
        LocalGemm( alpha, AR.Orient(orientA), BL.N(), Ring(1), CBL );
    else
        LocalGemm( alpha, AL.Orient(orientA), BR.N(), Ring(1), CTR );

    DTL.AlignWith( CTL );
    LocalGemm( alpha, AL.Orient(orientA), BL.N(), DTL );
    AxpyTrapezoid( uplo, Ring(1), DTL, CTL );

    DBR.AlignWith( CBR );
    LocalGemm( alpha, AR.Orient(orientA), BR.N(), DBR );
    AxpyTrapezoid( uplo, Ring(1), DBR, CBR );
}

// Local C := alpha A^{T/H} B^{T/H} + C
template<typename Ring>
inline void
TrrkTTKernel
( UpperOrLower uplo,
  Orientation orientA,
  Orientation orientB,
  Ring alpha, const Matrix<Ring>& A, 
              const Matrix<Ring>& B,
                    Matrix<Ring>& C )
{
    DEBUG_ONLY(
      CSE cse("TrrkTTKernel");
      CheckInputTT( orientA, orientB, A, B, C );
    )
    Matrix<Ring> AL, AR;
    Matrix<Ring> BT, BB;
    Matrix<Ring> CTL, CTR,
                 CBL, CBR;
    Matrix<Ring> DTL, DBR;

    const Int half = C.Height()/2;
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionDown( B, BT, BB, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    if( uplo == LOWER )
        Gemm( alpha, AR.Orient(orientA), BT.Orient(orientB), Ring(1), CBL );
    else
        Gemm( alpha, AL.Orient(orientA), BB.Orient(orientB), Ring(1), CTR );

    Gemm( alpha, AL.Orient(orientA), BT.Orient(orientB), DTL );
    AxpyTrapezoid( uplo, Ring(1), DTL, CTL );

    Gemm( alpha, AR.Orient(orientA), BB.Orient(orientB), DBR );
    AxpyTrapezoid( uplo, Ring(1), DBR, CBR );
}

// Distributed C := alpha A^{T/H} B^{T/H} + C
template<typename Ring>
inline void
LocalTrrkKernel
( UpperOrLower uplo,
  Orientation orientA,
  Orientation orientB,
  Ring alpha, const DistMatrix<Ring,STAR,MC  >& A,
              const DistMatrix<Ring,MR,  STAR>& B,
                    DistMatrix<Ring>& C )
{
    DEBUG_ONLY(
      CSE cse("LocalTrrkKernel");
      CheckInput( A, B, C );
    )
    const Grid& g = C.Grid();

    DistMatrix<Ring,STAR,MC> AL(g), AR(g);
    DistMatrix<Ring,MR,STAR> BT(g), BB(g);
    DistMatrix<Ring> CTL(g), CTR(g),
                     CBL(g), CBR(g);
    DistMatrix<Ring> DTL(g), DBR(g);

    const Int half = C.Height()/2;
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionDown( B, BT, BB, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    if( uplo == LOWER )
        LocalGemm
        ( alpha, AR.Orient(orientA), BT.Orient(orientB), Ring(1), CBL );
    else
        LocalGemm
        ( alpha, AL.Orient(orientA), BB.Orient(orientB), Ring(1), CTR );

    DTL.AlignWith( CTL );
    LocalGemm( alpha, AL.Orient(orientA), BT.Orient(orientB), DTL );
    AxpyTrapezoid( uplo, Ring(1), DTL, CTL );

    DBR.AlignWith( CBR );
    LocalGemm( alpha, AR.Orient(orientA), BB.Orient(orientB), DBR );
    AxpyTrapezoid( uplo, Ring(1), DBR, CBR );
}

// Local C := alpha A B + C
template<typename Ring>
void TrrkNN
( UpperOrLower uplo,
  Ring alpha, const Matrix<Ring>& A, 
              const Matrix<Ring>& B,
                    Matrix<Ring>& C )
{
    using namespace trrk;
    DEBUG_ONLY(
      CSE cse("trrk::TrrkNN");
      CheckInputNN( A, B, C );
    )
    if( C.Height() < LocalTrrkBlocksize<Ring>() )
    {
        TrrkNNKernel( uplo, alpha, A, B, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        Matrix<Ring> AT, AB;
        Matrix<Ring> BL, BR;
        Matrix<Ring> CTL, CTR,
                     CBL, CBR;

        const Int half = C.Height() / 2;
        LockedPartitionDown( A, AT, AB, half );
        LockedPartitionRight( B, BL, BR, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            Gemm( alpha, AB.N(), BL.N(), Ring(1), CBL );
        else
            Gemm( alpha, AT.N(), BR.N(), Ring(1), CTR );

        // Recurse
        TrrkNN( uplo, alpha, AT, BL, CTL );
        TrrkNN( uplo, alpha, AB, BR, CBR );
    }
}

// Local C := alpha A B^{T/H} + C
template<typename Ring>
void TrrkNT
( UpperOrLower uplo,
  Orientation orientB,
  Ring alpha, const Matrix<Ring>& A, 
              const Matrix<Ring>& B,
                    Matrix<Ring>& C )
{
    using namespace trrk;
    DEBUG_ONLY(
      CSE cse("trrk::TrrkNT");
      CheckInputNT( orientB, A, B, C );
    )
    if( C.Height() < LocalTrrkBlocksize<Ring>() )
    {
        TrrkNTKernel( uplo, orientB, alpha, A, B, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        Matrix<Ring> AT, AB;
        Matrix<Ring> BT, BB;
        Matrix<Ring> CTL, CTR,
                     CBL, CBR;

        const Int half = C.Height() / 2;
        LockedPartitionDown( A, AT, AB, half );
        LockedPartitionDown( B, BT, BB, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            Gemm( alpha, AB.N(), BT.Orient(orientB), Ring(1), CBL );
        else
            Gemm( alpha, AT.N(), BB.Orient(orientB), Ring(1), CTR );

        // Recurse
        TrrkNT( uplo, orientB, alpha, AT, BT, CTL );
        TrrkNT( uplo, orientB, alpha, AB, BB, CBR );
    }
}

// Local C := alpha A^{T/H} B + C
template<typename Ring>
void TrrkTN
( UpperOrLower uplo,
  Orientation orientA,
  Ring alpha, const Matrix<Ring>& A, 
              const Matrix<Ring>& B,
                    Matrix<Ring>& C )
{
    using namespace trrk;
    DEBUG_ONLY(
      CSE cse("trrk::TrrkTN");
      CheckInputTN( orientA, A, B, C );
    )
    if( C.Height() < LocalTrrkBlocksize<Ring>() )
    {
        TrrkTNKernel( uplo, orientA, alpha, A, B, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        Matrix<Ring> AL, AR;
        Matrix<Ring> BL, BR;
        Matrix<Ring> CTL, CTR,
                     CBL, CBR;

        const Int half = C.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionRight( B, BL, BR, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            Gemm( alpha, AR.Orient(orientA), BL.N(), Ring(1), CBL );
        else
            Gemm( alpha, AL.Orient(orientA), BR.N(), Ring(1), CTR );

        // Recurse
        TrrkTN( uplo, orientA, alpha, AL, BL, CTL );
        TrrkTN( uplo, orientA, alpha, AR, BR, CBR );
    }
}

// Local C := alpha A^{T/H} B^{T/H} + C
template<typename Ring>
void TrrkTT
( UpperOrLower uplo,
  Orientation orientA, Orientation orientB,
  Ring alpha, const Matrix<Ring>& A, 
              const Matrix<Ring>& B,
                    Matrix<Ring>& C )
{
    using namespace trrk;
    DEBUG_ONLY(
      CSE cse("trrk::TrrkTT");
      CheckInputTT( orientA, orientB, A, B, C );
    )
    if( C.Height() < LocalTrrkBlocksize<Ring>() )
    {
        TrrkTTKernel( uplo, orientA, orientB, alpha, A, B, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        Matrix<Ring> AL, AR;
        Matrix<Ring> BT, BB;
        Matrix<Ring> CTL, CTR,
                     CBL, CBR;

        const Int half = C.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionDown( B, BT, BB, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            Gemm( alpha, AR.Orient(orientA), BT.Orient(orientB), Ring(1), CBL );
        else
            Gemm( alpha, AL.Orient(orientA), BB.Orient(orientB), Ring(1), CTR );

        // Recurse
        TrrkTT( uplo, orientA, orientB, alpha, AL, BT, CTL );
        TrrkTT( uplo, orientA, orientB, alpha, AR, BB, CBR );
    }
}

} // namespace trrk

// Distributed C := alpha A B + beta C
template<typename Ring>
void LocalTrrk
( UpperOrLower uplo,
  Ring alpha, const DistMatrix<Ring,MC,  STAR>& A,
              const DistMatrix<Ring,STAR,MR  >& B,
  Ring beta,        DistMatrix<Ring>& C )
{
    using namespace trrk;
    DEBUG_ONLY(
      CSE cse("LocalTrrk");
      CheckInput( A, B, C );
    )
    const Grid& g = C.Grid();
    ScaleTrapezoid( beta, uplo, C );

    if( C.Height() < g.Width()*LocalTrrkBlocksize<Ring>() )
    {
        LocalTrrkKernel( uplo, alpha, A, B, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        DistMatrix<Ring,MC,STAR> AT(g), AB(g);
        DistMatrix<Ring,STAR,MR> BL(g), BR(g);
        DistMatrix<Ring> CTL(g), CTR(g),
                         CBL(g), CBR(g);

        const Int half = C.Height() / 2;
        LockedPartitionDown( A, AT, AB, half );
        LockedPartitionRight( B, BL, BR, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            LocalGemm( alpha, AB.N(), BL.N(), Ring(1), CBL );
        else
            LocalGemm( alpha, AT.N(), BR.N(), Ring(1), CTR );

        // Recurse
        LocalTrrk( uplo, alpha, AT, BL, Ring(1), CTL );
        LocalTrrk( uplo, alpha, AB, BR, Ring(1), CBR );
    }
}

// Distributed C := alpha A B^{T/H} + beta C
template<typename Ring>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientB,
  Ring alpha, const DistMatrix<Ring,MC,STAR>& A,
              const DistMatrix<Ring,MR,STAR>& B,
  Ring beta,        DistMatrix<Ring>& C )
{
    using namespace trrk;
    DEBUG_ONLY(
      CSE cse("LocalTrrk");
      CheckInput( A, B, C );
    )
    const Grid& g = C.Grid();
    ScaleTrapezoid( beta, uplo, C );

    if( C.Height() < g.Width()*LocalTrrkBlocksize<Ring>() )
    {
        LocalTrrkKernel( uplo, orientB, alpha, A, B, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        DistMatrix<Ring,MC,STAR> AT(g), AB(g);
        DistMatrix<Ring,MR,STAR> BT(g), BB(g);
        DistMatrix<Ring> CTL(g), CTR(g),
                         CBL(g), CBR(g);

        const Int half = C.Height() / 2;
        LockedPartitionDown( A, AT, AB, half );
        LockedPartitionDown( B, BT, BB, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            LocalGemm( alpha, AB.N(), BT.Orient(orientB), Ring(1), CBL );
        else
            LocalGemm( alpha, AT.N(), BB.Orient(orientB), Ring(1), CTR );

        // Recurse
        LocalTrrk( uplo, orientB, alpha, AT, BT, Ring(1), CTL );
        LocalTrrk( uplo, orientB, alpha, AB, BB, Ring(1), CBR );
    }
}

// Distributed C := alpha A^{T/H} B + beta C
template<typename Ring>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientA,
  Ring alpha, const DistMatrix<Ring,STAR,MC>& A,
              const DistMatrix<Ring,STAR,MR>& B,
  Ring beta,        DistMatrix<Ring>& C )
{
    using namespace trrk;
    DEBUG_ONLY(
      CSE cse("LocalTrrk");
      CheckInput( A, B, C );
    )
    const Grid& g = C.Grid();
    ScaleTrapezoid( beta, uplo, C );

    if( C.Height() < g.Width()*LocalTrrkBlocksize<Ring>() )
    {
        LocalTrrkKernel( uplo, orientA, alpha, A, B, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        DistMatrix<Ring,STAR,MC> AL(g), AR(g);
        DistMatrix<Ring,STAR,MR> BL(g), BR(g);
        DistMatrix<Ring> CTL(g), CTR(g),
                         CBL(g), CBR(g);

        const Int half = C.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionRight( B, BL, BR, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            LocalGemm( alpha, AR.Orient(orientA), BL.N(), Ring(1), CBL );
        else
            LocalGemm( alpha, AL.Orient(orientA), BR.N(), Ring(1), CTR );

        // Recurse
        LocalTrrk( uplo, orientA, alpha, AL, BL, Ring(1), CTL );
        LocalTrrk( uplo, orientA, alpha, AR, BR, Ring(1), CBR );
    }
}

// Distributed C := alpha A^{T/H} B^{T/H} + beta C
template<typename Ring>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientA, Orientation orientB,
  Ring alpha, const DistMatrix<Ring,STAR,MC  >& A,
              const DistMatrix<Ring,MR,  STAR>& B,
  Ring beta,        DistMatrix<Ring>& C )
{
    using namespace trrk;
    DEBUG_ONLY(
      CSE cse("LocalTrrk");
      CheckInput( A, B, C );
    )
    const Grid& g = C.Grid();
    ScaleTrapezoid( beta, uplo, C );

    if( C.Height() < g.Width()*LocalTrrkBlocksize<Ring>() )
    {
        LocalTrrkKernel( uplo, orientA, orientB, alpha, A, B, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        DistMatrix<Ring,STAR,MC> AL(g), AR(g);
        DistMatrix<Ring,MR,STAR> BT(g), BB(g);
        DistMatrix<Ring> CTL(g), CTR(g),
                         CBL(g), CBR(g);

        const Int half = C.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionDown( B, BT, BB, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            LocalGemm
            ( alpha, AR.Orient(orientA), BT.Orient(orientB), Ring(1), CBL );
        else
            LocalGemm
            ( alpha, AL.Orient(orientA), BB.Orient(orientB), Ring(1), CTR );

        // Recurse
        LocalTrrk
        ( uplo, orientA, orientB, alpha, AL, BT, Ring(1), CTL );
        LocalTrrk
        ( uplo, orientA, orientB, alpha, AR, BB, Ring(1), CBR );
    }
}

} // namespace El

#endif // ifndef EL_TRRK_LOCAL_HPP
