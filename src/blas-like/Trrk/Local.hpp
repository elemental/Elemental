/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef BLAS_TRRK_LOCAL_HPP
#define BLAS_TRRK_LOCAL_HPP

#include "elemental/blas-like/level1/AxpyTriangle.hpp"
#include "elemental/blas-like/level1/ScaleTrapezoid.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {

namespace trrk {

#ifndef RELEASE
// Local C := alpha A B + beta C
template<typename T>
void CheckInputNN( const Matrix<T>& A, const Matrix<T>& B, const Matrix<T>& C )
{
    if( A.Height() != C.Height() || B.Width()  != C.Width() ||
        A.Width()  != B.Height() || A.Height() != B.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrrk: \n"
            << "  A ~ " << A.Height() << " x "
                        << A.Width()  << "\n"
            << "  B ~ " << B.Height() << " x "
                        << B.Width()  << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// Distributed C := alpha A B + beta C
template<typename T>
void CheckInput
( const DistMatrix<T,MC,  STAR>& A, 
  const DistMatrix<T,STAR,MR  >& B,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A.Height() != C.Height() || B.Width()  != C.Width() ||
        A.Width()  != B.Height() || A.Height() != B.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrrk: \n"
            << "  A[MC,* ] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[* ,MR] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.ColAlignment() != C.ColAlignment() ||
        B.RowAlignment() != C.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrrk: \n"
            << "  A[MC,* ] ~ " << A.ColAlignment() << "\n"
            << "  B[* ,MR] ~ " << B.RowAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// Local C := alpha A B^{T/H} + beta C
template<typename T>
void CheckInputNT
( Orientation orientationOfB,
  const Matrix<T>& A, const Matrix<T>& B, const Matrix<T>& C )
{
    if( orientationOfB == NORMAL )
        throw std::logic_error("B must be (Conjugate)Transpose'd");
    if( A.Height() != C.Height() || B.Height() != C.Width() ||
        A.Width()  != B.Width()  || A.Height() != B.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrrk: \n"
            << "  A ~ " << A.Height() << " x "
                        << A.Width()  << "\n"
            << "  B ~ " << B.Height() << " x "
                        << B.Width()  << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// Distributed C := alpha A B^{T/H} + beta C
template<typename T>
void CheckInput
( Orientation orientationOfB,
  const DistMatrix<T,MC,STAR>& A,
  const DistMatrix<T,MR,STAR>& B,
  const DistMatrix<T>& C )
{
    if( orientationOfB == NORMAL )
        throw std::logic_error("B[MR,* ] must be (Conjugate)Transpose'd");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A.Height() != C.Height() || B.Height() != C.Width() ||
        A.Width()  != B.Width()  || A.Height() != B.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrrk: \n"
            << "  A[MC,* ] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[MR,* ] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.ColAlignment() != C.ColAlignment() ||
        B.ColAlignment() != C.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrrk: \n"
            << "  A[MC,* ] ~ " << A.ColAlignment() << "\n"
            << "  B[MR,* ] ~ " << B.ColAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// Local C := alpha A^{T/H} B + beta C
template<typename T>
void CheckInputTN
( Orientation orientationOfA,
  const Matrix<T>& A, const Matrix<T>& B, const Matrix<T>& C )
{
    if( orientationOfA == NORMAL )
        throw std::logic_error("A must be (Conjugate)Transpose'd");
    if( A.Width() != C.Height() || B.Width() != C.Width() ||
        A.Height() != B.Height() || A.Width() != B.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrrk: \n"
            << "  A ~ " << A.Height() << " x "
                        << A.Width()  << "\n"
            << "  B ~ " << B.Height() << " x "
                        << B.Width()  << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// Distributed C := alpha A^{T/H} B + beta C
template<typename T>
void CheckInput
( Orientation orientationOfA,
  const DistMatrix<T,STAR,MC>& A,
  const DistMatrix<T,STAR,MR>& B,
  const DistMatrix<T,MC,  MR>& C )
{
    if( orientationOfA == NORMAL )
        throw std::logic_error("A[* ,MC] must be (Conjugate)Transpose'd");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A.Width() != C.Height() || B.Width() != C.Width() ||
        A.Height() != B.Height() || A.Width() != B.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrrk: \n"
            << "  A[* ,MC] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[* ,MR] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.RowAlignment() != C.ColAlignment() ||
        B.RowAlignment() != C.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrrk: \n"
            << "  A[* ,MC] ~ " << A.RowAlignment() << "\n"
            << "  B[* ,MR] ~ " << B.RowAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// Local C := alpha A^{T/H} B^{T/H} + beta C
template<typename T>
void CheckInputTT
( Orientation orientationOfA,
  Orientation orientationOfB,
  const Matrix<T>& A, const Matrix<T>& B, const Matrix<T>& C )
{
    if( orientationOfA == NORMAL )
        throw std::logic_error("A must be (Conjugate)Transpose'd");
    if( orientationOfB == NORMAL )
        throw std::logic_error("B must be (Conjugate)Transpose'd");
    if( A.Width() != C.Height() || B.Height() != C.Width() ||
        A.Height() != B.Width() || A.Width() != B.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrrk: \n"
            << "  A ~ " << A.Height() << " x "
                        << A.Width()  << "\n"
            << "  B ~ " << B.Height() << " x "
                        << B.Width()  << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// Distributed C := alpha A^{T/H} B^{T/H} + beta C
template<typename T>
void CheckInput
( Orientation orientationOfA,
  Orientation orientationOfB,
  const DistMatrix<T,STAR,MC  >& A,
  const DistMatrix<T,MR,  STAR>& B,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA == NORMAL )
        throw std::logic_error("A[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfB == NORMAL )
        throw std::logic_error("B[MR,* ] must be (Conjugate)Transpose'd");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A.Width() != C.Height() || B.Height() != C.Width() ||
        A.Height() != B.Width() || A.Width() != B.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrrk: \n"
            << "  A[* ,MC] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[MR,* ] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.RowAlignment() != C.ColAlignment() ||
        B.ColAlignment() != C.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrrk: \n"
            << "  A[* ,MC] ~ " << A.RowAlignment() << "\n"
            << "  B[MR,* ] ~ " << B.ColAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}
#endif // ifndef RELEASE

// Local C := alpha A B + beta C
template<typename T>
inline void
TrrkNNKernel
( UpperOrLower uplo, 
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("TrrkNNKernel");
    CheckInputNN( A, B, C );
#endif
    Matrix<T> AT,
              AB;
    Matrix<T> BL, BR;
    Matrix<T> CTL, CTR,
              CBL, CBR;
    Matrix<T> DTL, DBR;

    const int half = C.Height()/2;
    ScaleTrapezoid( beta, LEFT, uplo, 0, C );
    LockedPartitionDown
    ( A, AT,
         AB, half );
    LockedPartitionRight( B, BL, BR, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    //------------------------------------------------------------------------//
    if( uplo == LOWER )
        Gemm( NORMAL, NORMAL, alpha, AB, BL, T(1), CBL );
    else
        Gemm( NORMAL, NORMAL, alpha, AT, BR, T(1), CTR );

    Zeros( CTL.Height(), CTL.Width(), DTL );
    Gemm( NORMAL, NORMAL, alpha, AT, BL, T(0), DTL );
    AxpyTriangle( uplo, T(1), DTL, CTL );

    Zeros( CBR.Height(), CBR.Width(), DBR );
    Gemm( NORMAL, NORMAL, alpha, AB, BR, T(0), DBR );
    AxpyTriangle( uplo, T(1), DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// Distributed C := alpha A B + beta C
template<typename T>
inline void
LocalTrrkKernel
( UpperOrLower uplo, 
  T alpha, const DistMatrix<T,MC,  STAR>& A,
           const DistMatrix<T,STAR,MR  >& B,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTrrkKernel");
    CheckInput( A, B, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,MC,STAR> AT(g), 
                          AB(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g);
    DistMatrix<T> CTL(g), CTR(g),
                  CBL(g), CBR(g);
    DistMatrix<T> DTL(g), DBR(g);

    const int half = C.Height()/2;
    ScaleTrapezoid( beta, LEFT, uplo, 0, C );
    LockedPartitionDown
    ( A, AT,
         AB, half );
    LockedPartitionRight( B, BL, BR, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    DTL.AlignWith( CTL );
    DBR.AlignWith( CBR );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
        LocalGemm( NORMAL, NORMAL, alpha, AB, BL, T(1), CBL );
    else
        LocalGemm( NORMAL, NORMAL, alpha, AT, BR, T(1), CTR );

    Zeros( CTL.Height(), CTL.Width(), DTL );
    LocalGemm( NORMAL, NORMAL, alpha, AT, BL, T(0), DTL );
    AxpyTriangle( uplo, T(1), DTL, CTL );

    Zeros( CBR.Height(), CBR.Width(), DBR );
    LocalGemm( NORMAL, NORMAL, alpha, AB, BR, T(0), DBR );
    AxpyTriangle( uplo, T(1), DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// Local C := alpha A B^{T/H} + beta C
template<typename T>
inline void
TrrkNTKernel
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("TrrkNTKernel");
    CheckInputNT( orientationOfB, A, B, C );
#endif
    Matrix<T> AT,
              AB;
    Matrix<T> BT,
              BB;
    Matrix<T> CTL, CTR,
              CBL, CBR;
    Matrix<T> DTL, DBR;

    const int half = C.Height()/2;
    ScaleTrapezoid( beta, LEFT, uplo, 0, C );
    LockedPartitionDown
    ( A, AT,
         AB, half );
    LockedPartitionDown
    ( B, BT, 
         BB, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    //------------------------------------------------------------------------//
    if( uplo == LOWER )
        Gemm( NORMAL, orientationOfB, alpha, AB, BT, T(1), CBL );
    else
        Gemm( NORMAL, orientationOfB, alpha, AT, BB, T(1), CTR );

    Zeros( CTL.Height(), CTL.Width(), DTL );
    Gemm( NORMAL, orientationOfB, alpha, AT, BT, T(0), DTL );
    AxpyTriangle( uplo, T(1), DTL, CTL );

    Zeros( CBR.Height(), CBR.Width(), DBR );
    Gemm( NORMAL, orientationOfB, alpha, AB, BB, T(0), DBR );
    AxpyTriangle( uplo, T(1), DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// Distributed C := alpha A B^{T/H} + beta C
template<typename T>
inline void
LocalTrrkKernel
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,STAR>& A,
           const DistMatrix<T,MR,STAR>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("LocalTrrkKernel");
    CheckInput( orientationOfB, A, B, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,MC,STAR> AT(g),
                          AB(g);
    DistMatrix<T,MR,STAR> BT(g), 
                          BB(g);
    DistMatrix<T> CTL(g), CTR(g),
                  CBL(g), CBR(g);
    DistMatrix<T> DTL(g), DBR(g);

    const int half = C.Height()/2;
    ScaleTrapezoid( beta, LEFT, uplo, 0, C );
    LockedPartitionDown
    ( A, AT,
         AB, half );
    LockedPartitionDown
    ( B, BT, 
         BB, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    DTL.AlignWith( CTL );
    DBR.AlignWith( CBR );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
        LocalGemm( NORMAL, orientationOfB, alpha, AB, BT, T(1), CBL );
    else
        LocalGemm( NORMAL, orientationOfB, alpha, AT, BB, T(1), CTR );

    Zeros( CTL.Height(), CTL.Width(), DTL );
    LocalGemm( NORMAL, orientationOfB, alpha, AT, BT, T(0), DTL );
    AxpyTriangle( uplo, T(1), DTL, CTL );

    Zeros( CBR.Height(), CBR.Width(), DBR );
    LocalGemm( NORMAL, orientationOfB, alpha, AB, BB, T(0), DBR );
    AxpyTriangle( uplo, T(1), DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// Local C := alpha A^{T/H} B + beta C
template<typename T>
inline void
TrrkTNKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("TrrkTNKernel");
    CheckInputTN( orientationOfA, A, B, C );
#endif
    Matrix<T> AL, AR;
    Matrix<T> BL, BR;
    Matrix<T> CTL, CTR,
              CBL, CBR;
    Matrix<T> DTL, DBR;

    const int half = C.Height()/2;
    ScaleTrapezoid( beta, LEFT, uplo, 0, C );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionRight( B, BL, BR, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    //------------------------------------------------------------------------//
    if( uplo == LOWER )
        Gemm( orientationOfA, NORMAL, alpha, AR, BL, T(1), CBL );
    else
        Gemm( orientationOfA, NORMAL, alpha, AL, BR, T(1), CTR );

    Zeros( CTL.Height(), CTL.Width(), DTL );
    Gemm( orientationOfA, NORMAL, alpha, AL, BL, T(0), DTL );
    AxpyTriangle( uplo, T(1), DTL, CTL );

    Zeros( CBR.Height(), CBR.Width(), DBR );
    Gemm( orientationOfA, NORMAL, alpha, AR, BR, T(0), DBR );
    AxpyTriangle( uplo, T(1), DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// Distributed C := alpha A^{T/H} B + beta C
template<typename T>
inline void
LocalTrrkKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T,STAR,MC>& A,
           const DistMatrix<T,STAR,MR>& B,
  T beta,        DistMatrix<T,MC,  MR>& C )
{
#ifndef RELEASE
    PushCallStack("LocalTrrkKernel");
    CheckInput( orientationOfA, A, B, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g);
    DistMatrix<T> CTL(g), CTR(g),
                  CBL(g), CBR(g);
    DistMatrix<T> DTL(g), DBR(g);

    const int half = C.Height()/2;
    ScaleTrapezoid( beta, LEFT, uplo, 0, C );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionRight( B, BL, BR, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    DTL.AlignWith( CTL );
    DBR.AlignWith( CBR );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
        LocalGemm( orientationOfA, NORMAL, alpha, AR, BL, T(1), CBL );
    else
        LocalGemm( orientationOfA, NORMAL, alpha, AL, BR, T(1), CTR );

    Zeros( CTL.Height(), CTL.Width(), DTL );
    LocalGemm( orientationOfA, NORMAL, alpha, AL, BL, T(0), DTL );
    AxpyTriangle( uplo, T(1), DTL, CTL );

    Zeros( CBR.Height(), CBR.Width(), DBR );
    LocalGemm( orientationOfA, NORMAL, alpha, AR, BR, T(0), DBR );
    AxpyTriangle( uplo, T(1), DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// Local C := alpha A^{T/H} B^{T/H} + beta C
template<typename T>
inline void
TrrkTTKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("TrrkTTKernel");
    CheckInputTT( orientationOfA, orientationOfB, A, B, C );
#endif
    Matrix<T> AL, AR;
    Matrix<T> BT,
              BB;
    Matrix<T> CTL, CTR,
              CBL, CBR;
    Matrix<T> DTL, DBR;

    const int half = C.Height()/2;
    ScaleTrapezoid( beta, LEFT, uplo, 0, C );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionDown
    ( B, BT, 
         BB, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    //------------------------------------------------------------------------//
    if( uplo == LOWER )
        Gemm( orientationOfA, orientationOfB, alpha, AR, BT, T(1), CBL );
    else
        Gemm( orientationOfA, orientationOfB, alpha, AL, BB, T(1), CTR );

    Zeros( CTL.Height(), CTL.Width(), DTL );
    Gemm( orientationOfA, orientationOfB, alpha, AL, BT, T(0), DTL );
    AxpyTriangle( uplo, T(1), DTL, CTL );

    Zeros( CBR.Height(), CBR.Width(), DBR );
    Gemm( orientationOfA, orientationOfB, alpha, AR, BB, T(0), DBR );
    AxpyTriangle( uplo, T(1), DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// Distributed C := alpha A^{T/H} B^{T/H} + beta C
template<typename T>
inline void
LocalTrrkKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,STAR,MC  >& A,
           const DistMatrix<T,MR,  STAR>& B,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTrrkKernel");
    CheckInput( orientationOfA, orientationOfB, A, B, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g);
    DistMatrix<T,MR,STAR> BT(g), 
                          BB(g);
    DistMatrix<T> CTL(g), CTR(g),
                  CBL(g), CBR(g);
    DistMatrix<T> DTL(g), DBR(g);

    const int half = C.Height()/2;
    ScaleTrapezoid( beta, LEFT, uplo, 0, C );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionDown
    ( B, BT, 
         BB, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    DTL.AlignWith( CTL );
    DBR.AlignWith( CBR );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
        LocalGemm( orientationOfA, orientationOfB, alpha, AR, BT, T(1), CBL );
    else
        LocalGemm( orientationOfA, orientationOfB, alpha, AL, BB, T(1), CTR );

    Zeros( CTL.Height(), CTL.Width(), DTL );
    LocalGemm( orientationOfA, orientationOfB, alpha, AL, BT, T(0), DTL );
    AxpyTriangle( uplo, T(1), DTL, CTL );

    Zeros( CBR.Height(), CBR.Width(), DBR );
    LocalGemm( orientationOfA, orientationOfB, alpha, AR, BB, T(0), DBR );
    AxpyTriangle( uplo, T(1), DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace trrk

namespace internal {

// Local C := alpha A B + beta C
template<typename T>
void TrrkNN
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
    using namespace trrk;
#ifndef RELEASE
    PushCallStack("internal::TrrkNN");
    CheckInputNN( A, B, C );
#endif
    if( C.Height() < LocalTrrkBlocksize<T>() )
    {
        TrrkNNKernel( uplo, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        Matrix<T> AT,
                  AB;
        Matrix<T> BL, BR;
        Matrix<T> CTL, CTR,
                  CBL, CBR;

        const int half = C.Height() / 2;
        LockedPartitionDown
        ( A, AT,
             AB, half );
        LockedPartitionRight( B, BL, BR, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            Gemm( NORMAL, NORMAL, alpha, AB, BL, beta, CBL );
        else
            Gemm( NORMAL, NORMAL, alpha, AT, BR, beta, CTR );

        // Recurse
        TrrkNN( uplo, alpha, AT, BL, beta, CTL );
        TrrkNN( uplo, alpha, AB, BR, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Local C := alpha A B^{T/H} + beta C
template<typename T>
void TrrkNT
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
    using namespace trrk;
#ifndef RELEASE
    PushCallStack("internal::TrrkNT");
    CheckInputNT( orientationOfB, A, B, C );
#endif
    if( C.Height() < LocalTrrkBlocksize<T>() )
    {
        TrrkNTKernel( uplo, orientationOfB, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        Matrix<T> AT,
                  AB;
        Matrix<T> BT,
                  BB;
        Matrix<T> CTL, CTR,
                  CBL, CBR;

        const int half = C.Height() / 2;
        LockedPartitionDown
        ( A, AT,
             AB, half );
        LockedPartitionDown
        ( B, BT, 
             BB, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            Gemm( NORMAL, orientationOfB, alpha, AB, BT, beta, CBL );
        else
            Gemm( NORMAL, orientationOfB, alpha, AT, BB, beta, CTR );

        // Recurse
        TrrkNT( uplo, orientationOfB, alpha, AT, BT, beta, CTL );
        TrrkNT( uplo, orientationOfB, alpha, AB, BB, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Local C := alpha A^{T/H} B + beta C
template<typename T>
void TrrkTN
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
    using namespace trrk;
#ifndef RELEASE
    PushCallStack("internal::TrrkTN");
    CheckInputTN( orientationOfA, A, B, C );
#endif
    if( C.Height() < LocalTrrkBlocksize<T>() )
    {
        TrrkTNKernel( uplo, orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        Matrix<T> AL, AR;
        Matrix<T> BL, BR;
        Matrix<T> CTL, CTR,
                  CBL, CBR;

        const int half = C.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionRight( B, BL, BR, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            Gemm( orientationOfA, NORMAL, alpha, AR, BL, beta, CBL );
        else
            Gemm( orientationOfA, NORMAL, alpha, AL, BR, beta, CTR );

        // Recurse
        TrrkTN( uplo, orientationOfA, alpha, AL, BL, beta, CTL );
        TrrkTN( uplo, orientationOfA, alpha, AR, BR, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Local C := alpha A^{T/H} B^{T/H} + beta C
template<typename T>
void TrrkTT
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
    using namespace trrk;
#ifndef RELEASE
    PushCallStack("internal::TrrkTT");
    CheckInputTT( orientationOfA, orientationOfB, A, B, C );
#endif
    if( C.Height() < LocalTrrkBlocksize<T>() )
    {
        TrrkTTKernel
        ( uplo, orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        Matrix<T> AL, AR;
        Matrix<T> BT,
                  BB;
        Matrix<T> CTL, CTR,
                  CBL, CBR;

        const int half = C.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionDown
        ( B, BT, 
             BB, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            Gemm( orientationOfA, orientationOfB, alpha, AR, BT, beta, CBL );
        else
            Gemm( orientationOfA, orientationOfB, alpha, AL, BB, beta, CTR );

        // Recurse
        TrrkTT
        ( uplo, orientationOfA, orientationOfB, alpha, AL, BT, beta, CTL );
        TrrkTT
        ( uplo, orientationOfA, orientationOfB, alpha, AR, BB, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal

// Distributed C := alpha A B + beta C
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,  STAR>& A,
           const DistMatrix<T,STAR,MR  >& B,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
    using namespace trrk;
#ifndef RELEASE
    PushCallStack("LocalTrrk");
    CheckInput( A, B, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTrrkBlocksize<T>() )
    {
        LocalTrrkKernel( uplo, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        DistMatrix<T,MC,STAR> AT(g),
                              AB(g);
        DistMatrix<T,STAR,MR> BL(g), BR(g);
        DistMatrix<T> CTL(g), CTR(g),
                      CBL(g), CBR(g);

        const int half = C.Height() / 2;
        LockedPartitionDown
        ( A, AT,
             AB, half );
        LockedPartitionRight( B, BL, BR, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            LocalGemm( NORMAL, NORMAL, alpha, AB, BL, beta, CBL );
        else
            LocalGemm( NORMAL, NORMAL, alpha, AT, BR, beta, CTR );

        // Recurse
        LocalTrrk( uplo, alpha, AT, BL, beta, CTL );
        LocalTrrk( uplo, alpha, AB, BR, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
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
    using namespace trrk;
#ifndef RELEASE
    PushCallStack("LocalTrrk");
    CheckInput( orientationOfB, A, B, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTrrkBlocksize<T>() )
    {
        LocalTrrkKernel( uplo, orientationOfB, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        DistMatrix<T,MC,STAR> AT(g),
                              AB(g);
        DistMatrix<T,MR,STAR> BT(g), 
                              BB(g);
        DistMatrix<T> CTL(g), CTR(g),
                      CBL(g), CBR(g);

        const int half = C.Height() / 2;
        LockedPartitionDown
        ( A, AT,
             AB, half );
        LockedPartitionDown
        ( B, BT, 
             BB, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            LocalGemm( NORMAL, orientationOfB, alpha, AB, BT, beta, CBL );
        else
            LocalGemm( NORMAL, orientationOfB, alpha, AT, BB, beta, CTR );

        // Recurse
        LocalTrrk( uplo, orientationOfB, alpha, AT, BT, beta, CTL );
        LocalTrrk( uplo, orientationOfB, alpha, AB, BB, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Distributed C := alpha A^{T/H} B + beta C
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T,STAR,MC>& A,
           const DistMatrix<T,STAR,MR>& B,
  T beta,        DistMatrix<T,MC,  MR>& C )
{
    using namespace trrk;
#ifndef RELEASE
    PushCallStack("LocalTrrk");
    CheckInput( orientationOfA, A, B, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTrrkBlocksize<T>() )
    {
        LocalTrrkKernel( uplo, orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g);
        DistMatrix<T,STAR,MR> BL(g), BR(g);
        DistMatrix<T> CTL(g), CTR(g),
                      CBL(g), CBR(g);

        const int half = C.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionRight( B, BL, BR, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            LocalGemm( orientationOfA, NORMAL, alpha, AR, BL, beta, CBL );
        else
            LocalGemm( orientationOfA, NORMAL, alpha, AL, BR, beta, CTR );

        // Recurse
        LocalTrrk( uplo, orientationOfA, alpha, AL, BL, beta, CTL );
        LocalTrrk( uplo, orientationOfA, alpha, AR, BR, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Distributed C := alpha A^{T/H} B^{T/H} + beta C
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,STAR,MC  >& A,
           const DistMatrix<T,MR,  STAR>& B,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
    using namespace trrk;
#ifndef RELEASE
    PushCallStack("LocalTrrk");
    CheckInput( orientationOfA, orientationOfB, A, B, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTrrkBlocksize<T>() )
    {
        LocalTrrkKernel
        ( uplo, orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g);
        DistMatrix<T,MR,STAR> BT(g), 
                              BB(g);
        DistMatrix<T> CTL(g), CTR(g),
                      CBL(g), CBR(g);

        const int half = C.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionDown
        ( B, BT, 
             BB, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
        { 
            LocalGemm
            ( orientationOfA, orientationOfB, alpha, AR, BT, beta, CBL );
        }
        else
        {
            LocalGemm
            ( orientationOfA, orientationOfB, alpha, AL, BB, beta, CTR );
        }

        // Recurse
        LocalTrrk
        ( uplo, orientationOfA, orientationOfB, alpha, AL, BT, beta, CTL );
        LocalTrrk
        ( uplo, orientationOfA, orientationOfB, alpha, AR, BB, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef BLAS_TRRK_LOCAL_HPP
