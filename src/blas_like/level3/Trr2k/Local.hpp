/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRR2K_LOCAL_HPP
#define EL_TRR2K_LOCAL_HPP

namespace El {

namespace trr2k {

// TODO: Fuse pairs of Gemms

// E := alpha op(A) op(B) + beta op(C) op(D) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientA, Orientation orientB,
  Orientation orientC, Orientation orientD,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  T beta,  const AbstractDistMatrix<T>& C, const AbstractDistMatrix<T>& D,
  T gamma,       AbstractDistMatrix<T>& E )
{
    DEBUG_ONLY(CallStackEntry cse("LocalTrr2kKernel"))

    const bool transA = orientA != NORMAL;
    const bool transB = orientB != NORMAL;
    const bool transC = orientC != NORMAL;
    const bool transD = orientD != NORMAL;
    // TODO: Stringent distribution and alignment checks

    typedef AbstractDistMatrix<T> ADM;
    auto A0 = std::unique_ptr<ADM>( A.Construct(A.Grid(),A.Root()) );
    auto A1 = std::unique_ptr<ADM>( A.Construct(A.Grid(),A.Root()) );
    auto B0 = std::unique_ptr<ADM>( B.Construct(B.Grid(),B.Root()) );
    auto B1 = std::unique_ptr<ADM>( B.Construct(B.Grid(),B.Root()) );
    auto C0 = std::unique_ptr<ADM>( C.Construct(C.Grid(),C.Root()) );
    auto C1 = std::unique_ptr<ADM>( C.Construct(C.Grid(),C.Root()) );
    auto D0 = std::unique_ptr<ADM>( D.Construct(D.Grid(),D.Root()) );
    auto D1 = std::unique_ptr<ADM>( D.Construct(D.Grid(),D.Root()) );
    auto ETL = std::unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );
    auto ETR = std::unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );
    auto EBL = std::unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );
    auto EBR = std::unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );
    auto FTL = std::unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );
    auto FBR = std::unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );

    const Int half = E.Height() / 2;
    if( transA )
        LockedPartitionRight( A, *A0, *A1, half );
    else
        LockedPartitionDown( A, *A0, *A1, half );
    if( transB )
        LockedPartitionDown( B, *B0, *B1, half );
    else
        LockedPartitionRight( B, *B0, *B1, half );
    if( transC )
        LockedPartitionRight( C, *C0, *C1, half );
    else
        LockedPartitionDown( C, *C0, *C1, half );
    if( transD )
        LockedPartitionDown( D, *D0, *D1, half );
    else
        LockedPartitionRight( D, *D0, *D1, half );
    PartitionDownDiagonal( E, *ETL, *ETR, *EBL, *EBR, half );

    ScaleTrapezoid( gamma, uplo, E );
    if( uplo == LOWER )
    {
        Gemm
        ( orientA, orientB, 
          alpha, A1->LockedMatrix(), B0->LockedMatrix(), 
          T(1), EBL->Matrix() );
        Gemm
        ( orientC, orientD, 
          beta, C1->LockedMatrix(), D0->LockedMatrix(), 
          T(1), EBL->Matrix() );
    }
    else
    {
        Gemm
        ( orientA, orientB, 
          alpha, A0->LockedMatrix(), B1->LockedMatrix(), 
          T(1), ETR->Matrix() );
        Gemm
        ( orientC, orientD, 
          beta, C0->LockedMatrix(), D1->LockedMatrix(), 
          T(1), ETR->Matrix() );
    }

    FTL->AlignWith( *ETL );
    FTL->Resize( ETL->Height(), ETL->Width() );
    Gemm
    ( orientA, orientB, 
      alpha, A0->LockedMatrix(), B0->LockedMatrix(),
      T(0), FTL->Matrix() );
    Gemm
    ( orientC, orientD,
      beta, C0->LockedMatrix(), D0->LockedMatrix(),
      T(1), FTL->Matrix() );
    AxpyTrapezoid( uplo, T(1), *FTL, *ETL );

    FBR->AlignWith( *EBR );
    FBR->Resize( EBR->Height(), EBR->Width() );
    Gemm
    ( orientA, orientB, 
      alpha, A1->LockedMatrix(), B1->LockedMatrix(),
      T(0), FBR->Matrix() );
    Gemm
    ( orientC, orientD,
      beta, C1->LockedMatrix(), D1->LockedMatrix(),
      T(1), FBR->Matrix() );
    AxpyTrapezoid( uplo, T(1), *FBR, *EBR );
}

} // namespace trr2k

// E := alpha op(A) op(B) + beta op(C) op(D) + beta E
template<typename T>
void LocalTrr2k
( UpperOrLower uplo, 
  Orientation orientA, Orientation orientB,
  Orientation orientC, Orientation orientD,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  T beta,  const AbstractDistMatrix<T>& C, const AbstractDistMatrix<T>& D,
  T gamma,       AbstractDistMatrix<T>& E )
{
    using namespace trr2k;
    DEBUG_ONLY(CallStackEntry cse("LocalTrr2k"))

    const bool transA = orientA != NORMAL;
    const bool transB = orientB != NORMAL;
    const bool transC = orientC != NORMAL;
    const bool transD = orientD != NORMAL;
    // TODO: Stringent distribution and alignment checks

    if( E.Height() < E.Grid().Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientA, orientB, orientC, orientD, 
          alpha, A, B, beta, C, D, gamma, E );
    }
    else
    {
        typedef AbstractDistMatrix<T> ADM;
        auto A0 = std::unique_ptr<ADM>( A.Construct(A.Grid(),A.Root()) );
        auto A1 = std::unique_ptr<ADM>( A.Construct(A.Grid(),A.Root()) );
        auto B0 = std::unique_ptr<ADM>( B.Construct(B.Grid(),B.Root()) );
        auto B1 = std::unique_ptr<ADM>( B.Construct(B.Grid(),B.Root()) );
        auto C0 = std::unique_ptr<ADM>( C.Construct(C.Grid(),C.Root()) );
        auto C1 = std::unique_ptr<ADM>( C.Construct(C.Grid(),C.Root()) );
        auto D0 = std::unique_ptr<ADM>( D.Construct(D.Grid(),D.Root()) );
        auto D1 = std::unique_ptr<ADM>( D.Construct(D.Grid(),D.Root()) );
        auto ETL = std::unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );
        auto ETR = std::unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );
        auto EBL = std::unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );
        auto EBR = std::unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );

        const Int half = E.Height() / 2;
        if( transA )
            LockedPartitionRight( A, *A0, *A1, half );
        else
            LockedPartitionDown( A, *A0, *A1, half );
        if( transB )
            LockedPartitionDown( B, *B0, *B1, half );
        else
            LockedPartitionRight( B, *B0, *B1, half );
        if( transC )
            LockedPartitionRight( C, *C0, *C1, half );
        else
            LockedPartitionDown( C, *C0, *C1, half );
        if( transD )
            LockedPartitionDown( D, *D0, *D1, half );
        else
            LockedPartitionRight( D, *D0, *D1, half );
        PartitionDownDiagonal( E, *ETL, *ETR, *EBL, *EBR, half );

        if( uplo == LOWER )
        { 
            Gemm
            ( orientA, orientB, 
              alpha, A1->LockedMatrix(), B0->LockedMatrix(), 
              gamma, EBL->Matrix() );
            Gemm
            ( orientC, orientD, 
              beta,  C1->LockedMatrix(), D0->LockedMatrix(), 
              T(1), EBL->Matrix() );
        }
        else
        {
            Gemm
            ( orientA, orientB,
              alpha, A0->LockedMatrix(), B1->LockedMatrix(), 
              gamma, ETR->Matrix() );
            Gemm
            ( orientC, orientD,
              beta,  C0->LockedMatrix(), D1->LockedMatrix(), 
              T(1), ETR->Matrix() );
        }

        // Recurse
        LocalTrr2k
        ( uplo, orientA, orientB, orientC, orientD, 
          alpha, *A0, *B0, beta, *C0, *D0, gamma, *ETL );
        LocalTrr2k
        ( uplo, orientA, orientB, orientC, orientD, 
          alpha, *A1, *B1, beta, *C1, *D1, gamma, *EBR );
    }
}

} // namespace El

#endif // ifndef EL_TRR2K_LOCAL_HPP
