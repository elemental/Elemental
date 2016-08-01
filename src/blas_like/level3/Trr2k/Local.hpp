/*
   Copyright (c) 2009-2016, Jack Poulson
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

// E := alpha op(A) op(B) + beta op(C) op(D) + E
template<typename T>
void LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientA, Orientation orientB,
  Orientation orientC, Orientation orientD,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  T beta,  const AbstractDistMatrix<T>& C, const AbstractDistMatrix<T>& D,
                 AbstractDistMatrix<T>& E )
{
    DEBUG_CSE

    const bool transA = orientA != NORMAL;
    const bool transB = orientB != NORMAL;
    const bool transC = orientC != NORMAL;
    const bool transD = orientD != NORMAL;
    // TODO: Stringent distribution and alignment checks

    typedef AbstractDistMatrix<T> ADM;
    auto A0 = unique_ptr<ADM>( A.Construct(A.Grid(),A.Root()) );
    auto A1 = unique_ptr<ADM>( A.Construct(A.Grid(),A.Root()) );
    auto B0 = unique_ptr<ADM>( B.Construct(B.Grid(),B.Root()) );
    auto B1 = unique_ptr<ADM>( B.Construct(B.Grid(),B.Root()) );
    auto C0 = unique_ptr<ADM>( C.Construct(C.Grid(),C.Root()) );
    auto C1 = unique_ptr<ADM>( C.Construct(C.Grid(),C.Root()) );
    auto D0 = unique_ptr<ADM>( D.Construct(D.Grid(),D.Root()) );
    auto D1 = unique_ptr<ADM>( D.Construct(D.Grid(),D.Root()) );
    auto ETL = unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );
    auto ETR = unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );
    auto EBL = unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );
    auto EBR = unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );

    const Int half = E.Height() / 2;
    const auto indTL = IR(0,half);
    const auto indBR = IR(half,END); 
    if( transA )
    {
        LockedView( *A0, A, ALL, indTL );
        LockedView( *A1, A, ALL, indBR );
    }
    else
    {
        LockedView( *A0, A, indTL, ALL );
        LockedView( *A1, A, indBR, ALL );
    }
    if( transB )
    {
        LockedView( *B0, B, indTL, ALL );
        LockedView( *B1, B, indBR, ALL );
    }
    else
    {
        LockedView( *B0, B, ALL, indTL );
        LockedView( *B1, B, ALL, indBR );
    }
    if( transC )
    {
        LockedView( *C0, C, ALL, indTL );
        LockedView( *C1, C, ALL, indBR );
    }
    else
    {
        LockedView( *C0, C, indTL, ALL );
        LockedView( *C1, C, indBR, ALL );
    }
    if( transD )
    {
        LockedView( *D0, D, indTL, ALL );
        LockedView( *D1, D, indBR, ALL );
    }
    else
    {
        LockedView( *D0, D, ALL, indTL );
        LockedView( *D1, D, ALL, indBR );
    }
    View( *ETL, E, indTL, indTL );
    View( *EBR, E, indBR, indBR );

    if( uplo == LOWER )
    {
        View( *EBL, E, indBR, indTL ); 
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
        View( *ETR, E, indTL, indBR );
        Gemm
        ( orientA, orientB, 
          alpha, A0->LockedMatrix(), B1->LockedMatrix(), 
          T(1), ETR->Matrix() );
        Gemm
        ( orientC, orientD, 
          beta, C0->LockedMatrix(), D1->LockedMatrix(), 
          T(1), ETR->Matrix() );
    }

    auto FTL = unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );
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
    LocalAxpyTrapezoid( uplo, T(1), *FTL, *ETL );

    auto FBR = unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );
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
    LocalAxpyTrapezoid( uplo, T(1), *FBR, *EBR );
}

} // namespace trr2k

// E := alpha op(A) op(B) + beta op(C) op(D) + gamma E
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
    DEBUG_CSE

    const bool transA = orientA != NORMAL;
    const bool transB = orientB != NORMAL;
    const bool transC = orientC != NORMAL;
    const bool transD = orientD != NORMAL;
    // TODO: Stringent distribution and alignment checks

    ScaleTrapezoid( gamma, uplo, E );
    if( E.Height() < E.Grid().Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientA, orientB, orientC, orientD, 
          alpha, A, B, beta, C, D, E );
    }
    else
    {
        typedef AbstractDistMatrix<T> ADM;
        // Ugh. This is likely too much overhead. It should be measured.
        auto A0 = unique_ptr<ADM>( A.Construct(A.Grid(),A.Root()) );
        auto A1 = unique_ptr<ADM>( A.Construct(A.Grid(),A.Root()) );
        auto B0 = unique_ptr<ADM>( B.Construct(B.Grid(),B.Root()) );
        auto B1 = unique_ptr<ADM>( B.Construct(B.Grid(),B.Root()) );
        auto C0 = unique_ptr<ADM>( C.Construct(C.Grid(),C.Root()) );
        auto C1 = unique_ptr<ADM>( C.Construct(C.Grid(),C.Root()) );
        auto D0 = unique_ptr<ADM>( D.Construct(D.Grid(),D.Root()) );
        auto D1 = unique_ptr<ADM>( D.Construct(D.Grid(),D.Root()) );
        auto ETL = unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );
        auto ETR = unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );
        auto EBL = unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );
        auto EBR = unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );

        const Int half = E.Height() / 2;
        const auto indTL = IR(0,half);
        const auto indBR = IR(half,END); 
        if( transA )
        {
            LockedView( *A0, A, ALL, indTL );
            LockedView( *A1, A, ALL, indBR );
        }
        else
        {
            LockedView( *A0, A, indTL, ALL );
            LockedView( *A1, A, indBR, ALL );
        }
        if( transB )
        {
            LockedView( *B0, B, indTL, ALL );
            LockedView( *B1, B, indBR, ALL );
        }
        else
        {
            LockedView( *B0, B, ALL, indTL );
            LockedView( *B1, B, ALL, indBR );
        }
        if( transC )
        {
            LockedView( *C0, C, ALL, indTL );
            LockedView( *C1, C, ALL, indBR );
        }
        else
        {
            LockedView( *C0, C, indTL, ALL );
            LockedView( *C1, C, indBR, ALL );
        }
        if( transD )
        {
            LockedView( *D0, D, indTL, ALL );
            LockedView( *D1, D, indBR, ALL );
        }
        else
        {
            LockedView( *D0, D, ALL, indTL );
            LockedView( *D1, D, ALL, indBR );
        }
        View( *ETL, E, indTL, indTL );
        View( *EBR, E, indBR, indBR );

        if( uplo == LOWER )
        { 
            View( *EBL, E, indBR, indTL ); 
            Gemm
            ( orientA, orientB, 
              alpha, A1->LockedMatrix(), B0->LockedMatrix(), 
              T(1), EBL->Matrix() );
            Gemm
            ( orientC, orientD, 
              beta,  C1->LockedMatrix(), D0->LockedMatrix(), 
              T(1), EBL->Matrix() );
        }
        else
        {
            View( *ETR, E, indTL, indBR );
            Gemm
            ( orientA, orientB,
              alpha, A0->LockedMatrix(), B1->LockedMatrix(), 
              T(1), ETR->Matrix() );
            Gemm
            ( orientC, orientD,
              beta,  C0->LockedMatrix(), D1->LockedMatrix(), 
              T(1), ETR->Matrix() );
        }

        // Recurse
        LocalTrr2k
        ( uplo, orientA, orientB, orientC, orientD, 
          alpha, *A0, *B0, beta, *C0, *D0, T(1), *ETL );
        LocalTrr2k
        ( uplo, orientA, orientB, orientC, orientD, 
          alpha, *A1, *B1, beta, *C1, *D1, T(1), *EBR );
    }
}

} // namespace El

#endif // ifndef EL_TRR2K_LOCAL_HPP
