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

// E := alpha op(A) op(B) + beta op(C) op(D) + E
template<typename Ring>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientA, Orientation orientB,
  Orientation orientC, Orientation orientD,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
  Ring beta,  const AbstractDistMatrix<Ring>& C, 
              const AbstractDistMatrix<Ring>& D,
                    AbstractDistMatrix<Ring>& E )
{
    DEBUG_ONLY(CSE cse("LocalTrr2kKernel"))

    const bool transA = orientA != NORMAL;
    const bool transB = orientB != NORMAL;
    const bool transC = orientC != NORMAL;
    const bool transD = orientD != NORMAL;
    // TODO: Stringent distribution and alignment checks

    typedef AbstractDistMatrix<Ring> ADM;
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
    auto FTL = unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );
    auto FBR = unique_ptr<ADM>( E.Construct(E.Grid(),E.Root()) );

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
        ( alpha, A1->LockedMatrix().Orient(orientA), 
                 B0->LockedMatrix().Orient(orientB), 
          Ring(1), EBL->Matrix() );
        Gemm
        ( beta, C1->LockedMatrix().Orient(orientC), 
                D0->LockedMatrix().Orient(orientD), 
          Ring(1), EBL->Matrix() );
    }
    else
    {
        Gemm
        ( alpha, A0->LockedMatrix().Orient(orientA), 
                 B1->LockedMatrix().Orient(orientB), 
          Ring(1), ETR->Matrix() );
        Gemm
        ( beta, C0->LockedMatrix().Orient(orientC), 
                D1->LockedMatrix().Orient(orientD), 
          Ring(1), ETR->Matrix() );
    }

    FTL->AlignWith( *ETL );
    FTL->Resize( ETL->Height(), ETL->Width() );
    Gemm
    ( alpha, A0->LockedMatrix().Orient(orientA), 
             B0->LockedMatrix().Orient(orientB),
      Ring(0), FTL->Matrix() );
    Gemm
    ( beta, C0->LockedMatrix().Orient(orientC), 
            D0->LockedMatrix().Orient(orientD),
      Ring(1), FTL->Matrix() );
    AxpyTrapezoid( uplo, Ring(1), *FTL, *ETL );

    FBR->AlignWith( *EBR );
    FBR->Resize( EBR->Height(), EBR->Width() );
    Gemm
    ( alpha, A1->LockedMatrix().Orient(orientA), 
             B1->LockedMatrix().Orient(orientB),
      Ring(0), FBR->Matrix() );
    Gemm
    ( beta, C1->LockedMatrix().Orient(orientC), 
            D1->LockedMatrix().Orient(orientD),
      Ring(1), FBR->Matrix() );
    AxpyTrapezoid( uplo, Ring(1), *FBR, *EBR );
}

} // namespace trr2k

// E := alpha op(A) op(B) + beta op(C) op(D) + gamma E
template<typename Ring>
void LocalTrr2k
( UpperOrLower uplo, 
  Orientation orientA, Orientation orientB,
  Orientation orientC, Orientation orientD,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
  Ring beta,  const AbstractDistMatrix<Ring>& C, 
              const AbstractDistMatrix<Ring>& D,
  Ring gamma,       AbstractDistMatrix<Ring>& E )
{
    using namespace trr2k;
    DEBUG_ONLY(CSE cse("LocalTrr2k"))

    const bool transA = orientA != NORMAL;
    const bool transB = orientB != NORMAL;
    const bool transC = orientC != NORMAL;
    const bool transD = orientD != NORMAL;
    // TODO: Stringent distribution and alignment checks

    ScaleTrapezoid( gamma, uplo, E );
    if( E.Height() < E.Grid().Width()*LocalTrr2kBlocksize<Ring>() )
    {
        LocalTrr2kKernel
        ( uplo, orientA, orientB, orientC, orientD, 
          alpha, A, B, beta, C, D, E );
    }
    else
    {
        typedef AbstractDistMatrix<Ring> ADM;
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
            ( alpha, A1->LockedMatrix().Orient(orientA), 
                     B0->LockedMatrix().Orient(orientB), 
              Ring(1), EBL->Matrix() );
            Gemm
            ( beta,  C1->LockedMatrix().Orient(orientC), 
                     D0->LockedMatrix().Orient(orientD), 
              Ring(1), EBL->Matrix() );
        }
        else
        {
            Gemm
            ( alpha, A0->LockedMatrix().Orient(orientA), 
                     B1->LockedMatrix().Orient(orientB), 
              Ring(1), ETR->Matrix() );
            Gemm
            ( beta,  C0->LockedMatrix().Orient(orientC), 
                     D1->LockedMatrix().Orient(orientD), 
              Ring(1), ETR->Matrix() );
        }

        // Recurse
        LocalTrr2k
        ( uplo, orientA, orientB, orientC, orientD, 
          alpha, *A0, *B0, beta, *C0, *D0, Ring(1), *ETL );
        LocalTrr2k
        ( uplo, orientA, orientB, orientC, orientD, 
          alpha, *A1, *B1, beta, *C1, *D1, Ring(1), *EBR );
    }
}

} // namespace El

#endif // ifndef EL_TRR2K_LOCAL_HPP
