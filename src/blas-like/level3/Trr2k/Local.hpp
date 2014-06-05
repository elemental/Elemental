/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRR2K_LOCAL_HPP
#define EL_TRR2K_LOCAL_HPP

namespace El {

namespace trr2k {

#ifndef EL_RELEASE

void EnsureSame
( const Grid& gA, const Grid& gB, const Grid& gC, 
  const Grid& gD, const Grid& gE )
{
    if( gA != gB || gB != gC || gC != gD || gD != gE )
        LogicError("Grids must be the same");
}

template<typename T>
void EnsureConformal
( const DistMatrix<T,MC,STAR>& A, const DistMatrix<T>& E, std::string name )
{
    if( A.Height() != E.Height() || A.ColAlign() != E.ColAlign() )
        LogicError(name," not conformal with E");
}

template<typename T>
void EnsureConformal
( const DistMatrix<T,STAR,MC>& A, const DistMatrix<T>& E, std::string name )
{
    if( A.Width() != E.Height() || A.RowAlign() != E.ColAlign() )
        LogicError(name," not conformal with E");
}

template<typename T>
void EnsureConformal
( const DistMatrix<T,MR,STAR>& A, const DistMatrix<T>& E, std::string name )
{
    if( A.Height() != E.Width() || A.ColAlign() != E.RowAlign() )
        LogicError(name," not conformal with E");
}

template<typename T>
void EnsureConformal
( const DistMatrix<T,STAR,MR>& A, const DistMatrix<T>& E, std::string name )
{
    if( A.Width() != E.Width() || A.RowAlign() != E.RowAlign() )
        LogicError(name," not conformal with E");
}

template<typename T,Distribution UA,Distribution VA,
                    Distribution UB,Distribution VB,
                    Distribution UC,Distribution VC,
                    Distribution UD,Distribution VD>
void CheckInput
( const DistMatrix<T,UA,VA>& A, const DistMatrix<T,UB,VB>& B, 
  const DistMatrix<T,UC,VC>& C, const DistMatrix<T,UD,VD>& D,
  const DistMatrix<T>& E )
{
    EnsureSame( A.Grid(), B.Grid(), C.Grid(), D.Grid(), E.Grid() );
    EnsureConformal( A, E, "A" );
    EnsureConformal( B, E, "B" );
    EnsureConformal( C, E, "C" );
    EnsureConformal( D, E, "D" );
}

#endif // ifndef EL_RELEASE

// E := alpha (A B + C D) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2kKernel");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    DistMatrix<T,MC,STAR> AT(g),  CT(g),
                          AB(g),  CB(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g),
                          DL(g), DR(g);
    DistMatrix<T> ETL(g), ETR(g),
                  EBL(g), EBR(g);
    DistMatrix<T> FTL(g), FBR(g);

    const Int half = E.Height()/2;
    ScaleTrapezoid( beta, uplo, E );
    LockedPartitionDown( A, AT, AB, half );
    LockedPartitionRight( B, BL, BR, half );
    LockedPartitionDown( C, CT, CB, half );
    LockedPartitionRight( D, DL, DR, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    if( uplo == LOWER )
    {
        LocalGemm( NORMAL, NORMAL, alpha, AB, BL, T(1), EBL );
        LocalGemm( NORMAL, NORMAL, alpha, CB, DL, T(1), EBL );
    }
    else
    {
        LocalGemm( NORMAL, NORMAL, alpha, AT, BR, T(1), ETR );
        LocalGemm( NORMAL, NORMAL, alpha, CT, DR, T(1), ETR );
    }

    FTL.AlignWith( ETL );
    LocalGemm( NORMAL, NORMAL, alpha, AT, BL, FTL );
    LocalGemm( NORMAL, NORMAL, alpha, CT, DL, T(1), FTL );
    AxpyTriangle( uplo, T(1), FTL, ETL );

    FBR.AlignWith( EBR );
    LocalGemm( NORMAL, NORMAL, alpha, AB, BR, FBR );
    LocalGemm( NORMAL, NORMAL, alpha, CB, DR, T(1), FBR );
    AxpyTriangle( uplo, T(1), FBR, EBR );
}

// E := alpha (A B + C D^{T/H}) + beta C
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo, Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2kKernel");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    DistMatrix<T,MC,STAR> AT(g),  CT(g),
                          AB(g),  CB(g);
    DistMatrix<T,MR,STAR> DT(g), 
                          DB(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g);
    DistMatrix<T> ETL(g), ETR(g),
                  EBL(g), EBR(g);
    DistMatrix<T> FTL(g), FBR(g);

    const Int half = E.Height()/2;
    ScaleTrapezoid( beta, uplo, E );
    LockedPartitionDown( A, AT, AB, half );
    LockedPartitionRight( B, BL, BR, half );
    LockedPartitionDown( C, CT, CB, half );
    LockedPartitionDown( D, DT, DB, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    if( uplo == LOWER )
    {
        LocalGemm( NORMAL, NORMAL, alpha, AB, BL, T(1), EBL );
        LocalGemm( NORMAL, orientationOfD, alpha, CB, DT, T(1), EBL );
    }
    else
    {
        LocalGemm( NORMAL, NORMAL, alpha, AT, BR, T(1), ETR );
        LocalGemm( NORMAL, orientationOfD, alpha, CT, DB, T(1), ETR );
    }

    FTL.AlignWith( ETL );
    LocalGemm( NORMAL, NORMAL, alpha, AT, BL, FTL );
    LocalGemm( NORMAL, orientationOfD, alpha, CT, DT, T(1), FTL );
    AxpyTriangle( uplo, T(1), FTL, ETL );

    FBR.AlignWith( EBR );
    LocalGemm( NORMAL, NORMAL, alpha, AB, BR, FBR );
    LocalGemm( NORMAL, orientationOfD, alpha, CB, DB, T(1), FBR );
    AxpyTriangle( uplo, T(1), FBR, EBR );
}

// E := alpha (A B + C^{T/H} D) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo, Orientation orientationOfC,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2kKernel");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    DistMatrix<T,MC,STAR> AT(g), AB(g);
    DistMatrix<T,STAR,MC> CL(g), CR(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g),
                          DL(g), DR(g);
    DistMatrix<T> ETL(g), ETR(g),
                  EBL(g), EBR(g);
    DistMatrix<T> FTL(g), FBR(g);

    const Int half = E.Height()/2;
    ScaleTrapezoid( beta, uplo, E );
    LockedPartitionDown( A, AT, AB, half );
    LockedPartitionRight( B, BL, BR, half );
    LockedPartitionRight( C, CL, CR, half );
    LockedPartitionRight( D, DL, DR, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    if( uplo == LOWER )
    {
        LocalGemm( NORMAL, NORMAL, alpha, AB, BL, T(1), EBL );
        LocalGemm( orientationOfC, NORMAL, alpha, CR, DL, T(1), EBL );
    }
    else
    {
        LocalGemm( NORMAL, NORMAL, alpha, AT, BR, T(1), ETR );
        LocalGemm( orientationOfC, NORMAL, alpha, CL, DR, T(1), ETR );
    }

    FTL.AlignWith( ETL );
    LocalGemm( NORMAL, NORMAL, alpha, AT, BL, FTL );
    LocalGemm( orientationOfC, NORMAL, alpha, CL, DL, T(1), FTL );
    AxpyTriangle( uplo, T(1), FTL, ETL );

    FBR.AlignWith( EBR );
    LocalGemm( NORMAL, NORMAL, alpha, AB, BR, FBR );
    LocalGemm( orientationOfC, NORMAL, alpha, CR, DR, T(1), FBR );
    AxpyTriangle( uplo, T(1), FBR, EBR );
}

// E := alpha (A B + C^{T/H} D^{T/H}) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo, Orientation orientationOfC, Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2kKernel");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    DistMatrix<T,MC,STAR> AT(g), AB(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g);
    DistMatrix<T,STAR,MC> CL(g), CR(g);
    DistMatrix<T,MR,STAR> DT(g), DB(g);
    DistMatrix<T> ETL(g), ETR(g),
                  EBL(g), EBR(g);
    DistMatrix<T> FTL(g), FBR(g);

    const Int half = E.Height()/2;
    ScaleTrapezoid( beta, uplo, E );
    LockedPartitionDown( A, AT, AB, half );
    LockedPartitionRight( B, BL, BR, half );
    LockedPartitionRight( C, CL, CR, half );
    LockedPartitionDown( D, DT, DB, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    if( uplo == LOWER )
    {
        LocalGemm( NORMAL, NORMAL, alpha, AB, BL, T(1), EBL );
        LocalGemm( orientationOfC, orientationOfD, alpha, CR, DT, T(1), EBL );
    }
    else
    {
        LocalGemm( NORMAL, NORMAL, alpha, AT, BR, T(1), ETR );
        LocalGemm( orientationOfC, orientationOfD, alpha, CL, DB, T(1), ETR );
    }

    FTL.AlignWith( ETL );
    LocalGemm( NORMAL, NORMAL, alpha, AT, BL, FTL );
    LocalGemm( orientationOfC, orientationOfD, alpha, CL, DT, T(1), FTL );
    AxpyTriangle( uplo, T(1), FTL, ETL );

    FBR.AlignWith( EBR );
    LocalGemm( NORMAL, NORMAL, alpha, AB, BR, FBR );
    LocalGemm( orientationOfC, orientationOfD, alpha, CR, DB, T(1), FBR );
    AxpyTriangle( uplo, T(1), FBR, EBR );
}

// E := alpha (A B^{T/H} + C D) + beta C
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo, Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2kKernel");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    DistMatrix<T,MC,STAR> AT(g),  CT(g),
                          AB(g),  CB(g);
    DistMatrix<T,MR,STAR> BT(g), BB(g);
    DistMatrix<T,STAR,MR> DL(g), DR(g);
    DistMatrix<T> ETL(g), ETR(g),
                  EBL(g), EBR(g);
    DistMatrix<T> FTL(g), FBR(g);

    const Int half = E.Height()/2;
    ScaleTrapezoid( beta, uplo, E );
    LockedPartitionDown( A, AT, AB, half );
    LockedPartitionDown( B, BT, BB, half );
    LockedPartitionDown( C, CT, CB, half );
    LockedPartitionRight( D, DL, DR, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    if( uplo == LOWER )
    {
        LocalGemm( NORMAL, orientationOfB, alpha, AB, BT, T(1), EBL );
        LocalGemm( NORMAL, NORMAL, alpha, CB, DL, T(1), EBL );
    }
    else
    {
        LocalGemm( NORMAL, orientationOfB, alpha, AT, BB, T(1), ETR );
        LocalGemm( NORMAL, NORMAL, alpha, CT, DR, T(1), ETR );
    }

    FTL.AlignWith( ETL );
    LocalGemm( NORMAL, orientationOfB, alpha, AT, BT, FTL );
    LocalGemm( NORMAL, NORMAL, alpha, CT, DL, T(1), FTL );
    AxpyTriangle( uplo, T(1), FTL, ETL );

    FBR.AlignWith( EBR );
    LocalGemm( NORMAL, orientationOfB, alpha, AB, BB, FBR );
    LocalGemm( NORMAL, NORMAL, alpha, CB, DR, T(1), FBR );
    AxpyTriangle( uplo, T(1), FBR, EBR );
}

// E := alpha (A B^{T/H} + C D^{T/H}) + beta C
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo, Orientation orientationOfB, Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2kKernel");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    DistMatrix<T,MC,STAR> AT(g),  CT(g),
                          AB(g),  CB(g);
    DistMatrix<T,MR,STAR> BT(g),  DT(g),
                          BB(g),  DB(g);
    DistMatrix<T> ETL(g), ETR(g),
                  EBL(g), EBR(g);
    DistMatrix<T> FTL(g), FBR(g);

    const Int half = E.Height()/2;
    ScaleTrapezoid( beta, uplo, E );
    LockedPartitionDown( A, AT, AB, half );
    LockedPartitionDown( B, BT, BB, half );
    LockedPartitionDown( C, CT, CB, half );
    LockedPartitionDown( D, DT, DB, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    if( uplo == LOWER )
    {
        LocalGemm( NORMAL, orientationOfB, alpha, AB, BT, T(1), EBL );
        LocalGemm( NORMAL, orientationOfD, alpha, CB, DT, T(1), EBL );
    }
    else
    {
        LocalGemm( NORMAL, orientationOfB, alpha, AT, BB, T(1), ETR );
        LocalGemm( NORMAL, orientationOfD, alpha, CT, DB, T(1), ETR );
    }

    FTL.AlignWith( ETL );
    LocalGemm( NORMAL, orientationOfB, alpha, AT, BT, FTL );
    LocalGemm( NORMAL, orientationOfD, alpha, CT, DT, T(1), FTL );
    AxpyTriangle( uplo, T(1), FTL, ETL );

    FBR.AlignWith( EBR );
    LocalGemm( NORMAL, orientationOfB, alpha, AB, BB, FBR );
    LocalGemm( NORMAL, orientationOfD, alpha, CB, DB, T(1), FBR );
    AxpyTriangle( uplo, T(1), FBR, EBR );
}

// E := alpha (A B^{T/H} + C^{T/H} D) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo, Orientation orientationOfB, Orientation orientationOfC,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2kKernel");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    DistMatrix<T,MC,STAR> AT(g), AB(g);
    DistMatrix<T,MR,STAR> BT(g), BB(g);
    DistMatrix<T,STAR,MC> CL(g), CR(g);
    DistMatrix<T,STAR,MR> DL(g), DR(g);
    DistMatrix<T> ETL(g), ETR(g),
                  EBL(g), EBR(g);
    DistMatrix<T> FTL(g), FBR(g);

    const Int half = E.Height()/2;
    ScaleTrapezoid( beta, uplo, E );
    LockedPartitionDown( A, AT, AB, half );
    LockedPartitionDown( B, BT, BB, half );
    LockedPartitionRight( C, CL, CR, half );
    LockedPartitionRight( D, DL, DR, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    if( uplo == LOWER )
    {
        LocalGemm( NORMAL, orientationOfB, alpha, AB, BT, T(1), EBL );
        LocalGemm( orientationOfC, NORMAL, alpha, CR, DL, T(1), EBL );
    }
    else
    {
        LocalGemm( NORMAL, orientationOfB, alpha, AT, BB, T(1), ETR );
        LocalGemm( orientationOfC, NORMAL, alpha, CL, DR, T(1), ETR );
    }

    FTL.AlignWith( ETL );
    LocalGemm( NORMAL, orientationOfB, alpha, AT, BT, FTL );
    LocalGemm( orientationOfC, NORMAL, alpha, CL, DL, T(1), FTL );
    AxpyTriangle( uplo, T(1), FTL, ETL );

    FBR.AlignWith( EBR );
    LocalGemm( NORMAL, orientationOfB, alpha, AB, BB, FBR );
    LocalGemm( orientationOfC, NORMAL, alpha, CR, DR, T(1), FBR );
    AxpyTriangle( uplo, T(1), FBR, EBR );
}

// E := alpha (A B^{T/H} + C^{T/H} D^{T/H}) + beta C
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2kKernel");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    DistMatrix<T,MC,STAR> AT(g), AB(g);
    DistMatrix<T,MR,STAR> BT(g), BB(g);
    DistMatrix<T,STAR,MC> CL(g), CR(g);
    DistMatrix<T,MR,STAR> DT(g), DB(g);
    DistMatrix<T> ETL(g), ETR(g),
                  EBL(g), EBR(g);
    DistMatrix<T> FTL(g), FBR(g);

    const Int half = E.Height()/2;
    ScaleTrapezoid( beta, uplo, E );
    LockedPartitionDown( A, AT, AB, half );
    LockedPartitionDown( B, BT, BB, half );
    LockedPartitionRight( C, CL, CR, half );
    LockedPartitionDown( D, DT, DB, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    if( uplo == LOWER )
    {
        LocalGemm( NORMAL, orientationOfB, alpha, AB, BT, T(1), EBL );
        LocalGemm( orientationOfC, orientationOfD, alpha, CR, DT, T(1), EBL );
    }
    else
    {
        LocalGemm( NORMAL, orientationOfB, alpha, AT, BB, T(1), ETR );
        LocalGemm( orientationOfC, orientationOfD, alpha, CL, DB, T(1), ETR );
    }

    FTL.AlignWith( ETL );
    LocalGemm( NORMAL, orientationOfB, alpha, AT, BT, FTL );
    LocalGemm( orientationOfC, orientationOfD, alpha, CL, DT, T(1), FTL );
    AxpyTriangle( uplo, T(1), FTL, ETL );

    FBR.AlignWith( EBR );
    LocalGemm( NORMAL, orientationOfB, alpha, AB, BB, FBR );
    LocalGemm( orientationOfC, orientationOfD, alpha, CR, DB, T(1), FBR );
    AxpyTriangle( uplo, T(1), FBR, EBR );
}

// E := alpha (A^{T/H} B + C D) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo, Orientation orientationOfA,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2kKernel");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g);
    DistMatrix<T,MC,STAR> CT(g), CB(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g),
                          DL(g), DR(g);
    DistMatrix<T> ETL(g), ETR(g),
                  EBL(g), EBR(g);
    DistMatrix<T> FTL(g), FBR(g);

    const Int half = E.Height()/2;
    ScaleTrapezoid( beta, uplo, E );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionRight( B, BL, BR, half );
    LockedPartitionDown( C, CT, CB, half );
    LockedPartitionRight( D, DL, DR, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    if( uplo == LOWER )
    {
        LocalGemm( orientationOfA, NORMAL, alpha, AR, BL, T(1), EBL );
        LocalGemm( NORMAL, NORMAL, alpha, CB, DL, T(1), EBL );
    }
    else
    {
        LocalGemm( orientationOfA, NORMAL, alpha, AL, BR, T(1), ETR );
        LocalGemm( NORMAL, NORMAL, alpha, CT, DR, T(1), ETR );
    }

    FTL.AlignWith( ETL );
    LocalGemm( orientationOfA, NORMAL, alpha, AL, BL, FTL );
    LocalGemm( NORMAL, NORMAL, alpha, CT, DL, T(1), FTL );
    AxpyTriangle( uplo, T(1), FTL, ETL );

    FBR.AlignWith( EBR );
    LocalGemm( orientationOfA, NORMAL, alpha, AR, BR, FBR );
    LocalGemm( NORMAL, NORMAL, alpha, CB, DR, T(1), FBR );
    AxpyTriangle( uplo, T(1), FBR, EBR );
}

// E := alpha (A^{T/H} B + C D^{T/H}) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo, Orientation orientationOfA, Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2kKernel");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g);
    DistMatrix<T,MC,STAR> CT(g), CB(g);
    DistMatrix<T,MR,STAR> DT(g), DB(g);
    DistMatrix<T> ETL(g), ETR(g),
                  EBL(g), EBR(g);
    DistMatrix<T> FTL(g), FBR(g);

    const Int half = E.Height()/2;
    ScaleTrapezoid( beta, uplo, E );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionRight( B, BL, BR, half );
    LockedPartitionDown( C, CT, CB, half );
    LockedPartitionDown( D, DT, DB, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    if( uplo == LOWER )
    {
        LocalGemm( orientationOfA, NORMAL, alpha, AR, BL, T(1), EBL );
        LocalGemm( NORMAL, orientationOfD, alpha, CB, DT, T(1), EBL );
    }
    else
    {
        LocalGemm( orientationOfA, NORMAL, alpha, AL, BR, T(1), ETR );
        LocalGemm( NORMAL, orientationOfD, alpha, CT, DB, T(1), ETR );
    }

    FTL.AlignWith( ETL );
    LocalGemm( orientationOfA, NORMAL, alpha, AL, BL, FTL );
    LocalGemm( NORMAL, orientationOfD, alpha, CT, DT, T(1), FTL );
    AxpyTriangle( uplo, T(1), FTL, ETL );

    FBR.AlignWith( EBR );
    LocalGemm( orientationOfA, NORMAL, alpha, AR, BR, FBR );
    LocalGemm( NORMAL, orientationOfD, alpha, CB, DB, T(1), FBR );
    AxpyTriangle( uplo, T(1), FBR, EBR );
}

// E := alpha (A^{T/H} B + C^{T/H} D) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo, Orientation orientationOfA, Orientation orientationOfC,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2kKernel");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g),
                          CL(g), CR(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g),
                          DL(g), DR(g);
    DistMatrix<T> ETL(g), ETR(g),
                  EBL(g), EBR(g);
    DistMatrix<T> FTL(g), FBR(g);

    const Int half = E.Height()/2;
    ScaleTrapezoid( beta, uplo, E );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionRight( B, BL, BR, half );
    LockedPartitionRight( C, CL, CR, half );
    LockedPartitionRight( D, DL, DR, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    if( uplo == LOWER )
    {
        LocalGemm( orientationOfA, NORMAL, alpha, AR, BL, T(1), EBL );
        LocalGemm( orientationOfC, NORMAL, alpha, CR, DL, T(1), EBL );
    }
    else
    {
        LocalGemm( orientationOfA, NORMAL, alpha, AL, BR, T(1), ETR );
        LocalGemm( orientationOfC, NORMAL, alpha, CL, DR, T(1), ETR );
    }

    FTL.AlignWith( ETL );
    LocalGemm( orientationOfA, NORMAL, alpha, AL, BL, FTL );
    LocalGemm( orientationOfC, NORMAL, alpha, CL, DL, T(1), FTL );
    AxpyTriangle( uplo, T(1), FTL, ETL );

    FBR.AlignWith( EBR );
    LocalGemm( orientationOfA, NORMAL, alpha, AR, BR, FBR );
    LocalGemm( orientationOfC, NORMAL, alpha, CR, DR, T(1), FBR );
    AxpyTriangle( uplo, T(1), FBR, EBR );
}

// E := alpha (A^{T/H} B + C^{T/H} D^{T/H}) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2kKernel");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g),
                          CL(g), CR(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g);
    DistMatrix<T,MR,STAR> DT(g), DB(g);
    DistMatrix<T> ETL(g), ETR(g),
                  EBL(g), EBR(g);
    DistMatrix<T> FTL(g), FBR(g);

    const Int half = E.Height()/2;
    ScaleTrapezoid( beta, uplo, E );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionRight( B, BL, BR, half );
    LockedPartitionRight( C, CL, CR, half );
    LockedPartitionDown( D, DT, DB, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    if( uplo == LOWER )
    {
        LocalGemm( orientationOfA, NORMAL, alpha, AR, BL, T(1), EBL );
        LocalGemm( orientationOfC, orientationOfD, alpha, CR, DT, T(1), EBL );
    }
    else
    {
        LocalGemm( orientationOfA, NORMAL, alpha, AL, BR, T(1), ETR );
        LocalGemm( orientationOfC, orientationOfD, alpha, CL, DB, T(1), ETR );
    }

    FTL.AlignWith( ETL );
    LocalGemm( orientationOfA, NORMAL, alpha, AL, BL, FTL );
    LocalGemm( orientationOfC, orientationOfD, alpha, CL, DT, T(1), FTL );
    AxpyTriangle( uplo, T(1), FTL, ETL );

    FBR.AlignWith( EBR );
    LocalGemm( orientationOfA, NORMAL, alpha, AR, BR, FBR );
    LocalGemm( orientationOfC, orientationOfD, alpha, CR, DB, T(1), FBR );
    AxpyTriangle( uplo, T(1), FBR, EBR );
}

// E := alpha (A^{T/H} B^{T/H} + C D) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo, Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2kKernel");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g);
    DistMatrix<T,MR,STAR> BT(g), BB(g);
    DistMatrix<T,MC,STAR> CT(g), CB(g);
    DistMatrix<T,STAR,MR> DL(g), DR(g);
    DistMatrix<T> ETL(g), ETR(g),
                  EBL(g), EBR(g);
    DistMatrix<T> FTL(g), FBR(g);

    const Int half = E.Height()/2;
    ScaleTrapezoid( beta, uplo, E );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionDown( B, BT, BB, half );
    LockedPartitionDown( C, CT, CB, half );
    LockedPartitionRight( D, DL, DR, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    if( uplo == LOWER )
    {
        LocalGemm( orientationOfA, orientationOfB, alpha, AR, BT, T(1), EBL );
        LocalGemm( NORMAL, NORMAL, alpha, CB, DL, T(1), EBL );
    }
    else
    {
        LocalGemm( orientationOfA, orientationOfB, alpha, AL, BB, T(1), ETR );
        LocalGemm( NORMAL, NORMAL, alpha, CT, DR, T(1), ETR );
    }

    FTL.AlignWith( ETL );
    LocalGemm( orientationOfA, orientationOfB, alpha, AL, BT, FTL );
    LocalGemm( NORMAL, NORMAL, alpha, CT, DL, T(1), FTL );
    AxpyTriangle( uplo, T(1), FTL, ETL );

    FBR.AlignWith( EBR );
    LocalGemm( orientationOfA, orientationOfB, alpha, AR, BB, FBR );
    LocalGemm( NORMAL, NORMAL, alpha, CB, DR, T(1), FBR );
    AxpyTriangle( uplo, T(1), FBR, EBR );
}

// E := alpha (A^{T/H} B^{T/H} + C D^{T/H}) + beta C
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2kKernel");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g);
    DistMatrix<T,MR,STAR> BT(g), BB(g);
    DistMatrix<T,MC,STAR> CT(g), CB(g);
    DistMatrix<T,MR,STAR> DT(g), DB(g);
    DistMatrix<T> ETL(g), ETR(g),
                  EBL(g), EBR(g);
    DistMatrix<T> FTL(g), FBR(g);

    const Int half = E.Height()/2;
    ScaleTrapezoid( beta, uplo, E );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionDown( B, BT, BB, half );
    LockedPartitionDown( C, CT, CB, half );
    LockedPartitionDown( D, DT, DB, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    if( uplo == LOWER )
    {
        LocalGemm( orientationOfA, orientationOfB, alpha, AR, BT, T(1), EBL );
        LocalGemm( NORMAL, orientationOfD, alpha, CB, DT, T(1), EBL );
    }
    else
    {
        LocalGemm( orientationOfA, orientationOfB, alpha, AL, BB, T(1), ETR );
        LocalGemm( NORMAL, orientationOfD, alpha, CT, DB, T(1), ETR );
    }

    FTL.AlignWith( ETL );
    LocalGemm( orientationOfA, orientationOfB, alpha, AL, BT, FTL );
    LocalGemm( NORMAL, orientationOfD, alpha, CT, DT, T(1), FTL );
    AxpyTriangle( uplo, T(1), FTL, ETL );

    FBR.AlignWith( EBR );
    LocalGemm( orientationOfA, orientationOfB, alpha, AR, BB, FBR );
    LocalGemm( NORMAL, orientationOfD, alpha, CB, DB, T(1), FBR );
    AxpyTriangle( uplo, T(1), FBR, EBR );
}

// E := alpha (A^{T/H} B^{T/H} + C^{T/H} D) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2kKernel");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g),
                          CL(g), CR(g);
    DistMatrix<T,MR,STAR> BT(g), BB(g);
    DistMatrix<T,STAR,MR> DL(g), DR(g);
    DistMatrix<T> ETL(g), ETR(g),
                  EBL(g), EBR(g);
    DistMatrix<T> FTL(g), FBR(g);

    const Int half = E.Height()/2;
    ScaleTrapezoid( beta, uplo, E );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionDown( B, BT, BB, half );
    LockedPartitionRight( C, CL, CR, half );
    LockedPartitionRight( D, DL, DR, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    if( uplo == LOWER )
    {
        LocalGemm( orientationOfA, orientationOfB, alpha, AR, BT, T(1), EBL );
        LocalGemm( orientationOfC, NORMAL, alpha, CR, DL, T(1), EBL );
    }
    else
    {
        LocalGemm( orientationOfA, orientationOfB, alpha, AL, BB, T(1), ETR );
        LocalGemm( orientationOfC, NORMAL, alpha, CL, DR, T(1), ETR );
    }

    FTL.AlignWith( ETL );
    LocalGemm( orientationOfA, orientationOfB, alpha, AL, BT, FTL );
    LocalGemm( orientationOfC, NORMAL, alpha, CL, DL, T(1), FTL );
    AxpyTriangle( uplo, T(1), FTL, ETL );

    FBR.AlignWith( EBR );
    LocalGemm( orientationOfA, orientationOfB, alpha, AR, BB, FBR );
    LocalGemm( orientationOfC, NORMAL, alpha, CR, DR, T(1), FBR );
    AxpyTriangle( uplo, T(1), FBR, EBR );
}

// E := alpha (A^{T/H} B^{T/H} + C^{T/H} D^{T/H}) + beta C
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  Orientation orientationOfC, Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2kKernel");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g),
                          CL(g), CR(g);
    DistMatrix<T,MR,STAR> BT(g),  DT(g),
                          BB(g),  DB(g);
    DistMatrix<T> ETL(g), ETR(g),
                  EBL(g), EBR(g);
    DistMatrix<T> FTL(g), FBR(g);

    const Int half = E.Height()/2;
    ScaleTrapezoid( beta, uplo, E );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionDown( B, BT, BB, half );
    LockedPartitionRight( C, CL, CR, half );
    LockedPartitionDown( D, DT, DB, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    if( uplo == LOWER )
    {
        LocalGemm( orientationOfA, orientationOfB, alpha, AR, BT, T(1), EBL );
        LocalGemm( orientationOfC, orientationOfD, alpha, CR, DT, T(1), EBL );
    }
    else
    {
        LocalGemm( orientationOfA, orientationOfB, alpha, AL, BB, T(1), ETR );
        LocalGemm( orientationOfC, orientationOfD, alpha, CL, DB, T(1), ETR );
    }

    FTL.AlignWith( ETL );
    LocalGemm( orientationOfA, orientationOfB, alpha, AL, BT, FTL );
    LocalGemm( orientationOfC, orientationOfD, alpha, CL, DT, T(1), FTL );
    AxpyTriangle( uplo, T(1), FTL, ETL );

    FBR.AlignWith( EBR );
    LocalGemm( orientationOfA, orientationOfB, alpha, AR, BB, FBR );
    LocalGemm( orientationOfC, orientationOfD, alpha, CR, DB, T(1), FBR );
    AxpyTriangle( uplo, T(1), FBR, EBR );
}

} // namespace trr2k

// E := alpha (A B + C D) + beta E
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E )
{
    using namespace trr2k;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2k");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel( uplo, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,MC,STAR> AT(g),  CT(g),
                              AB(g),  CB(g);
        DistMatrix<T,STAR,MR> BL(g), BR(g),
                              DL(g), DR(g);
        DistMatrix<T> ETL(g), ETR(g),
                      EBL(g), EBR(g);

        const Int half = E.Height() / 2;
        LockedPartitionDown( A, AT, AB, half );
        LockedPartitionRight( B, BL, BR, half );
        LockedPartitionDown( C, CT, CB, half );
        LockedPartitionRight( D, DL, DR, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            LocalGemm( NORMAL, NORMAL, alpha, AB, BL, beta, EBL );
            LocalGemm( NORMAL, NORMAL, alpha, CB, DL, T(1), EBL );
        }
        else
        {
            LocalGemm( NORMAL, NORMAL, alpha, AT, BR, beta, ETR );
            LocalGemm( NORMAL, NORMAL, alpha, CT, DR, T(1), ETR );
        }

        // Recurse
        LocalTrr2k( uplo, alpha, AT, BL, CT, DL, beta, ETL );
        LocalTrr2k( uplo, alpha, AB, BR, CB, DR, beta, EBR );
    }
}

// E := alpha (A B + C D^{T/H}) + beta E
template<typename T>
void LocalTrr2k
( UpperOrLower uplo, Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E  )
{
    using namespace trr2k;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2k");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel( uplo, orientationOfD, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,MC,STAR> AT(g),  CT(g),
                              AB(g),  CB(g);
        DistMatrix<T,STAR,MR> BL(g), BR(g);
        DistMatrix<T,MR,STAR> DT(g), DB(g);
        DistMatrix<T> ETL(g), ETR(g),
                      EBL(g), EBR(g);

        const Int half = E.Height() / 2;
        LockedPartitionDown( A, AT, AB, half );
        LockedPartitionRight( B, BL, BR, half );
        LockedPartitionDown( C, CT, CB, half );
        LockedPartitionDown( D, DT, DB, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            LocalGemm( NORMAL, NORMAL, alpha, AB, BL, T(1), EBL );
            LocalGemm( NORMAL, orientationOfD, alpha, CB, DT, beta, EBL );
        }
        else
        {
            LocalGemm( NORMAL, NORMAL, alpha, AT, BR, T(1), ETR );
            LocalGemm( NORMAL, orientationOfD, alpha, CT, DB, beta, ETR );
        }

        // Recurse
        LocalTrr2k( uplo, orientationOfD, alpha, AT, BL, CT, DT, beta, ETL );
        LocalTrr2k( uplo, orientationOfD, alpha, AB, BR, CB, DB, beta, EBR );
    }
}

// E := alpha (A B + C^{T/H} D) + beta E
template<typename T>
void LocalTrr2k
( UpperOrLower uplo, Orientation orientationOfC,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E  )
{
    using namespace trr2k;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2k");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel( uplo, orientationOfC, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,MC,STAR> AT(g), AB(g);
        DistMatrix<T,STAR,MR> BL(g), BR(g),
                              DL(g), DR(g);
        DistMatrix<T,STAR,MC> CL(g), CR(g);
        DistMatrix<T> ETL(g), ETR(g),
                      EBL(g), EBR(g);

        const Int half = E.Height() / 2;
        LockedPartitionDown( A, AT, AB, half );
        LockedPartitionRight( B, BL, BR, half );
        LockedPartitionRight( C, CL, CR, half );
        LockedPartitionRight( D, DL, DR, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            LocalGemm( NORMAL, NORMAL, alpha, AB, BL, beta, EBL );
            LocalGemm( orientationOfC, NORMAL, alpha, CR, DL, T(1), EBL );
        }
        else
        {
            LocalGemm( NORMAL, NORMAL, alpha, AT, BR, beta, ETR );
            LocalGemm( orientationOfC, NORMAL, alpha, CL, DR, T(1), ETR );
        }

        // Recurse
        LocalTrr2k( uplo, orientationOfC, alpha, AT, BL, CL, DL, beta, ETL );
        LocalTrr2k( uplo, orientationOfC, alpha, AB, BR, CR, DR, beta, EBR );
    }
}

// E := alpha (A B + C^{T/H} D^{T/H}) + beta E
template<typename T>
void LocalTrr2k
( UpperOrLower uplo, Orientation orientationOfC, Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E  )
{
    using namespace trr2k;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2k");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfC, orientationOfD, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,MC,STAR> AT(g), AB(g); 
        DistMatrix<T,STAR,MR> BL(g), BR(g);
        DistMatrix<T,STAR,MC> CL(g), CR(g);
        DistMatrix<T,MR,STAR> DT(g), DB(g);
        DistMatrix<T> ETL(g), ETR(g),
                      EBL(g), EBR(g);

        const Int half = E.Height() / 2;
        LockedPartitionDown( A, AT, AB, half );
        LockedPartitionRight( B, BL, BR, half );
        LockedPartitionRight( C, CL, CR, half );
        LockedPartitionDown( D, DT, DB, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            LocalGemm( NORMAL, NORMAL, alpha, AB, BL, beta, EBL );
            LocalGemm
            ( orientationOfC, orientationOfD, alpha, CR, DT, T(1), EBL );
        }
        else
        {
            LocalGemm( NORMAL, NORMAL, alpha, AT, BR, beta, ETR );
            LocalGemm
            ( orientationOfC, orientationOfD, alpha, CL, DB, T(1), ETR );
        }

        // Recurse
        LocalTrr2k
        ( uplo, orientationOfC, orientationOfD, 
          alpha, AT, BL, CL, DT, beta, ETL );
        LocalTrr2k
        ( uplo, orientationOfC, orientationOfD,
          alpha, AB, BR, CR, DB, beta, EBR );
    }
}

// E := alpha (A B^{T/H} + C D) + beta E
template<typename T>
void LocalTrr2k
( UpperOrLower uplo, Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E  )
{
    using namespace trr2k;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2k");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel( uplo, orientationOfB, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,MC,STAR> AT(g),  CT(g),
                              AB(g),  CB(g);
        DistMatrix<T,MR,STAR> BT(g), BB(g);
        DistMatrix<T,STAR,MR> DL(g), DR(g);
        DistMatrix<T> ETL(g), ETR(g),
                      EBL(g), EBR(g);

        const Int half = E.Height() / 2;
        LockedPartitionDown( A, AT, AB, half );
        LockedPartitionDown( B, BT, BB, half );
        LockedPartitionDown( C, CT, CB, half );
        LockedPartitionRight( D, DL, DR, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            LocalGemm( NORMAL, orientationOfB, alpha, AB, BT, T(1), EBL );
            LocalGemm( NORMAL, NORMAL, alpha, CB, DL, beta, EBL );
        }
        else
        {
            LocalGemm( NORMAL, orientationOfB, alpha, AT, BB, T(1), ETR );
            LocalGemm( NORMAL, NORMAL, alpha, CT, DR, beta, ETR );
        }

        // Recurse
        LocalTrr2k( uplo, orientationOfB, alpha, AT, BT, CT, DL, beta, ETL );
        LocalTrr2k( uplo, orientationOfB, alpha, AB, BB, CB, DR, beta, EBR );
    }
}

// E := alpha (A B^{T/H} + C D^{T/H}) + beta E
template<typename T>
void LocalTrr2k
( UpperOrLower uplo, Orientation orientationOfB, Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E )
{
    using namespace trr2k;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2k");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfB, orientationOfD, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,MC,STAR> AT(g),  CT(g),
                              AB(g),  CB(g);
        DistMatrix<T,MR,STAR> BT(g),  DT(g),
                              BB(g),  DB(g);
        DistMatrix<T> ETL(g), ETR(g),
                      EBL(g), EBR(g);

        const Int half = E.Height() / 2;
        LockedPartitionDown( A, AT, AB, half );
        LockedPartitionDown( B, BT, BB, half );
        LockedPartitionDown( C, CT, CB, half );
        LockedPartitionDown( D, DT, DB, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            LocalGemm( NORMAL, orientationOfB, alpha, AB, BT, beta, EBL );
            LocalGemm( NORMAL, orientationOfD, alpha, CB, DT, T(1), EBL );
        }
        else
        {
            LocalGemm( NORMAL, orientationOfB, alpha, AT, BB, beta, ETR );
            LocalGemm( NORMAL, orientationOfD, alpha, CT, DB, T(1), ETR );
        }

        // Recurse
        LocalTrr2k
        ( uplo, orientationOfB, orientationOfD,
          alpha, AT, BT, CT, DT, beta, ETL );
        LocalTrr2k
        ( uplo, orientationOfB, orientationOfD,
          alpha, AB, BB, CB, DB, beta, EBR );
    }
}

// E := alpha (A B^{T/H} + C^{T/H} D) + beta E
template<typename T>
void LocalTrr2k
( UpperOrLower uplo, Orientation orientationOfB, Orientation orientationOfC,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E  )
{
    using namespace trr2k;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2k");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfB, orientationOfC, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,MC,STAR> AT(g), AB(g); 
        DistMatrix<T,MR,STAR> BT(g), BB(g);
        DistMatrix<T,STAR,MC> CL(g), CR(g);
        DistMatrix<T,STAR,MR> DL(g), DR(g);
        DistMatrix<T> ETL(g), ETR(g),
                      EBL(g), EBR(g);

        const Int half = E.Height() / 2;
        LockedPartitionDown( A, AT, AB, half );
        LockedPartitionDown( B, BT, BB, half );
        LockedPartitionRight( C, CL, CR, half );
        LockedPartitionRight( D, DL, DR, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            LocalGemm( NORMAL, orientationOfB, alpha, AB, BT, beta, EBL );
            LocalGemm( orientationOfC, NORMAL, alpha, CR, DL, T(1), EBL );
        }
        else
        {
            LocalGemm( NORMAL, orientationOfB, alpha, AT, BB, beta, ETR );
            LocalGemm( orientationOfC, NORMAL, alpha, CL, DR, T(1), ETR );
        }

        // Recurse
        LocalTrr2k
        ( uplo, orientationOfB, orientationOfC,
          alpha, AT, BT, CL, DL, beta, ETL );
        LocalTrr2k
        ( uplo, orientationOfB, orientationOfC,
          alpha, AB, BB, CR, DR, beta, EBR );
    }
}

// E := alpha (A B^{T/H} + C^{T/H} D^{T/H}) + beta E
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,MR,STAR>& B, 
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E  )
{
    using namespace trr2k;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2k");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfB, orientationOfC, orientationOfD, 
          alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,MC,STAR> AT(g), AB(g); 
        DistMatrix<T,MR,STAR> BT(g),  DT(g),
                              BB(g),  DB(g);
        DistMatrix<T,STAR,MC> CL(g), CR(g);
        DistMatrix<T> ETL(g), ETR(g),
                      EBL(g), EBR(g);

        const Int half = E.Height() / 2;
        LockedPartitionDown( A, AT, AB, half );
        LockedPartitionDown( B, BT, BB, half );
        LockedPartitionRight( C, CL, CR, half );
        LockedPartitionDown( D, DT, DB, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            LocalGemm( NORMAL, orientationOfB, alpha, AB, BT, beta, EBL );
            LocalGemm
            ( orientationOfC, orientationOfD, alpha, CR, DT, T(1), EBL );
        }
        else
        {
            LocalGemm( NORMAL, orientationOfB, alpha, AT, BB, beta, ETR );
            LocalGemm
            ( orientationOfC, orientationOfD, alpha, CL, DB, T(1), ETR );
        }

        // Recurse
        LocalTrr2k
        ( uplo, orientationOfB, orientationOfC, orientationOfD,
          alpha, AT, BT, CL, DT, beta, ETL );
        LocalTrr2k
        ( uplo, orientationOfB, orientationOfC, orientationOfD,
          alpha, AB, BB, CR, DB, beta, EBR );
    }
}

// E := alpha (A^{T/H} B + C D) + beta E
template<typename T>
void LocalTrr2k
( UpperOrLower uplo, Orientation orientationOfA,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E  )
{
    using namespace trr2k;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2k");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel( uplo, orientationOfA, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g);
        DistMatrix<T,STAR,MR> BL(g), BR(g),
                              DL(g), DR(g);
        DistMatrix<T,MC,STAR> CT(g), CB(g);
        DistMatrix<T> ETL(g), ETR(g),
                      EBL(g), EBR(g);

        const Int half = E.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionRight( B, BL, BR, half );
        LockedPartitionDown( C, CT, CB, half );
        LockedPartitionRight( D, DL, DR, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            LocalGemm( orientationOfA, NORMAL, alpha, AR, BL, beta, EBL );
            LocalGemm( NORMAL, NORMAL, alpha, CB, DL, T(1), EBL );
        }
        else
        {
            LocalGemm( orientationOfA, NORMAL, alpha, AL, BR, beta, ETR );
            LocalGemm( NORMAL, NORMAL, alpha, CT, DR, T(1), ETR );
        }

        // Recurse
        LocalTrr2k( uplo, orientationOfA, alpha, AL, BL, CT, DL, beta, ETL );
        LocalTrr2k( uplo, orientationOfA, alpha, AR, BR, CB, DR, beta, EBR );
    }
}

// E := alpha (A^{T/H} B + C D^{T/H}) + beta E
template<typename T>
void LocalTrr2k
( UpperOrLower uplo, Orientation orientationOfA, Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E  )
{
    using namespace trr2k;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2k");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfA, orientationOfD, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g);
        DistMatrix<T,STAR,MR> BL(g), BR(g);
        DistMatrix<T,MC,STAR> CT(g), CB(g); 
        DistMatrix<T,MR,STAR> DT(g), DB(g);
        DistMatrix<T> ETL(g), ETR(g),
                      EBL(g), EBR(g);

        const Int half = E.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionRight( B, BL, BR, half );
        LockedPartitionDown( C, CT, CB, half );
        LockedPartitionDown( D, DT, DB, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            LocalGemm( orientationOfA, NORMAL, alpha, AR, BL, beta, EBL );
            LocalGemm( NORMAL, orientationOfD, alpha, CB, DT, T(1), EBL );
        }
        else
        {
            LocalGemm( orientationOfA, NORMAL, alpha, AL, BR, beta, ETR );
            LocalGemm( NORMAL, orientationOfD, alpha, CT, DB, T(1), ETR );
        }

        // Recurse
        LocalTrr2k
        ( uplo, orientationOfA, orientationOfD,
          alpha, AL, BL, CT, DT, beta, ETL );
        LocalTrr2k
        ( uplo, orientationOfA, orientationOfD,
          alpha, AR, BR, CB, DB, beta, EBR );
    }
}

// E := alpha (A^{T/H} B + C^{T/H} D) + beta E
template<typename T>
void LocalTrr2k
( UpperOrLower uplo, Orientation orientationOfA, Orientation orientationOfC,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E  )
{
    using namespace trr2k;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2k");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfA, orientationOfC, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g),
                              CL(g), CR(g);
        DistMatrix<T,STAR,MR> BL(g), BR(g),
                              DL(g), DR(g);
        DistMatrix<T> ETL(g), ETR(g),
                      EBL(g), EBR(g);

        const Int half = E.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionRight( B, BL, BR, half );
        LockedPartitionRight( C, CL, CR, half );
        LockedPartitionRight( D, DL, DR, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            LocalGemm( orientationOfA, NORMAL, alpha, AR, BL, beta, EBL );
            LocalGemm( orientationOfC, NORMAL, alpha, CR, DL, T(1), EBL );
        }
        else
        {
            LocalGemm( orientationOfA, NORMAL, alpha, AL, BR, beta, ETR );
            LocalGemm( orientationOfC, NORMAL, alpha, CL, DR, T(1), ETR );
        }

        // Recurse
        LocalTrr2k
        ( uplo, orientationOfA, orientationOfC,
          alpha, AL, BL, CL, DL, beta, ETL );
        LocalTrr2k
        ( uplo, orientationOfA, orientationOfC, 
          alpha, AR, BR, CR, DR, beta, EBR );
    }
}

// E := alpha (A^{T/H} B + C^{T/H} D^{T/H}) + beta E
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E  )
{
    using namespace trr2k;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2k");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfA, orientationOfC, orientationOfD,
          alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g),
                              CL(g), CR(g);
        DistMatrix<T,STAR,MR> BL(g), BR(g);
        DistMatrix<T,MR,STAR> DT(g), DB(g);
        DistMatrix<T> ETL(g), ETR(g),
                      EBL(g), EBR(g);

        const Int half = E.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionRight( B, BL, BR, half );
        LockedPartitionRight( C, CL, CR, half );
        LockedPartitionDown( D, DT, DB, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            LocalGemm( orientationOfA, NORMAL, alpha, AR, BL, beta, EBL );
            LocalGemm
            ( orientationOfC, orientationOfD, alpha, CR, DT, T(1), EBL );
        }
        else
        {
            LocalGemm( orientationOfA, NORMAL, alpha, AL, BR, beta, ETR );
            LocalGemm
            ( orientationOfC, orientationOfD, alpha, CL, DB, T(1), ETR );
        }

        // Recurse
        LocalTrr2k
        ( uplo, orientationOfA, orientationOfC, orientationOfD,
          alpha, AL, BL, CL, DT, beta, ETL );
        LocalTrr2k
        ( uplo, orientationOfA, orientationOfC, orientationOfD,
          alpha, AR, BR, CR, DB, beta, EBR );
    }
}

// E := alpha (A^{T/H} B^{T/H} + C D) + beta E
template<typename T>
void LocalTrr2k
( UpperOrLower uplo, Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E  )
{
    using namespace trr2k;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2k");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfA, orientationOfB, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g);
        DistMatrix<T,MR,STAR> BT(g), BB(g);
        DistMatrix<T,MC,STAR> CT(g), CB(g); 
        DistMatrix<T,STAR,MR> DL(g), DR(g);
        DistMatrix<T> ETL(g), ETR(g),
                      EBL(g), EBR(g);

        const Int half = E.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionDown( B, BT, BB, half );
        LockedPartitionDown( C, CT, CB, half );
        LockedPartitionRight( D, DL, DR, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            LocalGemm
            ( orientationOfA, orientationOfB, alpha, AR, BT, beta, EBL );
            LocalGemm( NORMAL, NORMAL, alpha, CB, DL, T(1), EBL );
        }
        else
        {
            LocalGemm
            ( orientationOfA, orientationOfB, alpha, AL, BB, beta, ETR );
            LocalGemm( NORMAL, NORMAL, alpha, CT, DR, T(1), ETR );
        }

        // Recurse
        LocalTrr2k
        ( uplo, orientationOfA, orientationOfB,
          alpha, AL, BT, CT, DL, beta, ETL );
        LocalTrr2k
        ( uplo, orientationOfA, orientationOfB,
          alpha, AR, BB, CB, DR, beta, EBR );
    }
}

// E := alpha (A^{T/H} B^{T/H} + C D^{T/H}) + beta E
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E  )
{
    using namespace trr2k;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2k");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfA, orientationOfB, orientationOfD, 
          alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g);
        DistMatrix<T,MR,STAR> BT(g),  DT(g),
                              BB(g),  DB(g);
        DistMatrix<T,MC,STAR> CT(g), CB(g); 
        DistMatrix<T> ETL(g), ETR(g),
                      EBL(g), EBR(g);

        const Int half = E.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionDown( B, BT, BB, half );
        LockedPartitionDown( C, CT, CB, half );
        LockedPartitionDown( D, DT, DB, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            LocalGemm
            ( orientationOfA, orientationOfB, alpha, AR, BT, beta, EBL );
            LocalGemm( NORMAL, orientationOfD, alpha, CB, DT, T(1), EBL );
        }
        else
        {
            LocalGemm
            ( orientationOfA, orientationOfB, alpha, AL, BB, beta, ETR );
            LocalGemm( NORMAL, orientationOfD, alpha, CT, DB, T(1), ETR );
        }

        // Recurse
        LocalTrr2k
        ( uplo, orientationOfA, orientationOfB, orientationOfD,
          alpha, AL, BT, CT, DT, beta, ETL );
        LocalTrr2k
        ( uplo, orientationOfA, orientationOfB, orientationOfD,
          alpha, AR, BB, CB, DB, beta, EBR );
    }
}

// E := alpha (A^{T/H} B^{T/H} + C^{T/H} D) + beta E
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E  )
{
    using namespace trr2k;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2k");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfA, orientationOfB, orientationOfC,
          alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g),
                              CL(g), CR(g);
        DistMatrix<T,MR,STAR> BT(g), BB(g);
        DistMatrix<T,STAR,MR> DL(g), DR(g);
        DistMatrix<T> ETL(g), ETR(g),
                      EBL(g), EBR(g);

        const Int half = E.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionDown( B, BT, BB, half );
        LockedPartitionRight( C, CL, CR, half );
        LockedPartitionRight( D, DL, DR, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            LocalGemm
            ( orientationOfA, orientationOfB, alpha, AR, BT, beta, EBL );
            LocalGemm( orientationOfC, NORMAL, alpha, CR, DL, T(1), EBL );
        }
        else
        {
            LocalGemm
            ( orientationOfA, orientationOfB, alpha, AL, BB, beta, ETR );
            LocalGemm( orientationOfC, NORMAL, alpha, CL, DR, T(1), ETR );
        }

        // Recurse
        LocalTrr2k
        ( uplo, orientationOfA, orientationOfB, orientationOfC,
          alpha, AL, BT, CL, DL, beta, ETL );
        LocalTrr2k
        ( uplo, orientationOfA, orientationOfB, orientationOfC,
          alpha, AR, BB, CR, DR, beta, EBR );
    }
}

// E := alpha (A^{T/H} B^{T/H} + C^{T/H} D^{T/H}) + beta E
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E  )
{
    using namespace trr2k;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrr2k");
        CheckInput( A, B, C, D, E );
    )
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfA, orientationOfB, orientationOfC, orientationOfD,
          alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g),
                              CL(g), CR(g);
        DistMatrix<T,MR,STAR> BT(g),  DT(g),
                              BB(g),  DB(g);
        DistMatrix<T> ETL(g), ETR(g),
                      EBL(g), EBR(g);

        const Int half = E.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionDown( B, BT, BB, half );
        LockedPartitionRight( C, CL, CR, half );
        LockedPartitionDown( D, DT, DB, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            LocalGemm
            ( orientationOfA, orientationOfB, alpha, AR, BT, beta, EBL );
            LocalGemm
            ( orientationOfC, orientationOfD, alpha, CR, DT, T(1), EBL );
        }
        else
        {
            LocalGemm
            ( orientationOfA, orientationOfB, alpha, AL, BB, beta, ETR );
            LocalGemm
            ( orientationOfC, orientationOfD, alpha, CL, DB, T(1), ETR );
        }

        // Recurse
        LocalTrr2k
        ( uplo, 
          orientationOfA, orientationOfB, orientationOfC, orientationOfD, 
          alpha, AL, BT, CL, DT, beta, ETL );

        LocalTrr2k
        ( uplo, 
          orientationOfA, orientationOfB, orientationOfC, orientationOfD,
          alpha, AR, BB, CR, DB, beta, EBR );
    }
}

} // namespace El

#endif // ifndef EL_TRR2K_LOCAL_HPP
