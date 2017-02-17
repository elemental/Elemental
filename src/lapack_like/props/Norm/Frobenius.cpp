/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename Field>
Base<Field> FrobeniusNorm( const Matrix<Field>& A )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    Real scale = 0;
    Real scaledSquare = 1;
    const Int width = A.Width();
    const Int height = A.Height();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            UpdateScaledSquare( A(i,j), scale, scaledSquare );
    return scale*Sqrt(scaledSquare);
}

template<typename Field>
Base<Field> FrobeniusNorm( const SparseMatrix<Field>& A )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    Real scale = 0;
    Real scaledSquare = 1;
    const Int numEntries = A.NumEntries();
    const Field* valBuf = A.LockedValueBuffer();
    for( Int k=0; k<numEntries; ++k )
        UpdateScaledSquare( valBuf[k], scale, scaledSquare );
    return scale*Sqrt(scaledSquare);
}

template<typename Field>
Base<Field> HermitianFrobeniusNorm( UpperOrLower uplo, const Matrix<Field>& A )
{
    EL_DEBUG_CSE
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    typedef Base<Field> Real;
    Real scale = 0;
    Real scaledSquare = 1;
    const Int height = A.Height();
    const Int width = A.Width();
    if( uplo == UPPER )
    {
        for( Int j=0; j<width; ++j )
        {
            for( Int i=0; i<j; ++i )
            {
                UpdateScaledSquare( A(i,j), scale, scaledSquare );
                UpdateScaledSquare( A(i,j), scale, scaledSquare );
            }
            UpdateScaledSquare( A(j,j), scale, scaledSquare );
        }
    }
    else
    {
        for( Int j=0; j<width; ++j )
        {
            for( Int i=j+1; i<height; ++i )
            {
                UpdateScaledSquare( A(i,j), scale, scaledSquare );
                UpdateScaledSquare( A(i,j), scale, scaledSquare );
            }
            UpdateScaledSquare( A(j,j), scale, scaledSquare );
        }
    }
    return scale*Sqrt(scaledSquare);
}

template<typename Field>
Base<Field>
HermitianFrobeniusNorm( UpperOrLower uplo, const SparseMatrix<Field>& A )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    Real scale = 0;
    Real scaledSquare = 1;
    const Int numEntries = A.NumEntries();
    const Field* valBuf = A.LockedValueBuffer();
    for( Int k=0; k<numEntries; ++k )
    {
        const Int i = A.Row(k);
        const Int j = A.Col(k);
        if( (uplo==LOWER && i>j) || (uplo==UPPER && i<j) )
        {
            UpdateScaledSquare( valBuf[k], scale, scaledSquare );
            UpdateScaledSquare( valBuf[k], scale, scaledSquare );
        }
        else if( i == j )
        {
            UpdateScaledSquare( valBuf[k], scale, scaledSquare );
        }
    }
    return scale*Sqrt(scaledSquare);
}

template<typename Field>
Base<Field> SymmetricFrobeniusNorm( UpperOrLower uplo, const Matrix<Field>& A )
{
    EL_DEBUG_CSE
    return HermitianFrobeniusNorm( uplo, A );
}

template<typename Field>
Base<Field>
SymmetricFrobeniusNorm( UpperOrLower uplo, const SparseMatrix<Field>& A )
{
    EL_DEBUG_CSE
    return HermitianFrobeniusNorm( uplo, A );
}

template<typename Real>
Real NormFromScaledSquare
( Real localScale, Real localScaledSquare, mpi::Comm comm )
{
    // Find the maximum relative scale
    const Real scale = mpi::AllReduce( localScale, mpi::MAX, comm );

    if( scale != Real(0) )
    {
        // Equilibrate our local scaled sum to the maximum scale
        Real relScale = localScale/scale;
        localScaledSquare *= relScale*relScale;

        // The scaled square is now the sum of the local contributions
        const Real scaledSquare = mpi::AllReduce( localScaledSquare, comm );
        return scale*Sqrt(scaledSquare);
    }
    else
      return 0;
}

template<typename Field>
Base<Field> FrobeniusNorm( const AbstractDistMatrix<Field>& A )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    Real norm;
    if( A.Participating() )
    {
        Real localScale=0, localScaledSquare=1;
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        const Matrix<Field>& ALoc = A.LockedMatrix();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                UpdateScaledSquare
                ( ALoc(iLoc,jLoc), localScale, localScaledSquare );

        norm = NormFromScaledSquare
          ( localScale, localScaledSquare, A.DistComm() );
    }
    mpi::Broadcast( norm, A.Root(), A.CrossComm() );
    return norm;
}

template<typename Field>
Base<Field> FrobeniusNorm( const DistSparseMatrix<Field>& A )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    Real localScale=0, localScaledSquare=1;
    const Int numLocalEntries = A.NumLocalEntries();
    const Field* valBuf = A.LockedValueBuffer();
    for( Int k=0; k<numLocalEntries; ++k )
        UpdateScaledSquare( valBuf[k], localScale, localScaledSquare );

    return NormFromScaledSquare
      ( localScale, localScaledSquare, A.Grid().Comm() );
}

template<typename Field>
Base<Field> HermitianFrobeniusNorm
( UpperOrLower uplo, const AbstractDistMatrix<Field>& A )
{
    EL_DEBUG_CSE
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    typedef Base<Field> Real;
    Real norm;
    if( A.Participating() )
    {
        Real localScale = 0;
        Real localScaledSquare = 1;
        const Int localWidth = A.LocalWidth();
        const Int localHeight = A.LocalHeight();
        const Matrix<Field>& ALoc = A.LockedMatrix();
        if( uplo == UPPER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int numUpperRows = A.LocalRowOffset(j+1);
                for( Int iLoc=0; iLoc<numUpperRows; ++iLoc )
                {
                    const Int i = A.GlobalRow(iLoc);
                    UpdateScaledSquare
                    ( ALoc(iLoc,jLoc), localScale, localScaledSquare );
                    if( i != j )
                        UpdateScaledSquare
                        ( ALoc(iLoc,jLoc), localScale, localScaledSquare );
                }
            }
        }
        else
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int numStrictlyUpperRows = A.LocalRowOffset(j);
                for( Int iLoc=numStrictlyUpperRows; iLoc<localHeight; ++iLoc )
                {
                    const Int i = A.GlobalRow(iLoc);
                    UpdateScaledSquare
                    ( ALoc(iLoc,jLoc), localScale, localScaledSquare );
                    if( i != j )
                        UpdateScaledSquare
                        ( ALoc(iLoc,jLoc), localScale, localScaledSquare );
                }
            }
        }

        norm = NormFromScaledSquare
          ( localScale, localScaledSquare, A.DistComm() );
    }
    mpi::Broadcast( norm, A.Root(), A.CrossComm() );
    return norm;
}

template<typename Field>
Base<Field> HermitianFrobeniusNorm
( UpperOrLower uplo, const DistSparseMatrix<Field>& A )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    Real localScale=0, localScaledSquare=1;
    const Int numLocalEntries = A.NumLocalEntries();
    const Int* rowBuf = A.LockedSourceBuffer();
    const Int* colBuf = A.LockedTargetBuffer();
    const Field* valBuf = A.LockedValueBuffer();
    for( Int k=0; k<numLocalEntries; ++k )
    {
        const Int i = rowBuf[k];
        const Int j = colBuf[k];
        const Field value = valBuf[k];
        if( (uplo==UPPER && i<j) || (uplo==LOWER && i>j) )
        {
            UpdateScaledSquare( value, localScale, localScaledSquare );
            UpdateScaledSquare( value, localScale, localScaledSquare );
        }
        else if( i ==j )
            UpdateScaledSquare( value, localScale, localScaledSquare );
    }

    return NormFromScaledSquare
      ( localScale, localScaledSquare, A.Grid().Comm() );
}

template<typename Field>
Base<Field> SymmetricFrobeniusNorm
( UpperOrLower uplo, const AbstractDistMatrix<Field>& A )
{
    EL_DEBUG_CSE
    return HermitianFrobeniusNorm( uplo, A );
}

template<typename Field>
Base<Field> SymmetricFrobeniusNorm
( UpperOrLower uplo, const DistSparseMatrix<Field>& A )
{
    EL_DEBUG_CSE
    return HermitianFrobeniusNorm( uplo, A );
}

template<typename Field>
Base<Field> FrobeniusNorm( const DistMultiVec<Field>& A )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    Real localScale=0, localScaledSquare=1;
    const Int localHeight = A.LocalHeight();
    const Int width = A.Width();
    const Matrix<Field>& ALoc = A.LockedMatrix();
    for( Int j=0; j<width; ++j )
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            UpdateScaledSquare( ALoc(iLoc,j), localScale, localScaledSquare );

    return NormFromScaledSquare
      ( localScale, localScaledSquare, A.Grid().Comm() );
}

#define PROTO(Field) \
  template Base<Field> FrobeniusNorm( const Matrix<Field>& A ); \
  template Base<Field> FrobeniusNorm ( const AbstractDistMatrix<Field>& A ); \
  template Base<Field> FrobeniusNorm( const SparseMatrix<Field>& A ); \
  template Base<Field> FrobeniusNorm( const DistSparseMatrix<Field>& A ); \
  template Base<Field> FrobeniusNorm ( const DistMultiVec<Field>& A ); \
  template Base<Field> HermitianFrobeniusNorm \
  ( UpperOrLower uplo, const Matrix<Field>& A ); \
  template Base<Field> HermitianFrobeniusNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<Field>& A ); \
  template Base<Field> HermitianFrobeniusNorm \
  ( UpperOrLower uplo, const SparseMatrix<Field>& A ); \
  template Base<Field> HermitianFrobeniusNorm \
  ( UpperOrLower uplo, const DistSparseMatrix<Field>& A ); \
  template Base<Field> SymmetricFrobeniusNorm \
  ( UpperOrLower uplo, const Matrix<Field>& A ); \
  template Base<Field> SymmetricFrobeniusNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<Field>& A ); \
  template Base<Field> SymmetricFrobeniusNorm \
  ( UpperOrLower uplo, const SparseMatrix<Field>& A ); \
  template Base<Field> SymmetricFrobeniusNorm \
  ( UpperOrLower uplo, const DistSparseMatrix<Field>& A );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
