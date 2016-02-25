/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F> 
Base<F> FrobeniusNorm( const Matrix<F>& A )
{
    DEBUG_ONLY(CSE cse("FrobeniusNorm"))
    typedef Base<F> Real;
    Real scale = 0;
    Real scaledSquare = 1;
    const Int width = A.Width();
    const Int height = A.Height();
    const F* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            UpdateScaledSquare( ABuf[i+j*ALDim], scale, scaledSquare );
    return scale*Sqrt(scaledSquare);
}

template<typename F> 
Base<F> FrobeniusNorm( const SparseMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("FrobeniusNorm"))
    typedef Base<F> Real;
    Real scale = 0;
    Real scaledSquare = 1;
    const Int numEntries = A.NumEntries();
    const F* valBuf = A.LockedValueBuffer();
    for( Int k=0; k<numEntries; ++k )
        UpdateScaledSquare( valBuf[k], scale, scaledSquare );
    return scale*Sqrt(scaledSquare);
}

template<typename F>
Base<F> HermitianFrobeniusNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CSE cse("HermitianFrobeniusNorm"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    typedef Base<F> Real;
    Real scale = 0;
    Real scaledSquare = 1;
    const Int height = A.Height();
    const Int width = A.Width();
    const F* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    if( uplo == UPPER )
    {
        for( Int j=0; j<width; ++j )
        {
            for( Int i=0; i<j; ++i )
            {
                UpdateScaledSquare( ABuf[i+j*ALDim], scale, scaledSquare );
                UpdateScaledSquare( ABuf[i+j*ALDim], scale, scaledSquare );
            }
            UpdateScaledSquare( ABuf[j+j*ALDim], scale, scaledSquare );
        }
    }
    else
    {
        for( Int j=0; j<width; ++j )
        {
            for( Int i=j+1; i<height; ++i )
            {
                UpdateScaledSquare( ABuf[i+j*ALDim], scale, scaledSquare );
                UpdateScaledSquare( ABuf[i+j*ALDim], scale, scaledSquare );
            }
            UpdateScaledSquare( ABuf[j+j*ALDim], scale, scaledSquare );
        }
    }
    return scale*Sqrt(scaledSquare);
}

template<typename F> 
Base<F> HermitianFrobeniusNorm( UpperOrLower uplo, const SparseMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("HermitianFrobeniusNorm"))
    typedef Base<F> Real;
    Real scale = 0;
    Real scaledSquare = 1;
    const Int numEntries = A.NumEntries();
    const F* valBuf = A.LockedValueBuffer();
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

template<typename F>
Base<F> SymmetricFrobeniusNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CSE cse("SymmetricFrobeniusNorm"))
    return HermitianFrobeniusNorm( uplo, A );
}

template<typename F>
Base<F> SymmetricFrobeniusNorm( UpperOrLower uplo, const SparseMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("SymmetricFrobeniusNorm"))
    return HermitianFrobeniusNorm( uplo, A );
}

template<typename F> 
Base<F> FrobeniusNorm( const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("FrobeniusNorm"))
    typedef Base<F> Real;
    Real norm;
    if( A.Participating() )
    {
        Real locScale=0, locScaledSquare=1;
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        const F* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                UpdateScaledSquare
                ( ABuf[iLoc+jLoc*ALDim], locScale, locScaledSquare );

        // Find the maximum relative scale
        mpi::Comm comm = A.DistComm();
        const Real scale = mpi::AllReduce( locScale, mpi::MAX, comm );

        norm = 0;
        if( scale != Real(0) )
        {
            // Equilibrate our local scaled sum to the maximum scale
            Real relScale = locScale/scale;
            locScaledSquare *= relScale*relScale;

            // The scaled square is now the sum of the local contributions
            const Real scaledSquare = mpi::AllReduce( locScaledSquare, comm );
            norm = scale*Sqrt(scaledSquare);
        }
    }
    mpi::Broadcast( norm, A.Root(), A.CrossComm() );
    return norm;
}

template<typename F> 
Base<F> FrobeniusNorm( const DistSparseMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("FrobeniusNorm"))
    typedef Base<F> Real;
    Real norm;

    Real locScale=0, locScaledSquare=1;
    const Int numLocalEntries = A.NumLocalEntries();
    const F* valBuf = A.LockedValueBuffer();
    for( Int k=0; k<numLocalEntries; ++k )
        UpdateScaledSquare( valBuf[k], locScale, locScaledSquare );

    // Find the maximum relative scale
    mpi::Comm comm = A.Comm();
    const Real scale = mpi::AllReduce( locScale, mpi::MAX, comm );

    norm = 0;
    if( scale != Real(0) )
    {
        // Equilibrate our local scaled sum to the maximum scale
        Real relScale = locScale/scale;
        locScaledSquare *= relScale*relScale;

        // The scaled square is now the sum of the local contributions
        const Real scaledSquare = mpi::AllReduce( locScaledSquare, comm );
        norm = scale*Sqrt(scaledSquare);
    }

    return norm;
}

template<typename F>
Base<F> HermitianFrobeniusNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("HermitianFrobeniusNorm"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    typedef Base<F> Real;
    Real norm;
    if( A.Participating() )
    {
        Real locScale = 0;
        Real locScaledSquare = 1;
        const Int localWidth = A.LocalWidth();
        const Int localHeight = A.LocalHeight();
        const F* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
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
                    ( ABuf[iLoc+jLoc*ALDim], locScale, locScaledSquare );
                    if( i != j )
                        UpdateScaledSquare
                        ( ABuf[iLoc+jLoc*ALDim], locScale, locScaledSquare );
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
                    ( ABuf[iLoc+jLoc*ALDim], locScale, locScaledSquare );
                    if( i != j )
                        UpdateScaledSquare
                        ( ABuf[iLoc+jLoc*ALDim], locScale, locScaledSquare );
                }
            }
        }

        // Find the maximum relative scale
        const Real scale = mpi::AllReduce( locScale, mpi::MAX, A.DistComm() );

        norm = 0;
        if( scale != Real(0) )
        {
            // Equilibrate our local scaled sum to the maximum scale
            Real relScale = locScale/scale;
            locScaledSquare *= relScale*relScale;

            // The scaled square is now the sum of the local contributions
            const Real scaledSquare = 
                mpi::AllReduce( locScaledSquare, A.DistComm() );
            norm = scale*Sqrt(scaledSquare);
        }
    }
    mpi::Broadcast( norm, A.Root(), A.CrossComm() );
    return norm;
}

template<typename F> 
Base<F> HermitianFrobeniusNorm
( UpperOrLower uplo, const DistSparseMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("HermitianFrobeniusNorm"))
    typedef Base<F> Real;
    Real norm;

    Real locScale=0, locScaledSquare=1;
    const Int numLocalEntries = A.NumLocalEntries();
    const Int* rowBuf = A.LockedSourceBuffer();
    const Int* colBuf = A.LockedTargetBuffer();
    const F* valBuf = A.LockedValueBuffer();
    for( Int k=0; k<numLocalEntries; ++k )
    {
        const Int i = rowBuf[k];
        const Int j = colBuf[k];
        const F value = valBuf[k];
        if( (uplo==UPPER && i<j) || (uplo==LOWER && i>j) )
        {
            UpdateScaledSquare( value, locScale, locScaledSquare );
            UpdateScaledSquare( value, locScale, locScaledSquare );
        }
        else if( i ==j )
            UpdateScaledSquare( value, locScale, locScaledSquare );
    }

    // Find the maximum relative scale
    mpi::Comm comm = A.Comm();
    const Real scale = mpi::AllReduce( locScale, mpi::MAX, comm );

    norm = 0;
    if( scale != Real(0) )
    {
        // Equilibrate our local scaled sum to the maximum scale
        Real relScale = locScale/scale;
        locScaledSquare *= relScale*relScale;

        // The scaled square is now the sum of the local contributions
        const Real scaledSquare = mpi::AllReduce( locScaledSquare, comm );
        norm = scale*Sqrt(scaledSquare);
    }

    return norm;
}

template<typename F>
Base<F> SymmetricFrobeniusNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("SymmetricFrobeniusNorm"))
    return HermitianFrobeniusNorm( uplo, A );
}

template<typename F>
Base<F> SymmetricFrobeniusNorm
( UpperOrLower uplo, const DistSparseMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("SymmetricFrobeniusNorm"))
    return HermitianFrobeniusNorm( uplo, A );
}

template<typename F> 
Base<F> FrobeniusNorm( const DistMultiVec<F>& A )
{
    DEBUG_ONLY(CSE cse("FrobeniusNorm"))
    typedef Base<F> Real;
    Real norm;
    Real locScale=0, locScaledSquare=1;
    const Int localHeight = A.LocalHeight();
    const Int width = A.Width();
    const F* ABuf = A.LockedMatrix().LockedBuffer();
    const Int ALDim = A.LockedMatrix().LDim();
    for( Int j=0; j<width; ++j )
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            UpdateScaledSquare( ABuf[iLoc+j*ALDim], locScale, locScaledSquare );

    // Find the maximum relative scale
    mpi::Comm comm = A.Comm();
    const Real scale = mpi::AllReduce( locScale, mpi::MAX, comm );

    norm = 0;
    if( scale != Real(0) )
    {
        // Equilibrate our local scaled sum to the maximum scale
        Real relScale = locScale/scale;
        locScaledSquare *= relScale*relScale;

        // The scaled square is now the sum of the local contributions
        const Real scaledSquare = mpi::AllReduce( locScaledSquare, comm );
        norm = scale*Sqrt(scaledSquare);
    }
    return norm;
}

#define PROTO(F) \
  template Base<F> FrobeniusNorm( const Matrix<F>& A ); \
  template Base<F> FrobeniusNorm ( const AbstractDistMatrix<F>& A ); \
  template Base<F> FrobeniusNorm( const SparseMatrix<F>& A ); \
  template Base<F> FrobeniusNorm( const DistSparseMatrix<F>& A ); \
  template Base<F> FrobeniusNorm ( const DistMultiVec<F>& A ); \
  template Base<F> HermitianFrobeniusNorm \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> HermitianFrobeniusNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A ); \
  template Base<F> HermitianFrobeniusNorm \
  ( UpperOrLower uplo, const SparseMatrix<F>& A ); \
  template Base<F> HermitianFrobeniusNorm \
  ( UpperOrLower uplo, const DistSparseMatrix<F>& A ); \
  template Base<F> SymmetricFrobeniusNorm \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> SymmetricFrobeniusNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A ); \
  template Base<F> SymmetricFrobeniusNorm \
  ( UpperOrLower uplo, const SparseMatrix<F>& A ); \
  template Base<F> SymmetricFrobeniusNorm \
  ( UpperOrLower uplo, const DistSparseMatrix<F>& A );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
