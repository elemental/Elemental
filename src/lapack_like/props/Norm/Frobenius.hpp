/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_NORM_FROBENIUS_HPP
#define EL_NORM_FROBENIUS_HPP

namespace El {

template<typename F> 
Base<F> FrobeniusNorm( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("FrobeniusNorm"))
    typedef Base<F> Real;
    Real scale = 0;
    Real scaledSquare = 1;
    const Int width = A.Width();
    const Int height = A.Height();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            UpdateScaledSquare( A.Get(i,j), scale, scaledSquare );
    return scale*Sqrt(scaledSquare);
}

template<typename F> 
Base<F> FrobeniusNorm( const SparseMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("FrobeniusNorm"))
    typedef Base<F> Real;
    Real scale = 0;
    Real scaledSquare = 1;
    const Int numEntries = A.NumEntries();
    for( Int k=0; k<numEntries; ++k )
        UpdateScaledSquare( A.Value(k), scale, scaledSquare );
    return scale*Sqrt(scaledSquare);
}

template<typename F>
Base<F> HermitianFrobeniusNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianFrobeniusNorm"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    typedef Base<F> Real;
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
                UpdateScaledSquare( A.Get(i,j), scale, scaledSquare );
                UpdateScaledSquare( A.Get(i,j), scale, scaledSquare );
            }
            UpdateScaledSquare( A.Get(j,j), scale, scaledSquare );
        }
    }
    else
    {
        for( Int j=0; j<width; ++j )
        {
            for( Int i=j+1; i<height; ++i )
            {
                UpdateScaledSquare( A.Get(i,j), scale, scaledSquare );
                UpdateScaledSquare( A.Get(i,j), scale, scaledSquare );
            }
            UpdateScaledSquare( A.Get(j,j), scale, scaledSquare );
        }
    }
    return scale*Sqrt(scaledSquare);
}

template<typename F> 
Base<F> HermitianFrobeniusNorm( UpperOrLower uplo, const SparseMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianFrobeniusNorm"))
    typedef Base<F> Real;
    Real scale = 0;
    Real scaledSquare = 1;
    const Int numEntries = A.NumEntries();
    for( Int k=0; k<numEntries; ++k )
    {
        const Int i = A.Row(k);
        const Int j = A.Col(k);
        if( (uplo==LOWER && i>j) || (uplo==UPPER && i<j) )
        { 
            UpdateScaledSquare( A.Value(k), scale, scaledSquare );
            UpdateScaledSquare( A.Value(k), scale, scaledSquare );
        }
        else if( i == j )
        {
            UpdateScaledSquare( A.Value(k), scale, scaledSquare );
        }
    }
    return scale*Sqrt(scaledSquare);
}

template<typename F>
Base<F> SymmetricFrobeniusNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricFrobeniusNorm"))
    return HermitianFrobeniusNorm( uplo, A );
}

template<typename F>
Base<F> SymmetricFrobeniusNorm( UpperOrLower uplo, const SparseMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricFrobeniusNorm"))
    return HermitianFrobeniusNorm( uplo, A );
}

template<typename F> 
Base<F> FrobeniusNorm( const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("FrobeniusNorm"))
    typedef Base<F> Real;
    Real norm;
    if( A.Participating() )
    {
        Real locScale=0, locScaledSquare=1;
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                UpdateScaledSquare
                ( A.GetLocal(iLoc,jLoc), locScale, locScaledSquare );

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
    DEBUG_ONLY(CallStackEntry cse("FrobeniusNorm"))
    typedef Base<F> Real;
    Real norm;

    Real locScale=0, locScaledSquare=1;
    const Int numLocalEntries = A.NumLocalEntries();
    for( Int k=0; k<numLocalEntries; ++k )
        UpdateScaledSquare( A.Value(k), locScale, locScaledSquare );

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
    DEBUG_ONLY(CallStackEntry cse("HermitianFrobeniusNorm"))
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
                    ( A.GetLocal(iLoc,jLoc), locScale, locScaledSquare );
                    if( i != j )
                        UpdateScaledSquare
                        ( A.GetLocal(iLoc,jLoc), locScale, locScaledSquare );
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
                    ( A.GetLocal(iLoc,jLoc), locScale, locScaledSquare );
                    if( i != j )
                        UpdateScaledSquare
                        ( A.GetLocal(iLoc,jLoc), locScale, locScaledSquare );
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
    DEBUG_ONLY(CallStackEntry cse("HermitianFrobeniusNorm"))
    typedef Base<F> Real;
    Real norm;

    Real locScale=0, locScaledSquare=1;
    const Int numLocalEntries = A.NumLocalEntries();
    for( Int k=0; k<numLocalEntries; ++k )
    {
        const Int i = A.Row(k);
        const Int j = A.Col(k);
        const F value = A.Value(k);
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
    DEBUG_ONLY(CallStackEntry cse("SymmetricFrobeniusNorm"))
    return HermitianFrobeniusNorm( uplo, A );
}

template<typename F>
Base<F> SymmetricFrobeniusNorm
( UpperOrLower uplo, const DistSparseMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricFrobeniusNorm"))
    return HermitianFrobeniusNorm( uplo, A );
}

template<typename F> 
Base<F> FrobeniusNorm( const DistMultiVec<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("FrobeniusNorm"))
    typedef Base<F> Real;
    Real norm;
    Real locScale=0, locScaledSquare=1;
    const Int localHeight = A.LocalHeight();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            UpdateScaledSquare( A.GetLocal(iLoc,j), locScale, locScaledSquare );

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

} // namespace El

#endif // ifndef EL_NORM_FROBENIUS_HPP
