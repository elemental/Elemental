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
Base<F> EntrywiseNorm( const Matrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CSE cse("EntrywiseNorm"))
    // TODO: Make this more numerically stable
    typedef Base<F> Real;
    Real sum = 0;
    const Int width = A.Width();
    const Int height = A.Height();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            sum += Pow( Abs(A.Get(i,j)), p );
    return Pow( sum, 1/p );
}

template<typename F> 
Base<F> EntrywiseNorm( const SparseMatrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CSE cse("EntrywiseNorm"))
    // TODO: Make this more numerically stable
    typedef Base<F> Real;
    Real sum = 0;
    const Int numEntries = A.NumEntries();
    for( Int k=0; k<numEntries; ++k )
        sum += Pow( Abs(A.Value(k)), p );
    return Pow( sum, 1/p );
}

template<typename F>
Base<F> HermitianEntrywiseNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CSE cse("HermitianEntrywiseNorm"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    // TODO: make this more numerically stable
    typedef Base<F> Real;
    Real sum = 0;
    const Int height = A.Height();
    const Int width = A.Width();
    if( uplo == UPPER )
    {
        for( Int j=0; j<width; ++j )
        {
            for( Int i=0; i<j; ++i )
            {
                const Real term = Pow( Abs(A.Get(i,j)), p );
                if( i ==j )
                    sum += term;
                else
                    sum += 2*term;
            }
        }
    }
    else
    {
        for( Int j=0; j<width; ++j )
        {
            for( Int i=j+1; i<height; ++i )
            {
                const Real term = Pow( Abs(A.Get(i,j)), p );
                if( i ==j )
                    sum += term;
                else
                    sum += 2*term;
            }
        }
    }
    return Pow( sum, 1/p );
}

template<typename F> 
Base<F> HermitianEntrywiseNorm
( UpperOrLower uplo, const SparseMatrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CSE cse("HermitianEntrywiseNorm"))
    // TODO: Make this more numerically stable
    typedef Base<F> Real;
    Real sum = 0;
    const Int numEntries = A.NumEntries();
    for( Int k=0; k<numEntries; ++k )
    {
        const Int i = A.Row(k);
        const Int j = A.Col(k);
        if( (uplo==UPPER && i<j) || (uplo==LOWER && i>j) )
            sum += 2*Pow( Abs(A.Value(k)), p );
        else if( i == j )
            sum += Pow( Abs(A.Value(k)), p );
    }
    return Pow( sum, 1/p );
}

template<typename F>
Base<F> SymmetricEntrywiseNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CSE cse("SymmetricEntrywiseNorm"))
    return HermitianEntrywiseNorm( uplo, A, p );
}

template<typename F>
Base<F> SymmetricEntrywiseNorm
( UpperOrLower uplo, const SparseMatrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CSE cse("SymmetricEntrywiseNorm"))
    return HermitianEntrywiseNorm( uplo, A, p );
}

template<typename F> 
Base<F> EntrywiseNorm( const AbstractDistMatrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CSE cse("EntrywiseNorm"))
    typedef Base<F> Real;
    Real norm;
    if( A.Participating() )
    {
        Real localSum = 0;
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                localSum += Pow( Abs(A.GetLocal(iLoc,jLoc)), p ); 
        const Real sum = mpi::AllReduce( localSum, A.DistComm() );
        norm = Pow( sum, 1/p );
    }
    mpi::Broadcast( norm, A.Root(), A.CrossComm() );
    return norm;
}

template<typename F> 
Base<F> EntrywiseNorm( const DistSparseMatrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CSE cse("EntrywiseNorm"))
    typedef Base<F> Real;

    Real localSum = 0;
    const Int numLocalEntries = A.NumLocalEntries();
    for( Int k=0; k<numLocalEntries; ++k )
        localSum += Pow( Abs(A.Value(k)), p ); 

    const Real sum = mpi::AllReduce( localSum, A.Comm() );
    return Pow( sum, Real(1)/p );
}

template<typename F> 
Base<F> EntrywiseNorm( const DistMultiVec<F>& A, Base<F> p )
{
    DEBUG_ONLY(CSE cse("EntrywiseNorm"))
    typedef Base<F> Real;

    Real localSum = 0;
    for( Int j=0; j<A.Width(); ++j )
        for( Int iLoc=0; iLoc<A.LocalHeight(); ++iLoc )
            localSum += Pow( Abs(A.GetLocal(iLoc,j)), p ); 

    const Real sum = mpi::AllReduce( localSum, A.Comm() );
    return Pow( sum, Real(1)/p );
}

template<typename F>
Base<F> HermitianEntrywiseNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CSE cse("HermitianEntrywiseNorm"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    typedef Base<F> Real;
    Real sum;
    if( A.Participating() )
    {
        Real localSum = 0;
        const Int localWidth = A.LocalWidth();
        if( uplo == UPPER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int numUpperRows = A.LocalRowOffset(j+1);
                for( Int iLoc=0; iLoc<numUpperRows; ++iLoc )
                {
                    const Int i = A.GlobalRow(iLoc);
                    const Real term = Pow( Abs(A.GetLocal(iLoc,jLoc)), p );
                    if( i ==j )
                        localSum += term;
                    else
                        localSum += 2*term;
                }
            }
        }
        else
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int numStrictlyUpperRows = A.LocalRowOffset(j);
                for( Int iLoc=numStrictlyUpperRows;
                     iLoc<A.LocalHeight(); ++iLoc )
                {
                    const Int i = A.GlobalRow(iLoc);
                    const Real term = Pow( Abs(A.GetLocal(iLoc,jLoc)), p );
                    if( i ==j )
                        localSum += term;
                    else
                        localSum += 2*term;
                }
            }
        }
        sum = mpi::AllReduce( localSum, A.DistComm() );
    }
    mpi::Broadcast( sum, A.Root(), A.CrossComm() );
    return Pow( sum, 1/p );
}

template<typename F> 
Base<F> HermitianEntrywiseNorm
( UpperOrLower uplo, const DistSparseMatrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CSE cse("HermitianEntrywiseNorm"))
    typedef Base<F> Real;

    Real localSum = 0;
    const Int numLocalEntries = A.NumLocalEntries();
    for( Int k=0; k<numLocalEntries; ++k )
    {
        const Int i = A.Row(k);
        const Int j = A.Col(k); 
        if( (uplo==UPPER && i<j) || (uplo==LOWER && i>j) )
            localSum += 2*Pow( Abs(A.Value(k)), p ); 
        else if( i == j )
            localSum += Pow( Abs(A.Value(k)), p ); 
    }

    const Real sum = mpi::AllReduce( localSum, A.Comm() );
    return Pow( sum, Real(1)/p );
}

template<typename F>
Base<F> SymmetricEntrywiseNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CSE cse("SymmetricEntrywiseNorm"))
    return HermitianEntrywiseNorm( uplo, A, p );
}

template<typename F>
Base<F> SymmetricEntrywiseNorm
( UpperOrLower uplo, const DistSparseMatrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CSE cse("SymmetricEntrywiseNorm"))
    return HermitianEntrywiseNorm( uplo, A, p );
}

#define PROTO(F) \
  template Base<F> EntrywiseNorm( const Matrix<F>& A, Base<F> p ); \
  template Base<F> EntrywiseNorm( const AbstractDistMatrix<F>& A, Base<F> p ); \
  template Base<F> EntrywiseNorm( const SparseMatrix<F>& A, Base<F> p ); \
  template Base<F> EntrywiseNorm( const DistSparseMatrix<F>& A, Base<F> p ); \
  template Base<F> EntrywiseNorm( const DistMultiVec<F>& A, Base<F> p ); \
  template Base<F> HermitianEntrywiseNorm \
  ( UpperOrLower uplo, const Matrix<F>& A, Base<F> p ); \
  template Base<F> HermitianEntrywiseNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A, Base<F> p ); \
  template Base<F> HermitianEntrywiseNorm \
  ( UpperOrLower uplo, const SparseMatrix<F>& A, Base<F> p ); \
  template Base<F> HermitianEntrywiseNorm \
  ( UpperOrLower uplo, const DistSparseMatrix<F>& A, Base<F> p ); \
  template Base<F> SymmetricEntrywiseNorm \
  ( UpperOrLower uplo, const Matrix<F>& A, Base<F> p ); \
  template Base<F> SymmetricEntrywiseNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A, Base<F> p ); \
  template Base<F> SymmetricEntrywiseNorm \
  ( UpperOrLower uplo, const SparseMatrix<F>& A, Base<F> p ); \
  template Base<F> SymmetricEntrywiseNorm \
  ( UpperOrLower uplo, const DistSparseMatrix<F>& A, Base<F> p );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
