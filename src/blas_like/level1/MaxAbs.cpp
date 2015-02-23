/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// TODO: Add options for FastAbs instead of Abs

template<typename T>
ValueInt<Base<T>> VectorMaxAbs( const Matrix<T>& x )
{
    DEBUG_ONLY(CallStackEntry cse("VectorMaxAbs"))
    typedef Base<T> Real;
    const Int m = x.Height();
    const Int n = x.Width();
    DEBUG_ONLY(
        if( m != 1 && n != 1 )
            LogicError("Input should have been a vector");
    )

    ValueInt<Real> pivot;
    if( Min(m,n) == 0 )
    {
        pivot.value = 0;
        pivot.index = -1;
        return pivot;
    }

    pivot.value = 0;
    pivot.index = 0;
    if( n == 1 )
    {
        for( Int i=0; i<m; ++i )
        {
            const Real absVal = Abs(x.Get(i,0));
            if( absVal > pivot.value )
            {
                pivot.index = i;
                pivot.value = absVal;
            }
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            const Real absVal = Abs(x.Get(0,j));
            if( absVal > pivot.value )
            {
                pivot.index = j;
                pivot.value = absVal;
            }
        }
    }
    return pivot;
}

template<typename T>
ValueInt<Base<T>> VectorMaxAbs( const AbstractDistMatrix<T>& x )
{
    DEBUG_ONLY(CallStackEntry cse("VectorMaxAbs"))
    typedef Base<T> Real;
    const Int m = x.Height();
    const Int n = x.Width();
    DEBUG_ONLY(
        if( m != 1 && n != 1 )
            LogicError("Input should have been a vector");
        if( !x.Grid().InGrid() )
            LogicError("Viewing processes are not allowed");
    )
    ValueInt<Real> pivot;
    if( Min(m,n) == 0 )
    {
        pivot.index = -1;
        pivot.value = 0;
        return pivot;
    }

    if( x.Participating() )
    {
        ValueInt<Real> localPivot;
        localPivot.index = 0;
        localPivot.value = 0;
        if( n == 1 )
        {
            if( x.RowRank() == x.RowAlign() )
            {
                const Int mLocal = x.LocalHeight();
                for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                {
                    const Real absVal = Abs(x.GetLocal(iLoc,0));
                    if( absVal > localPivot.value )
                    {
                        localPivot.index = x.GlobalRow(iLoc);
                        localPivot.value = absVal;
                    }
                }
            }
        }
        else
        {
            if( x.ColRank() == x.ColAlign() )
            {
                const Int nLocal = x.LocalWidth();
                for( Int jLoc=0; jLoc<nLocal; ++jLoc )
                {
                    const Real absVal = Abs(x.GetLocal(0,jLoc));
                    if( absVal > localPivot.value )
                    {
                        localPivot.index = x.GlobalCol(jLoc);
                        localPivot.value = absVal;
                    }
                }
            }
        }
        pivot = mpi::AllReduce
                ( localPivot, mpi::MaxLocOp<Real>(), x.DistComm() );
    }
    mpi::Broadcast( pivot.index, x.Root(), x.CrossComm() );
    mpi::Broadcast( pivot.value, x.Root(), x.CrossComm() );
    return pivot;
}

template<typename T>
ValueIntPair<Base<T>> MaxAbs( const Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MaxAbs"))
    typedef Base<T> Real;
    const Int m = A.Height();
    const Int n = A.Width();

    ValueIntPair<Real> pivot;
    if( Min(m,n) == 0 )
    {
        pivot.value = 0;
        pivot.indices[0] = -1;
        pivot.indices[1] = -1;
        return pivot;
    }

    pivot.value = 0;
    pivot.indices[0] = 0;
    pivot.indices[1] = 0;
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            const Real absVal = Abs(A.Get(i,j));
            if( absVal > pivot.value )
            {
                pivot.value = absVal;
                pivot.indices[0] = i;
                pivot.indices[1] = j;
            }
        }
    }
    return pivot;
}

template<typename T>
ValueIntPair<Base<T>> MaxAbs( const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("MaxAbs");
        if( !A.Grid().InGrid() )
            LogicError("Viewing processes are not allowed");
    )
    typedef Base<T> Real;
    ValueIntPair<Real> pivot;
    if( A.Height() == 0 )
    {
        pivot.value = 0;
        pivot.indices[0] = -1;
        pivot.indices[1] = -1;
        return pivot;
    }

    if( A.Participating() )
    {
        // Store the index/value of the local pivot candidate
        ValueIntPair<Real> localPivot;
        localPivot.value = 0;
        localPivot.indices[0] = 0;
        localPivot.indices[1] = 0;
        const Int mLocal = A.LocalHeight();
        const Int nLocal = A.LocalWidth();
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            {
                const Real value = Abs(A.GetLocal(iLoc,jLoc));
                if( value > localPivot.value )
                {
                    const Int i = A.GlobalRow(iLoc);
                    localPivot.value = value;
                    localPivot.indices[0] = i;
                    localPivot.indices[1] = j;
                }
            }
        }

        // Compute and store the location of the new pivot
        pivot = mpi::AllReduce
                ( localPivot, mpi::MaxLocPairOp<Real>(), A.DistComm() );
    }
    mpi::Broadcast( pivot.indices, 2, A.Root(), A.CrossComm() );
    mpi::Broadcast( pivot.value, A.Root(), A.CrossComm() );
    return pivot;
}

template<typename T>
ValueIntPair<Base<T>> MaxAbs( const SparseMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MaxAbs"))
    typedef Base<T> Real;

    ValueIntPair<Real> pivot;
    if( A.Height() == 0 || A.Width() == 0 )
    {
        pivot.value = 0;
        pivot.indices[0] = -1;
        pivot.indices[1] = -1;
        return pivot;
    }

    pivot.value = 0;
    pivot.indices[0] = 0;
    pivot.indices[1] = 0;
    const Int numEntries = A.NumEntries();
    for( Int e=0; e<numEntries; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        const Real absVal = Abs(A.Value(e));
        if( absVal > pivot.value )
        {
            pivot.value = absVal;
            pivot.indices[0] = i;
            pivot.indices[1] = j;
        }
    }
    return pivot;
}

template<typename T>
ValueIntPair<Base<T>> MaxAbs( const DistSparseMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MaxAbs"))
    typedef Base<T> Real;

    ValueIntPair<Real> pivot;
    if( A.Height() == 0 || A.Width() == 0 )
    {
        pivot.value = 0;
        pivot.indices[0] = -1;
        pivot.indices[1] = -1;
        return pivot;
    }

    pivot.value = 0;
    pivot.indices[0] = 0;
    pivot.indices[1] = 0;
    const Int numLocalEntries = A.NumLocalEntries();
    for( Int e=0; e<numLocalEntries; ++e )
    {
        const Int i = A.Row( e );
        const Int j = A.Col( e );
        const Real absVal = Abs(A.Value( e ));
        if( absVal > pivot.value )
        {
            pivot.value = absVal;
            pivot.indices[0] = i;
            pivot.indices[1] = j;
        }
    }
    return mpi::AllReduce( pivot, mpi::MaxLocPairOp<Real>(), A.Comm() );
}

template<typename T>
ValueIntPair<Base<T>> SymmetricMaxAbs( UpperOrLower uplo, const Matrix<T>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("SymmetricMaxAbs");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    typedef Base<T> Real;
    const Int n = A.Width();

    ValueIntPair<Real> pivot;
    if( n == 0 )
    {
        pivot.value = 0;
        pivot.indices[0] = -1;
        pivot.indices[1] = -1;
        return pivot;
    }

    pivot.value = 0;
    pivot.indices[0] = 0;
    pivot.indices[1] = 0;
    for( Int j=0; j<n; ++j )
    {
        if( uplo == LOWER )
        {
            for( Int i=j; i<n; ++i )
            {
                const Real absVal = Abs(A.Get(i,j));
                if( absVal > pivot.value )
                {
                    pivot.value = absVal;
                    pivot.indices[0] = i;
                    pivot.indices[1] = j;
                }
            }
        }
        else
        {
            for( Int i=0; i<=j; ++i )
            {
                const Real absVal = Abs(A.Get(i,j));
                if( absVal > pivot.value )
                {
                    pivot.value = absVal;
                    pivot.indices[0] = i;
                    pivot.indices[1] = j;
                }
            }
        }
    }
    return pivot;
}

template<typename T>
ValueIntPair<Base<T>> SymmetricMaxAbs
( UpperOrLower uplo, const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("SymmetricMaxAbs");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( !A.Grid().InGrid() )
            LogicError("Viewing processes are not allowed");
    )
    typedef Base<T> Real;
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();

    ValueIntPair<Real> pivot;
    if( A.Height() == 0 )
    {
        pivot.value = 0;
        pivot.indices[0] = -1;
        pivot.indices[1] = -1;
        return pivot;
    }

    if( A.Participating() )
    {
        ValueIntPair<Real> localPivot;
        localPivot.value = 0;
        localPivot.indices[0] = 0;
        localPivot.indices[1] = 0;

        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            if( uplo == LOWER )
            {
                const Int mLocBefore = A.LocalRowOffset(j);
                for( Int iLoc=mLocBefore; iLoc<mLocal; ++iLoc )
                {
                    const Real absVal = Abs(A.GetLocal(iLoc,jLoc));
                    if( absVal > localPivot.value )
                    {
                        const Int i = A.GlobalRow(iLoc);
                        localPivot.value = absVal;
                        localPivot.indices[0] = i;
                        localPivot.indices[1] = j;
                    }
                }
            }
            else
            {
                const Int mLocBefore = A.LocalRowOffset(j+1);
                for( Int iLoc=0; iLoc<mLocBefore; ++iLoc )
                {
                    const Real absVal = Abs(A.GetLocal(iLoc,jLoc));
                    if( absVal > localPivot.value )
                    {
                        const Int i = A.GlobalRow(iLoc);
                        localPivot.value = absVal;
                        localPivot.indices[0] = i;
                        localPivot.indices[1] = j;
                    }
                }
            }
        }

        // Compute and store the location of the new pivot
        pivot = mpi::AllReduce
                ( localPivot, mpi::MaxLocPairOp<Real>(), A.DistComm() );
    }
    mpi::Broadcast( pivot.indices, 2, A.Root(), A.CrossComm() );
    mpi::Broadcast( pivot.value, A.Root(), A.CrossComm() );
    return pivot;
}

template<typename T>
ValueIntPair<Base<T>> SymmetricMaxAbs
( UpperOrLower uplo, const SparseMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricMaxAbs"))
    typedef Base<T> Real;

    ValueIntPair<Real> pivot;
    if( A.Height() == 0 || A.Width() == 0 )
    {
        pivot.value = 0;
        pivot.indices[0] = -1;
        pivot.indices[1] = -1;
        return pivot;
    }

    pivot.value = 0;
    pivot.indices[0] = 0;
    pivot.indices[1] = 0;
    const Int numEntries = A.NumEntries();
    for( Int e=0; e<numEntries; ++e )
    {
        const Int i = A.Row( e );
        const Int j = A.Col( e );
        const Real absVal = Abs(A.Value( e ));
        const bool valid = ( uplo==LOWER ? i>=j : i<=j );
        if( valid && absVal > pivot.value )
        {
            pivot.value = absVal;
            pivot.indices[0] = i;
            pivot.indices[1] = j;
        }
    }
    return pivot;
}

template<typename T>
ValueIntPair<Base<T>> SymmetricMaxAbs
( UpperOrLower uplo, const DistSparseMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricMaxAbs"))
    typedef Base<T> Real;

    ValueIntPair<Real> pivot;
    if( A.Height() == 0 || A.Width() == 0 )
    {
        pivot.value = 0;
        pivot.indices[0] = -1;
        pivot.indices[1] = -1;
        return pivot;
    }

    pivot.value = 0;
    pivot.indices[0] = 0;
    pivot.indices[1] = 0;
    const Int numLocalEntries = A.NumLocalEntries();
    for( Int e=0; e<numLocalEntries; ++e )
    {
        const Int i = A.Row( e );
        const Int j = A.Col( e );
        const Real absVal = Abs(A.Value( e ));
        const bool valid = ( uplo==LOWER ? i>=j : i<=j );
        if( valid && absVal > pivot.value )
        {
            pivot.value = absVal;
            pivot.indices[0] = i;
            pivot.indices[1] = j;
        }
    }
    return mpi::AllReduce( pivot, mpi::MaxLocPairOp<Real>(), A.Comm() );
}

#define PROTO(T) \
  template ValueInt<Base<T>> VectorMaxAbs( const Matrix<T>& x ); \
  template ValueInt<Base<T>> VectorMaxAbs( const AbstractDistMatrix<T>& x ); \
  template ValueIntPair<Base<T>> MaxAbs( const Matrix<T>& x ); \
  template ValueIntPair<Base<T>> MaxAbs( const AbstractDistMatrix<T>& x ); \
  template ValueIntPair<Base<T>> MaxAbs( const SparseMatrix<T>& x ); \
  template ValueIntPair<Base<T>> MaxAbs( const DistSparseMatrix<T>& x ); \
  template ValueIntPair<Base<T>> SymmetricMaxAbs \
  ( UpperOrLower uplo, const Matrix<T>& A ); \
  template ValueIntPair<Base<T>> SymmetricMaxAbs \
  ( UpperOrLower uplo, const AbstractDistMatrix<T>& A ); \
  template ValueIntPair<Base<T>> SymmetricMaxAbs \
  ( UpperOrLower uplo, const SparseMatrix<T>& x ); \
  template ValueIntPair<Base<T>> SymmetricMaxAbs \
  ( UpperOrLower uplo, const DistSparseMatrix<T>& x );

#include "El/macros/Instantiate.h"

} // namespace El
