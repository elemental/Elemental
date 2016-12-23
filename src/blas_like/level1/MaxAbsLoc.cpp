/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>

namespace El {

// TODO(poulson): Add options for OneAbs instead of Abs

template<typename Ring>
ValueInt<Base<Ring>> VectorMaxAbsLoc( const Matrix<Ring>& x )
{
    EL_DEBUG_CSE
    typedef Base<Ring> RealRing;
    const Int m = x.Height();
    const Int n = x.Width();
    EL_DEBUG_ONLY(
      if( m != 1 && n != 1 )
          LogicError("Input should have been a vector");
    )

    ValueInt<RealRing> pivot;
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
            const RealRing absVal = Abs(x(i));
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
            const RealRing absVal = Abs(x(0,j));
            if( absVal > pivot.value )
            {
                pivot.index = j;
                pivot.value = absVal;
            }
        }
    }
    return pivot;
}

template<typename Ring>
ValueInt<Base<Ring>> VectorMaxAbsLoc( const AbstractDistMatrix<Ring>& x )
{
    EL_DEBUG_CSE
    typedef Base<Ring> RealRing;
    const Int m = x.Height();
    const Int n = x.Width();
    EL_DEBUG_ONLY(
        if( m != 1 && n != 1 )
            LogicError("Input should have been a vector");
        if( !x.Grid().InGrid() )
            LogicError("Viewing processes are not allowed");
    )
    ValueInt<RealRing> pivot;
    if( Min(m,n) == 0 )
    {
        pivot.index = -1;
        pivot.value = 0;
        return pivot;
    }

    if( x.Participating() )
    {
        ValueInt<RealRing> localPivot;
        localPivot.index = 0;
        localPivot.value = 0;
        if( n == 1 )
        {
            if( x.RowRank() == x.RowAlign() )
            {
                const Int mLocal = x.LocalHeight();
                for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                {
                    const RealRing absVal = Abs(x.GetLocal(iLoc,0));
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
                    const RealRing absVal = Abs(x.GetLocal(0,jLoc));
                    if( absVal > localPivot.value )
                    {
                        localPivot.index = x.GlobalCol(jLoc);
                        localPivot.value = absVal;
                    }
                }
            }
        }
        pivot = mpi::AllReduce
                ( localPivot, mpi::MaxLocOp<RealRing>(), x.DistComm() );
    }
    mpi::Broadcast( pivot, x.Root(), x.CrossComm() );
    return pivot;
}

template<typename Ring>
Entry<Base<Ring>> MaxAbsLoc( const Matrix<Ring>& A )
{
    EL_DEBUG_CSE
    typedef Base<Ring> RealRing;
    const Int m = A.Height();
    const Int n = A.Width();

    Entry<RealRing> pivot;
    if( Min(m,n) == 0 )
    {
        pivot.i = -1;
        pivot.j = -1;
        pivot.value = 0;
        return pivot;
    }

    pivot.i = 0;
    pivot.j = 0;
    pivot.value = 0;
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            const RealRing absVal = Abs(A(i,j));
            if( absVal > pivot.value )
            {
                pivot.i = i;
                pivot.j = j;
                pivot.value = absVal;
            }
        }
    }
    return pivot;
}

template<typename Ring>
Entry<Base<Ring>> MaxAbsLoc( const AbstractDistMatrix<Ring>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( !A.Grid().InGrid() )
          LogicError("Viewing processes are not allowed");
    )
    typedef Base<Ring> RealRing;
    Entry<RealRing> pivot;
    if( A.Height() == 0 )
    {
        pivot.i = -1;
        pivot.j = -1;
        pivot.value = 0;
        return pivot;
    }

    if( A.Participating() )
    {
        // Store the index/value of the local pivot candidate
        Entry<RealRing> localPivot;
        localPivot.i = 0;
        localPivot.j = 0;
        localPivot.value = 0;
        const Int mLocal = A.LocalHeight();
        const Int nLocal = A.LocalWidth();
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            {
                const RealRing value = Abs(A.GetLocal(iLoc,jLoc));
                if( value > localPivot.value )
                {
                    const Int i = A.GlobalRow(iLoc);
                    localPivot.i = i;
                    localPivot.j = j;
                    localPivot.value = value;
                }
            }
        }

        // Compute and store the location of the new pivot
        pivot = mpi::AllReduce
                ( localPivot, mpi::MaxLocPairOp<RealRing>(), A.DistComm() );
    }
    mpi::Broadcast( pivot, A.Root(), A.CrossComm() );
    return pivot;
}

template<typename Ring>
Entry<Base<Ring>> MaxAbsLoc( const SparseMatrix<Ring>& A )
{
    EL_DEBUG_CSE
    typedef Base<Ring> RealRing;

    Entry<RealRing> pivot;
    if( A.Height() == 0 || A.Width() == 0 )
    {
        pivot.i = -1;
        pivot.j = -1;
        pivot.value = 0;
        return pivot;
    }

    pivot.i = 0;
    pivot.j = 0;
    pivot.value = 0;
    const Int numEntries = A.NumEntries();
    for( Int e=0; e<numEntries; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        const RealRing absVal = Abs(A.Value(e));
        if( absVal > pivot.value )
        {
            pivot.i = i;
            pivot.j = j;
            pivot.value = absVal;
        }
    }
    return pivot;
}

template<typename Ring>
Entry<Base<Ring>> MaxAbsLoc( const DistSparseMatrix<Ring>& A )
{
    EL_DEBUG_CSE
    typedef Base<Ring> RealRing;

    Entry<RealRing> pivot;
    if( A.Height() == 0 || A.Width() == 0 )
    {
        pivot.i = -1;
        pivot.j = -1;
        pivot.value = 0;
        return pivot;
    }

    pivot.i = 0;
    pivot.j = 0;
    pivot.value = 0;
    const Int numLocalEntries = A.NumLocalEntries();
    for( Int e=0; e<numLocalEntries; ++e )
    {
        const Int i = A.Row( e );
        const Int j = A.Col( e );
        const RealRing absVal = Abs(A.Value( e ));
        if( absVal > pivot.value )
        {
            pivot.i = i;
            pivot.j = j;
            pivot.value = absVal;
        }
    }
    return mpi::AllReduce
      ( pivot, mpi::MaxLocPairOp<RealRing>(), A.Grid().Comm() );
}

template<typename Ring>
Entry<Base<Ring>> SymmetricMaxAbsLoc( UpperOrLower uplo, const Matrix<Ring>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
    )
    typedef Base<Ring> RealRing;
    const Int n = A.Width();

    Entry<RealRing> pivot;
    if( n == 0 )
    {
        pivot.i = -1;
        pivot.j = -1;
        pivot.value = 0;
        return pivot;
    }

    pivot.i = 0;
    pivot.j = 0;
    pivot.value = 0;
    for( Int j=0; j<n; ++j )
    {
        if( uplo == LOWER )
        {
            for( Int i=j; i<n; ++i )
            {
                const RealRing absVal = Abs(A(i,j));
                if( absVal > pivot.value )
                {
                    pivot.i = i;
                    pivot.j = j;
                    pivot.value = absVal;
                }
            }
        }
        else
        {
            for( Int i=0; i<=j; ++i )
            {
                const RealRing absVal = Abs(A(i,j));
                if( absVal > pivot.value )
                {
                    pivot.i = i;
                    pivot.j = j;
                    pivot.value = absVal;
                }
            }
        }
    }
    return pivot;
}

template<typename Ring>
Entry<Base<Ring>> SymmetricMaxAbsLoc
( UpperOrLower uplo, const AbstractDistMatrix<Ring>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( !A.Grid().InGrid() )
          LogicError("Viewing processes are not allowed");
    )
    typedef Base<Ring> RealRing;
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();

    Entry<RealRing> pivot;
    if( A.Height() == 0 )
    {
        pivot.i = -1;
        pivot.j = -1;
        pivot.value = 0;
        return pivot;
    }

    if( A.Participating() )
    {
        Entry<RealRing> localPivot;
        localPivot.i = 0;
        localPivot.j = 0;
        localPivot.value = 0;

        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            if( uplo == LOWER )
            {
                const Int mLocBefore = A.LocalRowOffset(j);
                for( Int iLoc=mLocBefore; iLoc<mLocal; ++iLoc )
                {
                    const RealRing absVal = Abs(A.GetLocal(iLoc,jLoc));
                    if( absVal > localPivot.value )
                    {
                        const Int i = A.GlobalRow(iLoc);
                        localPivot.i = i;
                        localPivot.j = j;
                        localPivot.value = absVal;
                    }
                }
            }
            else
            {
                const Int mLocBefore = A.LocalRowOffset(j+1);
                for( Int iLoc=0; iLoc<mLocBefore; ++iLoc )
                {
                    const RealRing absVal = Abs(A.GetLocal(iLoc,jLoc));
                    if( absVal > localPivot.value )
                    {
                        const Int i = A.GlobalRow(iLoc);
                        localPivot.i = i;
                        localPivot.j = j;
                        localPivot.value = absVal;
                    }
                }
            }
        }

        // Compute and store the location of the new pivot
        pivot = mpi::AllReduce
                ( localPivot, mpi::MaxLocPairOp<RealRing>(), A.DistComm() );
    }
    mpi::Broadcast( pivot, A.Root(), A.CrossComm() );
    return pivot;
}

template<typename Ring>
Entry<Base<Ring>> SymmetricMaxAbsLoc
( UpperOrLower uplo, const SparseMatrix<Ring>& A )
{
    EL_DEBUG_CSE
    typedef Base<Ring> RealRing;

    Entry<RealRing> pivot;
    if( A.Height() == 0 || A.Width() == 0 )
    {
        pivot.i = -1;
        pivot.j = -1;
        pivot.value = 0;
        return pivot;
    }

    pivot.i = 0;
    pivot.j = 0;
    pivot.value = 0;
    const Int numEntries = A.NumEntries();
    for( Int e=0; e<numEntries; ++e )
    {
        const Int i = A.Row( e );
        const Int j = A.Col( e );
        const RealRing absVal = Abs(A.Value( e ));
        const bool valid = ( uplo==LOWER ? i>=j : i<=j );
        if( valid && absVal > pivot.value )
        {
            pivot.i = i;
            pivot.j = j;
            pivot.value = absVal;
        }
    }
    return pivot;
}

template<typename Ring>
Entry<Base<Ring>> SymmetricMaxAbsLoc
( UpperOrLower uplo, const DistSparseMatrix<Ring>& A )
{
    EL_DEBUG_CSE
    typedef Base<Ring> RealRing;

    Entry<RealRing> pivot;
    if( A.Height() == 0 || A.Width() == 0 )
    {
        pivot.i = -1;
        pivot.j = -1;
        pivot.value = 0;
        return pivot;
    }

    pivot.i = 0;
    pivot.j = 0;
    pivot.value = 0;
    const Int numLocalEntries = A.NumLocalEntries();
    for( Int e=0; e<numLocalEntries; ++e )
    {
        const Int i = A.Row( e );
        const Int j = A.Col( e );
        const RealRing absVal = Abs(A.Value( e ));
        const bool valid = ( uplo==LOWER ? i>=j : i<=j );
        if( valid && absVal > pivot.value )
        {
            pivot.i = i;
            pivot.j = j;
            pivot.value = absVal;
        }
    }
    return mpi::AllReduce
      ( pivot, mpi::MaxLocPairOp<RealRing>(), A.Grid().Comm() );
}

#define PROTO(Ring) \
  template ValueInt<Base<Ring>> VectorMaxAbsLoc \
  ( const Matrix<Ring>& x ); \
  template ValueInt<Base<Ring>> VectorMaxAbsLoc \
  ( const AbstractDistMatrix<Ring>& x ); \
  template Entry<Base<Ring>> MaxAbsLoc( const Matrix<Ring>& x ); \
  template Entry<Base<Ring>> MaxAbsLoc( const AbstractDistMatrix<Ring>& x ); \
  template Entry<Base<Ring>> MaxAbsLoc( const SparseMatrix<Ring>& x ); \
  template Entry<Base<Ring>> MaxAbsLoc( const DistSparseMatrix<Ring>& x ); \
  template Entry<Base<Ring>> SymmetricMaxAbsLoc \
  ( UpperOrLower uplo, const Matrix<Ring>& A ); \
  template Entry<Base<Ring>> SymmetricMaxAbsLoc \
  ( UpperOrLower uplo, const AbstractDistMatrix<Ring>& A ); \
  template Entry<Base<Ring>> SymmetricMaxAbsLoc \
  ( UpperOrLower uplo, const SparseMatrix<Ring>& x ); \
  template Entry<Base<Ring>> SymmetricMaxAbsLoc \
  ( UpperOrLower uplo, const DistSparseMatrix<Ring>& x );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
