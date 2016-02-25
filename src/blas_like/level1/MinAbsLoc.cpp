/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// TODO: Add options for FastAbs instead of Abs

template<typename F>
ValueInt<Base<F>> VectorMinAbsLoc( const Matrix<F>& x )
{
    DEBUG_ONLY(CSE cse("VectorMinAbsLoc"))
    typedef Base<F> Real;
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

    pivot.value = Abs(x.Get(0,0));
    pivot.index = 0;
    if( n == 1 )
    {
        for( Int i=1; i<m; ++i )
        {
            const Real absVal = Abs(x.Get(i,0));
            if( absVal < pivot.value )
            {
                pivot.index = i;
                pivot.value = absVal;
            }
        }
    }
    else
    {
        for( Int j=1; j<n; ++j )
        {
            const Real absVal = Abs(x.Get(0,j));
            if( absVal < pivot.value )
            {
                pivot.index = j;
                pivot.value = absVal;
            }
        }
    }
    return pivot;
}

template<typename F>
ValueInt<Base<F>> VectorMinAbsLoc( const AbstractDistMatrix<F>& x )
{
    DEBUG_ONLY(CSE cse("VectorMinAbsLoc"))
    typedef Base<F> Real;
    const Int m = x.Height();
    const Int n = x.Width();
    DEBUG_ONLY(
      if( m != 1 && n != 1 )
          LogicError("Input should have been a vector");
      if( !x.Grid().InGrid() )
          LogicError("viewing processes are not allowed");
    )
    ValueInt<Real> pivot;
    if( Min(m,n) == 0 )
    {
        pivot.value = 0;
        pivot.index = -1;
        return pivot;
    }
 
    ValueInt<Real> localPivot;
    localPivot.value = Abs(x.Get(0,0));
    localPivot.index = 0;
    if( x.Participating() )
    {
        if( n == 1 )
        {
            if( x.RowRank() == x.RowAlign() )
            {
                const Int mLocal = x.LocalHeight();
                for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                {
                    const Real absVal = Abs(x.GetLocal(iLoc,0));
                    if( absVal < localPivot.value )
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
                    if( absVal < localPivot.value )
                    {
                        localPivot.index = x.GlobalCol(jLoc);
                        localPivot.value = absVal;
                    }
                }
            }
        }
        pivot = mpi::AllReduce
                ( localPivot, mpi::MinLocOp<Real>(), x.DistComm() );
    }
    mpi::Broadcast( pivot, x.Root(), x.CrossComm() );
    return pivot;
}

template<typename F>
Entry<Base<F>> MinAbsLoc( const Matrix<F>& A )
{
    DEBUG_ONLY(CSE cse("MinAbsLoc"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();

    Entry<Real> pivot;
    if( Min(m,n) == 0 )
    {
        pivot.i = -1;
        pivot.j = -1;
        pivot.value = 0;
        return pivot;
    }

    pivot.i = 0;
    pivot.j = 0;
    pivot.value = Abs(A.Get(0,0));
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            const Real absVal = Abs(A.Get(i,j));
            if( absVal < pivot.value )
            {
                pivot.i = i;
                pivot.j = j;
                pivot.value = absVal;
            }
        }
    }
    return pivot;
}

template<typename F>
Entry<Base<F>> MinAbsLoc( const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(
      CSE cse("MinAbs");
      if( !A.Grid().InGrid() )
          LogicError("Viewing processes are not allowed");
    )
    typedef Base<F> Real;
    Entry<Real> pivot;
    if( Min(A.Height(),A.Width()) == 0 )
    {
        pivot.i = -1;
        pivot.j = -1;
        pivot.value = 0;
        return pivot;
    }

    Entry<Real> localPivot;
    localPivot.i = 0;
    localPivot.j = 0;
    localPivot.value = Abs(A.Get(0,0));
    if( A.Participating() )
    {
        // Store the index/value of the local pivot candidate
        const Int mLocal = A.LocalHeight();
        const Int nLocal = A.LocalWidth();
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            {
                const Real value = Abs(A.GetLocal(iLoc,jLoc));
                if( value < localPivot.value )
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
                ( localPivot, mpi::MinLocPairOp<Real>(), A.DistComm() );
    }
    mpi::Broadcast( pivot, A.Root(), A.CrossComm() );
    return pivot;
}

template<typename F>
Entry<Base<F>> SymmetricMinAbsLoc( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(
      CSE cse("SymmetricMinAbs");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
    )
    typedef Base<F> Real;
    const Int n = A.Width();
    Entry<Real> pivot;
    if( n == 0 )
    {
        pivot.i = -1;
        pivot.j = -1;
        pivot.value = 0;
        return pivot;
    }

    pivot.i = 0;
    pivot.j = 0;
    pivot.value = Abs(A.Get(0,0));
    if( uplo == LOWER )
    {
        for( Int j=0; j<n; ++j )
        {
            for( Int i=j; i<n; ++i )
            {
                const Real absVal = Abs(A.Get(i,j));
                if( absVal < pivot.value )
                {
                    pivot.i = i;
                    pivot.j = j;
                    pivot.value = absVal;
                }
            }
        }
    }
    else
    {
        for( Int j=0; j<n; ++j ) 
        { 
            for( Int i=0; i<=j; ++i )
            {
                const Real absVal = Abs(A.Get(i,j));
                if( absVal < pivot.value )
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

template<typename F>
Entry<Base<F>>
SymmetricMinAbsLoc( UpperOrLower uplo, const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(
      CSE cse("SymmetricMinAbs");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( !A.Grid().InGrid() )
          LogicError("Viewing processes are not allowed");
    )
    typedef Base<F> Real;
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();

    Entry<Real> pivot;
    if( A.Height() == 0 )
    {
        pivot.i = -1;
        pivot.j = -1;
        pivot.value = 0;
        return pivot;
    }

    Entry<Real> localPivot;
    localPivot.i = 0;
    localPivot.j = 0;
    localPivot.value = Abs(A.Get(0,0));
    if( A.Participating() )
    {
        if( uplo == LOWER )
        {
            for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int mLocBefore = A.LocalRowOffset(j);
                for( Int iLoc=mLocBefore; iLoc<mLocal; ++iLoc )
                {
                    const Real absVal = Abs(A.GetLocal(iLoc,jLoc));
                    if( absVal < localPivot.value )
                    {
                        const Int i = A.GlobalRow(iLoc);
                        localPivot.i = i;
                        localPivot.j = j;
                        localPivot.value = absVal;
                    }
                }
            }
        }
        else
        {
            for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int mLocBefore = A.LocalRowOffset(j+1);
                for( Int iLoc=0; iLoc<mLocBefore; ++iLoc )
                {
                    const Real absVal = Abs(A.GetLocal(iLoc,jLoc));
                    if( absVal < localPivot.value )
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
                ( localPivot, mpi::MinLocPairOp<Real>(), A.DistComm() );
    }
    mpi::Broadcast( pivot, A.Root(), A.CrossComm() );
    return pivot;
}

#define PROTO(F) \
  template ValueInt<Base<F>> VectorMinAbsLoc( const Matrix<F>& x ); \
  template ValueInt<Base<F>> VectorMinAbsLoc \
  ( const AbstractDistMatrix<F>& x ); \
  template Entry<Base<F>> MinAbsLoc( const Matrix<F>& x ); \
  template Entry<Base<F>> MinAbsLoc( const AbstractDistMatrix<F>& x ); \
  template Entry<Base<F>> SymmetricMinAbsLoc \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Entry<Base<F>> SymmetricMinAbsLoc \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
