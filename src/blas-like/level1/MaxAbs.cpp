/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

// TODO: Switch from F to T
// TODO: Add options for FastAbs instead of Abs

template<typename F>
ValueInt<Base<F>> VectorMaxAbs( const Matrix<F>& x )
{
    DEBUG_ONLY(CallStackEntry cse("VectorMaxAbs"))
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

    pivot.value = 0;
    pivot.index = 0;
    if( n == 1 )
    {
        for( Int i=0; i<m; ++i )
        {
            const Real abs = Abs(x.Get(i,0));
            if( abs > pivot.value )
            {
                pivot.index = i;
                pivot.value = abs;
            }
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            const Real abs = Abs(x.Get(0,j));
            if( abs > pivot.value )
            {
                pivot.index = j;
                pivot.value = abs;
            }
        }
    }
    return pivot;
}

template<typename F>
ValueInt<Base<F>> VectorMaxAbs( const AbstractDistMatrix<F>& x )
{
    DEBUG_ONLY(CallStackEntry cse("VectorMaxAbs"))
    typedef Base<F> Real;
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
                    const Real abs = Abs(x.GetLocal(iLoc,0));
                    if( abs > localPivot.value )
                    {
                        localPivot.index = x.GlobalRow(iLoc);
                        localPivot.value = abs;
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
                    const Real abs = Abs(x.GetLocal(0,jLoc));
                    if( abs > localPivot.value )
                    {
                        localPivot.index = x.GlobalCol(jLoc);
                        localPivot.value = abs;
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

template<typename F>
ValueIntPair<Base<F>> MaxAbs( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MaxAbs"))
    typedef Base<F> Real;
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
            const Real abs = Abs(A.Get(i,j));
            if( abs > pivot.value )
            {
                pivot.value = abs;
                pivot.indices[0] = i;
                pivot.indices[1] = j;
            }
        }
    }
    return pivot;
}

template<typename F>
ValueIntPair<Base<F>> MaxAbs( const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("MaxAbs");
        if( !A.Grid().InGrid() )
            LogicError("Viewing processes are not allowed");
    )
    typedef Base<F> Real;
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

template<typename F>
ValueIntPair<Base<F>> SymmetricMaxAbs( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("SymmetricMaxAbs");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    typedef Base<F> Real;
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
    if( uplo == LOWER )
    {
        for( Int j=0; j<n; ++j )
        {
            for( Int i=j; i<n; ++i )
            {
                const Real abs = Abs(A.Get(i,j));
                if( abs > pivot.value )
                {
                    pivot.value = abs;
                    pivot.indices[0] = i;
                    pivot.indices[1] = j;
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
                const Real abs = Abs(A.Get(i,j));
                if( abs > pivot.value )
                {
                    pivot.value = abs;
                    pivot.indices[0] = i;
                    pivot.indices[1] = j;
                }
            }
        }
    }
    return pivot;
}

template<typename F>
ValueIntPair<Base<F>> SymmetricMaxAbs
( UpperOrLower uplo, const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("SymmetricMaxAbs");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( !A.Grid().InGrid() )
            LogicError("Viewing processes are not allowed");
    )
    typedef Base<F> Real;
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

        if( uplo == LOWER )
        {
            for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int mLocBefore = A.LocalRowOffset(j);
                for( Int iLoc=mLocBefore; iLoc<mLocal; ++iLoc )
                {
                    const Real abs = Abs(A.GetLocal(iLoc,jLoc));
                    if( abs > localPivot.value )
                    {
                        const Int i = A.GlobalRow(iLoc);
                        localPivot.value = abs;
                        localPivot.indices[0] = i;
                        localPivot.indices[1] = j;
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
                    const Real abs = Abs(A.GetLocal(iLoc,jLoc));
                    if( abs > localPivot.value )
                    {
                        const Int i = A.GlobalRow(iLoc);
                        localPivot.value = abs;
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

template<typename F>
ValueInt<Base<F>> DiagonalMaxAbs( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalMaxAbs"))
    return VectorMaxAbs( A.GetDiagonal() );
}

template<typename F,Dist U,Dist V>
ValueInt<Base<F>> DiagonalMaxAbs( const DistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalMaxAbs"))
    return VectorMaxAbs( A.GetDiagonal() );
}

#define DIST_PROTO(F,U,V) \
  template ValueInt<Base<F>> DiagonalMaxAbs( const DistMatrix<F,U,V>& A );

#define PROTO(F) \
  template ValueInt<Base<F>> VectorMaxAbs( const Matrix<F>& x ); \
  template ValueInt<Base<F>> VectorMaxAbs( const AbstractDistMatrix<F>& x ); \
  template ValueIntPair<Base<F>> MaxAbs( const Matrix<F>& x ); \
  template ValueIntPair<Base<F>> MaxAbs( const AbstractDistMatrix<F>& x ); \
  template ValueIntPair<Base<F>> SymmetricMaxAbs \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template ValueIntPair<Base<F>> SymmetricMaxAbs \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A ); \
  template ValueInt<Base<F>> DiagonalMaxAbs( const Matrix<F>& A ); \
  DIST_PROTO(F,CIRC,CIRC); \
  DIST_PROTO(F,MC,  MR  ); \
  DIST_PROTO(F,MC,  STAR); \
  DIST_PROTO(F,MD,  STAR); \
  DIST_PROTO(F,MR,  MC  ); \
  DIST_PROTO(F,MR,  STAR); \
  DIST_PROTO(F,STAR,MC  ); \
  DIST_PROTO(F,STAR,MD  ); \
  DIST_PROTO(F,STAR,MR  ); \
  DIST_PROTO(F,STAR,STAR); \
  DIST_PROTO(F,STAR,VC  ); \
  DIST_PROTO(F,STAR,VR  ); \
  DIST_PROTO(F,VC,  STAR); \
  DIST_PROTO(F,VR,  STAR); 

PROTO(Int);
PROTO(float);
PROTO(double);
PROTO(Complex<float>);
PROTO(Complex<double>);

} // namespace El
