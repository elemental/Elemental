/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_MAX_HPP
#define ELEM_BLAS_MAX_HPP

namespace elem {

// TODO: Add options for FastAbs instead of Abs

template<typename F>
inline ValueInt<BASE(F)>
VectorMax( const Matrix<F>& x )
{
#ifndef RELEASE
    CallStackEntry cse("VectorMax");
#endif
    typedef BASE(F) Real;
    const Int m = x.Height();
    const Int n = x.Width();
#ifndef RELEASE
    if( m != 1 && n != 1 )
        LogicError("Input should have been a vector");
#endif

    ValueInt<Real> pivot;
    pivot.index = 0;
    pivot.value = 0;
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
inline ValueIntPair<BASE(F)>
Max( const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("Max");
#endif
    typedef BASE(F) Real;
    const Int m = A.Height();
    const Int n = A.Width();

    ValueIntPair<Real> pivot;
    pivot.value = 0;
    pivot.indices[0] = 0;
    pivot.indices[1] = 0;
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<n; ++i )
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

template<typename F,Distribution U,Distribution V>
inline ValueIntPair<BASE(F)>
Max( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry cse("Max");
#endif
    typedef BASE(F) Real;
    ValueIntPair<Real> pivot;
    if( A.Participating() )
    {
        // Store the index/value of the local pivot candidate
        ValueIntPair<Real> localPivot;
        localPivot.value = 0;
        localPivot.indices[0] = 0;
        localPivot.indices[1] = 0;
        const Int mLocal = A.LocalHeight();
        const Int nLocal = A.LocalWidth();
        const Int colShift = A.ColShift();
        const Int rowShift = A.RowShift();
        const Int colStride = A.ColStride();
        const Int rowStride = A.RowStride();
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            const Int j = rowShift + jLoc*rowStride;
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            {
                const Int i = colShift + iLoc*colStride;
                const Real value = Abs(A.GetLocal(iLoc,jLoc));
                if( value > localPivot.value )
                {
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
inline ValueIntPair<BASE(F)>
SymmetricMax( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("SymmetricMax");
#endif
    typedef BASE(F) Real;
    const Int m = A.Height();
    const Int n = A.Width();

    ValueIntPair<Real> pivot;
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

} // namespace elem

#endif // ifndef ELEM_BLAS_MAX_HPP
