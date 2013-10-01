/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_NORM_ENTRYWISE_HPP
#define ELEM_LAPACK_NORM_ENTRYWISE_HPP

namespace elem {

template<typename F> 
inline Base<F>
EntrywiseNorm( const Matrix<F>& A, Base<F> p )
{
#ifndef RELEASE
    CallStackEntry entry("EntrywiseNorm");
#endif
    // TODO: Make this more numerically stable
    typedef Base<F> R;
    R sum = 0;
    const Int width = A.Width();
    const Int height = A.Height();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            sum += Pow( Abs(A.Get(i,j)), p );
    return Pow( sum, 1/p );
}

template<typename F>
inline Base<F>
HermitianEntrywiseNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEntrywiseNorm");
#endif
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    // TODO: make this more numerically stable
    typedef Base<F> R;
    R sum = 0;
    const Int height = A.Height();
    const Int width = A.Width();
    if( uplo == UPPER )
    {
        for( Int j=0; j<width; ++j )
        {
            for( Int i=0; i<j; ++i )
            {
                const R term = Pow( Abs(A.Get(i,j)), p );
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
                const R term = Pow( Abs(A.Get(i,j)), p );
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
inline Base<F>
SymmetricEntrywiseNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p )
{
#ifndef RELEASE
    CallStackEntry entry("SymmetricEntrywiseNorm");
#endif
    return HermitianEntrywiseNorm( uplo, A, p );
}

template<typename F,Distribution U,Distribution V> 
inline Base<F>
EntrywiseNorm( const DistMatrix<F,U,V>& A, Base<F> p )
{
#ifndef RELEASE
    CallStackEntry entry("EntrywiseNorm");
#endif
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
inline Base<F>
HermitianEntrywiseNorm( UpperOrLower uplo, const DistMatrix<F>& A, Base<F> p )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEntrywiseNorm");
#endif
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    const Int r = A.Grid().Height();
    const Int c = A.Grid().Width();
    const Int colShift = A.ColShift();
    const Int rowShift = A.RowShift();

    typedef Base<F> R;
    R localSum = 0;
    const Int localWidth = A.LocalWidth();
    if( uplo == UPPER )
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            Int j = rowShift + jLoc*c;
            Int numUpperRows = Length(j+1,colShift,r);
            for( Int iLoc=0; iLoc<numUpperRows; ++iLoc )
            {
                Int i = colShift + iLoc*r;
                const R term = Pow( Abs(A.GetLocal(iLoc,jLoc)), p );
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
            Int j = rowShift + jLoc*c;
            Int numStrictlyUpperRows = Length(j,colShift,r);
            for( Int iLoc=numStrictlyUpperRows;
                 iLoc<A.LocalHeight(); ++iLoc )
            {
                Int i = colShift + iLoc*r;
                const R term = Pow( Abs(A.GetLocal(iLoc,jLoc)), p );
                if( i ==j )
                    localSum += term;
                else
                    localSum += 2*term;
            }
        }
    }

    const R sum = mpi::AllReduce( localSum, A.Grid().VCComm() );
    return Pow( sum, 1/p );
}

template<typename F,Distribution U,Distribution V>
inline Base<F>
SymmetricEntrywiseNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Base<F> p )
{
#ifndef RELEASE
    CallStackEntry entry("SymmetricEntrywiseNorm");
#endif
    return HermitianEntrywiseNorm( uplo, A, p );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_NORM_ENTRYWISE_HPP
