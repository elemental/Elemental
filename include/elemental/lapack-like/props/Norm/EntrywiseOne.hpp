/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_NORM_ENTRYWISEONE_HPP
#define ELEM_NORM_ENTRYWISEONE_HPP

namespace elem {

template<typename F> 
inline BASE(F)
EntrywiseOneNorm( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("EntrywiseOneNorm"))
    typedef Base<F> R;
    R norm = 0;
    const Int width = A.Width();
    const Int height = A.Height();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            norm += Abs(A.Get(i,j));
    return norm;
}

template<typename F>
inline BASE(F)
HermitianEntrywiseOneNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEntrywiseOneNorm"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    typedef Base<F> R;
    R norm = 0;
    const Int height = A.Height();
    const Int width = A.Width();
    if( uplo == UPPER )
    {
        for( Int j=0; j<width; ++j )
        {
            for( Int i=0; i<j; ++i )
            {
                const R alpha = Abs(A.Get(i,j));
                if( i ==j )
                    norm += alpha;
                else
                    norm += 2*alpha;
            }
        }
    }
    else
    {
        for( Int j=0; j<width; ++j )
        {
            for( Int i=j+1; i<height; ++i )
            {
                const R alpha = Abs(A.Get(i,j));
                if( i ==j )
                    norm += alpha;
                else
                    norm += 2*alpha;
            }
        }
    }
    return norm;
}

template<typename F>
inline BASE(F)
SymmetricEntrywiseOneNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricEntrywiseOneNorm"))
    return HermitianEntrywiseOneNorm( uplo, A );
}

template<typename F,Dist U,Dist V> 
inline BASE(F)
EntrywiseOneNorm( const DistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("EntrywiseOneNorm"))
    typedef Base<F> Real;
    Real norm;
    if( A.Participating() )
    {
        Real localSum = 0;
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                localSum += Abs(A.GetLocal(iLoc,jLoc)); 
        norm = mpi::AllReduce( localSum, A.DistComm() );
    }
    mpi::Broadcast( norm, A.Root(), A.CrossComm() );
    return norm;
}

template<typename F>
inline BASE(F)
HermitianEntrywiseOneNorm( UpperOrLower uplo, const DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEntrywiseOneNorm"))
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
                const R alpha = Abs(A.GetLocal(iLoc,jLoc));
                if( i ==j )
                    localSum += alpha;
                else
                    localSum += 2*alpha;
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
                const R alpha = Abs(A.GetLocal(iLoc,jLoc));
                if( i ==j )
                    localSum += alpha;
                else
                    localSum += 2*alpha;
            }
        }
    }

    return mpi::AllReduce( localSum, A.Grid().VCComm() );
}

template<typename F,Dist U,Dist V>
inline BASE(F)
SymmetricEntrywiseOneNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricEntrywiseOneNorm"))
    return HermitianEntrywiseOneNorm( uplo, A );
}

} // namespace elem

#endif // ifndef ELEM_NORM_ENTRYWISEONE_HPP
