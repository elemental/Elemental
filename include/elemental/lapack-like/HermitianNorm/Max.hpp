/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

template<typename F>
inline typename Base<F>::type
HermitianMaxNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::HermitianMaxNorm");
#endif
    typedef typename Base<F>::type R;

    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    R maxAbs = 0;
    const int height = A.Height();
    const int width = A.Width();
    if( uplo == UPPER )
    {
        for( int j=0; j<width; ++j )
        {
            for( int i=0; i<=j; ++i )
            {
                const R thisAbs = Abs(A.Get(i,j));
                maxAbs = std::max( maxAbs, thisAbs );
            }
        }
    }
    else
    {
        for( int j=0; j<width; ++j )
        {
            for( int i=j; i<height; ++i )
            {
                const R thisAbs = Abs(A.Get(i,j));
                maxAbs = std::max( maxAbs, thisAbs );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return maxAbs;
}

template<typename F>
inline typename Base<F>::type
HermitianMaxNorm( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::HermitianMaxNorm");
#endif
    typedef typename Base<F>::type R;

    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    const int r = A.Grid().Height();
    const int c = A.Grid().Width();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();

    R localMaxAbs = 0;
    const int localWidth = A.LocalWidth();
    if( uplo == UPPER )
    {
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int numUpperRows = LocalLength(j+1,colShift,r);
            for( int iLocal=0; iLocal<numUpperRows; ++iLocal )
            {
                const R thisAbs = Abs(A.GetLocal(iLocal,jLocal));
                localMaxAbs = std::max( localMaxAbs, thisAbs );
            }
        }
    }
    else
    {
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int numStrictlyUpperRows = LocalLength(j,colShift,r);
            for( int iLocal=numStrictlyUpperRows; 
                 iLocal<A.LocalHeight(); ++iLocal )
            {
                const R thisAbs = Abs(A.GetLocal(iLocal,jLocal));
                localMaxAbs = std::max( localMaxAbs, thisAbs );
            }
        }
    }

    R maxAbs;
    mpi::AllReduce( &localMaxAbs, &maxAbs, 1, mpi::MAX, A.Grid().VCComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return maxAbs;
}

} // namespace internal
} // namespace elem
