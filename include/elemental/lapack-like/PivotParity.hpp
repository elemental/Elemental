/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

// A permutation is even if and only if it is the product of an even number 
// of transpositions, so we can decide this by simply checking how many 
// nontrivial pivots were performed.

inline bool
PivotParity( const Matrix<int>& p, int pivotOffset )
{
#ifndef RELEASE
    PushCallStack("PivotParity");
    if( p.Width() != 1 )
        throw std::logic_error("p must be a column vector");
    if( pivotOffset < 0 )
        throw std::logic_error("pivot offset cannot be negative");
#endif
    const int n = p.Height();
    bool isOdd = false;
    for( int i=0; i<n; ++i )
        if( p.Get(i,0) != i+pivotOffset )
            isOdd = !isOdd;
#ifndef RELEASE
    PopCallStack();
#endif
    return isOdd;
}

inline bool
PivotParity( const DistMatrix<int,VC,STAR>& p, int pivotOffset ) 
{
#ifndef RELEASE
    PushCallStack("PivotParity");
    if( p.Width() != 1 )
        throw std::logic_error("p must be a column vector");
    if( pivotOffset < 0 )
        throw std::logic_error("pivot offset cannot be negative");
#endif
    const int mLocal = p.LocalHeight();
    const Grid& g = p.Grid();

    bool isLocallyOdd = false;
    const int colShift = p.ColShift();
    const int gridSize = g.Size();
    for( int iLocal=0; iLocal<mLocal; ++iLocal )
    {
        const int i = colShift + iLocal*gridSize;
        if( p.GetLocal(iLocal,0) != i+pivotOffset )
            isLocallyOdd = !isLocallyOdd;
    }

    int localContribution = isLocallyOdd;
    int globalContribution;
    mpi::AllReduce
    ( &localContribution, &globalContribution, 1, MPI_SUM, g.VCComm() );
    globalContribution = globalContribution % 2;

    bool isOdd = globalContribution;
#ifndef RELEASE
    PopCallStack();
#endif
    return isOdd;
}

} // namespace elem
