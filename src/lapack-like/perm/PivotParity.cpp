/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

// A permutation is even if and only if it is the product of an even number 
// of transpositions, so we can decide this by simply checking how many 
// nontrivial pivots were performed.

bool PivotParity( const Matrix<Int>& p, Int pivotOffset )
{
    DEBUG_ONLY(
        CallStackEntry cse("PivotParity");
        if( p.Width() != 1 )
            LogicError("p must be a column vector");
        if( pivotOffset < 0 )
            LogicError("pivot offset cannot be negative");
    )
    const Int n = p.Height();
    bool isOdd = false;
    for( Int i=0; i<n; ++i )
        if( p.Get(i,0) != i+pivotOffset )
            isOdd = !isOdd;
    return isOdd;
}

bool PivotParity( const DistMatrix<Int,VC,STAR>& p, Int pivotOffset ) 
{
    DEBUG_ONLY(
        CallStackEntry cse("PivotParity");
        if( p.Width() != 1 )
            LogicError("p must be a column vector");
        if( pivotOffset < 0 )
            LogicError("pivot offset cannot be negative");
    )
    bool isLocallyOdd = false;
    const Int mLocal = p.LocalHeight();
    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
    {
        const Int i = p.GlobalRow(iLoc);
        if( p.GetLocal(iLoc,0) != i+pivotOffset )
            isLocallyOdd = !isLocallyOdd;
    }

    Int localContribution = isLocallyOdd;
    Int globalContribution;
    mpi::AllReduce
    ( &localContribution, &globalContribution, 1, MPI_SUM, p.ColComm() );
    return globalContribution % 2;
}

} // namespace El
