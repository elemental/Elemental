/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

Int SOCDegree( const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCDegree"))
    const Int height = firstInds.Height();
    Int degree = 0;
    for( Int i=0; i<height; ++i )
        if( i == firstInds.Get(i,0) )
            ++degree;
    return degree;
}

Int SOCDegree( const AbstractDistMatrix<Int>& firstIndsPre )
{
    DEBUG_ONLY(CSE cse("SOCDegree"))
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre);
    auto& firstInds = *firstIndsPtr;

    Int localDegree = 0;
    const Int localHeight = firstInds.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = firstInds.GlobalRow(iLoc);
        if( i == firstInds.GetLocal(iLoc,0) )
            ++localDegree;
    }
    return mpi::AllReduce( localDegree, firstInds.DistComm() );
}

Int SOCDegree( const DistMultiVec<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCDegree"))
    Int localDegree = 0;
    const Int localHeight = firstInds.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = firstInds.GlobalRow(iLoc);
        if( i == firstInds.GetLocal(iLoc,0) )
            ++localDegree;
    }
    return mpi::AllReduce( localDegree, firstInds.Comm() );
}

} // namespace El
