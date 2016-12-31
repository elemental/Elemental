/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace soc {

// TODO(poulson): Lower-level access

Int Degree( const Matrix<Int>& firstInds )
{
    EL_DEBUG_CSE
    const Int height = firstInds.Height();
    Int degree = 0;
    for( Int i=0; i<height; ++i )
        if( i == firstInds(i) )
            ++degree;
    return degree;
}

Int Degree( const AbstractDistMatrix<Int>& firstIndsPre )
{
    EL_DEBUG_CSE

    DistMatrixReadProxy<Int,Int,VC,STAR> firstIndsProx( firstIndsPre );
    auto& firstInds = firstIndsProx.GetLocked();

    Int localDegree = 0;
    const Int localHeight = firstInds.LocalHeight();
    auto& firstIndsLoc = firstInds.LockedMatrix();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = firstInds.GlobalRow(iLoc);
        if( i == firstIndsLoc(iLoc) )
            ++localDegree;
    }
    return mpi::AllReduce( localDegree, firstInds.DistComm() );
}

Int Degree( const DistMultiVec<Int>& firstInds )
{
    EL_DEBUG_CSE
    Int localDegree = 0;
    const Int localHeight = firstInds.LocalHeight();
    auto& firstIndsLoc = firstInds.LockedMatrix();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = firstInds.GlobalRow(iLoc);
        if( i == firstIndsLoc(iLoc) )
            ++localDegree;
    }
    return mpi::AllReduce( localDegree, firstInds.Grid().Comm() );
}

} // namespace soc
} // namespace El
