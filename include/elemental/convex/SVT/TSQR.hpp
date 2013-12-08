/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CONVEX_SVT_TSQR_HPP
#define ELEM_CONVEX_SVT_TSQR_HPP

#include "elemental/lapack-like/Norm/Zero.hpp"
#include "elemental/lapack-like/QR.hpp"

namespace elem {
namespace svt {

// Singular-value soft-thresholding based on TSQR

template<typename F,Distribution U>
inline Int
TSQR( DistMatrix<F,U,STAR>& A, BASE(F) tau, bool relative=false )
{
#ifndef RELEASE
    CallStackEntry cse("SVT");
#endif
    const Int p = mpi::CommSize( A.ColComm() );
    if( p == 1 )
        return SVT( A.Matrix(), tau, relative );

    Int zeroNorm;
    qr::TreeData<F> treeData;
    treeData.QR0 = A.LockedMatrix();
    QR( treeData.QR0, treeData.t0 );
    qr::ts::Reduce( A, treeData );
    if( A.ColRank() == 0 )
        zeroNorm = SVT( qr::ts::RootQR(A,treeData), tau, relative );
    qr::ts::Scatter( A, treeData );

    mpi::Broadcast( zeroNorm, 0, A.ColComm() );
    return zeroNorm;
}

} // namespace svt
} // namespace elem

#endif // ifndef ELEM_CONVEX_SVT_TSQR_HPP
