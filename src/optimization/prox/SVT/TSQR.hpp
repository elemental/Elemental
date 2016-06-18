/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SVT_TSQR_HPP
#define EL_SVT_TSQR_HPP

namespace El {
namespace svt {

// Singular-value soft-thresholding based on TSQR

template<typename F>
Int TSQR( ElementalMatrix<F>& APre, Base<F> tau, bool relative )
{
    DEBUG_CSE

    DistMatrixReadWriteProxy<F,F,VC,STAR> AProx( APre );
    auto& A = AProx.Get();

    const Int p = mpi::Size( A.ColComm() );
    if( p == 1 )
        return SVT( A.Matrix(), tau, relative );

    Int zeroNorm;
    qr::TreeData<F> treeData;
    treeData.QR0 = A.LockedMatrix();
    QR( treeData.QR0, treeData.phase0, treeData.signature0 );
    qr::ts::Reduce( A, treeData );
    if( A.ColRank() == 0 )
        zeroNorm = SVT( qr::ts::RootQR(A,treeData), tau, relative );
    qr::ts::Scatter( A, treeData );

    mpi::Broadcast( zeroNorm, 0, A.ColComm() );
    return zeroNorm;
}

} // namespace svt
} // namespace El

#endif // ifndef EL_SVT_TSQR_HPP
