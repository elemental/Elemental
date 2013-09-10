/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_CHOLESKY_UVAR3SQUARE_HPP
#define ELEM_LAPACK_CHOLESKY_UVAR3SQUARE_HPP

#include "elemental/blas-like/level3/Herk.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"

// TODO: Reverse version

namespace elem {
namespace cholesky {

template<typename F>
inline void
UVar3Square( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("cholesky::UVar3Square");
    if( A.Height() != A.Width() )
        LogicError("Can only compute Cholesky factor of square matrices.");
    if( A.Grid().Height() != A.Grid().Width() )
        LogicError("CholeskyUVar3Square assumes a square process grid.");
#endif
    const Grid& g = A.Grid();

    // Find the process holding our transposed data
    const Int r = g.Height();
    Int transposeRank;
    {
        const Int colAlignment = A.ColAlignment();
        const Int rowAlignment = A.RowAlignment();
        const Int colShift = A.ColShift();
        const Int rowShift = A.RowShift();

        const Int transposeRow = (colAlignment+rowShift) % r;
        const Int transposeCol = (rowAlignment+colShift) % r;
        transposeRank = transposeRow + r*transposeCol;
    }
    const bool onDiagonal = ( transposeRank == g.VCRank() );

    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,STAR,MC  > A12_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A12_STAR_MR(g);

    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = std::min(bsize,n-k);
        auto A11 = View( A, k,    k,    nb,       nb       );
        auto A12 = View( A, k,    k+nb, nb,       n-(k+nb) );
        auto A22 = View( A, k+nb, k+nb, n-(k+nb), n-(k+nb) );

        A11_STAR_STAR = A11;
        LocalCholesky( UPPER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A12_STAR_VR.AlignWith( A22 );
        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A11_STAR_STAR, A12_STAR_VR );

        A12_STAR_MR.AlignWith( A22 );
        A12_STAR_MR = A12_STAR_VR;
        // SendRecv to form A12[* ,MC] from A12[* ,MR]
        A12_STAR_MC.AlignWith( A22 );
        A12_STAR_MC.ResizeTo( A12.Height(), A12.Width() );
        {
            if( onDiagonal )
            {
                const Int size = A11.Height()*A22.LocalWidth();
                MemCopy
                ( A12_STAR_MC.Buffer(), 
                  A12_STAR_MR.Buffer(), size );
            }
            else
            {
                const Int sendSize = A11.Height()*A22.LocalWidth();
                const Int recvSize = A11.Width()*A22.LocalHeight();
                // We know that the ldim is the height since we have manually
                // created both temporary matrices.
                mpi::SendRecv
                ( A12_STAR_MR.Buffer(), sendSize, transposeRank,
                  A12_STAR_MC.Buffer(), recvSize, transposeRank, g.VCComm() );
            }
        }
        LocalTrrk
        ( UPPER, ADJOINT, F(-1), A12_STAR_MC, A12_STAR_MR, F(1), A22 );
        A12 = A12_STAR_MR;
    }
}

} // namespace cholesky
} // namespace elem

#endif // ifndef ELEM_LAPACK_CHOLESKY_UVAR3SQUARE_HPP
