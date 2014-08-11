/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CHOLESKY_LVAR3SQUARE_HPP
#define EL_CHOLESKY_LVAR3SQUARE_HPP

// TODO: Reverse version

namespace El {
namespace cholesky {

template<typename F>
inline void
LVar3Square( AbstractDistMatrix<F>& APre )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::LVar3Square");
        if( APre.Height() != APre.Width() )
            LogicError("Can only compute Cholesky factor of square matrices");
        if( APre.Grid().Height() != APre.Grid().Width() )
            LogicError("CholeskyLVar3Square requires a square process grid");
    )
    const Grid& g = APre.Grid();

    DistMatrix<F> A(g);
    Copy( APre, A, READ_WRITE_PROXY );

    // Find the process holding our transposed data
    const Int transposeRank = 
        A.RowOwner(A.RowShift()) + A.ColStride()*A.ColOwner(A.ColShift());
    const bool onDiagonal = ( transposeRank == g.VCRank() );

    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,STAR,MC  > A21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A21Adj_STAR_MR(g);

    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const IndexRange ind1( k,    k+nb );
        const IndexRange ind2( k+nb, n    );

        auto A11 = View( A, ind1, ind1 );
        auto A21 = View( A, ind2, ind1 );
        auto A22 = View( A, ind2, ind2 );

        A11_STAR_STAR = A11;
        LocalCholesky( LOWER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A21_VC_STAR.AlignWith( A22 );
        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), A11_STAR_STAR, A21_VC_STAR );

        A21Trans_STAR_MC.AlignWith( A22 );
        A21_VC_STAR.TransposePartialColAllGather( A21Trans_STAR_MC );
        // SendRecv to form A21^T[* ,MR] from A21^T[* ,MC], then conjugate
        // the buffer to form A21^H[* ,MR]
        A21Adj_STAR_MR.AlignWith( A22 );
        A21Adj_STAR_MR.Resize( A21.Width(), A21.Height() ); 
        {
            if( onDiagonal )
            { 
                const Int size = A11.Height()*A22.LocalWidth();
                MemCopy
                ( A21Adj_STAR_MR.Buffer(), 
                  A21Trans_STAR_MC.Buffer(), size );
            }
            else
            {
                const Int sendSize = A21.LocalHeight()*A21.Width();
                const Int recvSize = A11.Height()*A22.LocalWidth();
                // We know that the ldim is the height since we have manually 
                // created both temporary matrices.
                mpi::SendRecv 
                ( A21Trans_STAR_MC.Buffer(), sendSize, transposeRank,
                  A21Adj_STAR_MR.Buffer(),   recvSize, transposeRank,
                  g.VCComm() );
            }
            Conjugate( A21Adj_STAR_MR );
        }

        // (A21^T[* ,MC])^T A21^H[* ,MR] = A21[MC,* ] A21^H[* ,MR]
        //                               = (A21 A21^H)[MC,MR]
        LocalTrrk
        ( LOWER, TRANSPOSE, 
          F(-1), A21Trans_STAR_MC, A21Adj_STAR_MR, F(1), A22 );

        A21.TransposeRowFilterFrom( A21Trans_STAR_MC );
    }
    Copy( A, APre, RESTORE_READ_WRITE_PROXY );
}

} // namespace cholesky
} // namespace El

#endif // ifndef EL_CHOLESKY_LVAR3SQUARE_HPP
