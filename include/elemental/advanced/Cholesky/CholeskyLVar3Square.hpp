/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

namespace elem {

/*
   Parallelization of Variant 3 Lower Cholesky factorization for
   square process grids.

   Original serial update:
   ------------------------
   A11 := Cholesky(A11) 
   A21 := A21 tril(A11)^-H
   A22 := A22 - A21 A21^H
   ------------------------

   Corresponding parallel update:
   -----------------------------------------------------
   A11[* ,* ] <- A11[MC,MR] 
   A11[* ,* ] := Cholesky(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]
   
   A21[VC,* ] <- A21[MC,MR]
   A21[VC,* ] := A21[VC,* ] tril(A11[* ,* ])^-H
   
   A21^T[* ,MC] <- A21[VC,* ]
   A21^H[* ,MR] <- A21^T[* ,MC]
   A22[MC,MR] := A22[MC,MR] - (A21^T[* ,MC])^T A21^H[* ,MR]
   A21[MC,MR] <- A21^T[* ,MC]
   -----------------------------------------------------
*/
template<typename F>
inline void
internal::CholeskyLVar3Square( DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("internal::CholeskyLVar3Square");
    if( A.Height() != A.Width() )
        throw std::logic_error
        ("Can only compute Cholesky factor of square matrices");
    if( A.Grid().Height() != A.Grid().Width() )
        throw std::logic_error
        ("CholeskyLVar3Square requires a square process grid");
#endif
    const Grid& g = A.Grid();

    // Find the process holding our transposed data
    const int r = g.Height();
    int transposeRank;
    {
        const int colAlignment = A.ColAlignment();
        const int rowAlignment = A.RowAlignment();
        const int colShift = A.ColShift();
        const int rowShift = A.RowShift();

        const int transposeRow = (colAlignment+rowShift) % r;
        const int transposeCol = (rowAlignment+colShift) % r;
        transposeRank = transposeRow + r*transposeCol;
    }
    const bool onDiagonal = ( transposeRank == g.VCRank() );

    // Matrix views
    DistMatrix<F,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary matrices
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,STAR,MC  > A21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A21Adj_STAR_MR(g);

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ABR.Height() > 0 )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/   
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        A21_VC_STAR.AlignWith( A22 );
        A21Trans_STAR_MC.AlignWith( A22 );
        A21Adj_STAR_MR.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_STAR_STAR = A11;
        internal::LocalCholesky( LOWER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A21_VC_STAR = A21;
        internal::LocalTrsm
        ( RIGHT, LOWER, ADJOINT, NON_UNIT, (F)1, A11_STAR_STAR, A21_VC_STAR );

        A21Trans_STAR_MC.TransposeFrom( A21_VC_STAR );
        // SendRecv to form A21^T[* ,MR] from A21^T[* ,MC], then conjugate
        // the buffer to form A21^H[* ,MR]
        A21Adj_STAR_MR.ResizeTo( A21.Width(), A21.Height() ); 
        {
            if( onDiagonal )
            { 
                const int size = A11.Height()*A22.LocalWidth();
                std::memcpy
                ( A21Adj_STAR_MR.LocalBuffer(), A21Trans_STAR_MC.LocalBuffer(),
                  size*sizeof(F) );
            }
            else
            {
                const int sendSize = A22.LocalHeight()*A11.Width();
                const int recvSize = A22.LocalWidth()*A11.Height();
                // We know that the ldim is the height since we have manually 
                // created both temporary matrices.
                mpi::SendRecv 
                ( A21Trans_STAR_MC.LocalBuffer(), sendSize, transposeRank, 0,
                  A21Adj_STAR_MR.LocalBuffer(),  recvSize, transposeRank, 0,
                  g.VCComm() );
            }
            Conjugate( A21Adj_STAR_MR );
        }

        // (A21^T[* ,MC])^T A21^H[* ,MR] = A21[MC,* ] A21^H[* ,MR]
        //                               = (A21 A21^H)[MC,MR]
        internal::LocalTrrk
        ( LOWER, TRANSPOSE, 
          (F)-1, A21Trans_STAR_MC, A21Adj_STAR_MR, (F)1, A22 );

        A21.TransposeFrom( A21Trans_STAR_MC );
        //--------------------------------------------------------------------//
        A21_VC_STAR.FreeAlignments();
        A21Trans_STAR_MC.FreeAlignments();
        A21Adj_STAR_MR.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
} 

} // namespace elem
