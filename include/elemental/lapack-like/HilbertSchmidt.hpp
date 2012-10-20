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

namespace internal {
template<typename F,Distribution U,Distribution V>
mpi::Comm NormComm( const DistMatrix<F,U,V>& A );
}

// TODO: Think about using a more stable accumulation algorithm?

template<typename F> 
inline F
HilbertSchmidt( const Matrix<F>& A, const Matrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("HilbertSchmidt");
#endif
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        throw std::logic_error("Matrices must be the same size");
    F innerProd(0);
    const int width = A.Width();
    const int height = A.Height();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            innerProd += Conj(A.Get(i,j)*B.Get(i,j));
#ifndef RELEASE
    PopCallStack();
#endif
    return innerProd;
}

template<typename F,Distribution U,Distribution V> 
inline F
HilbertSchmidt( const DistMatrix<F,U,V>& A, const DistMatrix<F,U,V>& B )
{
#ifndef RELEASE
    PushCallStack("HilbertSchmidt");
#endif
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        throw std::logic_error("Matrices must be the same size");
    if( A.Grid() != B.Grid() )
        throw std::logic_error("Grids must match");
    if( A.ColAlignment() != B.ColAlignment() || 
        A.RowAlignment() != B.RowAlignment() )
        throw std::logic_error("Matrices must be aligned");

    F localInnerProd(0);
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            localInnerProd += Conj(A.GetLocal(iLocal,jLocal))*
                                   B.GetLocal(iLocal,jLocal);

    F innerProd;
    mpi::Comm comm = internal::NormComm( A );
    mpi::AllReduce( &localInnerProd, &innerProd, 1, mpi::SUM, comm );
#ifndef RELEASE
    PopCallStack();
#endif
    return innerProd;
}

} // namespace elem
