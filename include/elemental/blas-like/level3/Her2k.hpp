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

#include "./Her2k/LC.hpp"
#include "./Her2k/LN.hpp"
#include "./Her2k/UC.hpp"
#include "./Her2k/UN.hpp"

namespace elem {

template<typename T>
inline void
Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Her2k");
    if( orientation == NORMAL )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() ||
            B.Height() != C.Height() ||B.Height() != C.Width() )
            throw std::logic_error("Nonconformal Her2k");
    }
    else if( orientation == ADJOINT )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() ||
            B.Width() != C.Height() || B.Width() != C.Width() )
            throw std::logic_error("Nonconformal Her2k");
    }
    else
        throw std::logic_error
        ("Her2k only accepts NORMAL and ADJOINT options");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const int k = ( orientation == NORMAL ? A.Width() : A.Height() );
    blas::Her2k
    ( uploChar, transChar, C.Height(), k,
      alpha, A.LockedBuffer(), A.LDim(),
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Her2k
( UpperOrLower uplo, 
  Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("Her2k");
    if( orientation == TRANSPOSE )
        throw std::logic_error("Her2k accepts NORMAL and ADJOINT options");
#endif
    if( uplo == LOWER && orientation == NORMAL )
        internal::Her2kLN( alpha, A, B, beta, C );
    else if( uplo == LOWER )
        internal::Her2kLC( alpha, A, B, beta, C );
    else if( orientation == NORMAL )
        internal::Her2kUN( alpha, A, B, beta, C );
    else
        internal::Her2kUC( alpha, A, B, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
