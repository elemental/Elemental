/*
   Copyright (c) 2009-2010, Jack Poulson
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
#include "elemental/blas_internal.hpp"
using namespace elemental;

template<typename T>
void
elemental::blas::Trmm
( Side side, 
  Shape shape, 
  Orientation orientation, 
  Diagonal diagonal,
  T alpha, 
  const DistMatrix<T,MC,MR>& A,
        DistMatrix<T,MC,MR>& X   )
{
#ifndef RELEASE
    PushCallStack("blas::Trmm");
#endif
    if( side == Left && shape == Lower )
    {
        if( orientation == Normal )
            blas::internal::TrmmLLN( diagonal, alpha, A, X );
        else
            blas::internal::TrmmLLT( orientation, diagonal, alpha, A, X );
    }
    else if( side == Left && shape == Upper )
    {
        if( orientation == Normal )
            blas::internal::TrmmLUN( diagonal, alpha, A, X );
        else
            blas::internal::TrmmLUT( orientation, diagonal, alpha, A, X );
    }
    else if( side == Right && shape == Lower )
    {
        if( orientation == Normal )
            blas::internal::TrmmRLN( diagonal, alpha, A, X );
        else
            blas::internal::TrmmRLT( orientation, diagonal, alpha, A, X );
    }
    else if( side == Right && shape == Upper )
    {
        if( orientation == Normal )
            blas::internal::TrmmRUN( diagonal, alpha, A, X );
        else
            blas::internal::TrmmRUT( orientation, diagonal, alpha, A, X );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::Trmm
( Side side, 
  Shape shape, 
  Orientation orientation, 
  Diagonal diagonal,
  float alpha, 
  const DistMatrix<float,MC,MR>& A,
        DistMatrix<float,MC,MR>& X );

template void elemental::blas::Trmm
( Side side, 
  Shape shape, 
  Orientation orientation, 
  Diagonal diagonal,
  double alpha, 
  const DistMatrix<double,MC,MR>& A,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::Trmm
( Side side, 
  Shape shape, 
  Orientation orientation, 
  Diagonal diagonal,
  scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& A,
        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::Trmm
( Side side, 
  Shape shape, 
  Orientation orientation, 
  Diagonal diagonal,
  dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& A,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

