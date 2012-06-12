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

#include "./ExpandPackedReflectors/LV.hpp"

namespace elem {

template<typename R> 
inline void
ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir,
  int offset, Matrix<R>& H )
{
#ifndef RELEASE
    PushCallStack("ExpandPackedReflectors");
#endif
    // Since the complex version does not have the same argument list, there is
    // currently no good way to ensure that this version is not called with 
    // complex datatypes. Until C++11 compilers are commonplace, we cannot
    // use static_assert either.
    if( IsComplex<R>::val )
        throw std::logic_error("Called real routine with complex datatype");

    if( uplo == LOWER && dir == VERTICAL )
        internal::ExpandPackedReflectorsLV( offset, H );
    else
        throw std::logic_error("This option is not yet supported");
#ifndef RELEASE
    PopCallStack();
#endif
}

// None of the underlying routines are implemented yet
/*
template<typename R> 
inline void
ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir, 
  int offset, DistMatrix<R>& H )
{
#ifndef RELEASE
    PushCallStack("ExpandPackedReflectors");
#endif
    // Since the complex version does not have the same argument list, there is
    // currently no good way to ensure that this version is not called with 
    // complex datatypes. Until C++11 compilers are commonplace, we cannot
    // use static_assert either.
    if( IsComplex<R>::val )
        throw std::logic_error("Called real routine with complex datatype");

    if( uplo == LOWER && dir == VERTICAL )
        internal::ExpandPackedReflectorsLV( offset, H );
    else
        throw std::logic_error("This option is not yet supported");
#ifndef RELEASE
    PopCallStack();
#endif
}
*/

template<typename R> 
inline void
ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation,
  int offset, Matrix<Complex<R> >& H, const Matrix<Complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("ExpandPackedReflectors");
#endif
    if( uplo == LOWER && dir == VERTICAL )
        internal::ExpandPackedReflectorsLV( conjugation, offset, H, t );
    else
        throw std::logic_error("This option is not yet supported");
#ifndef RELEASE
    PopCallStack();
#endif
}

// None of the underlying routines are implemented yet
/*
template<typename R> 
inline void
ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation,
  int offset, 
        DistMatrix<Complex<R> >& H, 
  const DistMatrix<Complex<R>,MD,STAR>& t )
{
#ifndef RELEASE
    PushCallStack("ExpandPackedReflectors");
#endif
    if( uplo == LOWER && dir == VERTICAL )
        internal::ExpandPackedReflectorsLV( conjugation, offset, H, t );
    else
        throw std::logic_error("This option is not yet supported");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation,
  int offset,
        DistMatrix<Complex<R> >& H, 
  const DistMatrix<Complex<R>,STAR,STAR>& t )
{
#ifndef RELEASE
    PushCallStack("ExpandPackedReflectors");
#endif
    DistMatrix<Complex<R>,MD,STAR> tDiag(A.Grid());
    tDiag.AlignWithDiagonal( A, offset );
    tDiag = t;
    ExpandPackedReflectors( uplo, dir, conjugation, offset, H, tDiag );
#ifndef RELEASE
    PopCallStack();
#endif
}
*/

} // namespace elem
