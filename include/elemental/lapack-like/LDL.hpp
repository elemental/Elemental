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

#include "./LDL/Var3.hpp"

namespace elem {

template<typename F>
inline void
LDLH( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("LDLH");
#endif
    Matrix<F> d;
    LDLH( A, d );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void 
LDLH( DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("LDLH");
#endif
    DistMatrix<F,MC,STAR> d( A.Grid() );
    LDLH( A, d );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
LDLH( Matrix<F>& A, Matrix<F>& d )
{
#ifndef RELEASE
    PushCallStack("LDLH");
#endif
    internal::LDLVar3( ADJOINT, A, d );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void 
LDLH( DistMatrix<F>& A, DistMatrix<F,MC,STAR>& d )
{
#ifndef RELEASE
    PushCallStack("LDLH");
#endif
    internal::LDLVar3( ADJOINT, A, d );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
LDLT( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("LDLT");
#endif
    Matrix<F> d;
    LDLT( A, d );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void 
LDLT( DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("LDLT");
#endif
    DistMatrix<F,MC,STAR> d( A.Grid() );
    LDLT( A, d );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
LDLT( Matrix<F>& A, Matrix<F>& d )
{
#ifndef RELEASE
    PushCallStack("LDLT");
#endif
    internal::LDLVar3( TRANSPOSE, A, d );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void 
LDLT( DistMatrix<F>& A, DistMatrix<F,MC,STAR>& d )
{
#ifndef RELEASE
    PushCallStack("LDLT");
#endif
    internal::LDLVar3( TRANSPOSE, A, d );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
