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

template<typename T>
inline void
TrtrmmLUnblocked( Orientation orientation, Matrix<T>& L )
{
#ifndef RELEASE
    PushCallStack("internal::TrtrmmLUnblocked");
    if( L.Height() != L.Width() )
        throw std::logic_error("L must be square");
    if( orientation == NORMAL )
        throw std::logic_error("Trtrmm requires (conjugate-)transpose");
#endif
    const int n = L.Height();

    T* LBuffer = L.Buffer();
    const int ldim = L.LDim();
    for( int j=0; j<n; ++j )
    {
        T* RESTRICT l10 = &LBuffer[j];
        if( orientation == ADJOINT )
        {
            // L00 := L00 + l10^H l10
            for( int k=0; k<j; ++k )
            {
                const T gamma = l10[k*ldim];
                T* RESTRICT L00Col = &LBuffer[k*ldim];
                for( int i=k; i<j; ++i )
                    L00Col[i] += Conj(l10[i*ldim])*gamma;
            }
        }
        else
        {
            // L00 := L00 + l10^T l10
            for( int k=0; k<j; ++k )
            {
                const T gamma = l10[k*ldim];
                T* RESTRICT L00Col = &LBuffer[k*ldim];
                for( int i=k; i<j; ++i )
                    L00Col[i] += l10[i*ldim]*gamma;
            }
        }

        // l10 := l10 lambda11
        const T lambda11 = LBuffer[j+j*ldim];
        for( int k=0; k<j; ++k )
            l10[k*ldim] *= lambda11;

        // lambda11 := lambda11^2 or |lambda11|^2
        if( orientation == ADJOINT )
            LBuffer[j+j*ldim] = lambda11*Conj(lambda11);
        else
            LBuffer[j+j*ldim] = lambda11*lambda11;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
TrtrmmUUnblocked( Orientation orientation, Matrix<T>& U )
{
#ifndef RELEASE
    PushCallStack("internal::TrtrmmUUnblocked");
    if( U.Height() != U.Width() )
        throw std::logic_error("U must be square");
    if( orientation == NORMAL )
        throw std::logic_error("Trtrmm requires (conjugate-)transpose");
#endif
    const int n = U.Height();

    T* UBuffer = U.Buffer();
    const int ldim = U.LDim();
    for( int j=0; j<n; ++j )
    {
        T* RESTRICT u01 = &UBuffer[j*ldim];
        if( orientation == ADJOINT )
        {
            // U00 := U00 + u01 u01^H
            for( int k=0; k<j; ++k )
            {
                const T gamma = Conj(u01[k]);
                T* RESTRICT U00Col = &UBuffer[k*ldim];
                for( int i=0; i<k; ++i )
                    U00Col[i] += u01[i]*gamma;
            }
        }
        else
        {
            // U00 := U00 + u01 u01^T
            for( int k=0; k<j; ++k )
            {
                const T gamma = u01[k];
                T* RESTRICT U00Col = &UBuffer[k*ldim];
                for( int i=0; i<k; ++i )
                    U00Col[i] += u01[i]*gamma;
            }
        }

        // u01 := u01 upsilon11
        const T upsilon11 = UBuffer[j+j*ldim];
        for( int k=0; k<j; ++k )
            u01[k] *= upsilon11;

        // upsilon11 := upsilon11^2 or |upsilon11|^2
        if( orientation == ADJOINT )
            UBuffer[j+j*ldim] = upsilon11*Conj(upsilon11);
        else
            UBuffer[j+j*ldim] = upsilon11*upsilon11;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
