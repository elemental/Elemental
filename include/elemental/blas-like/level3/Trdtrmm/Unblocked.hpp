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

template<typename F>
inline void
TrdtrmmLUnblocked( Orientation orientation, Matrix<F>& L )
{
#ifndef RELEASE
    PushCallStack("internal::TrdtrmmLUnblocked");
    if( L.Height() != L.Width() )
        throw std::logic_error("L must be square");
    if( orientation == NORMAL )
        throw std::logic_error("Trdtrmm requires (conjugate-)transpose");
#endif
    const int n = L.Height();

    F* LBuffer = L.Buffer();
    const int ldim = L.LDim();
    for( int j=0; j<n; ++j )
    {
        const F delta11 = LBuffer[j+j*ldim];
        if( delta11 == (F)0 )
            throw SingularMatrixException();

        F* RESTRICT l10 = &LBuffer[j];
        if( orientation == ADJOINT )
        {
            // L00 := L00 + l10^H (l10 / delta11)
            for( int k=0; k<j; ++k )
            {
                const F gamma = l10[k*ldim] / delta11; 
                F* RESTRICT L00Col = &LBuffer[k*ldim];
                for( int i=k; i<j; ++i )
                    L00Col[i] += Conj(l10[i*ldim])*gamma;
            }
        }
        else
        {
            // L00 := L00 + l10^T (l10 / delta11)
            for( int k=0; k<j; ++k )
            {
                const F gamma = l10[k*ldim] / delta11;
                F* RESTRICT L00Col = &LBuffer[k*ldim];
                for( int i=k; i<j; ++i )
                    L00Col[i] += l10[i*ldim]*gamma;
            }
        }

        // l10 := l10 / delta11
        for( int k=0; k<j; ++k )
            l10[k*ldim] /= delta11;

        // lambda11 := 1 / delta11
        LBuffer[j+j*ldim] = 1 / delta11;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
TrdtrmmUUnblocked( Orientation orientation, Matrix<F>& U )
{
#ifndef RELEASE
    PushCallStack("internal::TrdtrmmUUnblocked");
    if( U.Height() != U.Width() )
        throw std::logic_error("U must be square");
    if( orientation == NORMAL )
        throw std::logic_error("Trdtrmm requires (conjugate-)transpose");
#endif
    const int n = U.Height();

    F* UBuffer = U.Buffer();
    const int ldim = U.LDim();
    for( int j=0; j<n; ++j )
    {
        const F delta11 = UBuffer[j+j*ldim];
        if( delta11 == (F)0 )
            throw SingularMatrixException();

        F* RESTRICT u01 = &UBuffer[j*ldim];
        if( orientation == ADJOINT )
        {
            // U00 := U00 + u01 (u01 / conj(delta11))^H
            for( int k=0; k<j; ++k )
            {
                const F gamma = Conj(u01[k]) / delta11;
                F* RESTRICT U00Col = &UBuffer[k*ldim];
                for( int i=0; i<k; ++i )
                    U00Col[i] += u01[i]*gamma;
            }
        }
        else
        {
            // U00 := U00 + u01 (u01 / delta11)^T
            for( int k=0; k<j; ++k )
            {
                const F gamma = u01[k] / delta11;
                F* RESTRICT U00Col = &UBuffer[k*ldim];
                for( int i=0; i<k; ++i )
                    U00Col[i] += u01[i]*gamma;
            }
        }

        // u01 := u01 / delta11
        for( int k=0; k<j; ++k )
            u01[k] /= delta11;

        // lambda11 := 1 / delta11
        UBuffer[j+j*ldim] = 1 / delta11;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
