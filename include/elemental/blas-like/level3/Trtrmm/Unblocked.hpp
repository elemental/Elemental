/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
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
