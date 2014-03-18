/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_KMS_HPP
#define ELEM_KMS_HPP

namespace elem {

template<typename T> 
inline void
KMS( Matrix<T>& K, Int n, T rho )
{
    DEBUG_ONLY(CallStackEntry cse("KMS"))
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<j; ++i )
            K.Set( i, j, Pow(rho,T(j-i)) );
        for( Int i=j; i<n; ++i )
            K.Set( i, j, Conj(Pow(rho,T(i-j))) );
    }
}

#ifndef SWIG
template<typename T> 
inline Matrix<T>
KMS( Int n, T rho )
{
    Matrix<T> K;
    KMS( K, n, rho );
    return K;
}
#endif

template<typename T,Dist U,Dist V>
inline void
KMS( DistMatrix<T,U,V>& K, Int n, T rho )
{
    DEBUG_ONLY(CallStackEntry cse("KMS"))
    const Int localHeight = K.LocalHeight();
    const Int localWidth = K.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = K.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = K.GlobalRow(iLoc);
            if( i < j )
                K.SetLocal( iLoc, jLoc, Pow(rho,T(j-i)) );
            else
                K.SetLocal( iLoc, jLoc, Conj(Pow(rho,T(i-j))) );
        }
    }
}

#ifndef SWIG
template<typename T,Dist U=MC,Dist V=MR>
inline DistMatrix<T,U,V>
KMS( const Grid& g, Int n, T rho )
{
    DistMatrix<T,U,V> K(g);
    KMS( K, n, rho );
    return K;
}
#endif

} // namespace elem

#endif // ifndef ELEM_KMS_HPP
