/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_KMS_HPP
#define ELEM_MATRICES_KMS_HPP

namespace elem {

template<typename T> 
inline void
KMS( Matrix<T>& K, Int n, T rho )
{
#ifndef RELEASE
    CallStackEntry cse("KMS");
#endif
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<j; ++i )
            K.Set( i, j, Pow(rho,T(j-i)) );
        for( Int i=j; i<n; ++i )
            K.Set( i, j, Conj(Pow(rho,T(i-j))) );
    }
}

template<typename T> 
inline Matrix<T>
KMS( Int n, T rho )
{
    Matrix<T> K;
    KMS( K, n, rho );
    return K;
}

template<typename T,Distribution U,Distribution V>
inline void
KMS( DistMatrix<T,U,V>& K, Int n, T rho )
{
#ifndef RELEASE
    CallStackEntry cse("KMS");
#endif
    const Int localHeight = K.LocalHeight();
    const Int localWidth = K.LocalWidth();
    const Int colShift = K.ColShift();
    const Int rowShift = K.RowShift();
    const Int colStride = K.ColStride();
    const Int rowStride = K.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
            if( i < j )
                K.SetLocal( iLoc, jLoc, Pow(rho,T(j-i)) );
            else
                K.SetLocal( iLoc, jLoc, Conj(Pow(rho,T(i-j))) );
        }
    }
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
KMS( const Grid& g, Int n, T rho )
{
    DistMatrix<T,U,V> K(g);
    KMS( K, n, rho );
    return K;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_KMS_HPP
