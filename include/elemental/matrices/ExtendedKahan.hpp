/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_EXTENDEDKAHAN_HPP
#define MATRICES_EXTENDEDKAHAN_HPP

// Generate a 3(2^k) x 3(2^k) Extended Kahan matrix, which has the form
// A = S R, where S = diag(1,zeta,...,zeta^(3 2^k - 1)), 
// 
//         | I -phi H_k    0     |
//     R = | 0    I      phi H_k |,
//         | 0    0        I     |
//
// 0 < mu << 1, and phi^2 + zeta^2 = 1.
//

namespace elem {

template<typename F>
inline void
ExtendedKahan( Matrix<F>& A, int k, BASE(F) phi, BASE(F) mu )
{
#ifndef RELEASE
    CallStackEntry entry("ExtendedKahan");
#endif
    const int n = 3*(1u<<k);
    A.ResizeTo( n, n );
    MakeExtendedKahan( A, phi, mu );
}

template<typename F,Distribution U,Distribution V>
inline void
ExtendedKahan( DistMatrix<F,U,V>& A, int k, BASE(F) phi, BASE(F) mu )
{
#ifndef RELEASE
    CallStackEntry entry("ExtendedKahan");
#endif
    const int n = 3*(1u<<k);
    A.ResizeTo( n, n );
    MakeExtendedKahan( A, phi, mu );
}

template<typename F> 
inline void
MakeExtendedKahan( Matrix<F>& A, BASE(F) phi, BASE(F) mu )
{
#ifndef RELEASE
    CallStackEntry entry("MakeExtendedKahan");
#endif
    typedef BASE(F) R;

    if( A.Height() != A.Width() )
        throw std::logic_error("Extended Kahan matrices must be square");
    const int n = A.Height();
    if( n % 3 != 0 )
        throw std::logic_error("Dimension must be an integer multiple of 3");
    const int l = n / 3;
    if( !l || (l & (l-1)) )
        throw std::logic_error("n/3 is not a power of two");
    int k=0;
    while( (1u<<k) < l )
        ++k;

    if( phi <= R(0) || phi >= R(1) )
        throw std::logic_error("phi must be in (0,1)");
    if( mu <= R(0) || mu >= R(1) )
        throw std::logic_error("mu must be in (0,1)");

    // Start by setting A to the identity, and then modify the necessary 
    // l x l blocks of its 3 x 3 partitioning.
    MakeIdentity( A );
    Matrix<F> ABlock;
    View( ABlock, A, 2*l, 2*l, l, l );
    Scale( mu, ABlock );
    View( ABlock, A, 0, l, l, l );
    MakeWalsh( ABlock, k );
    Scale( -phi, ABlock );
    View( ABlock, A, l, 2*l, l, l );
    MakeWalsh( ABlock, k );
    Scale( phi, ABlock );

    // Now scale R by S
    const R zeta = Sqrt(R(1)-phi*phi);
    Matrix<R> d( n, 1 );
    for( int i=0; i<n; ++i )
        d.Set( i, 0, Pow(zeta,i) );
    DiagonalScale( LEFT, NORMAL, d, A );
}

template<typename F,Distribution U,Distribution V>
inline void
MakeExtendedKahan( DistMatrix<F,U,V>& A, BASE(F) phi, BASE(F) mu )
{
#ifndef RELEASE
    CallStackEntry entry("MakeExtendedKahan");
#endif
    typedef BASE(F) R;

    if( A.Height() != A.Width() )
        throw std::logic_error("Extended Kahan matrices must be square");
    const int n = A.Height();
    if( n % 3 != 0 )
        throw std::logic_error("Dimension must be an integer multiple of 3");
    const int l = n / 3;
    if( !l || (l & (l-1)) )
        throw std::logic_error("n/3 is not a power of two");
    int k=0;
    while( (1u<<k) < l )
        ++k;

    if( phi <= R(0) || phi >= R(1) )
        throw std::logic_error("phi must be in (0,1)");
    if( mu <= R(0) || mu >= R(1) )
        throw std::logic_error("mu must be in (0,1)");

    // Start by setting A to the identity, and then modify the necessary 
    // l x l blocks of its 3 x 3 partitioning.
    MakeIdentity( A );
    DistMatrix<F,U,V> ABlock( A.Grid() );
    View( ABlock, A, 2*l, 2*l, l, l );
    Scale( mu, ABlock );
    View( ABlock, A, 0, l, l, l );
    MakeWalsh( ABlock, k );
    Scale( -phi, ABlock );
    View( ABlock, A, l, 2*l, l, l );
    MakeWalsh( ABlock, k );
    Scale( phi, ABlock );

    // Now scale R by S
    const R zeta = Sqrt(R(1)-phi*phi);
    DistMatrix<R,U,STAR> d( n, 1, A.Grid() );
    const int colShift = d.ColShift();
    const int colStride = d.ColStride();
    for( int iLocal=0; iLocal<d.LocalHeight(); ++iLocal )
    {
        const int i = colShift + iLocal*colStride;
        d.SetLocal( iLocal, 0, Pow(zeta,i) );
    }
    DiagonalScale( LEFT, NORMAL, d, A );
}

} // namespace elem

#endif // ifndef MATRICES_EXTENDEDKAHAN_HPP
