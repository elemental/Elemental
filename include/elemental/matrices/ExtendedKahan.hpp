/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_EXTENDEDKAHAN_HPP
#define ELEM_EXTENDEDKAHAN_HPP

#include "./Walsh.hpp"

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
MakeExtendedKahan( Matrix<F>& A, Base<F> phi, Base<F> mu )
{
    DEBUG_ONLY(CallStackEntry cse("MakeExtendedKahan"))
    typedef Base<F> R;

    if( A.Height() != A.Width() )
        LogicError("Extended Kahan matrices must be square");
    const Int n = A.Height();
    if( n % 3 != 0 )
        LogicError("Dimension must be an integer multiple of 3");
    const Int l = n / 3;
    if( !l || (l & (l-1)) )
        LogicError("n/3 is not a power of two");
    Int k=0;
    while( Int(1u<<k) < l )
        ++k;

    if( phi <= R(0) || phi >= R(1) )
        LogicError("phi must be in (0,1)");
    if( mu <= R(0) || mu >= R(1) )
        LogicError("mu must be in (0,1)");

    // Start by setting A to the identity, and then modify the necessary 
    // l x l blocks of its 3 x 3 partitioning.
    MakeIdentity( A );
    auto ABlock = View( A, 2*l, 2*l, l, l );
    Scale( mu, ABlock );
    ABlock = View( A, 0, l, l, l );
    MakeWalsh( ABlock, k );
    Scale( -phi, ABlock );
    ABlock = View( A, l, 2*l, l, l );
    MakeWalsh( ABlock, k );
    Scale( phi, ABlock );

    // Now scale R by S
    const R zeta = Sqrt(R(1)-phi*phi);
    Matrix<R> d( n, 1 );
    for( Int i=0; i<n; ++i )
        d.Set( i, 0, Pow(zeta,R(i)) );
    DiagonalScale( LEFT, NORMAL, d, A );
}

template<typename F,Dist U,Dist V>
inline void
MakeExtendedKahan( DistMatrix<F,U,V>& A, Base<F> phi, Base<F> mu )
{
    DEBUG_ONLY(CallStackEntry cse("MakeExtendedKahan"))
    typedef Base<F> R;

    if( A.Height() != A.Width() )
        LogicError("Extended Kahan matrices must be square");
    const Int n = A.Height();
    if( n % 3 != 0 )
        LogicError("Dimension must be an integer multiple of 3");
    const Int l = n / 3;
    if( !l || (l & (l-1)) )
        LogicError("n/3 is not a power of two");
    Int k=0;
    while( Int(1u<<k) < l )
        ++k;

    if( phi <= R(0) || phi >= R(1) )
        LogicError("phi must be in (0,1)");
    if( mu <= R(0) || mu >= R(1) )
        LogicError("mu must be in (0,1)");

    // Start by setting A to the identity, and then modify the necessary 
    // l x l blocks of its 3 x 3 partitioning.
    MakeIdentity( A );
    auto ABlock = View( A, 2*l, 2*l, l, l );
    Scale( mu, ABlock );
    ABlock = View( A, 0, l, l, l );
    MakeWalsh( ABlock, k );
    Scale( -phi, ABlock );
    ABlock = View( A, l, 2*l, l, l );
    MakeWalsh( ABlock, k );
    Scale( phi, ABlock );

    // Now scale R by S
    const R zeta = Sqrt(R(1)-phi*phi);
    DistMatrix<R,U,STAR> d( n, 1, A.Grid() );
    for( Int iLoc=0; iLoc<d.LocalHeight(); ++iLoc )
    {
        const Int i = d.GlobalRow(iLoc);
        d.SetLocal( iLoc, 0, Pow(zeta,R(i)) );
    }
    DiagonalScale( LEFT, NORMAL, d, A );
}

template<typename F,Dist U,Dist V>
inline void
MakeExtendedKahan( BlockDistMatrix<F,U,V>& A, Base<F> phi, Base<F> mu )
{
    DEBUG_ONLY(CallStackEntry cse("MakeExtendedKahan"))
    typedef Base<F> R;

    if( A.Height() != A.Width() )
        LogicError("Extended Kahan matrices must be square");
    const Int n = A.Height();
    if( n % 3 != 0 )
        LogicError("Dimension must be an integer multiple of 3");
    const Int l = n / 3;
    if( !l || (l & (l-1)) )
        LogicError("n/3 is not a power of two");
    Int k=0;
    while( Int(1u<<k) < l )
        ++k;

    if( phi <= R(0) || phi >= R(1) )
        LogicError("phi must be in (0,1)");
    if( mu <= R(0) || mu >= R(1) )
        LogicError("mu must be in (0,1)");

    // Start by setting A to the identity, and then modify the necessary 
    // l x l blocks of its 3 x 3 partitioning.
    MakeIdentity( A );
    auto ABlock = View( A, 2*l, 2*l, l, l );
    Scale( mu, ABlock );
    ABlock = View( A, 0, l, l, l );
    MakeWalsh( ABlock, k );
    Scale( -phi, ABlock );
    ABlock = View( A, l, 2*l, l, l );
    MakeWalsh( ABlock, k );
    Scale( phi, ABlock );

    // Now scale R by S
    const R zeta = Sqrt(R(1)-phi*phi);
    BlockDistMatrix<R,U,STAR> d( n, 1, A.Grid() );
    for( Int iLoc=0; iLoc<d.LocalHeight(); ++iLoc )
    {
        const Int i = d.GlobalRow(iLoc);
        d.SetLocal( iLoc, 0, Pow(zeta,R(i)) );
    }
    DiagonalScale( LEFT, NORMAL, d, A );
}

template<typename F>
inline void
ExtendedKahan( Matrix<F>& A, Int k, Base<F> phi, Base<F> mu )
{
    DEBUG_ONLY(CallStackEntry cse("ExtendedKahan"))
    const Int n = 3*(1u<<k);
    A.Resize( n, n );
    MakeExtendedKahan( A, phi, mu );
}

template<typename F,Dist U,Dist V>
inline void
ExtendedKahan( DistMatrix<F,U,V>& A, Int k, Base<F> phi, Base<F> mu )
{
    DEBUG_ONLY(CallStackEntry cse("ExtendedKahan"))
    const Int n = 3*(1u<<k);
    A.Resize( n, n );
    MakeExtendedKahan( A, phi, mu );
}

template<typename F,Dist U,Dist V>
inline void
ExtendedKahan( BlockDistMatrix<F,U,V>& A, Int k, Base<F> phi, Base<F> mu )
{
    DEBUG_ONLY(CallStackEntry cse("ExtendedKahan"))
    const Int n = 3*(1u<<k);
    A.Resize( n, n );
    MakeExtendedKahan( A, phi, mu );
}

} // namespace elem

#endif // ifndef ELEM_EXTENDEDKAHAN_HPP
