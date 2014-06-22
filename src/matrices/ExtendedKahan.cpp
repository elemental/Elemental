/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// Generate a 3(2^k) x 3(2^k) Extended Kahan matrix, which has the form
// A = S R, where S = diag(1,zeta,...,zeta^(3 2^k - 1)), 
// 
//         | I -phi H_k    0     |
//     R = | 0    I      phi H_k |,
//         | 0    0        I     |
//
// 0 < mu << 1, and phi^2 + zeta^2 = 1.
//

namespace El {

template<typename F> 
inline void MakeExtendedKahan
( Matrix<F>& A, Base<F> phi, Base<F> mu )
{
    DEBUG_ONLY(CallStackEntry cse("MakeExtendedKahan"))
    typedef Base<F> Real;

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

    if( phi <= Real(0) || phi >= Real(1) )
        LogicError("phi must be in (0,1)");
    if( mu <= Real(0) || mu >= Real(1) )
        LogicError("mu must be in (0,1)");

    // Start by setting A to the identity, and then modify the necessary 
    // l x l blocks of its 3 x 3 partitioning.
    MakeIdentity( A );
    auto ABlock = View( A, 2*l, 2*l, l, l );
    Scale( mu, ABlock );
    ABlock = View( A, 0, l, l, l );
    Walsh( ABlock, k );
    Scale( -phi, ABlock );
    ABlock = View( A, l, 2*l, l, l );
    Walsh( ABlock, k );
    Scale( phi, ABlock );

    // Now scale A by S
    const Real zeta = Sqrt(Real(1)-phi*phi);
    for( Int i=0; i<n; ++i )
    {
        const Real gamma = Pow(zeta,Real(i));
        for( Int j=0; j<n; ++j )
            A.Set( i, j, gamma*A.Get(i,j) );
    }
}

template<typename F,Dist U,Dist V>
inline void MakeExtendedKahan
( DistMatrix<F,U,V>& A, Base<F> phi, Base<F> mu )
{
    DEBUG_ONLY(CallStackEntry cse("MakeExtendedKahan"))
    typedef Base<F> Real;

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

    if( phi <= Real(0) || phi >= Real(1) )
        LogicError("phi must be in (0,1)");
    if( mu <= Real(0) || mu >= Real(1) )
        LogicError("mu must be in (0,1)");

    // Start by setting A to the identity, and then modify the necessary 
    // l x l blocks of its 3 x 3 partitioning.
    MakeIdentity( A );
    auto ABlock = View( A, 2*l, 2*l, l, l );
    Scale( mu, ABlock );
    ABlock = View( A, 0, l, l, l );
    Walsh( ABlock, k );
    Scale( -phi, ABlock );
    ABlock = View( A, l, 2*l, l, l );
    Walsh( ABlock, k );
    Scale( phi, ABlock );

    // Now scale A by S
    const Real zeta = Sqrt(Real(1)-phi*phi);
    for( Int iLoc=0; iLoc<A.LocalHeight(); ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        const Real gamma = Pow(zeta,Real(i));
        for( Int jLoc=0; jLoc<A.LocalWidth(); ++jLoc )
            A.SetLocal( iLoc, jLoc, gamma*A.GetLocal(iLoc,jLoc) );
    }
}

template<typename F>
void ExtendedKahan( Matrix<F>& A, Int k, Base<F> phi, Base<F> mu )
{
    DEBUG_ONLY(CallStackEntry cse("ExtendedKahan"))
    const Int n = 3*(1u<<k);
    A.Resize( n, n );
    MakeExtendedKahan( A, phi, mu );
}

template<typename F,Dist U,Dist V>
void ExtendedKahan( DistMatrix<F,U,V>& A, Int k, Base<F> phi, Base<F> mu )
{
    DEBUG_ONLY(CallStackEntry cse("ExtendedKahan"))
    const Int n = 3*(1u<<k);
    A.Resize( n, n );
    MakeExtendedKahan( A, phi, mu );
}

#define PROTO_DIST(F,U,V) \
  template void ExtendedKahan \
  ( DistMatrix<F,U,V>& A, Int k, Base<F> phi, Base<F> mu );

#define PROTO(F) \
  template void ExtendedKahan( Matrix<F>& A, Int k, Base<F> phi, Base<F> mu ); \
  PROTO_DIST(F,CIRC,CIRC) \
  PROTO_DIST(F,MC,  MR  ) \
  PROTO_DIST(F,MC,  STAR) \
  PROTO_DIST(F,MD,  STAR) \
  PROTO_DIST(F,MR,  MC  ) \
  PROTO_DIST(F,MR,  STAR) \
  PROTO_DIST(F,STAR,MC  ) \
  PROTO_DIST(F,STAR,MD  ) \
  PROTO_DIST(F,STAR,MR  ) \
  PROTO_DIST(F,STAR,STAR) \
  PROTO_DIST(F,STAR,VC  ) \
  PROTO_DIST(F,STAR,VR  ) \
  PROTO_DIST(F,VC,  STAR) \
  PROTO_DIST(F,VR,  STAR)

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
