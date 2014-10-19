/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// An example of Bunch-Kaufman A producing large element growth in 
// floating-point arithmetic. Please see Theorem 5 from:
//     http://www.alexdruinsky.com/pdfs/bkbound-revised.pdf

template<typename F> 
void DruinskyToledo( Matrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("DruinskyToledo"))
    if( n <= 2 )
    {
        Ones( A, n, n );
        return;
    }
    typedef Base<F> Real;
    const Real phi = Real(1) + 4*lapack::MachineEpsilon<Real>();
    const Real alphaPhi = LDLPivotConstant<Real>(BUNCH_KAUFMAN_A)*phi;
    std::vector<Real> d( n-2 );
    Real sigma(1);
    for( Int i=0; i<n-2; ++i )
    {
        d[i] = -alphaPhi/sigma;
        sigma -= 1/d[i];
    }

    A.Resize( n, n );

    auto ABL = A(IR(n-2,n),IR(0,n));
    Ones( ABL, ABL.Height(), ABL.Width() );

    auto ATR = A(IR(0,n-2),IR(n-2,n));
    Ones( ATR, ATR.Height(), ATR.Width() );

    auto ABR = A(IR(0,n-2),IR(0,n-2));
    Diagonal( ABR, d );
}

template<typename F> 
void DruinskyToledo( AbstractDistMatrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("DruinskyToledo"))
    if( n <= 2 )
    {
        Ones( A, n, n );
        return;
    }
    typedef Base<F> Real;
    const Real phi = Real(1) + 4*lapack::MachineEpsilon<Real>();
    const Real alphaPhi = LDLPivotConstant<Real>(BUNCH_KAUFMAN_A)*phi;
    std::vector<Real> d( n-2 );
    Real sigma(1);
    for( Int i=0; i<n-2; ++i )
    {
        d[i] = -alphaPhi/sigma;
        sigma -= 1/d[i];
    }

    A.Resize( n, n );

    auto ASub = std::unique_ptr<AbstractDistMatrix<F>>( A.Construct() );

    View( *ASub, A, IR(n-2,n), IR(0,  n) );
    Ones( *ASub, 2, n );

    View( *ASub, A, IR(0,n-2), IR(n-2,n) );
    Ones( *ASub, n-2, 2 );

    View( *ASub, A, IR(0,n-2), IR(0,n-2) );
    Diagonal( *ASub, d );
}

#define PROTO(F) \
  template void DruinskyToledo( Matrix<F>& A, Int n ); \
  template void DruinskyToledo( AbstractDistMatrix<F>& A, Int n );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
