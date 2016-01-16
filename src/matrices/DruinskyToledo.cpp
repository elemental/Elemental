/*
   Copyright (c) 2009-2016, Jack Poulson
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
//
// The modification used to make the matrix better-conditioned is to form
//   A = | G I |, where G is the matrix described in Section 4.
//       | I I |

template<typename F> 
void DruinskyToledo( Matrix<F>& A, Int k )
{
    DEBUG_ONLY(CSE cse("DruinskyToledo"))
    const Int n = 2*k;
    Zeros( A, n, n );
    if( k == 0 )
      return;
    if( k == 1 )
    {
        Ones( A, n, n );
        return;
    }
    typedef Base<F> Real;
    const Real phi = Real(1) + 4*limits::Epsilon<Real>();
    const Real alphaPhi = LDLPivotConstant<Real>(BUNCH_KAUFMAN_A)*phi;
    vector<Real> d( k-2 );
    Real sigma(1);
    for( Int i=0; i<k-2; ++i )
    {
        d[i] = -alphaPhi/sigma;
        sigma -= 1/d[i];
    }

    auto GBL = A(IR(k-2,k),IR(0,k));
    Ones( GBL, GBL.Height(), GBL.Width() );

    auto GTR = A(IR(0,k),IR(k-2,k));
    Ones( GTR, GTR.Height(), GTR.Width() );

    auto GTL = A(IR(0,k-2),IR(0,k-2));
    Diagonal( GTL, d );

    auto ABL = A(IR(k,n),IR(0,k));
    Identity( ABL, k, k );

    auto ABR = A(IR(k,n),IR(k,n));
    Identity( ABR, k, k );

    auto ATR = A(IR(0,k),IR(k,n));
    Identity( ATR, k, k );
}

template<typename F> 
void DruinskyToledo( ElementalMatrix<F>& A, Int k )
{
    DEBUG_ONLY(CSE cse("DruinskyToledo"))
    const Int n = 2*k;
    Zeros( A, n, n );
    if( k == 0 )
      return;
    if( k == 1 )
    {
        Ones( A, n, n );
        return;
    }
    typedef Base<F> Real;
    const Real phi = Real(1) + 4*limits::Epsilon<Real>();
    const Real alphaPhi = LDLPivotConstant<Real>(BUNCH_KAUFMAN_A)*phi;
    vector<Real> d( k-2 );
    Real sigma(1);
    for( Int i=0; i<k-2; ++i )
    {
        d[i] = -alphaPhi/sigma;
        sigma -= 1/d[i];
    }

    unique_ptr<ElementalMatrix<F>> ASub( A.Construct(A.Grid(),A.Root()) );

    View( *ASub, A, IR(k-2,k), IR(0,k) );
    Ones( *ASub, 2, k );

    View( *ASub, A, IR(0,k), IR(k-2,k) );
    Ones( *ASub, k, 2 );

    View( *ASub, A, IR(0,k-2), IR(0,k-2) );
    Diagonal( *ASub, d );

    View( *ASub, A, IR(k,n), IR(0,k) );
    Identity( *ASub, k, k );

    View( *ASub, A, IR(k,n), IR(k,n) );
    Identity( *ASub, k, k );

    View( *ASub, A, IR(0,k), IR(k,n) );
    Identity( *ASub, k, k );
}

#define PROTO(F) \
  template void DruinskyToledo( Matrix<F>& A, Int k ); \
  template void DruinskyToledo( ElementalMatrix<F>& A, Int k );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
