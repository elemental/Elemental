/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/matrices.hpp>

// The inverse of a scaled Jordan block

namespace El {

template<typename F> 
void Demmel( Matrix<F>& A, Int n )
{
    DEBUG_CSE
    typedef Base<F> Real;
    A.Resize( n, n );
    if( n == 1 )
    {
        A.Set( 0, 0, -Real(1) );
        return;
    }
    else if( n == 0 )
        return;

    const Real beta = Pow(Real(10),Real(4)/(n-1));

    const Int numDiags = 2*n-1;
    vector<F> a( numDiags, 0 );
    for( Int j=0; j<n; ++j )
        a[j] = -Pow(beta,Real(n-1-j));
    for( Int j=n; j<numDiags; ++j )
        a[j] = 0;
    Toeplitz( A, n, n, a );
}

template<typename F>
void Demmel( AbstractDistMatrix<F>& A, Int n )
{
    DEBUG_CSE
    typedef Base<F> Real;
    A.Resize( n, n );
    if( n == 1 )
    {
        A.Set( 0, 0, -Real(1) );
        return;
    }
    else if( n == 0 )
        return;
    
    const Real beta = Pow(Real(10),Real(4)/(n-1));

    const Int numDiags = 2*n-1;
    vector<F> a( numDiags, 0 );
    for( Int j=0; j<n; ++j )
        a[j] = -Pow(beta,Real(n-1-j));
    for( Int j=n; j<numDiags; ++j )
        a[j] = 0;
    Toeplitz( A, n, n, a );
}

#define PROTO(F) \
  template void Demmel( Matrix<F>& A, Int n ); \
  template void Demmel( AbstractDistMatrix<F>& A, Int n );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
