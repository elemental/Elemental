/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// Gram-Schmidt without normalization.
// TODO: Blocked algorithm
// TODO: Return R instead of M

namespace El {

template<typename F>
void LatticeGramSchmidt( const Matrix<F>& B, Matrix<F>& G, Matrix<F>& M )
{
    DEBUG_ONLY(CSE cse("LatticeGramSchmidt"))
    typedef Base<F> Real;
    const Int n = B.Width();
    G = B;
    Zeros( M, n, n );
    if( n == 0 )
        return;

    Matrix<Real> normsSquared;
    Zeros( normsSquared, n, 1 );
    normsSquared.Set( 0, 0, RealPart(Dot(G(ALL,IR(0)),G(ALL,IR(0)))) );
#ifdef NAIVE_LGS
    for( Int i=1; i<n; ++i )
    {
        auto bi = B( ALL, IR(i) );
        auto gi = G( ALL, IR(i) );
        for( Int j=0; j<i; ++j )
        {
            auto gj = G( ALL, IR(j) );
            const F mu = Dot(bi,gj) / normsSquared.Get(j,0);
            Axpy( -mu, gj, gi );
            M.Set( i, j, mu );
        }
        normsSquared.Set( i, 0, RealPart(Dot(gi,gi)) );
    }
#else
    Matrix<F> z;
    Zeros( z, n, 1 );
    for( Int i=1; i<n; ++i )
    {
        auto b1 = B( ALL, IR(i)   );
        auto g1 = G( ALL, IR(i)   );
        auto G0 = G( ALL, IR(0,i) );

        Gemv( ADJOINT, F(1), G0, b1, z );
        DiagonalSolve( LEFT, NORMAL, normsSquared(IR(0,i),ALL), z );
        Gemv( NORMAL, F(-1), G0, z, F(1), g1 );

        // Store the mu coefficients
        auto m10 = M( IR(i), IR(0,i) );
        Adjoint( z, m10 );

        normsSquared.Set( i, 0, RealPart(Dot(g1,g1)) );
    }
#endif
}

#define PROTO(F) \
  template void LatticeGramSchmidt \
  ( const Matrix<F>& B, Matrix<F>& G, Matrix<F>& M );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
