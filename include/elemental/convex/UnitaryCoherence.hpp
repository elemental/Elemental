/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CONVEX_UNITARYCOHERENCE_HPP
#define CONVEX_UNITARYCOHERENCE_HPP

// Definition taken from Eq. (1.8) from: 
//   Emmanuel J. Candes and Benjamin Recht, 
//   "Exact matrix completion via convex optimization",
//   http://www-stat.stanford.edu/~candes/papers/MatrixCompletion.pdf
//
// NOTE: U is assumed to have orthonormal columns, and that the return value
//       should be between 1 (when U is completely incoherent) and n/r, when 
//       the span of the columns of U contains a standard basis vector.

#include "elemental/blas-like/level1/MakeHermitian.hpp"
#include "elemental/blas-like/level3/Herk.hpp"
#include <algorithm>

namespace elem {

template<typename F>
inline BASE(F)
UnitaryCoherence( Matrix<F>& U )
{
#ifndef RELEASE
    CallStackEntry entry("UnitaryCoherence");
#endif
    typedef BASE(F) R;
    const int n = U.Height();
    const int r = U.Width();

    // Z := U U' in n^2 r work
    Matrix<F> Z;
    Herk( UPPER, NORMAL, F(1), U, Z );

    // Now make Z explicitly Hermitian so that our job is easier
    MakeHermitian( UPPER, Z );

    // Compute the maximum column two-norm
    R maxColNorm = 0;
    for( int j=0; j<n; ++j )
    {
        const R colNorm = blas::Nrm2( n, Z.LockedBuffer(0,j), 1 );
        maxColNorm = std::max( colNorm, maxColNorm );
    }
    return (n*maxColNorm*maxColNorm)/r;
}

template<typename F>
inline BASE(F)
UnitaryCoherence( DistMatrix<F>& U )
{
#ifndef RELEASE
    CallStackEntry entry("UnitaryCoherence");
#endif
    typedef BASE(F) R;
    const Grid& grid = U.Grid();
    const int n = U.Height();
    const int r = U.Width();

    // Z := U U' in n^2 r work
    DistMatrix<F> Z( grid );
    Herk( UPPER, NORMAL, F(1), U, Z );

    // Now make Z explicitly Hermitian so that our job is easier
    MakeHermitian( UPPER, Z );

    // Compute the maximum column two-norm squared
    const int localWidth = Z.LocalWidth();
    const int localHeight = Z.LocalHeight();
    std::vector<R> normsSquared( localWidth );
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const R localNorm = 
            blas::Nrm2( localHeight, Z.LockedBuffer(0,jLocal), 1 );
        normsSquared[jLocal] = localNorm*localNorm;
    }
    mpi::AllReduce( &normsSquared[0], localWidth, mpi::SUM, grid.ColComm() );
    R maxLocalNormSquared = 
        *std::max_element( normsSquared.begin(), normsSquared.end() );
    R maxNormSquared;
    mpi::AllReduce
    ( &maxLocalNormSquared, &maxNormSquared, 1, mpi::MAX, grid.RowComm() );

    return (n*maxNormSquared)/r;
}

} // namespace elem

#endif // ifndef CONVEX_UNITARYCOHERENCE_HPP
