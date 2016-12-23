/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>

namespace El {

template<typename Real>
void NormsFromScaledSquares
( const Matrix<Real>& localScales,
        Matrix<Real>& localScaledSquares,
        Matrix<Real>& normsLoc,
        mpi::Comm comm )
{
    EL_DEBUG_CSE
    const Int nLocal = localScales.Height();

    // Find the maximum relative scales
    Matrix<Real> scales( nLocal, 1 );
    mpi::AllReduce
    ( localScales.LockedBuffer(), scales.Buffer(), nLocal, mpi::MAX, comm );

    // Equilibrate the local scaled sums
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Real scale = scales(jLoc);
        if( scale != Real(0) )
        {
            // Equilibrate our local scaled sum to the maximum scale
            Real relScale = localScales(jLoc)/scale;
            localScaledSquares(jLoc) *= relScale*relScale;
        }
        else
            localScaledSquares(jLoc) = 0;
    }

    // Combine the local contributions
    Matrix<Real> scaledSquares( nLocal, 1 );
    mpi::AllReduce
    ( localScaledSquares.Buffer(),
      scaledSquares.Buffer(), nLocal, mpi::SUM, comm );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        normsLoc(jLoc) = scales(jLoc)*Sqrt(scaledSquares(jLoc));
}

} // namespace El
