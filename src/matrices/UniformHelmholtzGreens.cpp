/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Generate the "random Green's matrix" from 
//   A. Goetschy and S. E. Skipetrov's "Non-Hermitian Euclidean random matrix
//   theory".
// It is essentially the 3D Helmholtz Green's function restricted to a set
// of points chosen uniformly from the unit ball. The behaviour of the spectrum
// is known to change dramatically dependinging upon the number of points 
// sampled per wavelength-cubed.

template<typename Real>
void UniformHelmholtzGreens( Matrix<Complex<Real>>& A, Int n, Real lambda )
{
    DEBUG_ONLY(CSE cse("UniformHelmholtzGreens"))
    typedef Complex<Real> C;
    const Real pi = 4*Atan( Real(1) );
    const Real k0 = 2*pi/lambda;

    // Generate a list of n uniform samples from the 3D unit ball
    Matrix<Real> X(3,n);
    for( Int j=0; j<n; ++j )
    {
        Real x0, x1, x2;
        // Sample uniformly from [-1,+1]^3 until a point is drawn from the ball
        while( true )
        {
            x0 = SampleUniform( Real(-1), Real(1) ); 
            x1 = SampleUniform( Real(-1), Real(1) );
            x2 = SampleUniform( Real(-1), Real(1) );
            const Real radiusSq = x0*x0 + x1*x1 + x2*x2;
            if( radiusSq > 0 && radiusSq <= 1 )
                break;
        }
        X.Set( 0, j, x0 );
        X.Set( 1, j, x1 );
        X.Set( 2, j, x2 );
    }
    
    A.Resize( n, n );
    for( Int j=0; j<n; ++j )
    {
        const Real xj0 = X.Get(0,j);
        const Real xj1 = X.Get(1,j);
        const Real xj2 = X.Get(2,j);
        for( Int i=0; i<n; ++i )
        {
            if( i == j ) 
            {
                A.Set( i, j, 0 );
            }
            else
            {
                const Real d0 = X.Get(0,i)-xj0;
                const Real d1 = X.Get(1,i)-xj1;
                const Real d2 = X.Get(2,i)-xj2;
                const Real gamma = k0*Sqrt(d0*d0+d1*d1+d2*d2);
                const Real realPart = Cos(gamma)/gamma;
                const Real imagPart = Sin(gamma)/gamma;
                A.Set( i, j, C(realPart,imagPart) );
            }
        }
    }
}

template<typename Real>
void UniformHelmholtzGreens
( ElementalMatrix<Complex<Real>>& A, Int n, Real lambda )
{
    DEBUG_ONLY(CSE cse("UniformHelmholtzGreens"))
    typedef Complex<Real> C;
    const Real pi = 4*Atan( Real(1) );
    const Real k0 = 2*pi/lambda;
    const Grid& g = A.Grid();

    // Generate a list of n uniform samples from the 3D unit ball
    DistMatrix<Real,STAR,VR> X_STAR_VR(3,n,g);
    for( Int jLoc=0; jLoc<X_STAR_VR.LocalWidth(); ++jLoc )
    {
        Real x0, x1, x2;
        // Sample uniformly from [-1,+1]^3 until a point is drawn from the ball
        while( true )
        {
            x0 = SampleUniform( Real(-1), Real(1) );
            x1 = SampleUniform( Real(-1), Real(1) );
            x2 = SampleUniform( Real(-1), Real(1) );
            const Real radiusSq = x0*x0 + x1*x1 + x2*x2;
            if( radiusSq > 0 && radiusSq <= 1 )
                break;
        }
        X_STAR_VR.SetLocal( 0, jLoc, x0 );
        X_STAR_VR.SetLocal( 1, jLoc, x1 );
        X_STAR_VR.SetLocal( 2, jLoc, x2 );
    }
    DistMatrix<Real,STAR,STAR> X_STAR_STAR( X_STAR_VR );

    A.Resize( n, n );
    for( Int jLoc=0; jLoc<A.LocalWidth(); ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        const Real xj0 = X_STAR_STAR.GetLocal(0,j);
        const Real xj1 = X_STAR_STAR.GetLocal(1,j);
        const Real xj2 = X_STAR_STAR.GetLocal(2,j);
        for( Int iLoc=0; iLoc<A.LocalHeight(); ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            if( i == j )
            {
                A.SetLocal( iLoc, jLoc, 0 );
            }
            else
            {
                const Real d0 = X_STAR_STAR.GetLocal(0,i)-xj0;
                const Real d1 = X_STAR_STAR.GetLocal(1,i)-xj1;
                const Real d2 = X_STAR_STAR.GetLocal(2,i)-xj2;
                const Real gamma = k0*Sqrt(d0*d0+d1*d1+d2*d2);
                const Real realPart = Cos(gamma)/gamma;
                const Real imagPart = Sin(gamma)/gamma;
                A.SetLocal( iLoc, jLoc, C(realPart,imagPart) );
            }
        }
    }
}

#define PROTO(Real) \
  template void UniformHelmholtzGreens \
  ( Matrix<Complex<Real>>& A, Int n, Real lambda ); \
  template void UniformHelmholtzGreens \
  ( ElementalMatrix<Complex<Real>>& A, Int n, Real lambda );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
