/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_ID_HPP
#define LAPACK_ID_HPP

#include "elemental/lapack-like/QR/BusingerGolub.hpp"

// TODO: Add detailed references to Tygert et al.'s ID package and the papers
//       "Randomized algorithms for the low-rank approximation of matrices", 
//       "A randomized algorithm for principal component analysis", and 
//       "On the compression of low-rank matrices"

namespace elem {

namespace id {

template<typename F>
inline void
PseudoTrsm( const Matrix<F>& RL, Matrix<F>& RR, BASE(F) tol )
{
#ifndef RELEASE
    CallStackEntry entry("id::PseudoTrsm");
#endif
    typedef BASE(F) Real;
    const int m = RR.Height();
    const int n = RR.Width();

    // Compute the spectral radius of the triangular matrix
    Real maxAbsEig = 0;
    for( int i=0; i<m; ++i )
        maxAbsEig = std::max( maxAbsEig, Abs(RL.Get(i,i)) );

    for( int i=m-1; i>=0; --i )
    {
        // Apply the pseudo-inverse of the i'th diagonal value of RL 
        const F rho = RL.Get(i,i);
        const Real rhoAbs = Abs(rho);
        if( rhoAbs >= tol*maxAbsEig )
        {
            for( int j=0; j<n; ++j ) 
            {
                const F zeta = RR.Get(i,j);
                RR.Set(i,j,zeta/rho);
            }
        }
        else
        {
            for( int j=0; j<n; ++j )
                RR.Set(i,j,0);
        }

        // Now update RR using an outer-product of the column of RL above the 
        // i'th diagonal with the i'th row of RR
        blas::Geru
        ( i, n, 
          F(-1), RL.LockedBuffer(0,i), 1, RR.LockedBuffer(i,0), RR.LDim(), 
                 RR.Buffer(0,0), RR.LDim() );
    }
}

template<typename F>
inline void
PseudoTrsm( const DistMatrix<F>& RL, DistMatrix<F>& RR, BASE(F) tol )
{
#ifndef RELEASE
    CallStackEntry entry("id::PseudoTrsm");
#endif
    // For now, assume that RL is sufficiently small and give each process a 
    // full copy so that we may independently apply its pseudoinverse to each
    // column of RR
    DistMatrix<F,STAR,STAR> RL_STAR_STAR( RL );
    DistMatrix<F,STAR,VR> RR_STAR_VR( RR );
    PseudoTrsm( RL_STAR_STAR.Matrix(), RR_STAR_VR.Matrix(), tol );
    RR = RR_STAR_VR;
}

} // namespace id

// On output, the matrix Z contains the non-trivial portion of the interpolation
// matrix, and p contains the pivots used during the iterations of 
// pivoted QR. Either 'maxSteps' iterations are reached, or a pivot value less 
// than or equal to tol times the original pivot value was found. 
// The input matrix A is unchanged.

template<typename Real> 
inline void
ID( const Matrix<Real>& A, Matrix<int>& p, Matrix<Real>& Z, 
    int maxSteps, Real tol )
{
#ifndef RELEASE
    CallStackEntry entry("ID");
#endif
    const int n = A.Width();

    // Perform the pivoted QR factorization on a copy of A
    Matrix<Real> ACopy( A );
    qr::BusingerGolub( ACopy, p, maxSteps, tol );
    const int numSteps = p.Height();

    Real pinvTol;
    if( tol < Real(0) )
    {
        const Real epsilon = lapack::MachineEpsilon<Real>();
        pinvTol = numSteps*epsilon;
    }
    else
        pinvTol = tol;

    // Now form a minimizer of || RL Z - RR ||_2 via pseudo triangular solves
    Matrix<Real> RL, RR;
    LockedView( RL, ACopy, 0, 0, numSteps, numSteps );
    LockedView( RR, ACopy, 0, numSteps, numSteps, n-numSteps );
    Z = RR;
    id::PseudoTrsm( RL, Z, pinvTol );
}

template<typename Real> 
inline void
ID( const Matrix<Real>& A, Matrix<int>& p, Matrix<Real>& Z, int numSteps )
{
#ifndef RELEASE
    CallStackEntry entry("ID");
#endif
    // Use a negative tolerance to guarantee numSteps iterations of QR
    ID( A, p, Z, numSteps, Real(-1) );
}

// This only exists since complex QR has an extra return argument related to
// the phases of the pseudo-Householder transformations used
template<typename Real> 
inline void
ID( const Matrix<Complex<Real> >& A, Matrix<int>& p, Matrix<Complex<Real> >& Z, 
    int maxSteps, Real tol )
{
#ifndef RELEASE
    CallStackEntry entry("ID");
#endif
    typedef Complex<Real> C;
    const int n = A.Width();

    // Perform the pivoted QR factorization on a copy of A
    Matrix<C> ACopy( A ), t;
    qr::BusingerGolub( ACopy, t, p, maxSteps, tol );
    const int numSteps = p.Height();

    Real pinvTol;
    if( tol < Real(0) )
    {
        const Real epsilon = lapack::MachineEpsilon<Real>();
        pinvTol = numSteps*epsilon;
    }
    else
        pinvTol = tol;

    // Now form a minimizer of || RL Z - RR ||_2 via pseudo triangular solves
    Matrix<C> RL, RR;
    LockedView( RL, ACopy, 0, 0, numSteps, numSteps );
    LockedView( RR, ACopy, 0, numSteps, numSteps, n-numSteps );
    Z = RR;
    id::PseudoTrsm( RL, Z, pinvTol );
}

template<typename Real> 
inline void
ID( const Matrix<Complex<Real> >& A, Matrix<int>& p, Matrix<Complex<Real> >& Z, 
    int numSteps )
{
#ifndef RELEASE
    CallStackEntry entry("ID");
#endif
    ID( A, p, Z, numSteps, Real(-1) );
}

template<typename Real> 
inline void
ID
( const DistMatrix<Real>& A, DistMatrix<int,VR,STAR>& p, DistMatrix<Real>& Z, 
  int maxSteps, Real tol )
{
#ifndef RELEASE
    CallStackEntry entry("ID");
#endif
    const Grid& g = A.Grid();
    const int n = A.Width();

    // Perform the pivoted QR factorization on a copy of A
    DistMatrix<Real> ACopy( A );
    qr::BusingerGolub( ACopy, p, maxSteps, tol );
    const int numSteps = p.Height();

    Real pinvTol;
    if( tol < Real(0) )
    {
        const Real epsilon = lapack::MachineEpsilon<Real>();
        pinvTol = numSteps*epsilon;
    }
    else
        pinvTol = tol;

    // Now form a minimizer of || RL Z - RR ||_2 via pseudo triangular solves
    DistMatrix<Real> RL(g), RR(g);
    LockedView( RL, ACopy, 0, 0, numSteps, numSteps );
    LockedView( RR, ACopy, 0, numSteps, numSteps, n-numSteps );
    Z = RR;
    id::PseudoTrsm( RL, Z, pinvTol );
}

template<typename Real> 
inline void
ID
( const DistMatrix<Real>& A, DistMatrix<int,VR,STAR>& p, DistMatrix<Real>& Z, 
  int numSteps )
{
#ifndef RELEASE
    CallStackEntry entry("ID");
#endif
    ID( A, p, Z, numSteps, Real(-1) );
}

// This only exists since complex QR has an extra return argument related to
// the phases of the pseudo-Householder transformations used
template<typename Real> 
inline void
ID
( const DistMatrix<Complex<Real> >& A, DistMatrix<int,VR,STAR>& p, 
        DistMatrix<Complex<Real> >& Z, int maxSteps, Real tol )
{
#ifndef RELEASE
    CallStackEntry entry("ID");
#endif
    typedef Complex<Real> C;
    const Grid& g = A.Grid();
    const int n = A.Width();

    // Perform the pivoted QR factorization on a copy of A
    DistMatrix<C> ACopy( A );
    DistMatrix<C,MD,STAR> t(g);
    qr::BusingerGolub( ACopy, t, p, maxSteps, tol );
    const int numSteps = p.Height();

    Real pinvTol;
    if( tol < Real(0) )
    {
        const Real epsilon = lapack::MachineEpsilon<Real>();
        pinvTol = numSteps*epsilon;
    }
    else
        pinvTol = tol;

    // Now form a minimizer of || RL Z - RR ||_2 via pseudo triangular solves
    DistMatrix<C> RL(g), RR(g);
    LockedView( RL, ACopy, 0, 0, numSteps, numSteps );
    LockedView( RR, ACopy, 0, numSteps, numSteps, n-numSteps );
    Z = RR;
    id::PseudoTrsm( RL, Z, pinvTol );
}

template<typename Real> 
inline void
ID
( const DistMatrix<Complex<Real> >& A, DistMatrix<int,VR,STAR>& p, 
        DistMatrix<Complex<Real> >& Z, int numSteps )
{
#ifndef RELEASE
    CallStackEntry entry("ID");
#endif
    ID( A, p, Z, numSteps, Real(-1) );
}

} // namespace elem

#endif // ifndef LAPACK_ID_HPP
