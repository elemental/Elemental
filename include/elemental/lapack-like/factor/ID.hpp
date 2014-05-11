/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_ID_HPP
#define ELEM_ID_HPP

#include ELEM_QR_BUSINGERGOLUB_INC

// TODO: Add detailed references to Tygert et al.'s ID package and the papers
//       "Randomized algorithms for the low-rank approximation of matrices", 
//       "A randomized algorithm for principal component analysis", and 
//       "On the compression of low-rank matrices"

namespace elem {

namespace id {

template<typename F>
inline void
PseudoTrsm( const Matrix<F>& RL, Matrix<F>& RR, Base<F> tol )
{
    DEBUG_ONLY(CallStackEntry cse("id::PseudoTrsm"))
    typedef Base<F> Real;
    const Int m = RR.Height();
    const Int n = RR.Width();

    // Compute the spectral radius of the triangular matrix
    Real maxAbsEig = 0;
    for( Int i=0; i<m; ++i )
        maxAbsEig = std::max( maxAbsEig, Abs(RL.Get(i,i)) );

    for( Int i=m-1; i>=0; --i )
    {
        // Apply the pseudo-inverse of the i'th diagonal value of RL 
        const F rho = RL.Get(i,i);
        const Real rhoAbs = Abs(rho);
        if( rhoAbs >= tol*maxAbsEig )
        {
            for( Int j=0; j<n; ++j ) 
            {
                const F zeta = RR.Get(i,j);
                RR.Set(i,j,zeta/rho);
            }
        }
        else
        {
            for( Int j=0; j<n; ++j )
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

// For now, assume that RL is sufficiently small and give each process a full
// copy so that we may independently apply its pseudoinverse to each column of
// RR
template<typename F>
inline void
PseudoTrsm( const DistMatrix<F>& RL, DistMatrix<F,STAR,VR>& RR, Base<F> tol )
{
    DEBUG_ONLY(CallStackEntry cse("id::PseudoTrsm"))
    DistMatrix<F,STAR,STAR> RL_STAR_STAR( RL );
    PseudoTrsm( RL_STAR_STAR.Matrix(), RR.Matrix(), tol );
}

// On output, the matrix Z contains the non-trivial portion of the interpolation
// matrix, and p contains the pivots used during the iterations of 
// pivoted QR. Either 'maxSteps' iterations are reached, or a pivot value less 
// than or equal to tol times the original pivot value was found. 
// The input matrix A is unchanged.

template<typename F> 
inline void
BusingerGolub
( Matrix<F>& A, Matrix<Int>& pPerm, Matrix<F>& Z, Int maxSteps, Base<F> tol )
{
    DEBUG_ONLY(CallStackEntry cse("id::BusingerGolub"))
    typedef Base<F> Real;
    const Int n = A.Width();

    // Perform the pivoted QR factorization
    const Int numSteps = qr::BusingerGolub( A, pPerm, maxSteps, tol );

    Real pinvTol;
    if( tol < Real(0) )
    {
        const Real epsilon = lapack::MachineEpsilon<Real>();
        pinvTol = numSteps*epsilon;
    }
    else
        pinvTol = tol;

    // Now form a minimizer of || RL Z - RR ||_2 via pseudo triangular solves
    auto RL = LockedViewRange( A, 0, 0,        numSteps, numSteps );
    auto RR = LockedViewRange( A, 0, numSteps, numSteps, n        );
    Z = RR;
    PseudoTrsm( RL, Z, pinvTol );
}

template<typename F,Dist UPerm> 
inline void
BusingerGolub
( DistMatrix<F>& A, 
  DistMatrix<Int,UPerm,STAR>& pPerm, 
  DistMatrix<F,STAR,VR>& Z, 
  Int maxSteps, Base<F> tol )
{
    DEBUG_ONLY(CallStackEntry cse("id::BusingerGolub"))
    typedef Base<F> Real;
    const Int n = A.Width();

    // Perform the pivoted QR factorization on a copy of A
    const Int numSteps = qr::BusingerGolub( A, pPerm, maxSteps, tol );

    Real pinvTol;
    if( tol < Real(0) )
    {
        const Real epsilon = lapack::MachineEpsilon<Real>();
        pinvTol = numSteps*epsilon;
    }
    else
        pinvTol = tol;

    // Now form a minimizer of || RL Z - RR ||_2 via pseudo triangular solves
    auto RL = LockedViewRange( A, 0, 0,        numSteps, numSteps );
    auto RR = LockedViewRange( A, 0, numSteps, numSteps, n        );
    Z = RR;
    PseudoTrsm( RL, Z, pinvTol );
}

} // namespace id

template<typename F> 
inline void
ID
( const Matrix<F>& A, Matrix<Int>& pPerm, Matrix<F>& Z, 
  Int maxSteps, Base<F> tol )
{
    DEBUG_ONLY(CallStackEntry cse("ID"))
    Matrix<F> B( A );
    id::BusingerGolub( B, pPerm, Z, maxSteps, tol );
}

template<typename F> 
inline void
ID
( Matrix<F>& A, Matrix<Int>& pPerm, Matrix<F>& Z, 
  Int maxSteps, Base<F> tol, bool canOverwrite=false )
{
    DEBUG_ONLY(CallStackEntry cse("ID"))
    Matrix<F> B;
    if( canOverwrite )
        View( B, A );
    else
        B = A;
    id::BusingerGolub( B, pPerm, Z, maxSteps, tol );
}

template<typename F> 
inline void
ID( const Matrix<F>& A, Matrix<Int>& pPerm, Matrix<F>& Z, Int numSteps )
{
    DEBUG_ONLY(CallStackEntry cse("ID"))
    ID( A, pPerm, Z, numSteps, Base<F>(-1) );
}

template<typename F> 
inline void
ID
( Matrix<F>& A, Matrix<Int>& pPerm, Matrix<F>& Z, Int numSteps, 
  bool canOverwrite=false )
{
    DEBUG_ONLY(CallStackEntry cse("ID"))
    ID( A, pPerm, Z, numSteps, Base<F>(-1), canOverwrite );
}

template<typename F,Dist UPerm> 
inline void
ID
( const DistMatrix<F>& A, 
        DistMatrix<Int,UPerm,STAR>& pPerm, 
        DistMatrix<F,STAR,VR>& Z, 
  Int maxSteps, Base<F> tol )
{
    DEBUG_ONLY(CallStackEntry cse("ID"))
    DistMatrix<F> B( A );
    id::BusingerGolub( B, pPerm, Z, maxSteps, tol );
}

template<typename F,Dist UPerm> 
inline void
ID
( DistMatrix<F>& A, 
  DistMatrix<Int,UPerm,STAR>& pPerm, 
  DistMatrix<F,STAR,VR>& Z, 
  Int maxSteps, Base<F> tol, bool canOverwrite=false )
{
    DEBUG_ONLY(CallStackEntry cse("ID"))
    DistMatrix<F> B( A.Grid() );
    if( canOverwrite )
        View( B, A );
    else
        B = A;
    id::BusingerGolub( B, pPerm, Z, maxSteps, tol );
}

template<typename F,Dist UPerm> 
inline void
ID
( const DistMatrix<F>& A, 
        DistMatrix<Int,UPerm,STAR>& pPerm, 
        DistMatrix<F,STAR,VR>& Z, 
  Int numSteps )
{
    DEBUG_ONLY(CallStackEntry cse("ID"))
    ID( A, pPerm, Z, numSteps, Base<F>(-1) );
}

template<typename F,Dist UPerm> 
inline void
ID
( DistMatrix<F>& A, 
  DistMatrix<Int,UPerm,STAR>& pPerm, 
  DistMatrix<F,STAR,VR>& Z, 
  Int numSteps, bool canOverwrite=false )
{
    DEBUG_ONLY(CallStackEntry cse("ID"))
    ID( A, pPerm, Z, numSteps, Base<F>(-1), canOverwrite );
}

} // namespace elem

#endif // ifndef ELEM_ID_HPP
