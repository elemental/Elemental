/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HERMITIANGENDEFINITEEIG_HPP
#define LAPACK_HERMITIANGENDEFINITEEIG_HPP

#ifndef WITHOUT_PMRRR

#include "elemental/lapack-like/Cholesky.hpp"
#include "elemental/lapack-like/HermitianEig.hpp"

namespace elem {

// Grab the full set of eigenpairs.
template<typename R> 
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<R>& A,
  DistMatrix<R>& B,
  DistMatrix<R,VR,STAR>& w,
  DistMatrix<R>& X )
{
#ifndef RELEASE
    PushCallStack("HermitianGenDefiniteEig");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X );
    if( type == AXBX || type == ABX )
    {
        if( uplo == LOWER )
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, R(1), B, X );
        else
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, R(1), B, X );
    }
    else /* type == BAX */
    {
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, NORMAL, NON_UNIT, R(1), B, X );
        else
            Trmm( LEFT, UPPER, ADJOINT, NON_UNIT, R(1), B, X );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Grab a partial set of eigenpairs. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
template<typename R>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<R>& A,
  DistMatrix<R>& B,
  DistMatrix<R,VR,STAR>& w,
  DistMatrix<R>& X,
  int a, int b )
{
#ifndef RELEASE
    PushCallStack("HermitianGenDefiniteEig");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, a, b );
    if( type == AXBX || type == ABX )
    {
        if( uplo == LOWER )
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, R(1), B, X );
        else
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, R(1), B, X );
    }
    else /* type == BAX */
    {
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, NORMAL, NON_UNIT, R(1), B, X );
        else
            Trmm( LEFT, UPPER, ADJOINT, NON_UNIT, R(1), B, X );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Grab a partial set of eigenpairs.
// The partial set is determined by the half-open interval (a,b]
template<typename R>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<R>& A,
  DistMatrix<R>& B,
  DistMatrix<R,VR,STAR>& w,
  DistMatrix<R>& X,
  R a, R b )
{
#ifndef RELEASE
    PushCallStack("HermitianGenDefiniteEig");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, a, b );
    if( type == AXBX || type == ABX )
    {
        if( uplo == LOWER )
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, R(1), B, X );
        else
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, R(1), B, X );
    }
    else /* type == BAX */
    {
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, NORMAL, NON_UNIT, R(1), B, X );
        else
            Trmm( LEFT, UPPER, ADJOINT, NON_UNIT, R(1), B, X );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Grab the full set of eigenvalues
template<typename R>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<R>& A,
  DistMatrix<R>& B,
  DistMatrix<R,VR,STAR>& w )
{
#ifndef RELEASE
    PushCallStack("HermitianGenDefiniteEig");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w );
#ifndef RELEASE
    PopCallStack();
#endif
}

// Grab a partial set of eigenvalues. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
template<typename R>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<R>& A,
  DistMatrix<R>& B,
  DistMatrix<R,VR,STAR>& w,
  int a, int b )
{
#ifndef RELEASE
    PushCallStack("HermitianGenDefiniteEig");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}

// Grab a partial set of eigenvalues.
// The partial set is determined by the half-open interval (a,b]
template<typename R> 
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<R>& A,
  DistMatrix<R>& B,
  DistMatrix<R,VR,STAR>& w,
  R a, R b )
{
#ifndef RELEASE
    PushCallStack("HermitianGenDefiniteEig");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}

// Grab the full set of eigenpairs
template<typename R>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& B,
  DistMatrix<R,VR,STAR>& w,
  DistMatrix<Complex<R> >& X )
{
#ifndef RELEASE
    PushCallStack("HermitianGenDefiniteEig");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X );
    if( type == AXBX || type == ABX )
    {
        if( uplo == LOWER )
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, Complex<R>(1), B, X );
        else
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, Complex<R>(1), B, X );
    }
    else /* type == BAX */
    {
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, NORMAL, NON_UNIT, Complex<R>(1), B, X );
        else
            Trmm( LEFT, UPPER, ADJOINT, NON_UNIT, Complex<R>(1), B, X );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Grab a partial set of eigenpairs. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
template<typename R>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& B,
  DistMatrix<R,VR,STAR>& w,
  DistMatrix<Complex<R> >& X,
  int a, int b )
{
#ifndef RELEASE
    PushCallStack("HermitianGenDefiniteEig");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, a, b );
    if( type == AXBX || type == ABX )
    {
        if( uplo == LOWER )
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, Complex<R>(1), B, X );
        else
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, Complex<R>(1), B, X );
    }
    else /* type == BAX */
    {
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, NORMAL, NON_UNIT, Complex<R>(1), B, X );
        else
            Trmm( LEFT, UPPER, ADJOINT, NON_UNIT, Complex<R>(1), B, X );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Grab a partial set of eigenpairs.
// The partial set is determined by the half-open interval (a,b]
template<typename R> 
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& B,
  DistMatrix<R,VR,STAR>& w,
  DistMatrix<Complex<R> >& X,
  R a, R b )
{
#ifndef RELEASE
    PushCallStack("HermitianGenDefiniteEig");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, a, b );
    if( type == AXBX || type == ABX )
    {
        if( uplo == LOWER )
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, Complex<R>(1), B, X );
        else
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, Complex<R>(1), B, X );
    }
    else /* type == BAX */
    {
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, NORMAL, NON_UNIT, Complex<R>(1), B, X );
        else
            Trmm( LEFT, UPPER, ADJOINT, NON_UNIT, Complex<R>(1), B, X );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Grab a full set of eigenvalues
template<typename R>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& B,
  DistMatrix<        R, VR,STAR>& w )
{
#ifndef RELEASE
    PushCallStack("HermitianGenDefiniteEig");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w );
#ifndef RELEASE
    PopCallStack();
#endif
}

// Grab a partial set of eigenvalues. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
template<typename R>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& B,
  DistMatrix<R,VR,STAR>& w,
  int a, int b )
{
#ifndef RELEASE
    PushCallStack("HermitianGenDefiniteEig");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}

// Grab a partial set of eigenvalues.
// The partial set is determined by the half-open interval (a,b]
template<typename R>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& B,
  DistMatrix<R,VR,STAR>& w,
  R a, R b )
{
#ifndef RELEASE
    PushCallStack("HermitianGenDefiniteEig");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // WITHOUT_PMRRR

#endif // ifndef LAPACK_HERMITIANGENDEFINITEEIG_HPP
