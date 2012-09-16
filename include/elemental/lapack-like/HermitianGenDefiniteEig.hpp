/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef WITHOUT_PMRRR

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
        TwoSidedTrsm( uplo, A, B );
    else
        TwoSidedTrmm( uplo, A, B );
    HermitianEig( uplo, A, w, X );
    if( type == AXBX || type == ABX )
    {
        if( uplo == LOWER )
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, (R)1, B, X );
        else
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (R)1, B, X );
    }
    else /* type == BAX */
    {
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, NORMAL, NON_UNIT, (R)1, B, X );
        else
            Trmm( LEFT, UPPER, ADJOINT, NON_UNIT, (R)1, B, X );
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
        TwoSidedTrsm( uplo, A, B );
    else
        TwoSidedTrmm( uplo, A, B );
    HermitianEig( uplo, A, w, X, a, b );
    if( type == AXBX || type == ABX )
    {
        if( uplo == LOWER )
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, (R)1, B, X );
        else
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (R)1, B, X );
    }
    else /* type == BAX */
    {
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, NORMAL, NON_UNIT, (R)1, B, X );
        else
            Trmm( LEFT, UPPER, ADJOINT, NON_UNIT, (R)1, B, X );
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
        TwoSidedTrsm( uplo, A, B );
    else
        TwoSidedTrmm( uplo, A, B );
    HermitianEig( uplo, A, w, X, a, b );
    if( type == AXBX || type == ABX )
    {
        if( uplo == LOWER )
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, (R)1, B, X );
        else
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (R)1, B, X );
    }
    else /* type == BAX */
    {
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, NORMAL, NON_UNIT, (R)1, B, X );
        else
            Trmm( LEFT, UPPER, ADJOINT, NON_UNIT, (R)1, B, X );
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
        TwoSidedTrsm( uplo, A, B );
    else
        TwoSidedTrmm( uplo, A, B );
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
        TwoSidedTrsm( uplo, A, B );
    else
        TwoSidedTrmm( uplo, A, B );
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
        TwoSidedTrsm( uplo, A, B );
    else
        TwoSidedTrmm( uplo, A, B );
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
        TwoSidedTrsm( uplo, A, B );
    else
        TwoSidedTrmm( uplo, A, B );
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
        TwoSidedTrsm( uplo, A, B );
    else
        TwoSidedTrmm( uplo, A, B );
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
        TwoSidedTrsm( uplo, A, B );
    else
        TwoSidedTrmm( uplo, A, B );
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
        TwoSidedTrsm( uplo, A, B );
    else
        TwoSidedTrmm( uplo, A, B );
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
        TwoSidedTrsm( uplo, A, B );
    else
        TwoSidedTrmm( uplo, A, B );
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
        TwoSidedTrsm( uplo, A, B );
    else
        TwoSidedTrmm( uplo, A, B );
    HermitianEig( uplo, A, w, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // WITHOUT_PMRRR
