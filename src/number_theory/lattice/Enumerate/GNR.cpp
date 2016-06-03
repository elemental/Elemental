/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

namespace svp {

// See Algorithm 2 from:
//
//   Nicolas Gama, Phong Q. Nguyen, and Oded Regev,
//   "Lattice enumeration using extreme pruning", Eurocrypt 2010.
//
// Note that our algorithm uses 'R' to denote the upper triangular matrix from
// the QR factorization of B and 'u' to denote the sequence of upper-bounds
//
//     u(0)^2 <= u(1)^2 <= ... <= u(n-1)^2.
//
//
// TODO: Extend the following to handle complex lattices.

namespace gnr_enum {

template<typename Real,typename=EnableIf<IsReal<Real>>>
inline Real NextClosest( const Real& position, const Real& center )
{
    const Real dist = Abs(position-center);
    if( position > center )
    {
        const Real reflection = center - dist;
        return Floor(reflection);
    } 
    else
    {
        const Real reflection = center + dist;
        if( reflection == Ceil(reflection) )
            return reflection + 1;
        else
            return Ceil(reflection);
    }
}

template<typename F,typename=DisableIf<IsReal<F>>,typename=void>
inline F NextClosest( const F& position, const F& center )
{
    typedef Base<F> Real;
    const Real realPos = position.real();
    const Real imagPos = position.imag();
    const Real realCenter = center.real();
    const Real imagCenter = center.imag();

    const Real realDist = Abs(realPos-realCenter);
    const Real imagDist = Abs(imagPos-imagCenter);

    LogicError("This routine is not yet written");
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
inline Real WalkOutward( const Real& position )
{
    // The coset Z modulo multiplication by -1 can be represented by Z+
    return position + 1;
}

template<typename F,typename=DisableIf<IsReal<F>>,typename=void>
inline F WalkOutward( const F& position )
{
    typedef Base<F> Real;
    Real realPos = position.real();
    Real imagPos = position.imag();

    // The coset Z^2 modulo multiplication by i can be represented by the set
    //
    //    { (a,b) in Z^2 : a > 0, |b| < a } U {(0,0)}, 
    //
    // and the origin can be ignored for our purposes. We traverse this set by
    // searching each admissible value of b after incrementing a.

    if( imagPos > Real(0) )
    {
        if( imagPos == realPos-1 )
            imagPos = -1;
        else
            imagPos += 1;
    }
    else if( imagPos < Real(0) )
    {
        if( imagPos == -(realPos-1) )
        {
            realPos += 1;
            imagPos = 0;
        }
        else
        {
            imagPos -= 1;
        }
    }
    else
    {
        imagPos += 1;
    }
    return Complex<Real>(realPos,imagPos);
}

template<typename F>
Base<F> Helper
( const Matrix<Base<F>>& d,
  const Matrix<F>& N,
  const Matrix<Base<F>>& upperBounds,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("svp::gnr_enum::Helper"))
    typedef Base<F> Real;
    const Int n = N.Height();
    if( n != N.Width() )
        LogicError("Expected height(N) = width(N)");

    Matrix<F> partialSums;
    Zeros( partialSums, n+1, n );

    // The 'sumIndices' vector is indexed differently from the 'r' from
    // Gama/Nguyen/Regev, and (r_0,r_1,r_2,...,r_n)=(0,1,2,...,n) becomes
    // (-1,0,1,...,n-1).
    Matrix<Int> sumIndices;
    Zeros( sumIndices, n+1, 1 );
    for( Int j=0; j<=n; ++j )
        sumIndices(j) = j-1;

    // Note: We maintain the norms rather than their squares
    Matrix<Real> partialNorms;
    Zeros( partialNorms, n+1, 1 );

    Zeros( v, n, 1 );
    if( n == 0 )
        return Real(0);
    v(0) = F(1);
    Int lastNonzero = 0; // will be set to -1 if all sumIndices are zero

    Matrix<F> centers;
    Zeros( centers, n, 1 );

    Int k=0;
    F* vBuf = &v(0);
    while( true )
    {
        const F entry = d(k)*(vBuf[k] - centers(k));
        const Real partialNorm = lapack::SafeNorm( partialNorms(k+1), entry );
        partialNorms(k) = partialNorm;
        if( partialNorm < upperBounds((n-1)-k) )
        {
            if( k == 0 )
            {
                // Success
                return partialNorm;
            }
            else
            {
                // Move down the tree
                --k;
                sumIndices(k) = Max(sumIndices(k),sumIndices(k+1));

                F* s = &partialSums(0,k);
                for( Int i=sumIndices(k+1); i>=k+1; --i )
                    s[i] = s[i+1] + N(k,i)*vBuf[i];

                centers(k) = -partialSums(k+1,k);
                vBuf[k] = Round(centers(k));
            }
        }
        else
        {
            // Move up the tree
            ++k;
            if( k == n )
            {
                // Return an arbitrary value greater than upperBounds(n-1)
                return 2*upperBounds(n-1)+1;
            }
            sumIndices(k) = k; // indicate that (i,j) are not synchronized
            if( k == lastNonzero+1 )
            {
                // Seed the new nonzero entry
                vBuf[k] = F(1);
                lastNonzero = k;
                if( ctrl.innerProgress )
                    Output("lastNonzero: ",lastNonzero);
            }
            else if( k == lastNonzero )
            {
                vBuf[k] = WalkOutward( vBuf[k] );
            }
            else
            {
                vBuf[k] = NextClosest( vBuf[k], centers(k) );
            }
        }
    }
}

template<typename F>
Base<F> TransposedHelper
( const Matrix<Base<F>>& d,
  const Matrix<F>& NTrans,
  const Matrix<Base<F>>& upperBounds,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("svp::gnr_enum::TransposedHelper"))
    typedef Base<F> Real;
    const Int n = NTrans.Height();
    if( n != NTrans.Width() )
        LogicError("Expected height(N) = width(N)");

    Matrix<F> partialSums;
    Zeros( partialSums, n+1, n );

    // The 'sumIndices' vector is indexed differently from the 'r' from
    // Gama/Nguyen/Regev, and (r_0,r_1,r_2,...,r_n)=(0,1,2,...,n) becomes
    // (-1,0,1,...,n-1).
    Matrix<Int> sumIndices;
    Zeros( sumIndices, n+1, 1 );
    for( Int j=0; j<=n; ++j )
        sumIndices(j) = j-1;

    // Note: We maintain the norms rather than their squares
    Matrix<Real> partialNorms;
    Zeros( partialNorms, n+1, 1 );

    Zeros( v, n, 1 );
    if( n == 0 )
        return Real(0);
    v(0) = F(1);
    Int lastNonzero = 0; // -1 if all sumIndices are zero

    Matrix<F> centers;
    Zeros( centers, n, 1 );

    Int k=0;
    F* vBuf = &v(0);
    while( true )
    {
        const F entry = d(k)*(vBuf[k] - centers(k));
        const Real partialNorm = lapack::SafeNorm( partialNorms(k+1), entry );
        partialNorms(k) = partialNorm;
        if( partialNorm < upperBounds((n-1)-k) )
        {
            if( k == 0 )
            {
                // Success
                return partialNorm;
            }
            else
            {
                // Move down the tree
                --k;
                sumIndices(k) = Max(sumIndices(k),sumIndices(k+1));

                      F* s = &partialSums(0,k);
                const F* nBuf = &NTrans(0,k);
                for( Int i=sumIndices(k+1); i>=k+1; --i )
                    s[i] = s[i+1] + nBuf[i]*vBuf[i];

                centers(k) = -partialSums(k+1,k);
                vBuf[k] = Round(centers(k));
            }
        }
        else
        {
            // Move up the tree
            ++k;
            if( k == n )
            {
                // Return an arbitrary value > than upperBounds(n-1)
                return 2*upperBounds(n-1)+1;
            }
            sumIndices(k) = k; // indicate that (i,j) are not synchronized
            if( k == lastNonzero+1 )
            {
                // Seed the new nonzero entry
                vBuf[k] = F(1);
                lastNonzero = k;
                if( ctrl.innerProgress )
                    Output("lastNonzero: ",lastNonzero);
            }
            else if( k == lastNonzero )
            {
                vBuf[k] = WalkOutward( vBuf[k] );
            }
            else
            {
                vBuf[k] = NextClosest( vBuf[k], centers(k) );
            }
        }
    }
}

} // namespace gnr_enum

template<typename F>
Base<F> GNREnumeration
( const Matrix<Base<F>>& d,
  const Matrix<F>& N,
  const Matrix<Base<F>>& upperBounds,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("svp::GNREnumeration"))
    if( ctrl.explicitTranspose )
    {
        Matrix<F> NTrans;
        Transpose( N, NTrans );
        return gnr_enum::TransposedHelper( d, NTrans, upperBounds, v, ctrl );
    }
    else
    {
        return gnr_enum::Helper( d, N, upperBounds, v, ctrl );
    }
}

} // namespace svp

#define PROTO(F) \
  template Base<F> svp::GNREnumeration \
  ( const Matrix<Base<F>>& d, \
    const Matrix<F>& N, \
    const Matrix<Base<F>>& u, \
          Matrix<F>& v, \
    const EnumCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
