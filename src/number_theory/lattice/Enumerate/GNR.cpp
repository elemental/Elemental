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
// the QR factorization of B and 'upperBounds' to denote the sequence of
// upper-bounds
//
//     upperBounds(0) <= upperBounds(1) <= ... <= upperBounds(n-1)
//
// such that upperBounds(j) bounds the two-norm of the last j+1 entries of
// R v = (diag(d) N) v.
//
// The "GNR" enumeration algorithm has been extended to support complex
// arithmetic below by reformulating the original algorithm in terms of 
// a vector of states, each of which allows for a simple traversal of a 
// discrete spiral in the complex plane that is centered about a particular
// point (with the initial direction and orientation of the spiral chosen
// based upon the difference between the center and its rounding).
//
// While it is possible to easily analytically solve for the "next closest"
// integer to a given real number relative to a given integer with O(1) work,
// and such a  solution can easily drive a walk "outward" from a center point, 
// the author is not aware of an O(1) analogue in the complex plane.
//
// The easier part of generalizing GNR enumeration is in generalizing the notion
// of traversing Z modulo multiplication by -1 (which can be represented by Z+)
// to Z^2 modulo multiplication by i. In particular, the latter can be 
// represented by 
//
//     { (a,b) in Z^2 : a > 0, |b| < a } U {(0,0)},
//
// which can be easily traversed by trying each admissible b after incrementing
// a.

// An example of traversing a clockwise Z^2 spiral which begins rightward
// (which makes the most sense if the target center was rounded up and left)
// is given here:
//
//           20- ...
//           |
//           19  6 - 7 - 8 - 9
//           |   |           |
//           18  5   0 - 1   10
//           |   |       |   |
//           17  4 - 3 - 2   11
//           |               |
//           16- 15- 14- 13- 12
//
// where the root node, marked as 0, was the complex rounding of the
// (non-integer) center. Since the legs of the spiral are of length 
//
//   1, 1, 2, 2, 3, 3, etc.
//
// we can easily traverse the spiral with O(1) state. And, while such an
// ordering of Z^2 does not precisely equate with an ordering based upon the
// distance to the center node, it is a reasonable approximation and there are
// no tie-breaking subtleties.
//

namespace gnr_enum {

template<typename Real>
struct SpiralState
{
    bool forward;
    Real jump;
    Real position;

    void Initialize( const Real& center )
    {
        position = Round(center);
        forward = (position <= center);
        jump = 1;
    }

    Real Step()
    {
        if( forward ) 
            position += jump;
        else
            position -= jump;
        forward = !forward;
        jump += 1;
        return position;
    }
};

template<typename Real>
struct SpiralState<Complex<Real>>
{
    Int legLength;
    Int numSteps;
    Int direction; // 0=right, 1=down, 2=left, 3=up
    bool clockwise;
    bool firstLeg;

    Complex<Real> position;

    void Initialize( const Complex<Real>& center )
    {
        legLength = 1;
        numSteps = 0; 
        firstLeg = true;

        const Real realCenter = center.real();
        const Real imagCenter = center.imag();
        const Real realStart = Round(realCenter);
        const Real imagStart = Round(imagCenter);
        position = Complex<Real>(realStart,imagStart);
        const Real realDist = Abs(realCenter-realStart);
        const Real imagDist = Abs(imagCenter-imagStart);
        if( realDist < imagDist )
        {
            if( realCenter > realStart )
            {
                direction = 0; // right
                clockwise = (imagCenter < imagStart);
            }
            else
            {
                direction = 2; // left
                clockwise = (imagCenter > imagStart);
            }
        }
        else
        {
            if( imagCenter > imagStart )
            {
                direction = 3; // up
                clockwise = (realCenter > realStart);
            }
            else
            {
                direction = 2; // left
                clockwise = (realCenter < realStart);
            }
        }
    }
    
    Complex<Real> Step()
    {
        Real realPos = position.real();
        Real imagPos = position.imag();
        if( direction == 0 )
            realPos += 1;
        else if( direction == 1 )
            imagPos -= 1;
        else if( direction == 2 )
            realPos -= 1;
        else if( direction == 3 )
            imagPos += 1;
            
        ++numSteps;
        if( numSteps == legLength )
        {
            numSteps = 0;

            if( clockwise )
                direction = Mod( direction+1, 4 );
            else
                direction = Mod( direction-1, 4 );

            if( firstLeg )
            {
                firstLeg = false;
            }
            else
            {
                ++legLength;
                firstLeg = true;
            }
        }

        position = Complex<Real>(realPos,imagPos);
        return position;
    }
};

template<typename Real,typename=EnableIf<IsReal<Real>>>
inline Real ConstrainedSpiral( const Real& position )
{
    // The coset Z modulo multiplication by -1 can be represented by Z+
    return position + 1;
}

template<typename F,typename=DisableIf<IsReal<F>>,typename=void>
inline F ConstrainedSpiral( const F& position )
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
    const Int m = N.Height();
    const Int n = N.Width();
    if( n > m )
        LogicError("Expected height(N) >= width(N)");

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

    vector<SpiralState<F>> spiralStates(n);

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
                spiralStates[k].Initialize( centers(k) );
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
                vBuf[k] = ConstrainedSpiral( vBuf[k] );
            }
            else
            {
                vBuf[k] = spiralStates[k].Step();
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
    const Int m = NTrans.Width();
    const Int n = NTrans.Height();
    if( n > m )
        LogicError("Expected height(N) >= width(N)");

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

    vector<SpiralState<F>> spiralStates(n);

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
                spiralStates[k].Initialize( centers(k) );
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
                vBuf[k] = ConstrainedSpiral( vBuf[k] );
            }
            else
            {
                vBuf[k] = spiralStates[k].Step();
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
