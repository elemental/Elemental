/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

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

namespace pruned_enum {

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real Helper
( const Matrix<Real>& R,
  const Matrix<Real>& u,
        Matrix<Real>& v )
{
    DEBUG_ONLY(CSE cse("svp::pruned_enum::Helper"))
    const Int m = R.Height();
    const Int n = R.Width();
    if( m < n )
        LogicError("Expected height(R) >= width(R)");

    // TODO: Make these arguments
    const bool progress = false;
    const bool track = false;

    Matrix<Real> S;
    Zeros( S, n+1, n );

    // The 'indices' vector is indexed differently from the 'r' from
    // Gama/Nguyen/Regev, and (r_0,r_1,r_2,...,r_n)=(0,1,2,...,n) becomes
    // (-1,0,1,...,n-1).
    Matrix<Int> indices;
    Zeros( indices, n+1, 1 );
    for( Int j=0; j<=n; ++j )
        indices.Set( j, 0, j-1 );

    // Note: We maintain the norms rather than their squares
    Matrix<Real> partialNorms;
    Zeros( partialNorms, n+1, 1 );

    Zeros( v, n, 1 );
    if( n == 0 )
        return Real(0);
    v.Set( 0, 0, Real(1) );

    Matrix<Real> centers;
    Zeros( centers, n, 1 );

    Matrix<Real> jumps;
    Zeros( jumps, n, 1 );

    Int lastNonzero = 0; // -1 if all indices are zero

    const Real* uBuf = u.LockedBuffer();

    const Real* RBuf = R.LockedBuffer();
    const Int RLDim = R.LDim();

    Int* indexBuf = indices.Buffer();
    Real* partialNormBuf = partialNorms.Buffer();
    Real* vBuf = v.Buffer();
    Real* centerBuf = centers.Buffer();
    Real* jumpBuf = jumps.Buffer();

    Real* SBuf = S.Buffer();
    const Int SLDim = S.LDim();

    Int k=0;
    while( true )
    {
        Real diff = vBuf[k]-centerBuf[k];
        Real rho_k =
            lapack::SafeNorm( partialNormBuf[k+1], diff*RBuf[k+k*RLDim] );
        partialNormBuf[k] = rho_k;
        if( rho_k < uBuf[(n-1)-k] )
        {
            if( k == 0 )
            {
                // Success
                return rho_k;
            }
            else
            {
                // Move down the tree
                if( track )
                    Output("  Moving down to k=",k-1," since rho_k=",rho_k," < u(n-1-k)=",uBuf[(n-1)-k]);
                --k;
                indexBuf[k] = Max(indexBuf[k],indexBuf[k+1]);
                for( Int i=indexBuf[k+1]; i>=k+1; --i )
                    SBuf[i+k*SLDim] =
                      SBuf[(i+1)+k*SLDim] +
                      vBuf[i]*(RBuf[k+i*RLDim]/RBuf[k+k*RLDim]);
                centerBuf[k] = -SBuf[(k+1)+k*SLDim];
                vBuf[k] = Round(centerBuf[k]);
                jumpBuf[k] = Real(1);
            }
        }
        else
        {
            // Move up the tree
            if( track )
                Output("  Moving up to k=",k+1," since rho_k=",rho_k," >= u(n-1-k)=",uBuf[(n-1)-k]);
            ++k;
            if( k == n )
                return 2*u.Get(n-1,0)+1; // An arbitrary value > than u(n-1)
            indexBuf[k] = k; // indicate that (i,j) are not synchronized
            if( k >= lastNonzero )
            {
                if( progress )
                    Output("lastNonzero: ",k);
                lastNonzero = k;
                vBuf[k] += Real(1);
            }
            else
            {
                if( vBuf[k] > centerBuf[k] )
                    vBuf[k] -= jumpBuf[k];
                else
                    vBuf[k] += jumpBuf[k];
                jumpBuf[k] += Real(1);
            }
        }
    }
}

// TODO: Complex version here

} // namespace pruned_enum

template<typename F>
Base<F> BoundedEnumeration
( const Matrix<F>& R,
  const Matrix<Base<F>>& u,
        Matrix<F>& v )
{
    DEBUG_ONLY(CSE cse("svp::BoundedEnumeration"))
    return pruned_enum::Helper( R, u, v );
}

} // namespace svp

// NOTE: This norm upper bound is *non-inclusive*
template<typename F>
Base<F> ShortVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
        Base<F> normUpperBound,
        Matrix<F>& v,
  bool probabalistic )
{
    DEBUG_ONLY(CSE cse("ShortVectorEnumeration"))
    typedef Base<F> Real;
    const Int n = B.Width();
    v.Resize( n, 1 );

    if( normUpperBound <= Real(0) )
        return Real(1); // This is an impossible bound to meet
    if( n == 0 )
        return Real(0);

    const bool time = true;
    const bool progress = true;

    const Real BOneNorm = OneNorm( B );
    const Real fudge = 1.5; // TODO: Make tunable
    const unsigned neededPrec = unsigned(Ceil(Log2(BOneNorm)*fudge));
    if( MantissaIsLonger<Real,float>::value &&
        MantissaBits<float>::value >= neededPrec )
    {
        typedef float RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        // TODO: Switch to read/write proxies
        Matrix<FLower> BLower, RLower, vLower;
        Copy( B, BLower );
        Copy( R, RLower );
        RealLower result =
          ShortVectorEnumeration
          ( BLower, RLower, RealLower(normUpperBound),
            vLower, probabalistic );
        Copy( vLower, v );
        return Real(result);
    }
    if( MantissaIsLonger<Real,double>::value &&
        MantissaBits<double>::value >= neededPrec )
    {
        typedef double RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        // TODO: Switch to read/write proxies
        Matrix<FLower> BLower, RLower, vLower;
        Copy( B, BLower );
        Copy( R, RLower );
        RealLower result =
          ShortVectorEnumeration
          ( BLower, RLower, RealLower(normUpperBound),
            vLower, probabalistic );
        Copy( vLower, v );
        return Real(result);
    }
#ifdef EL_HAVE_QD
    if( MantissaIsLonger<Real,DoubleDouble>::value &&
        MantissaBits<DoubleDouble>::value >= neededPrec )
    {
        typedef DoubleDouble RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        // TODO: Switch to read/write proxies
        Matrix<FLower> BLower, RLower, vLower;
        Copy( B, BLower );
        Copy( R, RLower );
        RealLower result =
          ShortVectorEnumeration
          ( BLower, RLower, RealLower(normUpperBound),
            vLower, probabalistic );
        Copy( vLower, v );
        return Real(result);
    }
    if( MantissaIsLonger<Real,QuadDouble>::value &&
        MantissaBits<QuadDouble>::value >= neededPrec )
    {
        typedef QuadDouble RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        // TODO: Switch to read/write proxies
        Matrix<FLower> BLower, RLower, vLower;
        Copy( B, BLower );
        Copy( R, RLower );
        RealLower result =
          ShortVectorEnumeration
          ( BLower, RLower, RealLower(normUpperBound),
            vLower, probabalistic );
        Copy( vLower, v );
        return Real(result);
    }
#endif
#ifdef EL_HAVE_QUAD
    if( MantissaIsLonger<Real,Quad>::value &&
        MantissaBits<Quad>::value >= neededPrec )
    {
        typedef Quad RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        // TODO: Switch to read/write proxies
        Matrix<FLower> BLower, RLower, vLower;
        Copy( B, BLower );
        Copy( R, RLower );
        RealLower result =
          ShortVectorEnumeration
          ( BLower, RLower, RealLower(normUpperBound),
            vLower, probabalistic );
        Copy( vLower, v );
        return Real(result);
    }
#endif
    // TODO: Arbitrary-precision drop?

    auto b0 = B( ALL, IR(0) );
    const Real b0Norm = FrobeniusNorm( b0 );
    if( b0Norm < normUpperBound )
    {
        Zeros( v, n, 1 );
        v.Set( 0, 0, F(1) );
        return b0Norm; 
    }

    Timer timer;
    if( probabalistic )
    {
        // TODO: Add support for different bounding functions
        Matrix<Real> upperBounds(n,1);
        for( Int j=0; j<n; ++j )
            upperBounds.Set( j, 0, Sqrt(Real(j+1)/Real(n))*normUpperBound );

        // Since we will manually build up a (weakly) pseudorandom
        // unimodular matrix so that the probabalistic enumerations traverse
        // different paths, we must keep track of the unimodular matrix so that
        //  'v' can be returned relative to the original lattice basis
        auto BNew( B );
        auto RNew( R );
        Matrix<F> U;
        Identity( U, n, n );

        Int numTrials = 10*n; // TODO: Make this tunable; probability is 1/n
        for( Int trial=0; trial<numTrials; ++trial )
        {
            if( progress )
                Output("Starting trial ",trial);
            if( time )
                timer.Start();
            Real result = svp::BoundedEnumeration( RNew, upperBounds, v );
            if( time )
                Output("  Probabalistic enumeration: ",timer.Stop()," seconds");
            if( result < normUpperBound )
            {
                if( progress )
                    Output("Found lattice member with norm ",result);
                if( trial > 0 )
                {
                    auto vCopy( v );
                    Gemv( NORMAL, F(1), U, vCopy, F(0), v );
                }
                return result;
            }
            else
            {
                // Apply a small random unimodular transformation to B
                const Int numCombines = n;
                for( Int j=0; j<numCombines; ++j )
                {
                    const Int c = SampleUniform( Int(0), n );
                    const Int scale = SampleUniform( Int(-5), Int(5) );
                    if( c == j || scale == 0 )
                        continue; // if scale=-1, we could have singularity
                    if( progress )
                        Output("  B(:,",j,") += ",scale,"*B(:,",c,")");

                    auto bj = BNew( ALL, j );
                    auto bc = BNew( ALL, c );
                    Axpy( scale, bc, bj );

                    auto uj = U(ALL,j);
                    auto uc = U(ALL,c);
                    Axpy( scale, uc, uj );
                }

                // NOTE: The LLL does not need to be particularly powerful
                LLLCtrl<Real> ctrl;
                ctrl.jumpstart = true; 
                ctrl.startCol = 0;
                ctrl.recursive = false;
                if( time )
                    timer.Start();
                LLL( BNew, U, RNew, ctrl );
                if( time )
                    Output("  Fix-up LLL: ",timer.Stop()," seconds");
            }
        }
        return 2*normUpperBound+1; // return a value above the upper bound
    }
    else
    {
        Matrix<Real> upperBounds;
        Zeros( upperBounds, n, 1 );
        Fill( upperBounds, normUpperBound );
        return svp::BoundedEnumeration( R, upperBounds, v );
    }
}

template<typename F>
Base<F> ShortestVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
        Matrix<F>& v,
  bool probabalistic )
{
    DEBUG_ONLY(CSE cse("ShortestVectorEnumeration"))
    typedef Base<F> Real;
    const Int n = B.Width();
    v.Resize( n, 1 );
    if( n == 0 )
        return Real(0);

    const Real normUpperBound = R.Get(0,0);
    return ShortestVectorEnumeration( B, R, normUpperBound, v, probabalistic );
}

// NOTE: This norm upper bound is *inclusive* so that setting it to || b_0 ||_2
//       is always sufficient for guaranteeing a solution
template<typename F>
Base<F> ShortestVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
        Base<F> normUpperBound,
        Matrix<F>& v,
  bool probabalistic )
{
    DEBUG_ONLY(CSE cse("ShortestVectorEnumeration"))
    typedef Base<F> Real;
    const Int n = B.Width();
    v.Resize( n, 1 );
    if( n == 0 )
        return Real(0);

    const Real BOneNorm = OneNorm( B );
    const Real fudge = 1.5; // TODO: Make tunable
    const unsigned neededPrec = unsigned(Ceil(Log2(BOneNorm)*fudge)); 
    if( MantissaIsLonger<Real,float>::value &&
        MantissaBits<float>::value >= neededPrec )
    {
        typedef float RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        // TODO: Switch to read/write proxies
        Matrix<FLower> BLower, RLower, vLower;
        Copy( B, BLower );
        Copy( R, RLower );
        RealLower result =
          ShortestVectorEnumeration
          ( BLower, RLower, RealLower(normUpperBound), vLower, probabalistic );
        Copy( vLower, v );
        return Real(result);
    }
    if( MantissaIsLonger<Real,double>::value &&
        MantissaBits<double>::value >= neededPrec )
    {
        typedef double RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        // TODO: Switch to read/write proxies
        Matrix<FLower> BLower, RLower, vLower;
        Copy( B, BLower );
        Copy( R, RLower );
        RealLower result =
          ShortestVectorEnumeration
          ( BLower, RLower, RealLower(normUpperBound), vLower, probabalistic );
        Copy( vLower, v );
        return Real(result);
    }
#ifdef EL_HAVE_QD
    if( MantissaIsLonger<Real,DoubleDouble>::value &&
        MantissaBits<DoubleDouble>::value >= neededPrec )
    {
        typedef DoubleDouble RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        // TODO: Switch to read/write proxies
        Matrix<FLower> BLower, RLower, vLower;
        Copy( B, BLower );
        Copy( R, RLower );
        RealLower result =
          ShortestVectorEnumeration
          ( BLower, RLower, RealLower(normUpperBound), vLower, probabalistic );
        Copy( vLower, v );
        return Real(result);
    }
    if( MantissaIsLonger<Real,QuadDouble>::value &&
        MantissaBits<QuadDouble>::value >= neededPrec )
    {
        typedef QuadDouble RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        // TODO: Switch to read/write proxies
        Matrix<FLower> BLower, RLower, vLower;
        Copy( B, BLower );
        Copy( R, RLower );
        RealLower result =
          ShortestVectorEnumeration
          ( BLower, RLower, RealLower(normUpperBound), vLower, probabalistic );
        Copy( vLower, v );
        return Real(result);
    }
#endif
#ifdef EL_HAVE_QUAD
    if( MantissaIsLonger<Real,Quad>::value &&
        MantissaBits<Quad>::value >= neededPrec )
    {
        typedef Quad RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        // TODO: Switch to read/write proxies
        Matrix<FLower> BLower, RLower, vLower;
        Copy( B, BLower );
        Copy( R, RLower );
        RealLower result =
          ShortestVectorEnumeration
          ( BLower, RLower, RealLower(normUpperBound), vLower, probabalistic );
        Copy( vLower, v );
        return Real(result);
    }
#endif
    // TODO: Arbitrary-precision drop?

    const Real b0Norm = R.Get(0,0);
    Zeros( v, n, 1 );
    v.Set( 0, 0, F(1) );

    bool satisfiedBound = ( b0Norm <= normUpperBound ? true : false );
    Real targetNorm = Min(normUpperBound,b0Norm);

    while( true )
    {
        Matrix<F> vCand;
        Real result =
          ShortVectorEnumeration( B, R, targetNorm, vCand, probabalistic );
        if( result < targetNorm )
        {
            v = vCand;
            targetNorm = result;
            satisfiedBound = true;
        }
        else if( satisfiedBound )
            return targetNorm;
        else
            RuntimeError("Could not satisfy (inclusive) norm upper bound");
    }
}

#define PROTO(F) \
  template Base<F> svp::BoundedEnumeration \
  ( const Matrix<F>& R, \
    const Matrix<Base<F>>& u, \
          Matrix<F>& v ); \
  template Base<F> ShortVectorEnumeration \
  ( const Matrix<F>& B, \
    const Matrix<F>& R, \
          Base<F> normUpperBound, \
          Matrix<F>& v, \
    bool probabalistic ); \
  template Base<F> ShortestVectorEnumeration \
  ( const Matrix<F>& B, \
    const Matrix<F>& R, \
          Matrix<F>& v, \
    bool probabalistic ); \
  template Base<F> ShortestVectorEnumeration \
  ( const Matrix<F>& B, \
    const Matrix<F>& R, \
          Base<F> normUpperBound, \
          Matrix<F>& v, \
    bool probabalistic );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
