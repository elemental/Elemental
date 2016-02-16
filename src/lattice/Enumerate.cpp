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

template<typename F>
void CoordinatesToSparse( const Matrix<F>& R, const Matrix<F>& v, Matrix<F>& y )
{
    DEBUG_ONLY(CSE cse("svp::CoordinatesToSparse"))
    const Int n = R.Height();
    Zeros( y, n, 1 );
    Gemv( NORMAL, F(1), R, v, F(0), y );
    auto r = GetDiagonal( R );
    DiagonalSolve( LEFT, NORMAL, r, y );
    Round( y );
}

template<typename F>
void SparseToCoordinates( const Matrix<F>& R, const Matrix<F>& y, Matrix<F>& v )
{
    DEBUG_ONLY(CSE cse("svp::SparseToCoordinates"))
    typedef Base<F> Real;
    const Int n = R.Height();
    v.Resize( n, 1 );

    // A custom rounded analogue of an upper-triangular solve
    const F* RBuf = R.LockedBuffer();
    const F* yBuf = y.LockedBuffer();
          F* vBuf = v.Buffer();
    const Int RLDim = R.LDim();
    for( Int j=n-1; j>=0; --j )
    {
        Real tau = 0;
        for( Int k=j+1; k<n; ++k ) 
            tau += RBuf[j+k*RLDim]*vBuf[k];
        tau /= RBuf[j+j*RLDim];
        vBuf[j] = yBuf[j] - Round(tau);
    }
}

template<typename F>
Base<F> CoordinatesToNorm( const Matrix<F>& R, const Matrix<F>& v )
{
    DEBUG_ONLY(CSE cse("svp::CoordinatesToNorm"))
    Matrix<F> z( v );
    Trmv( UPPER, NORMAL, NON_UNIT, R, z );
    return FrobeniusNorm( z );
}

template<typename F>
Base<F> SparseToNormLowerBound( const Matrix<F>& R, const Matrix<F>& y )
{
    DEBUG_ONLY(CSE cse("svp::SparseToNormLowerBound"))
    typedef Base<F> Real;
    const Int n = R.Height();
    const F* RBuf = R.LockedBuffer();
    const F* yBuf = y.LockedBuffer();
    const Int RLDim = R.LDim();

    Real lowerBoundSquared = 0;
    for( Int j=0; j<n; ++j )
    {
        if( yBuf[j] != F(0) )
        {
            const Real arg = (Abs(yBuf[j])-Real(0.5))*RealPart(RBuf[j+j*RLDim]);
            lowerBoundSquared += Pow(arg,Real(2));
        }
    }
    return Sqrt(lowerBoundSquared);
}

template<typename F>
Base<F> SparseToNorm( const Matrix<F>& R, const Matrix<F>& y )
{
    DEBUG_ONLY(CSE cse("svp::SparseToNorm"))
    const Int n = R.Height();
    Matrix<F> v;
    SparseToCoordinates( R, y, v );
    return CoordinatesToNorm( R, v );
}

template<typename Real>
void PhaseEnumerationBottom
( const Matrix<Real>& B,
  const Matrix<Real>& R,
        Matrix<Real>& y,
        Matrix<Real>& v,
        Real& normUpperBound,
        bool& foundVector,
        Int lastPhaseBeg,
        Int beg,
        Int minInf,
        Int maxInf,
        Int minOne,
        Int maxOne,
        bool zeroSoFar )
{
    DEBUG_ONLY(CSE cse("svp::PhaseEnumerationBottom"))
    const Int m = B.Height();
    const Int n = B.Width();
          Real* yBuf = y.Buffer();

    if( beg == lastPhaseBeg && minOne == 0 && minInf == 0 )
    {
        SparseToCoordinates( R, y, v );
        const Real bNorm = CoordinatesToNorm( R, v );

        if( bNorm < normUpperBound && bNorm != Real(0) )
        {
            Output("normUpperBound=",normUpperBound,", bNorm=",bNorm);
            Print( y, "y" );
            Print( v, "v" );

            Matrix<Real> b;
            Zeros( b, m, 1 );
            Gemv( NORMAL, Real(1), B, v, Real(0), b );
            Print( b, "b" );

            normUpperBound = bNorm;
            foundVector = true;
        }
    }

    // TODO: Accept as an input?
    Real baseOneNorm = 0;
    for( Int j=lastPhaseBeg; j<beg; ++j )
        baseOneNorm += Abs(yBuf[j]);

    // TODO: Accept as an input?
    Real baseInfNorm = 0;
    for( Int j=lastPhaseBeg; j<beg; ++j )
        baseInfNorm = Max(Abs(yBuf[j]),baseInfNorm);

    for( Int i=beg; i<n; ++i )
    {
        for( Int beta=-maxInf; beta<=maxInf; ++beta )
        {
            // Save a factor of two by demanding that the first entry is +
            if( beta == 0 || (zeroSoFar && beta<0) )
                continue;

            Real newOneNorm = baseOneNorm + Abs(beta);
            Real newInfNorm = Max(baseInfNorm,Real(Abs(beta)));
            if( minOne <= newOneNorm && newOneNorm <= maxOne &&
                minInf <= newInfNorm && newInfNorm <= maxInf )
            {
                yBuf[i] = beta;

                SparseToCoordinates( R, y, v );
                const Real bNorm = CoordinatesToNorm( R, v );

                if( bNorm < normUpperBound && bNorm != Real(0) )
                {
                    Output("normUpperBound=",normUpperBound,", bNorm=",bNorm);
                    Print( y, "y" );
                    Print( v, "v" );

                    Matrix<Real> b;
                    Zeros( b, m, 1 );
                    Gemv( NORMAL, Real(1), B, v, Real(0), b );
                    Print( b, "b" );

                    normUpperBound = bNorm;
                    foundVector = true;
                }

                if( newOneNorm < maxOne && i != n-1 )
                    PhaseEnumerationBottom
                    ( B, R, y, v, normUpperBound, foundVector, lastPhaseBeg,
                      i+1, minInf, maxInf, minOne, maxOne, false );

                yBuf[i] = 0;
            }
        }
    }
}

template<typename Real>
void PhaseEnumerationInner
( const Matrix<Real>& B,
  const Matrix<Real>& R,
        Matrix<Real>& y,
        Matrix<Real>& v,
        Real& normUpperBound,
        bool& foundVector,
        Int phaseLength,
        Int phase,
        Int beg,
        Int end,
  const vector<Int>& minInfNorms,
  const vector<Int>& maxInfNorms,
  const vector<Int>& minOneNorms,
  const vector<Int>& maxOneNorms,
  bool zeroSoFar,
  Int progressLevel )
{
    DEBUG_ONLY(CSE cse("svp::PhaseEnumerationInner"))
    const Int n = R.Height();
    Real* yBuf = y.Buffer();

    const Int phaseBeg = Max(end-phaseLength,0);
    const Int nextPhaseEnd = Min(end+phaseLength,n);

    if( end >= n )
    {
        PhaseEnumerationBottom
        ( B, R, y, v, normUpperBound, foundVector, beg, beg,
          minInfNorms.back(), maxInfNorms.back(),
          minOneNorms.back(), maxOneNorms.back(),
          zeroSoFar );
        return;
    }

    if( beg == phaseBeg && minInfNorms[phase] == 0 && minOneNorms[phase] == 0 )
    {
        // This phase can be all zeroes, so move to the next phase
        PhaseEnumerationInner
        ( B, R, y, v, normUpperBound, foundVector, phaseLength,
          phase+1, end, nextPhaseEnd,
          minInfNorms, maxInfNorms, minOneNorms, maxOneNorms, zeroSoFar,
          progressLevel );
    }

    Int baseOneNorm = 0;
    for( Int j=phaseBeg; j<beg; ++j ) 
        baseOneNorm += Int(Abs(yBuf[j])); 

    Int baseInfNorm = 0;
    for( Int j=phaseBeg; j<beg; ++j ) 
        baseInfNorm = Max( baseInfNorm, Int(Abs(yBuf[j])) );

    for( Int i=beg; i<end; ++i )
    {
        for( Int beta=-maxInfNorms[phase]; beta<=maxInfNorms[phase]; ++beta )
        {
            // Save a factor of two by ensuring that the first entry is +
            if( beta == 0 || (zeroSoFar && beta<0) )
                continue;

            const Int phaseOneNorm = baseOneNorm + Int(Abs(beta));
            const Int phaseInfNorm = Max( baseInfNorm, Int(Abs(beta)) );
            if( phaseOneNorm <= maxOneNorms[phase] )
            {
                if( phase < progressLevel )
                    Output("phase ",phase,": y[",i,"]=",beta);
                
                yBuf[i] = beta;

                if( phaseOneNorm >= minOneNorms[phase] &&
                    phaseOneNorm <= maxOneNorms[phase] &&
                    phaseInfNorm >= minInfNorms[phase] )
                {
                    // Fix y[i] = beta and move to the next phase
                    PhaseEnumerationInner
                    ( B, R, y, v, normUpperBound, foundVector, phaseLength,
                      phase+1, end, nextPhaseEnd,
                      minInfNorms, maxInfNorms, minOneNorms, maxOneNorms,
                      false, progressLevel );
                }
                if( phaseOneNorm < maxOneNorms[phase] )
                {
                    // Fix y[i] = beta and move to y[i+1]
                    PhaseEnumerationInner
                    ( B, R, y, v, normUpperBound, foundVector, phaseLength,
                      phase, i+1, end,
                      minInfNorms, maxInfNorms, minOneNorms, maxOneNorms,
                      false, progressLevel );
                }

                yBuf[i] = 0;
            }
        }
    }
}

template<typename Real>
Real PhaseEnumeration
( const Matrix<Real>& B,
  const Matrix<Real>& R,
        Real normUpperBound,
        Int startIndex,
        Int phaseLength,
  const vector<Int>& maxInfNorms,
  const vector<Int>& maxOneNorms,
        Matrix<Real>& v,
        Int progressLevel )
{
    DEBUG_ONLY(CSE cse("svp::PhaseEnumeration"))
    const Int n = R.Height();
    if( n <= 1 )
        return 2*normUpperBound + 1;

    // TODO: Make starting index modifiable
    const Int numPhases = ((n-startIndex)+phaseLength-1)/phaseLength;
    if( numPhases != Int(maxInfNorms.size()) )
        LogicError("Invalid length of maxInfNorms");
    if( numPhases != Int(maxOneNorms.size()) )
        LogicError("Invalid length of maxOneNorms");

    Matrix<Real> y;
    Zeros( y, n, 1 );

    // TODO: Make these values modifiable
    vector<Int> minInfNorms(numPhases,0), minOneNorms(numPhases,0);

    // TODO: Loop and increase bands for min and max one and inf norms?

    Int phase=0;
    Int beg=startIndex;
    Int end=Min(beg+phaseLength,n);
    bool foundVector = false;
    bool zeroSoFar = true;
    PhaseEnumerationInner
    ( B, R, y, v, normUpperBound, foundVector, phaseLength, phase, beg, end,
      minInfNorms, maxInfNorms, minOneNorms, maxOneNorms, zeroSoFar,
      progressLevel );

    if( foundVector )
        return normUpperBound;
    else
        return 2*normUpperBound+1;
}

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
        Matrix<Real>& v,
  const EnumCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("svp::pruned_enum::Helper"))
    const Int m = R.Height();
    const Int n = R.Width();
    if( m < n )
        LogicError("Expected height(R) >= width(R)");

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
            ++k;
            if( k == n )
                return 2*u.Get(n-1,0)+1; // An arbitrary value > than u(n-1)
            indexBuf[k] = k; // indicate that (i,j) are not synchronized
            if( k >= lastNonzero )
            {
                if( ctrl.innerProgress )
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
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("svp::BoundedEnumeration"))
    return pruned_enum::Helper( R, u, v, ctrl );
}

} // namespace svp

// NOTE: This norm upper bound is *non-inclusive*
template<typename F>
Base<F> ShortVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
        Base<F> normUpperBound,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("ShortVectorEnumeration"))
    typedef Base<F> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    v.Resize( n, 1 );

    if( normUpperBound <= Real(0) )
        return Real(1); // This is an impossible bound to meet
    if( n == 0 )
        return Real(0);

    const Real BOneNorm = OneNorm( B );
    const unsigned neededPrec = unsigned(Ceil(Log2(BOneNorm)*ctrl.fudge));
    if( MantissaIsLonger<Real,float>::value &&
        MantissaBits<float>::value >= neededPrec )
    {
        typedef float RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        try
        {
            // TODO: Switch to read/write proxies
            Matrix<FLower> BLower, RLower, vLower;
            Copy( B, BLower );
            Copy( R, RLower );
            RealLower result =
              ShortVectorEnumeration
              ( BLower, RLower, RealLower(normUpperBound), vLower, ctrl );
            Copy( vLower, v );
            return Real(result);
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
    if( MantissaIsLonger<Real,double>::value &&
        MantissaBits<double>::value >= neededPrec )
    {
        typedef double RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        try
        {
            // TODO: Switch to read/write proxies
            Matrix<FLower> BLower, RLower, vLower;
            Copy( B, BLower );
            Copy( R, RLower );
            EnumCtrl<RealLower> ctrlLower = ctrl;
            RealLower result =
              ShortVectorEnumeration
              ( BLower, RLower, RealLower(normUpperBound), vLower,
                ctrlLower );
            Copy( vLower, v );
            return Real(result);
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
#ifdef EL_HAVE_QD
    if( MantissaIsLonger<Real,DoubleDouble>::value &&
        MantissaBits<DoubleDouble>::value >= neededPrec )
    {
        typedef DoubleDouble RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        try
        {
            // TODO: Switch to read/write proxies
            Matrix<FLower> BLower, RLower, vLower;
            Copy( B, BLower );
            Copy( R, RLower );
            EnumCtrl<RealLower> ctrlLower = ctrl;
            RealLower result =
              ShortVectorEnumeration
              ( BLower, RLower, RealLower(normUpperBound), vLower, 
                ctrlLower );
            Copy( vLower, v );
            return Real(result);
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
    if( MantissaIsLonger<Real,QuadDouble>::value &&
        MantissaBits<QuadDouble>::value >= neededPrec )
    {
        typedef QuadDouble RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        try
        {
            // TODO: Switch to read/write proxies
            Matrix<FLower> BLower, RLower, vLower;
            Copy( B, BLower );
            Copy( R, RLower );
            EnumCtrl<RealLower> ctrlLower = ctrl;
            RealLower result =
              ShortVectorEnumeration
              ( BLower, RLower, RealLower(normUpperBound), vLower,
                ctrlLower );
            Copy( vLower, v );
            return Real(result);
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
#endif
#ifdef EL_HAVE_QUAD
    if( MantissaIsLonger<Real,Quad>::value &&
        MantissaBits<Quad>::value >= neededPrec )
    {
        typedef Quad RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        try
        {
            // TODO: Switch to read/write proxies
            Matrix<FLower> BLower, RLower, vLower;
            Copy( B, BLower );
            Copy( R, RLower );
            EnumCtrl<RealLower> ctrlLower = ctrl;
            RealLower result =
              ShortVectorEnumeration
              ( BLower, RLower, RealLower(normUpperBound), vLower,
                ctrlLower );
            Copy( vLower, v );
            return Real(result);
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
#endif
    // TODO: Arbitrary-precision drop?

    const Real b0ProjNorm = R.Get(0,0);
    if( b0ProjNorm < normUpperBound )
    {
        Zeros( v, n, 1 );
        v.Set( 0, 0, F(1) );
        return b0ProjNorm; 
    }
    Timer timer;

    if( ctrl.enumType == GNR_ENUM )
    {
        Matrix<Real> upperBounds(n,1);
        if( ctrl.linearBounding )
        {
            for( Int j=0; j<n; ++j )
                upperBounds.Set( j, 0, Sqrt(Real(j+1)/Real(n))*normUpperBound );
        }
        else
        {
            // See the n=140 column of Table 1 from Yoshinori Aono's
            // "A Faster Method for Computing Gama-Nguyen-Regev's Extreme
            // Pruning Coefficients".
            Matrix<Real> s(17,1);
            s.Set( 0, 0, Real(0.1318) ); 
            s.Set( 1, 0, Real(0.1859) );
            s.Set( 2, 0, Real(0.2240) ); 
            s.Set( 3, 0, Real(0.2326) );
            s.Set( 4, 0, Real(0.2336) );
            s.Set( 5, 0, Real(0.2565) );
            s.Set( 6, 0, Real(0.2871) );
            s.Set( 7, 0, Real(0.3353) );
            s.Set( 8, 0, Real(0.3978) );
            s.Set( 9, 0, Real(0.4860) );
            s.Set( 10, 0, Real(0.5808) );
            s.Set( 11, 0, Real(0.6936) );
            s.Set( 12, 0, Real(0.8241) );
            s.Set( 13, 0, Real(0.9191) );
            s.Set( 14, 0, Real(1) );
            s.Set( 15, 0, Real(1) );
            s.Set( 16, 0, Real(1) );
            for( Int j=0; j<n; ++j )
            {
                const Real percent = Real(j+1)/Real(n);
                const Real realIndex = percent*16;
                const Int floorIndex = Int(Floor(realIndex));
                const Int ceilIndex = Int(Ceil(realIndex));
                const Real indexFrac = realIndex-floorIndex;
                DEBUG_ONLY(
                  if( ceilIndex > 16 )
                      LogicError("Invalid ceiling index of ",ceilIndex);
                )
                // TODO: Use spline instead of linear interpolation?
                const Real floorVal = s.Get(floorIndex,0);
                const Real ceilVal = s.Get(ceilIndex,0);
                const Real interp = ceilVal*indexFrac + floorVal*(1-indexFrac);
                upperBounds.Set( j, 0, Sqrt(interp)*normUpperBound );
            }
        }

        // Since we will manually build up a (weakly) pseudorandom
        // unimodular matrix so that the probabalistic enumerations traverse
        // different paths, we must keep track of the unimodular matrix so that
        //  'v' can be returned relative to the original lattice basis
        auto RNew( R );
        Matrix<F> U;

        for( Int trial=0; trial<ctrl.numTrials; ++trial )
        {
            Matrix<F> BNew, RNew, U;
            BNew = B;
            Identity( U, n, n );
            if( trial != 0 )
            {
                // Apply a small random unimodular transformation to B
                const Int numCombines = n;
                for( Int j=0; j<numCombines; ++j )
                {
                    const Int c = SampleUniform( Int(0), n );
                    const Int scale = SampleUniform( Int(-5), Int(5) );
                    if( c == j || scale == 0 )
                        continue; // if scale=-1, we could have singularity
                    if( ctrl.progress )
                        Output("  B(:,",j,") += ",scale,"*B(:,",c,")");

                    auto bj = BNew( ALL, j );
                    auto bc = BNew( ALL, c );
                    Axpy( scale, bc, bj );

                    auto uj = U(ALL,j);
                    auto uc = U(ALL,c);
                    Axpy( scale, uc, uj );
                }

                // The BKZ does not need to be particularly powerful
                BKZCtrl<Real> ctrl;
                ctrl.jumpstart = true; // accumulate into U
                ctrl.blocksize = 10;
                ctrl.recursive = false;
                ctrl.lllCtrl.recursive = false;
                if( ctrl.time )
                    timer.Start();
                BKZ( BNew, U, RNew, ctrl );
                if( ctrl.time )
                    Output("  Fix-up BKZ: ",timer.Stop()," seconds");
            }
            RNew = BNew;
            qr::ExplicitTriang( RNew ); 

            if( ctrl.progress )
                Output("Starting trial ",trial);
            if( ctrl.time )
                timer.Start();
            Real result = svp::BoundedEnumeration( RNew, upperBounds, v, ctrl );
            if( ctrl.time )
                Output("  Probabalistic enumeration: ",timer.Stop()," seconds");
            if( result < normUpperBound )
            {
                if( ctrl.progress )
                    Output("Found lattice member with norm ",result);
                if( trial > 0 )
                {
                    if( ctrl.progress )
                    {
                        Print( v, "vInner" );
                        Matrix<F> y;
                        Zeros( y, RNew.Height(), 1 );
                        Gemv( NORMAL, F(1), RNew, v, F(0), y );
                        auto rNew = GetDiagonal( RNew );
                        DiagonalSolve( LEFT, NORMAL, rNew, y );
                        Print( y, "yBef" );
                        Round( y );
                        Print( y, "y" );
                    }
                    auto vCopy( v );
                    Gemv( NORMAL, F(1), U, vCopy, F(0), v );
                }
                if( ctrl.progress )
                {
                    Matrix<F> b;
                    Zeros( b, m, 1 );
                    Gemv( NORMAL, F(1), B, v, F(0), b );
                    Print( v, "v" );
                    Print( b, "b" );
                }
                return result;
            }
        }
        return 2*normUpperBound+1; // return a value above the upper bound
    }
    else if( ctrl.enumType == YSPARSE_ENUM )
    {
        const Int phaseLength = ctrl.phaseLength;
        const Int startIndex =
          ( ctrl.customStartIndex ? ctrl.startIndex : Max(n/2-1,0) );

        const Int numPhases = ((n-startIndex)+phaseLength-1)/phaseLength;

        vector<Int> maxInfNorms(numPhases,1), maxOneNorms(numPhases,1);
        if( numPhases >= 1 ) maxOneNorms[numPhases-1] = 2;

        if( ctrl.customMaxInfNorms )
            maxInfNorms = ctrl.maxInfNorms;
        if( ctrl.customMaxOneNorms )
            maxOneNorms = ctrl.maxOneNorms;

        return svp::PhaseEnumeration
          ( B, R, normUpperBound, startIndex, phaseLength,
            maxInfNorms, maxOneNorms, v, ctrl.progressLevel );
    }
    else
    {
        Matrix<Real> upperBounds;
        Zeros( upperBounds, n, 1 );
        Fill( upperBounds, normUpperBound );
        return svp::BoundedEnumeration( R, upperBounds, v, ctrl );
    }
}

template<typename F>
Base<F> ShortestVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("ShortestVectorEnumeration"))
    typedef Base<F> Real;
    const Int n = B.Width();
    v.Resize( n, 1 );
    if( n == 0 )
        return Real(0);

    const Real normUpperBound = R.Get(0,0);
    return ShortestVectorEnumeration( B, R, normUpperBound, v, ctrl );
}

// NOTE: This norm upper bound is *inclusive* so that setting it to || b_0 ||_2
//       is always sufficient for guaranteeing a solution
template<typename F>
Base<F> ShortestVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
        Base<F> normUpperBound,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("ShortestVectorEnumeration"))
    typedef Base<F> Real;
    const Int n = B.Width();
    v.Resize( n, 1 );
    if( n == 0 )
        return Real(0);

    const Real BOneNorm = OneNorm( B );
    const unsigned neededPrec = unsigned(Ceil(Log2(BOneNorm)*ctrl.fudge)); 
    if( MantissaIsLonger<Real,float>::value &&
        MantissaBits<float>::value >= neededPrec )
    {
        typedef float RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        try
        {
            // TODO: Switch to read/write proxies
            Matrix<FLower> BLower, RLower, vLower;
            Copy( B, BLower );
            Copy( R, RLower );
            RealLower result =
              ShortestVectorEnumeration
              ( BLower, RLower, RealLower(normUpperBound), vLower, ctrl );
            Copy( vLower, v );
            return Real(result);
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
    if( MantissaIsLonger<Real,double>::value &&
        MantissaBits<double>::value >= neededPrec )
    {
        typedef double RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        try
        {
            // TODO: Switch to read/write proxies
            Matrix<FLower> BLower, RLower, vLower;
            Copy( B, BLower );
            Copy( R, RLower );
            RealLower result =
              ShortestVectorEnumeration
              ( BLower, RLower, RealLower(normUpperBound), vLower, ctrl );
            Copy( vLower, v );
            return Real(result);
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
#ifdef EL_HAVE_QD
    if( MantissaIsLonger<Real,DoubleDouble>::value &&
        MantissaBits<DoubleDouble>::value >= neededPrec )
    {
        typedef DoubleDouble RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        try
        {
            // TODO: Switch to read/write proxies
            Matrix<FLower> BLower, RLower, vLower;
            Copy( B, BLower );
            Copy( R, RLower );
            RealLower result =
              ShortestVectorEnumeration
              ( BLower, RLower, RealLower(normUpperBound), vLower, ctrl );
            Copy( vLower, v );
            return Real(result);
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
    if( MantissaIsLonger<Real,QuadDouble>::value &&
        MantissaBits<QuadDouble>::value >= neededPrec )
    {
        typedef QuadDouble RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        try
        {
            // TODO: Switch to read/write proxies
            Matrix<FLower> BLower, RLower, vLower;
            Copy( B, BLower );
            Copy( R, RLower );
            RealLower result =
              ShortestVectorEnumeration
              ( BLower, RLower, RealLower(normUpperBound), vLower, ctrl );
            Copy( vLower, v );
            return Real(result);
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
#endif
#ifdef EL_HAVE_QUAD
    if( MantissaIsLonger<Real,Quad>::value &&
        MantissaBits<Quad>::value >= neededPrec )
    {
        typedef Quad RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        try
        {
            // TODO: Switch to read/write proxies
            Matrix<FLower> BLower, RLower, vLower;
            Copy( B, BLower );
            Copy( R, RLower );
            RealLower result =
              ShortestVectorEnumeration
              ( BLower, RLower, RealLower(normUpperBound), vLower, ctrl );
            Copy( vLower, v );
            return Real(result);
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
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
          ShortVectorEnumeration( B, R, targetNorm, vCand, ctrl );
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
  template void svp::CoordinatesToSparse \
  ( const Matrix<F>& R, const Matrix<F>& v, Matrix<F>& y ); \
  template void svp::SparseToCoordinates \
  ( const Matrix<F>& R, const Matrix<F>& y, Matrix<F>& v ); \
  template Base<F> svp::CoordinatesToNorm \
  ( const Matrix<F>& R, const Matrix<F>& v ); \
  template Base<F> svp::SparseToNorm \
  ( const Matrix<F>& R, const Matrix<F>& y ); \
  template Base<F> svp::PhaseEnumeration \
  ( const Matrix<F>& B, \
    const Matrix<F>& R, \
          Base<F> normUpperBound, \
          Int startIndex, \
          Int phaseLength, \
    const vector<Int>& maxInfNorms, \
    const vector<Int>& maxOneNorms, \
          Matrix<F>& v, \
          Int progressLevel ); \
  template Base<F> svp::BoundedEnumeration \
  ( const Matrix<F>& R, \
    const Matrix<Base<F>>& u, \
          Matrix<F>& v, \
    const EnumCtrl<Base<F>>& ctrl ); \
  template Base<F> ShortVectorEnumeration \
  ( const Matrix<F>& B, \
    const Matrix<F>& R, \
          Base<F> normUpperBound, \
          Matrix<F>& v, \
    const EnumCtrl<Base<F>>& ctrl ); \
  template Base<F> ShortestVectorEnumeration \
  ( const Matrix<F>& B, \
    const Matrix<F>& R, \
          Matrix<F>& v, \
    const EnumCtrl<Base<F>>& ctrl ); \
  template Base<F> ShortestVectorEnumeration \
  ( const Matrix<F>& B, \
    const Matrix<F>& R, \
          Base<F> normUpperBound, \
          Matrix<F>& v, \
    const EnumCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
