/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

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
            Matrix<FLower> BLower, RLower, vLower;
            Copy( B, BLower );
            Copy( R, RLower );
            EnumCtrl<RealLower> ctrlLower = ctrl;
            RealLower result =
              ShortVectorEnumeration
              ( BLower, RLower, RealLower(normUpperBound), vLower, ctrlLower );
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

    auto d = GetDiagonal( R );
    auto N( R );
    DiagonalSolve( LEFT, NORMAL, d, N );

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
        Matrix<F> BNew, U;

        for( Int trial=0; trial<ctrl.numTrials; ++trial )
        {
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

            auto dNew = GetDiagonal( RNew );
            auto NNew( RNew );
            DiagonalSolve( LEFT, NORMAL, dNew, NNew );

            if( ctrl.progress )
                Output("Starting trial ",trial);
            if( ctrl.time )
                timer.Start();
            Real result =
              svp::GNREnumeration( dNew, NNew, upperBounds, v, ctrl );
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
                        svp::CoordinatesToSparse( NNew, v, y );
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

        if( ctrl.progress )
            Output("Starting YSPARSE_ENUM(",n,")");
        if( ctrl.time )
            timer.Start();
        Real result = svp::PhaseEnumeration
          ( B, d, N, normUpperBound, startIndex, phaseLength,
            maxInfNorms, maxOneNorms, v, ctrl.progressLevel );
        if( ctrl.time )
            Output("YSPARSE_ENUM(",n,"): ",timer.Stop()," seconds");
        return result;
    }
    else
    {
        Matrix<Real> upperBounds;
        Zeros( upperBounds, n, 1 );
        Fill( upperBounds, normUpperBound );
        if( ctrl.progress )
            Output("Starting FULL_ENUM(",n,")");
        if( ctrl.time )
            timer.Start();
        Real result = svp::GNREnumeration( d, N, upperBounds, v, ctrl );
        if( ctrl.time )
            Output("FULL_ENUM(",n,"): ",timer.Stop()," seconds");
        return result;
    }
}

template<typename F>
std::pair<Base<F>,Int>
MultiShortVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
  const Matrix<Base<F>>& normUpperBounds,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("MultiShortVectorEnumeration"))
    typedef Base<F> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    v.Resize( n, 1 );

    // TODO: Guarantee normUpperBounds are positive?

    if( n == 0 )
        return std::pair<Real,Int>(Real(0),0);

    const Real BOneNorm = OneNorm( B );
    const unsigned neededPrec = unsigned(Ceil(Log2(BOneNorm)*ctrl.fudge));
    if( MantissaIsLonger<Real,float>::value &&
        MantissaBits<float>::value >= neededPrec )
    {
        typedef float RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        try
        {
            Matrix<FLower> BLower, RLower, vLower;
            Matrix<RealLower> normUpperBoundsLower;
            Copy( B, BLower );
            Copy( R, RLower );
            Copy( normUpperBounds, normUpperBoundsLower );
            EnumCtrl<RealLower> ctrlLower = ctrl;
            auto result =
              MultiShortVectorEnumeration
              ( BLower, RLower, normUpperBoundsLower, vLower, ctrlLower );
            Copy( vLower, v );
            return std::pair<Real,Int>(Real(result.first),result.second);
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
            Matrix<FLower> BLower, RLower, vLower;
            Matrix<RealLower> normUpperBoundsLower;
            Copy( B, BLower );
            Copy( R, RLower );
            Copy( normUpperBounds, normUpperBoundsLower );
            EnumCtrl<RealLower> ctrlLower = ctrl;
            auto result =
              MultiShortVectorEnumeration
              ( BLower, RLower, normUpperBoundsLower, vLower, ctrlLower );
            Copy( vLower, v );
            return std::pair<Real,Int>(Real(result.first),result.second);
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
            Matrix<FLower> BLower, RLower, vLower;
            Matrix<RealLower> normUpperBoundsLower;
            Copy( B, BLower );
            Copy( R, RLower );
            Copy( normUpperBounds, normUpperBoundsLower );
            EnumCtrl<RealLower> ctrlLower = ctrl;
            auto result =
              MultiShortVectorEnumeration
              ( BLower, RLower, normUpperBoundsLower, vLower, ctrlLower );
            Copy( vLower, v );
            return std::pair<Real,Int>(Real(result.first),result.second);
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
            Matrix<FLower> BLower, RLower, vLower;
            Matrix<RealLower> normUpperBoundsLower;
            Copy( B, BLower );
            Copy( R, RLower );
            Copy( normUpperBounds, normUpperBoundsLower );
            EnumCtrl<RealLower> ctrlLower = ctrl;
            auto result =
              MultiShortVectorEnumeration
              ( BLower, RLower, normUpperBoundsLower, vLower, ctrlLower );
            Copy( vLower, v );
            return std::pair<Real,Int>(Real(result.first),result.second);
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
            Matrix<FLower> BLower, RLower, vLower;
            Matrix<RealLower> normUpperBoundsLower;
            Copy( B, BLower );
            Copy( R, RLower );
            Copy( normUpperBounds, normUpperBoundsLower );
            EnumCtrl<RealLower> ctrlLower = ctrl;
            auto result =
              MultiShortVectorEnumeration
              ( BLower, RLower, normUpperBoundsLower, vLower, ctrlLower );
            Copy( vLower, v );
            return std::pair<Real,Int>(Real(result.first),result.second);
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
#endif
    // TODO: Arbitrary-precision drop?

    const Int numNested = normUpperBounds.Height();
    auto modNormUpperBounds( normUpperBounds );
    for( Int j=0; j<numNested; ++j )
    {
        const Real bProjNorm = R.Get(j,j);
        if( bProjNorm < normUpperBounds.Get(j,0) )
            modNormUpperBounds.Set(j,0,bProjNorm);
    }

    Timer timer;

    auto d = GetDiagonal( R );
    auto N( R );
    DiagonalSolve( LEFT, NORMAL, d, N );

    if( ctrl.enumType == GNR_ENUM )
    {
        // GNR enumeration does not yet support multi-enumeration
        const Real normUpperBound = modNormUpperBounds.Get(0,0);

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
        Matrix<F> BNew, U;

        for( Int trial=0; trial<ctrl.numTrials; ++trial )
        {
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

            auto dNew = GetDiagonal( RNew );
            auto NNew( RNew );
            DiagonalSolve( LEFT, NORMAL, dNew, NNew );

            if( ctrl.progress )
                Output("Starting trial ",trial);
            if( ctrl.time )
                timer.Start();
            Real result =
              svp::GNREnumeration( dNew, NNew, upperBounds, v, ctrl );
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
                        svp::CoordinatesToSparse( NNew, v, y );
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
                return std::pair<Real,Int>(result,0);
            }
        }
        for( Int j=0; j<numNested; ++j )
        {
            if( modNormUpperBounds.Get(j,0) < normUpperBounds.Get(j,0) )
            {
                Zeros( v, n-j, 1 );
                v.Set( 0, 0, F(1) );
                return std::pair<Real,Int>(modNormUpperBounds.Get(j,0),j);
            }
        }
        return std::pair<Real,Int>(2*normUpperBound+1,0);
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

        if( ctrl.progress )
            Output("Starting YSPARSE_ENUM(",n,")");
        if( ctrl.time )
            timer.Start();
        auto result = svp::PhaseEnumeration
          ( B, d, N, modNormUpperBounds, startIndex, phaseLength,
            maxInfNorms, maxOneNorms, v, ctrl.progressLevel );
        if( ctrl.time )
            Output("YSPARSE_ENUM(",n,"): ",timer.Stop()," seconds");
        return result;
    }
    else
    {
        // Full enumeration does not (yet) support multi-enumeration
        const Real normUpperBound = modNormUpperBounds.Get(0,0);

        Matrix<Real> upperBounds;
        Zeros( upperBounds, n, 1 );
        Fill( upperBounds, normUpperBound );
        if( ctrl.progress )
            Output("Starting FULL_ENUM(",n,")");
        if( ctrl.time )
            timer.Start();
        Real result = svp::GNREnumeration( d, N, upperBounds, v, ctrl );
        if( ctrl.time )
            Output("FULL_ENUM(",n,"): ",timer.Stop()," seconds");

        if( result < normUpperBound )
        {
            return std::pair<Real,Int>(result,0);
        }
        else
        {
            for( Int j=0; j<numNested; ++j )
            {
                if( modNormUpperBounds.Get(j,0) < normUpperBounds.Get(j,0) )
                {
                    Zeros( v, n-j, 1 );
                    v.Set( 0, 0, F(1) );
                    return std::pair<Real,Int>(modNormUpperBounds.Get(j,0),j);
                }
            }
            return std::pair<Real,Int>(result,0);
        }
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

template<typename F>
Base<F> ShortestVectorEnrichment
(       Matrix<F>& B,
  const Matrix<F>& R,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("ShortestVectorEnrichment"))
    typedef Base<F> Real;
    const Int n = B.Width();
    v.Resize( n, 1 );
    if( n == 0 )
        return Real(0);

    const Real normUpperBound = R.Get(0,0);
    return ShortestVectorEnrichment( B, R, normUpperBound, v, ctrl );
}

template<typename F>
Base<F> ShortestVectorEnrichment
(       Matrix<F>& B,
        Matrix<F>& U,
  const Matrix<F>& R,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("ShortestVectorEnrichment"))
    typedef Base<F> Real;
    const Int n = B.Width();
    v.Resize( n, 1 );
    if( n == 0 )
        return Real(0);

    const Real normUpperBound = R.Get(0,0);
    return ShortestVectorEnrichment( B, U, R, normUpperBound, v, ctrl );
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
            Matrix<FLower> BLower, RLower, vLower;
            Copy( B, BLower );
            Copy( R, RLower );
            EnumCtrl<RealLower> ctrlLower = ctrl;
            RealLower result =
              ShortestVectorEnumeration
              ( BLower, RLower, RealLower(normUpperBound), vLower, ctrlLower );
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
            Matrix<FLower> BLower, RLower, vLower;
            Copy( B, BLower );
            Copy( R, RLower );
            EnumCtrl<RealLower> ctrlLower = ctrl;
            RealLower result =
              ShortestVectorEnumeration
              ( BLower, RLower, RealLower(normUpperBound), vLower, ctrlLower );
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
            Matrix<FLower> BLower, RLower, vLower;
            Copy( B, BLower );
            Copy( R, RLower );
            EnumCtrl<RealLower> ctrlLower = ctrl;
            RealLower result =
              ShortestVectorEnumeration
              ( BLower, RLower, RealLower(normUpperBound), vLower, ctrlLower );
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
            Matrix<FLower> BLower, RLower, vLower;
            Copy( B, BLower );
            Copy( R, RLower );
            EnumCtrl<RealLower> ctrlLower = ctrl;
            RealLower result =
              ShortestVectorEnumeration
              ( BLower, RLower, RealLower(normUpperBound), vLower, ctrlLower );
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
            Matrix<FLower> BLower, RLower, vLower;
            Copy( B, BLower );
            Copy( R, RLower );
            EnumCtrl<RealLower> ctrlLower = ctrl;
            RealLower result =
              ShortestVectorEnumeration
              ( BLower, RLower, RealLower(normUpperBound), vLower, ctrlLower );
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
            // Y-sparse enumeration does not benefit from repetition
            if( ctrl.enumType == YSPARSE_ENUM ) 
                return result;
        }
        else if( satisfiedBound )
            return targetNorm;
        else
        {
            Zeros( v, n, 1 );
            v.Set( 0, 0, F(1) );
            return b0Norm;
        }
    }
}

template<typename F>
std::pair<Base<F>,Int>
MultiShortestVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
  const Matrix<Base<F>>& normUpperBounds,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("MultiShortestVectorEnumeration"))
    typedef Base<F> Real;
    const Int n = B.Width();
    v.Resize( n, 1 );
    if( n == 0 )
        return std::pair<Real,Int>(Real(0),0);

    const Real BOneNorm = OneNorm( B );
    const unsigned neededPrec = unsigned(Ceil(Log2(BOneNorm)*ctrl.fudge)); 
    if( MantissaIsLonger<Real,float>::value &&
        MantissaBits<float>::value >= neededPrec )
    {
        typedef float RealLower;
        typedef ConvertBase<F,RealLower> FLower;
        try
        {
            Matrix<FLower> BLower, RLower, vLower;
            Matrix<RealLower> normUpperBoundsLower;
            Copy( B, BLower );
            Copy( R, RLower );
            Copy( normUpperBounds, normUpperBoundsLower );
            EnumCtrl<RealLower> ctrlLower( ctrl );
            auto result =
              MultiShortestVectorEnumeration
              ( BLower, RLower, normUpperBoundsLower, vLower, ctrlLower );
            Copy( vLower, v );
            return std::pair<Real,Int>(Real(result.first),result.second);
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
            Matrix<FLower> BLower, RLower, vLower;
            Matrix<RealLower> normUpperBoundsLower;
            Copy( B, BLower );
            Copy( R, RLower );
            Copy( normUpperBounds, normUpperBoundsLower );
            EnumCtrl<RealLower> ctrlLower( ctrl );
            auto result =
              MultiShortestVectorEnumeration
              ( BLower, RLower, normUpperBoundsLower, vLower, ctrlLower );
            Copy( vLower, v );
            return std::pair<Real,Int>(Real(result.first),result.second);
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
            Matrix<FLower> BLower, RLower, vLower;
            Matrix<RealLower> normUpperBoundsLower;
            Copy( B, BLower );
            Copy( R, RLower );
            Copy( normUpperBounds, normUpperBoundsLower );
            EnumCtrl<RealLower> ctrlLower( ctrl );
            auto result =
              MultiShortestVectorEnumeration
              ( BLower, RLower, normUpperBoundsLower, vLower, ctrlLower );
            Copy( vLower, v );
            return std::pair<Real,Int>(Real(result.first),result.second);
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
            Matrix<FLower> BLower, RLower, vLower;
            Matrix<RealLower> normUpperBoundsLower;
            Copy( B, BLower );
            Copy( R, RLower );
            Copy( normUpperBounds, normUpperBoundsLower );
            EnumCtrl<RealLower> ctrlLower( ctrl );
            auto result =
              MultiShortestVectorEnumeration
              ( BLower, RLower, normUpperBoundsLower, vLower, ctrlLower );
            Copy( vLower, v );
            return std::pair<Real,Int>(Real(result.first),result.second);
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
            Matrix<FLower> BLower, RLower, vLower;
            Matrix<RealLower> normUpperBoundsLower;
            Copy( B, BLower );
            Copy( R, RLower );
            Copy( normUpperBounds, normUpperBoundsLower );
            EnumCtrl<RealLower> ctrlLower( ctrl );
            auto result =
              MultiShortestVectorEnumeration
              ( BLower, RLower, normUpperBoundsLower, vLower, ctrlLower );
            Copy( vLower, v );
            return std::pair<Real,Int>(Real(result.first),result.second);
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
#endif
    // TODO: Arbitrary-precision drop?

    // TODO: Check if later vectors are shorter too
    const Real b0Norm = R.Get(0,0);
    Zeros( v, n, 1 );
    v.Set( 0, 0, F(1) );

    const Int numNested = normUpperBounds.Height();
    Matrix<Real> targetNorms( normUpperBounds );
  
    bool satisfiedBound = false;
    Int satisfiedIndex = -1;
    for( Int j=0; j<numNested; ++j )
    {
        if( R.Get(j,j) <= normUpperBounds.Get(j,0) )
        {
            satisfiedBound = true;
            satisfiedIndex = j;
            targetNorms.Set( j, 0, R.Get(j,j) );
        }
    }

    while( true )
    {
        Matrix<F> vCand;
        auto result =
          MultiShortVectorEnumeration( B, R, targetNorms, vCand, ctrl );
        const Real normCand = result.first;
        const Int indexCand = result.second;

        if( normCand < targetNorms.Get(indexCand,0) )
        {
            v = vCand;
            targetNorms.Set(indexCand,0,normCand);
            satisfiedBound = true;
            satisfiedIndex = indexCand;
            // Y-sparse enumeration does not benefit from repetition
            if( ctrl.enumType == YSPARSE_ENUM ) 
                return result;
        }
        else if( satisfiedBound )
            return std::pair<Real,Int>
                   (targetNorms.Get(satisfiedIndex,0),satisfiedIndex);
        else
        {
            Zeros( v, n, 1 );
            v.Set( 0, 0, F(1) );
            return std::pair<Real,Int>(b0Norm,0);
        }
    }
}

template<typename F>
Base<F> ShortestVectorEnrichment
(       Matrix<F>& B,
  const Matrix<F>& R,
        Base<F> normUpperBound,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("ShortestVectorEnrichment"))
    typedef Base<F> Real;
    const Int n = B.Width();
    v.Resize( n, 1 );
    if( n == 0 )
        return Real(0);

    Real retNorm = ShortestVectorEnumeration( B, R, normUpperBound, v, ctrl );
    if( retNorm < normUpperBound )
        EnrichLattice( B, v );
    return retNorm;
}

template<typename F>
std::pair<Base<F>,Int>
MultiShortestVectorEnrichment
(       Matrix<F>& B,
  const Matrix<F>& R,
  const Matrix<Base<F>>& normUpperBounds,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("MultiShortestVectorEnrichment"))
    typedef Base<F> Real;
    const Int n = B.Width();
    v.Resize( n, 1 );
    if( n == 0 )
        return std::pair<Real,Int>(Real(0),0);

    auto ret = MultiShortestVectorEnumeration( B, R, normUpperBounds, v, ctrl );
    const Real retNorm = ret.first;
    const Int retIndex = ret.second;

    if( retNorm < normUpperBounds.Get(retIndex,0) )
    {
        auto BSub = B(ALL,IR(retIndex,END));
        EnrichLattice( BSub, v );
    }

    return ret;
}

template<typename F>
Base<F> ShortestVectorEnrichment
(       Matrix<F>& B,
        Matrix<F>& U,
  const Matrix<F>& R,
        Base<F> normUpperBound,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("ShortestVectorEnrichment"))
    typedef Base<F> Real;
    const Int n = B.Width();
    v.Resize( n, 1 );
    if( n == 0 )
        return Real(0);

    Real retNorm = ShortestVectorEnumeration( B, R, normUpperBound, v, ctrl );
    if( retNorm < normUpperBound )
        EnrichLattice( B, U, v );
    return retNorm;
}

template<typename F>
std::pair<Base<F>,Int>
MultiShortestVectorEnrichment
(       Matrix<F>& B,
        Matrix<F>& U,
  const Matrix<F>& R,
  const Matrix<Base<F>>& normUpperBounds,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("MultiShortestVectorEnrichment"))
    typedef Base<F> Real;
    const Int n = B.Width();
    v.Resize( n, 1 );
    if( n == 0 )
        return std::pair<Real,Int>(Real(0),0);

    auto ret = MultiShortestVectorEnumeration( B, R, normUpperBounds, v, ctrl );
    const Real retNorm = ret.first;
    const Int retIndex = ret.second;

    if( retNorm < normUpperBounds.Get(retIndex,0) )
    {
        const Range<Int> subInd(retIndex,END);
        auto BSub = B( ALL, subInd );
        auto USub = U( ALL, subInd );
        EnrichLattice( BSub, USub, v );
    }

    return ret;
}

// TODO: Instantiate batched variants?
#define PROTO(F) \
  template Base<F> ShortVectorEnumeration \
  ( const Matrix<F>& B, \
    const Matrix<F>& R, \
          Base<F> normUpperBound, \
          Matrix<F>& v, \
    const EnumCtrl<Base<F>>& ctrl ); \
  template std::pair<Base<F>,Int> MultiShortVectorEnumeration \
  ( const Matrix<F>& B, \
    const Matrix<F>& R, \
    const Matrix<Base<F>>& normUpperBounds, \
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
    const EnumCtrl<Base<F>>& ctrl ); \
  template std::pair<Base<F>,Int> MultiShortestVectorEnumeration \
  ( const Matrix<F>& B, \
    const Matrix<F>& R, \
    const Matrix<Base<F>>& normUpperBounds, \
          Matrix<F>& v, \
    const EnumCtrl<Base<F>>& ctrl ); \
  template Base<F> ShortestVectorEnrichment \
  (       Matrix<F>& B, \
    const Matrix<F>& R, \
          Matrix<F>& v, \
    const EnumCtrl<Base<F>>& ctrl ); \
  template Base<F> ShortestVectorEnrichment \
  (       Matrix<F>& B, \
          Matrix<F>& U, \
    const Matrix<F>& R, \
          Matrix<F>& v, \
    const EnumCtrl<Base<F>>& ctrl ); \
  template Base<F> ShortestVectorEnrichment \
  (       Matrix<F>& B, \
    const Matrix<F>& R, \
          Base<F> normUpperBound, \
          Matrix<F>& v, \
    const EnumCtrl<Base<F>>& ctrl ); \
  template std::pair<Base<F>,Int> MultiShortestVectorEnrichment \
  (       Matrix<F>& B, \
    const Matrix<F>& R, \
    const Matrix<Base<F>>& normUpperBounds, \
          Matrix<F>& v, \
    const EnumCtrl<Base<F>>& ctrl ); \
  template Base<F> ShortestVectorEnrichment \
  (       Matrix<F>& B, \
          Matrix<F>& U, \
    const Matrix<F>& R, \
          Base<F> normUpperBound, \
          Matrix<F>& v, \
    const EnumCtrl<Base<F>>& ctrl ); \
  template std::pair<Base<F>,Int> MultiShortestVectorEnrichment \
  (       Matrix<F>& B, \
          Matrix<F>& U, \
    const Matrix<F>& R, \
    const Matrix<Base<F>>& normUpperBounds, \
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
