/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// The implementations of Householder-based LLL in this file are extensions of
// the algorithm discussed in:
//
//   Ivan Morel, Damien Stehle, and Gilles Villard,
//   "H-LLL: Using Householder inside LLL",
//   ISSAC '09, July 28--31, 2009, Seoul, Republic of Korea.
//
// Note that there is a subtle bug in Algorithm 4 of the paper, as step 7
// involves setting k := max(k-1,2), which, in the case of zero indexing,
// becomes k := max(k-1,1). The second branch of the 'max' is only activated
// in the case where k=1 and is kept at this value despite having just 
// swapped b_0 and b_1. But this would imply that the 0'th column of the 
// Householder QR factorization is now incorrect, and so it must be refreshed.
// This implementation fixes said issue via a call to lll::HouseholderStep
// with k=0.
//
// Furthermore, an analogue of the above algorithm which maintains an 
// accumulated set of Householder transformations, using a so-called 
// UT Transform,
//
//   Thierry Joffrain, Tze Meng Low, Enrique S. Quintana-Orti, and 
//   Robert van de Geijn,
//   "Accumulating Householder Transforms, Revisited",
//   ACM Transactions on Mathematical Software, Vol. 32, No. 2, pp. 169--179,
//   2006.
//
// However, it can be seen that the k'th iteration of the accumulated algorithm
// requires
//
//     6 k ( m - k/2 )
//
// operations, whereas the standard algorithm only requires
//
//     4 k ( m - k/2 )
// 
// operations. For this reason, there is not a clear benefit to using the
// accumulated algorithm for LLL in a sequential setting since the transition
// from level 1 to level 2 BLAS does not warrant the 50% increase in work.
// However, the situation might be significantly different in a distributed
// setting, where each inner product would involve a non-trivial amount of
// latency. Also, future work on 'greedy' variants of LLL may further swing
// the tide towards accumulated transforms.
//
// Also, the following paper was consulted for more general factorization
// viewpoints:
//
// Dirk Wubben, Ronald Bohnke, Volker Kuhn, and Karl-Dirk Kammeyer,
// "MMSE-Based Lattice-Reduction for Near-ML Detection of MIMO Systems",
// ITG Workshop on Smart Antennas, pp. 106--113, 2004
//
// Future work will involve investigating blocked algorithms and/or 
// distributed-memory and/or GPU implementations.
//
// The seminal work on distributed-memory implementations of LLL is
//
//   G. Villard, "Parallel lattice basis reduction",
//   Papers from the International Symposium on Symbolic and algebraic
//   computation, 1992, pp. 269--277.
//
// The key idea is that, because only |R(i,i+1)/R(i,i)| <= 1/2 is required
// in order to achieve the optimality bound of the first column of the 
// reduced basis (noticed by Lenstra et al. in "Factoring polynomials with
// rational coefficients"), one can transition to an "all swaps" algorithm
// which alternates between swapping all possible (b_i,b_{i+1}) columns with
// odd and even coefficients. The result is then only weakly LLL reduced,
// but the parallel algorithm is greatly simplified. Furthermore, the QR
// factorization can easily be patched up with an independent Givens rotation
// for each swap.
//
// Henri Cohen's "A course in computational algebraic number theory" was also
// heavily consulted for extending LLL to support linearly-dependent columns
// (and its extension to computing integer kernels).
//
// Lastly, insights from the excellent survey paper
//
//   Damien Stehle, "Floating-Point LLL: Theoretical and Practical Aspects" 
//
// are slowly being incorporated.
//

namespace El {

static Timer applyHouseTimer, roundTimer, formSInvTimer;

namespace lll {

// Return the achieved delta and eta reduction properties
template<typename F>
inline std::pair<Base<F>,Base<F>>
Achieved
( const Matrix<F>& R,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("lll::Achieved"))
    typedef Base<F> Real;
    const Int m = R.Height();
    const Int n = R.Width();
    const Int minDim = Min(m,n);

    // Find the maximum delta such that
    //
    //   delta R(k,k)^2 <= R(k+1,k+1)^2 + |R(k,k+1)|^2
    //
    // for 0 <= k < Min(m,n)-1.
    //
    // TODO: Decide if m < n requires checking the k=m-1 case.
    //
    Matrix<F> z;
    Real delta = limits::Max<Real>();
    for( Int i=0; i<minDim-1; ++i )
    {
        const Real rho_i_i = R.GetRealPart(i,i);
        if( Abs(rho_i_i) <= ctrl.zeroTol )
            continue;
        const Real rho_i_ip1 = Abs(R.Get(i,i+1));
        const Real rho_ip1_ip1 = R.GetRealPart(i+1,i+1);

        const Real deltaBound =
          (rho_ip1_ip1*rho_ip1_ip1+rho_i_ip1*rho_i_ip1)/(rho_i_i*rho_i_i);
 
        delta = Min(delta,deltaBound);
    }

    // Ensure that, for all 0 <= l < Min(m,n) and l < k < n,
    //
    //    | R(l,k) | <= phi(F) eta R(l,l),
    //
    // unless a weak reduction was requested in which case k=l+1.
    //
    // TODO: Decide if m < n requires checking the k=m-1 case.
    //
    // NOTE: phi(F) is 1 for real F and sqrt(2) for complex F.
    //
    Real eta = 0;
    if( ctrl.weak )
    {
        for( Int i=0; i<minDim-1; ++i )
        {
            const F rho_ii = R.Get(i,i);
            if( Abs(rho_ii) <= ctrl.zeroTol )
                continue;
            else
                eta = Max(eta,Abs(R.Get(i,i+1)/rho_ii));
        }
    }
    else
    {
        for( Int i=0; i<minDim-1; ++i )
        {
            const F rho_ii = R.Get(i,i);
            if( Abs(rho_ii) <= ctrl.zeroTol )
                continue;
            else
                for( Int j=i+1; j<n; ++j )
                    eta = Max(eta,Abs(R.Get(i,j)/rho_ii));
        }
    }
    eta /= ( IsComplex<F>::value ? Sqrt(Real(2)) : Real(1) );

    return std::make_pair(delta,eta);
}

// Return the log of the absolute value of the determinant of the lattice
// (the sum of the logs of the nonzero diagonal entries of R)
template<typename F>
inline Base<F> LogVolume( const Matrix<F>& R )
{
    DEBUG_ONLY(CSE cse("lll::LogVolume"))
    typedef Base<F> Real;
    const Int m = R.Height();
    const Int n = R.Width();
    const Int minDim = Min(m,n);

    Real logVol = 0;
    for( Int j=0; j<minDim; ++j )
    {
        Real rho_j_j = R.GetRealPart(j,j);
        if( rho_j_j > Real(0) )
            logVol += Log(rho_j_j);
    }
    return logVol;
}

} // namespace lll

} // namespace El

#include "./LLL/Unblocked.hpp"
#include "./LLL/Blocked.hpp"

namespace El {

template<typename F>
LLLInfo<Base<F>> LLL
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& UInv,
  Matrix<F>& R,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LLL"))
    typedef Base<F> Real;
    if( ctrl.delta < Real(1)/Real(2) )
        LogicError("delta is assumed to be at least 1/2");
    if( ctrl.eta <= Real(1)/Real(2) || ctrl.eta >= Pow(ctrl.delta,Real(0.5)) )
        LogicError("eta should be in (1/2,sqrt(delta))");

    const Int n = B.Width();
    Identity( U, n, n ); 
    Identity( UInv, n, n );

    if( ctrl.presort )
    {
        QRCtrl<Real> qrCtrl;
        qrCtrl.smallestFirst = ctrl.smallestFirst;

        auto BCopy = B;
        Matrix<F> t;
        Matrix<Real> d;
        Permutation Omega;
        // TODO: Add support for qr::ProxyHouseholder as well
        El::QR( BCopy, t, d, Omega, qrCtrl );
        Omega.PermuteCols( B );
        Omega.PermuteCols( U );
        Omega.PermuteRows( UInv );
    }

    const bool useBlocked = false;
    const bool formU = true;
    const bool formUInv = true;
    if( useBlocked )
    {
        return lll::BlockedAlg( B, U, UInv, R, formU, formUInv, ctrl );
    }
    else
    {
        if( ctrl.deep )
            return lll::UnblockedDeepAlg
                   ( B, U, UInv, R, formU, formUInv, ctrl );
        else
            return lll::UnblockedAlg( B, U, UInv, R, formU, formUInv, ctrl );
    }
}

template<typename F>
LLLInfo<Base<F>>
LLL
( Matrix<F>& B,
  Matrix<F>& R,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LLL"))
    typedef Base<F> Real;
    if( ctrl.delta < Real(1)/Real(2) )
        LogicError("delta is assumed to be at least 1/2");
    if( ctrl.eta <= Real(1)/Real(2) || ctrl.eta >= Pow(ctrl.delta,Real(0.5)) )
        LogicError("eta should be in (1/2,sqrt(delta))");

    if( ctrl.presort )
    {
        QRCtrl<Real> qrCtrl;
        qrCtrl.smallestFirst = ctrl.smallestFirst;

        auto BCopy = B;
        Matrix<F> t;
        Matrix<Real> d;
        Permutation Omega;
        // TODO: Add support for qr::ProxyHouseholder as well
        El::QR( BCopy, t, d, Omega, qrCtrl );
        Omega.PermuteCols( B );
    }

    const bool useBlocked = false;
    const bool formU = false;
    const bool formUInv = false;
    Matrix<F> U, UInv;
    if( useBlocked )
    {
        return lll::BlockedAlg( B, U, UInv, R, formU, formUInv, ctrl );
    }
    else
    {
        if( ctrl.deep )
            return lll::UnblockedDeepAlg
                   ( B, U, UInv, R, formU, formUInv, ctrl );
        else
            return lll::UnblockedAlg( B, U, UInv, R, formU, formUInv, ctrl );
    }
}

// Emulate the flavor of quicksort/mergesort by recursively splitting the 
// vectors in half, applying LLL to each half, and merging the halves
// by running LLL on the interwoven reduced basis vectors
// (notice that this should allow the highest levels to often run at a lower
//  precision since the reduced basis vectors are likely to be much smaller, 
//  especially with SVP challenge lattices).
//
// NOTE: Until Complex<BigFloat> exists, we must have different implementations
//       for real and complex F

// TODO: Provide a way to display when changing precision without showing
//       all of the backtracks and deep insertions.

namespace lll {

template<typename Real>
LLLInfo<Real>
RecursiveHelper
( Matrix<Real>& B,
  Int cutoff,
  const LLLCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("lll::RecursiveHelper"))
    const Int n = B.Width();
    if( n < cutoff )
        return LLL( B, ctrl );
    Timer timer;
   
    auto C( B ); 

    const Int firstHalf = n-(n/2);
    auto CL = C( ALL, IR(0,firstHalf) );
    auto CR = C( ALL, IR(firstHalf,n) );

    double leftTime, rightTime;
    if( ctrl.time )
        timer.Start();
    auto leftInfo = RecursiveLLL( CL, cutoff, ctrl ); 
    if( ctrl.time )
    {
        leftTime = timer.Stop(); 
        timer.Start();
    }
    auto rightInfo = RecursiveLLL( CR, cutoff, ctrl );
    if( ctrl.time )
        rightTime = timer.Stop();
    if( ctrl.progress )
    {
        Output("n=",n);
        Output("  left swaps=",leftInfo.numSwaps);
        Output("  right swaps=",rightInfo.numSwaps);
    }
    if( ctrl.time )
    {
        Output("  left time:  ",leftTime," seconds");
        Output("  right time: ",rightTime," seconds");
    }
    const Real CLOneNorm = OneNorm( CL );
    const Real CROneNorm = OneNorm( CR );
    const Real CLMaxNorm = MaxNorm( CL );
    const Real CRMaxNorm = MaxNorm( CR );
    if( ctrl.progress )
    {
        Output("  || C_L ||_1 = ",CLOneNorm);
        Output("  || C_R ||_1 = ",CROneNorm);
        Output("  || C_L ||_max = ",CLMaxNorm);
        Output("  || C_R ||_max = ",CRMaxNorm);
    }

    const Real maxOneNorm = Max(CLOneNorm,CROneNorm);
    const Real fudge = 1.5;
    const Int neededPrec = Int(Ceil(Log2(maxOneNorm)*fudge));

    // Attempt to move to a lower precision
#ifdef EL_HAVE_MPC
    mpfr_prec_t inputPrec = mpc::Precision();
#endif
    try {
    if( PrecisionIsGreater<Real,double>::value && neededPrec <= 53 )
    {
        if( ctrl.progress )
            Output("  Dropping to double-precision");
        Matrix<double> BLowerPrec;
        BLowerPrec.Resize( B.Height(), n );
        // Interleave CL and CR to reform B before running LLL again
        // NOTE: This does not seem to make a substantial difference
        for( Int jSub=0; jSub<n/2; ++jSub )
        {
            auto cl = CL( ALL, IR(jSub) );
            auto cr = CR( ALL, IR(jSub) ); 
            auto bl = BLowerPrec( ALL, IR(2*jSub) );
            auto br = BLowerPrec( ALL, IR(2*jSub+1) );
            Copy( cl, bl );
            Copy( cr, br );
        }
        if( firstHalf > n/2 )
        {
            auto cl = CL( ALL, IR(firstHalf-1) );
            auto bl = BLowerPrec( ALL, IR(n-1) ); 
            Copy( cl, bl );
        }

        LLLCtrl<double> ctrlLowerPrec;
        ctrlLowerPrec.delta = double(ctrl.delta);
        ctrlLowerPrec.weak = ctrl.weak;
        ctrlLowerPrec.deep = ctrl.deep;
        ctrlLowerPrec.presort = ctrl.presort;
        ctrlLowerPrec.smallestFirst = ctrl.smallestFirst;
        ctrlLowerPrec.reorthogTol = double(ctrl.reorthogTol);
        ctrlLowerPrec.numOrthog = ctrl.numOrthog;
        ctrlLowerPrec.progress = ctrl.progress;
        ctrlLowerPrec.time = ctrl.time;
        if( ctrl.time )
            timer.Start();
        auto parentInfoLowerPrec = LLL( BLowerPrec, ctrlLowerPrec );
        if( ctrl.time )
            Output("  double-precision LLL took ",timer.Stop()," seconds");
        Copy( BLowerPrec, B );
        LLLInfo<Real> parentInfo;
        parentInfo.delta = parentInfoLowerPrec.delta;
        parentInfo.eta = parentInfoLowerPrec.eta;
        parentInfo.numSwaps = parentInfoLowerPrec.numSwaps + 
                              leftInfo.numSwaps + rightInfo.numSwaps;
        parentInfo.rank = parentInfoLowerPrec.rank;
        parentInfo.nullity = parentInfoLowerPrec.nullity;
        parentInfo.logVol = parentInfoLowerPrec.logVol;
        return parentInfo;
    }
#ifdef EL_HAVE_QUAD
    else if( PrecisionIsGreater<Real,Quad>::value && neededPrec <= 113 )
    {
        if( ctrl.progress )
            Output("  Dropping to quad-precision");
        Matrix<Quad> BLowerPrec;
        BLowerPrec.Resize( B.Height(), n );
        // Interleave CL and CR to reform B before running LLL again
        // NOTE: This does not seem to make a substantial difference
        for( Int jSub=0; jSub<n/2; ++jSub )
        {   
            auto cl = CL( ALL, IR(jSub) ); 
            auto cr = CR( ALL, IR(jSub) ); 
            auto bl = BLowerPrec( ALL, IR(2*jSub) );
            auto br = BLowerPrec( ALL, IR(2*jSub+1) );
            Copy( cl, bl );
            Copy( cr, br );
        }
        if( firstHalf > n/2 )
        {   
            auto cl = CL( ALL, IR(firstHalf-1) );
            auto bl = BLowerPrec( ALL, IR(n-1) );
            Copy( cl, bl );
        }

        LLLCtrl<Quad> ctrlLowerPrec;
        ctrlLowerPrec.delta = Quad(ctrl.delta);
        ctrlLowerPrec.weak = ctrl.weak;
        ctrlLowerPrec.deep = ctrl.deep;
        ctrlLowerPrec.presort = ctrl.presort;
        ctrlLowerPrec.smallestFirst = ctrl.smallestFirst;
        ctrlLowerPrec.reorthogTol = Quad(ctrl.reorthogTol);
        ctrlLowerPrec.numOrthog = ctrl.numOrthog;
        ctrlLowerPrec.progress = ctrl.progress;
        ctrlLowerPrec.time = ctrl.time;
        if( ctrl.time )
            timer.Start();
        auto parentInfoLowerPrec = LLL( BLowerPrec, ctrlLowerPrec );
        if( ctrl.time )
            Output("  quad-precision LLL took ",timer.Stop()," seconds");
        Copy( BLowerPrec, B );
        LLLInfo<Real> parentInfo;
        parentInfo.delta = parentInfoLowerPrec.delta;
        parentInfo.eta = parentInfoLowerPrec.eta;
        parentInfo.numSwaps = parentInfoLowerPrec.numSwaps + 
                              leftInfo.numSwaps + rightInfo.numSwaps;
        parentInfo.rank = parentInfoLowerPrec.rank;
        parentInfo.nullity = parentInfoLowerPrec.nullity;
        parentInfo.logVol = Quad(parentInfoLowerPrec.logVol);
        return parentInfo;
    }
#endif
#ifdef EL_HAVE_MPC
    // TODO: Debug the problem with this process
    //else if( neededPrec <= inputPrec )
    else if( false )
    {
        if( ctrl.progress )
            Output("  Dropping to prec=",neededPrec);
        mpc::SetPrecision( neededPrec );
        Matrix<BigFloat> BLowerPrec;
        BLowerPrec.Resize( B.Height(), n );
        // Interleave CL and CR to reform B before running LLL again
        // NOTE: This does not seem to make a substantial difference
        for( Int jSub=0; jSub<n/2; ++jSub )
        {
            auto cl = CL( ALL, IR(jSub) );
            auto cr = CR( ALL, IR(jSub) ); 
            auto bl = BLowerPrec( ALL, IR(2*jSub) );
            auto br = BLowerPrec( ALL, IR(2*jSub+1) );
            Copy( cl, bl );
            Copy( cr, br );
        }
        if( firstHalf > n/2 )
        {
            auto cl = CL( ALL, IR(firstHalf-1) );
            auto bl = BLowerPrec( ALL, IR(n-1) ); 
            Copy( cl, bl );
        }

        LLLCtrl<BigFloat> ctrlLowerPrec;
        ctrlLowerPrec.delta = ctrl.delta;
        ctrlLowerPrec.weak = ctrl.weak;
        ctrlLowerPrec.deep = ctrl.deep;
        ctrlLowerPrec.presort = ctrl.presort;
        ctrlLowerPrec.smallestFirst = ctrl.smallestFirst;
        ctrlLowerPrec.reorthogTol = ctrl.reorthogTol;
        ctrlLowerPrec.numOrthog = ctrl.numOrthog;
        ctrlLowerPrec.progress = ctrl.progress;
        ctrlLowerPrec.time = ctrl.time;
        if( ctrl.time )
            timer.Start();
        auto parentInfoLowerPrec = LLL( BLowerPrec, ctrlLowerPrec );
        if( ctrl.time )
            Output("  reduced-precision LLL took ",timer.Stop()," seconds");
        if( ctrl.progress )
            Output("  reverting back to prec=",inputPrec);
        mpc::SetPrecision( inputPrec );
        Copy( BLowerPrec, B );
        LLLInfo<Real> parentInfo;
        parentInfo.delta = Real(parentInfoLowerPrec.delta);
        parentInfo.eta = Real(parentInfoLowerPrec.eta);
        parentInfo.numSwaps = parentInfoLowerPrec.numSwaps + 
                              leftInfo.numSwaps + rightInfo.numSwaps;
        parentInfo.rank = parentInfoLowerPrec.rank;
        parentInfo.nullity = parentInfoLowerPrec.nullity;
        parentInfo.logVol = Real(parentInfoLowerPrec.logVol);
        return parentInfo;
    }
#endif
    }
    catch( std::exception& e )
    {
#ifdef EL_HAVE_MPC
        mpc::SetPrecision( inputPrec );
#endif
        Output("e.what()=",e.what()); 
    }

    // Interleave CL and CR to reform B before running LLL again
    // NOTE: This does not seem to make a substantial difference
    for( Int jSub=0; jSub<n/2; ++jSub )
    {
        auto cl = CL( ALL, IR(jSub) );
        auto cr = CR( ALL, IR(jSub) ); 
        auto bl = B( ALL, IR(2*jSub) );
        auto br = B( ALL, IR(2*jSub+1) );
        bl = cl;
        br = cr;
    }
    if( firstHalf > n/2 )
    {
        auto cl = CL( ALL, IR(firstHalf-1) );
        auto bl = B( ALL, IR(n-1) ); 
        bl = cl;
    }

    auto parentInfo = LLL( B, ctrl );
    parentInfo.numSwaps += leftInfo.numSwaps + rightInfo.numSwaps;
    return parentInfo;
}

template<typename Real>
LLLInfo<Real>
RecursiveHelper
( Matrix<Complex<Real>>& B,
  Int cutoff,
  const LLLCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("lll::RecursiveHelper"))
    typedef Complex<Real> F;
    const Int n = B.Width();
    if( n < cutoff )
        return LLL( B, ctrl );
    Timer timer;
   
    auto C( B ); 

    const Int firstHalf = n-(n/2);
    auto CL = C( ALL, IR(0,firstHalf) );
    auto CR = C( ALL, IR(firstHalf,n) );

    double leftTime, rightTime;
    if( ctrl.time )
        timer.Start();
    auto leftInfo = RecursiveLLL( CL, cutoff, ctrl ); 
    if( ctrl.time )
    {
        leftTime = timer.Stop(); 
        timer.Start();
    }
    auto rightInfo = RecursiveLLL( CR, cutoff, ctrl );
    if( ctrl.time )
        rightTime = timer.Stop();
    if( ctrl.progress )
    {
        Output("n=",n);
        Output("  left swaps=",leftInfo.numSwaps);
        Output("  right swaps=",rightInfo.numSwaps);
    }
    if( ctrl.time )
    {
        Output("  left time:  ",leftTime," seconds");
        Output("  right time: ",rightTime," seconds");
    }
    const Real CLOneNorm = OneNorm( CL );
    const Real CROneNorm = OneNorm( CR );
    const Real CLMaxNorm = MaxNorm( CL );
    const Real CRMaxNorm = MaxNorm( CR );
    if( ctrl.progress )
    {
        Output("  || C_L ||_1 = ",CLOneNorm);
        Output("  || C_R ||_1 = ",CROneNorm);
        Output("  || C_L ||_max = ",CLMaxNorm);
        Output("  || C_R ||_max = ",CRMaxNorm);
    }

    const Real maxOneNorm = Max(CLOneNorm,CROneNorm);
    const Real fudge = 1.5;
    const Int neededPrec = Int(Ceil(Log2(maxOneNorm)*fudge));

    // Attempt to move to a lower precision
    try {
    if( PrecisionIsGreater<Real,double>::value && neededPrec <= 53 )
    {
        typedef ConvertBase<F,double> FDbl;
        if( ctrl.progress )
            Output("  Dropping to double-precision");
        Matrix<FDbl> BLowerPrec;
        BLowerPrec.Resize( B.Height(), n );
        // Interleave CL and CR to reform B before running LLL again
        // NOTE: This does not seem to make a substantial difference
        for( Int jSub=0; jSub<n/2; ++jSub )
        {
            auto cl = CL( ALL, IR(jSub) );
            auto cr = CR( ALL, IR(jSub) ); 
            auto bl = BLowerPrec( ALL, IR(2*jSub) );
            auto br = BLowerPrec( ALL, IR(2*jSub+1) );
            Copy( cl, bl );
            Copy( cr, br );
        }
        if( firstHalf > n/2 )
        {
            auto cl = CL( ALL, IR(firstHalf-1) );
            auto bl = BLowerPrec( ALL, IR(n-1) ); 
            Copy( cl, bl );
        }

        LLLCtrl<double> ctrlLowerPrec;
        ctrlLowerPrec.delta = double(ctrl.delta);
        ctrlLowerPrec.weak = ctrl.weak;
        ctrlLowerPrec.deep = ctrl.deep;
        ctrlLowerPrec.presort = ctrl.presort;
        ctrlLowerPrec.smallestFirst = ctrl.smallestFirst;
        ctrlLowerPrec.reorthogTol = double(ctrl.reorthogTol);
        ctrlLowerPrec.numOrthog = ctrl.numOrthog;
        ctrlLowerPrec.progress = ctrl.progress;
        ctrlLowerPrec.time = ctrl.time;
        if( ctrl.time )
            timer.Start();
        auto parentInfoLowerPrec = LLL( BLowerPrec, ctrlLowerPrec );
        if( ctrl.time )
            Output("  double-precision LLL took ",timer.Stop()," seconds");
        Copy( BLowerPrec, B );
        LLLInfo<Real> parentInfo;
        parentInfo.delta = parentInfoLowerPrec.delta;
        parentInfo.eta = parentInfoLowerPrec.eta;
        parentInfo.numSwaps = parentInfoLowerPrec.numSwaps + 
                              leftInfo.numSwaps + rightInfo.numSwaps;
        parentInfo.rank = parentInfoLowerPrec.rank;
        parentInfo.nullity = parentInfoLowerPrec.nullity;
        parentInfo.logVol = parentInfoLowerPrec.logVol;
        return parentInfo;
    }
#ifdef EL_HAVE_QUAD
    else if( PrecisionIsGreater<Real,Quad>::value && neededPrec <= 113 )
    {
        typedef ConvertBase<F,Quad> FQuad;
        if( ctrl.progress )
            Output("  Dropping to quad-precision");
        Matrix<FQuad> BLowerPrec;
        BLowerPrec.Resize( B.Height(), n );
        // Interleave CL and CR to reform B before running LLL again
        // NOTE: This does not seem to make a substantial difference
        for( Int jSub=0; jSub<n/2; ++jSub )
        {   
            auto cl = CL( ALL, IR(jSub) ); 
            auto cr = CR( ALL, IR(jSub) ); 
            auto bl = BLowerPrec( ALL, IR(2*jSub) );
            auto br = BLowerPrec( ALL, IR(2*jSub+1) );
            Copy( cl, bl );
            Copy( cr, br );
        }
        if( firstHalf > n/2 )
        {   
            auto cl = CL( ALL, IR(firstHalf-1) );
            auto bl = BLowerPrec( ALL, IR(n-1) );
            Copy( cl, bl );
        }

        LLLCtrl<Quad> ctrlLowerPrec;
        ctrlLowerPrec.delta = Quad(ctrl.delta);
        ctrlLowerPrec.weak = ctrl.weak;
        ctrlLowerPrec.deep = ctrl.deep;
        ctrlLowerPrec.presort = ctrl.presort;
        ctrlLowerPrec.smallestFirst = ctrl.smallestFirst;
        ctrlLowerPrec.reorthogTol = Quad(ctrl.reorthogTol);
        ctrlLowerPrec.numOrthog = ctrl.numOrthog;
        ctrlLowerPrec.progress = ctrl.progress;
        ctrlLowerPrec.time = ctrl.time;
        if( ctrl.time )
            timer.Start();
        auto parentInfoLowerPrec = LLL( BLowerPrec, ctrlLowerPrec );
        if( ctrl.time )
            Output("  quad-precision LLL took ",timer.Stop()," seconds");
        Copy( BLowerPrec, B );
        LLLInfo<Real> parentInfo;
        parentInfo.delta = parentInfoLowerPrec.delta;
        parentInfo.eta = parentInfoLowerPrec.eta;
        parentInfo.numSwaps = parentInfoLowerPrec.numSwaps + 
                              leftInfo.numSwaps + rightInfo.numSwaps;
        parentInfo.rank = parentInfoLowerPrec.rank;
        parentInfo.nullity = parentInfoLowerPrec.nullity;
        parentInfo.logVol = Quad(parentInfoLowerPrec.logVol);
        return parentInfo;
    }
#endif
    }
    catch( std::exception& e )
    {
        Output("e.what()=",e.what()); 
    }

    // Interleave CL and CR to reform B before running LLL again
    // NOTE: This does not seem to make a substantial difference
    for( Int jSub=0; jSub<n/2; ++jSub )
    {
        auto cl = CL( ALL, IR(jSub) );
        auto cr = CR( ALL, IR(jSub) ); 
        auto bl = B( ALL, IR(2*jSub) );
        auto br = B( ALL, IR(2*jSub+1) );
        bl = cl;
        br = cr;
    }
    if( firstHalf > n/2 )
    {
        auto cl = CL( ALL, IR(firstHalf-1) );
        auto bl = B( ALL, IR(n-1) ); 
        bl = cl;
    }

    auto parentInfo = LLL( B, ctrl );
    parentInfo.numSwaps += leftInfo.numSwaps + rightInfo.numSwaps;
    return parentInfo;
}

} // namespace lll

template<typename F>
LLLInfo<Base<F>>
RecursiveLLL
( Matrix<F>& B,
  Int cutoff,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("RecursiveLLL"))
    return lll::RecursiveHelper( B, cutoff, ctrl );
}

template<typename F>
LLLInfo<Base<F>>
LLL
( Matrix<F>& B,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LLL"))
    Matrix<F> R;
    return LLL( B, R, ctrl );
}

#define PROTO(F) \
  template LLLInfo<Base<F>> LLL \
  ( Matrix<F>& B, \
    const LLLCtrl<Base<F>>& ctrl ); \
  template LLLInfo<Base<F>> LLL \
  ( Matrix<F>& B, \
    Matrix<F>& R, \
    const LLLCtrl<Base<F>>& ctrl ); \
  template LLLInfo<Base<F>> LLL \
  ( Matrix<F>& B, \
    Matrix<F>& U, \
    Matrix<F>& UInv, \
    Matrix<F>& R, \
    const LLLCtrl<Base<F>>& ctrl ); \
  template LLLInfo<Base<F>> RecursiveLLL \
  ( Matrix<F>& B, \
    Int cutoff, \
    const LLLCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
