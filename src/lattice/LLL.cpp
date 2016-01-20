/*
   Copyright (c) 2009-2016, Jack Poulson
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

static Timer stepTimer, houseStepTimer,
       houseViewTimer, houseReflectTimer,
       applyHouseTimer, roundTimer,
       formSInvTimer;

namespace lll {

// Return the achieved delta and eta reduction properties
template<typename F>
std::pair<Base<F>,Base<F>>
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
    if( ctrl.variant == LLL_WEAK )
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
Base<F> LogVolume( const Matrix<F>& R )
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
LLLInfo<Base<F>> LLLWithQ
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LLLWithQ"))
    typedef Base<F> Real;
    const Int n = B.Width();
    if( ctrl.recursive && ctrl.cutoff < n )
        return RecursiveLLLWithQ( B, U, QR, t, d, ctrl );

    if( ctrl.delta < Real(1)/Real(2) )
        LogicError("delta is assumed to be at least 1/2");
    if( ctrl.eta <= Real(1)/Real(2) || ctrl.eta >= Sqrt(ctrl.delta) )
        LogicError
        ("eta=",ctrl.eta," should be in (1/2,sqrt(delta)=",
         Sqrt(ctrl.delta),")");


    if( ctrl.jumpstart )
    {
        if( U.Height() != n || U.Width() != n )
            LogicError("U should have been n x n on input");
    }
    else
    {
        Identity( U, n, n ); 
    }

    if( ctrl.presort )
    {
        if( ctrl.jumpstart )
            LogicError("Cannot combine jumpstarting with presorting");
        QRCtrl<Real> qrCtrl;
        qrCtrl.smallestFirst = ctrl.smallestFirst;

        auto BCopy = B;
        Matrix<F> tPre;
        Matrix<Real> dPre;
        Permutation Omega;
        // TODO: Add support for qr::ProxyHouseholder as well
        El::QR( BCopy, tPre, dPre, Omega, qrCtrl );
        Omega.PermuteCols( B );
        Omega.PermuteCols( U );
    }

    const bool useBlocked = false;
    const bool formU = true;
    if( useBlocked )
    {
        return lll::BlockedAlg( B, U, QR, t, d, formU, ctrl );
    }
    else
    {
        if( ctrl.variant == LLL_DEEP_REDUCE )
            return lll::UnblockedDeepReduceAlg( B, U, QR, t, d, formU, ctrl );
        else if( ctrl.variant == LLL_DEEP )
            return lll::UnblockedDeepAlg( B, U, QR, t, d, formU, ctrl );
        else
            return lll::UnblockedAlg( B, U, QR, t, d, formU, ctrl );
    }
}

template<typename F>
LLLInfo<Base<F>> LLL
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& R,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LLL"))
    if( ctrl.jumpstart && ctrl.startCol > 0 )
        LogicError("Cannot jumpstart from this interface");
    typedef Base<F> Real;
    Matrix<F> t;
    Matrix<Real> d;
    auto info = LLLWithQ( B, U, R, t, d, ctrl );
    MakeTrapezoidal( UPPER, R );
    return info;
}

template<typename F>
LLLInfo<Base<F>>
LLLWithQ
( Matrix<F>& B,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LLLWithQ"))
    typedef Base<F> Real;
    const Int n = B.Width();
    if( ctrl.recursive && ctrl.cutoff < n )
        return RecursiveLLLWithQ( B, QR, t, d, ctrl );

    if( ctrl.delta < Real(1)/Real(2) )
        LogicError("delta is assumed to be at least 1/2");
    if( ctrl.eta <= Real(1)/Real(2) || ctrl.eta >= Sqrt(ctrl.delta) )
        LogicError
        ("eta=",ctrl.eta," should be in (1/2,sqrt(delta)=",
         Sqrt(ctrl.delta),")");

    if( ctrl.presort )
    {
        if( ctrl.jumpstart )
            LogicError("Cannot combine jumpstarting with presorting");

        QRCtrl<Real> qrCtrl;
        qrCtrl.smallestFirst = ctrl.smallestFirst;

        auto BCopy = B;
        Matrix<F> tPre;
        Matrix<Real> dPre;
        Permutation Omega;
        // TODO: Add support for qr::ProxyHouseholder as well
        El::QR( BCopy, tPre, dPre, Omega, qrCtrl );
        Omega.PermuteCols( B );
    }

    const bool useBlocked = false;
    const bool formU = false;
    Matrix<F> U;
    if( useBlocked )
    {
        return lll::BlockedAlg( B, U, QR, t, d, formU, ctrl );
    }
    else
    {
        if( ctrl.variant == LLL_DEEP_REDUCE )
            return lll::UnblockedDeepReduceAlg( B, U, QR, t, d, formU, ctrl );
        else if( ctrl.variant == LLL_DEEP )
            return lll::UnblockedDeepAlg( B, U, QR, t, d, formU, ctrl );
        else
            return lll::UnblockedAlg( B, U, QR, t, d, formU, ctrl );
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
    if( ctrl.jumpstart && ctrl.startCol > 0 )
        LogicError("Cannot jumpstart from this interface");
    typedef Base<F> Real;
    Matrix<F> t;
    Matrix<Real> d;
    auto info = LLLWithQ( B, R, t, d, ctrl );
    MakeTrapezoidal( UPPER, R );
    return info;
}

// Emulate the flavor of quicksort/mergesort by recursively splitting the 
// vectors in half, applying LLL to each half, and merging the halves
// by running LLL on the interwoven reduced basis vectors
// (notice that this should allow the highest levels to often run at a lower
//  precision since the reduced basis vectors are likely to be much smaller, 
//  especially with SVP challenge lattices).
//
// C.f. The analogue of Lehmer's version of Euclid's algorithm that Schnorr
// mentions at the end of "Progress on LLL and Lattice Reduction".
//
// NOTE: Until Complex<BigFloat> exists, we must have different implementations
//       for real and complex F

// TODO: Provide a way to display when changing precision without showing
//       all of the backtracks and deep insertions.

namespace lll {

template<typename F,typename RealLower>
LLLInfo<RealLower>
LowerPrecisionMerge
( const Matrix<F>& CL,
  const Matrix<F>& CR,
        Matrix<F>& B,
        Matrix<F>& U,
        Matrix<F>& QR,
        Matrix<F>& t,
        Matrix<Base<F>>& d,
        bool maintainU,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("lll::LowerPrecisionMerge"))
    typedef ConvertBase<F,RealLower> FLower;
    const string typeString = TypeName<RealLower>();

    const Int n = B.Width();
    const Int firstHalf = n-(n/2);

    if( ctrl.progress || ctrl.time )
        Output("  Dropping to " + typeString);
    Matrix<FLower> BLower;
    BLower.Resize( B.Height(), n );
    // Interleave CL and CR to reform B before running LLL again
    // NOTE: This does not seem to make a substantial difference
    for( Int jSub=0; jSub<n/2; ++jSub )
    {
        auto cl = CL( ALL, IR(jSub) );
        auto cr = CR( ALL, IR(jSub) );
        auto bl = BLower( ALL, IR(2*jSub) );
        auto br = BLower( ALL, IR(2*jSub+1) );
        Copy( cl, bl );
        Copy( cr, br );
    }
    if( firstHalf > n/2 )
    {
        auto cl = CL( ALL, IR(firstHalf-1) );
        auto bl = BLower( ALL, IR(n-1) );
        Copy( cl, bl );
    }

    LLLCtrl<RealLower> ctrlLower( ctrl );
    ctrlLower.recursive = false;
    RealLower eps = limits::Epsilon<RealLower>();
    RealLower minEta = RealLower(1)/RealLower(2)+Pow(eps,RealLower(0.9));
    ctrlLower.eta = Max(minEta,ctrlLower.eta);
    Timer timer;
    Matrix<FLower> QRLower;
    if( ctrl.time )
        timer.Start();
    LLLInfo<RealLower> infoLower;
    Matrix<FLower> UNewLower, tLower;
    Matrix<RealLower> dLower;
    if( maintainU )
        infoLower =
          LLLWithQ( BLower, UNewLower, QRLower, tLower, dLower, ctrlLower );
    else
        infoLower =
          LLLWithQ( BLower, QRLower, tLower, dLower, ctrlLower );
    if( ctrl.time )
        Output("  " + typeString + " LLL took ",timer.Stop()," seconds");
    Copy( BLower, B );
    Copy( QRLower, QR );
    Copy( tLower, t );
    Copy( dLower, d );

    if( maintainU )
    {
        Matrix<F> UNew;
        Copy( UNewLower, UNew );
        auto UCopy( U );
        Gemm( NORMAL, NORMAL, F(1), UCopy, UNew, F(0), U );
    }

    return infoLower;
}

template<typename Real>
LLLInfo<Real>
RecursiveHelper
( Matrix<Real>& B,
  Matrix<Real>& U,
  Matrix<Real>& QR,
  Matrix<Real>& t,
  Matrix<Real>& d,
  Int numShuffles,
  bool maintainU,
  const LLLCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("lll::RecursiveHelper"))

    typedef Real F;
    const Int n = B.Width();
    if( n < ctrl.cutoff )
    {
        auto ctrlMod( ctrl );
        ctrlMod.recursive = false;
        if( maintainU )
            return LLLWithQ( B, U, QR, t, d, ctrlMod );
        else
            return LLLWithQ( B, QR, t, d, ctrlMod );
    }
    Timer timer;

    LLLInfo<Real> info;
    info.numSwaps = 0;
    for( Int shuffle=0; shuffle<=numShuffles; ++shuffle )
    {
        if( ctrl.progress || ctrl.time )
            Output("Shuffle=",shuffle);
        auto C( B ); 

        const Int firstHalf = n-(n/2);
        Range<Int> indL(0,firstHalf), indR(firstHalf,n);
        auto CL = C( ALL, indL );
        auto CR = C( ALL, indR );

        double leftTime;
        if( ctrl.time )
            timer.Start();
        LLLInfo<Real> leftInfo;
        if( maintainU )
        {
            Matrix<F> ULNew, QRL, tL;
            Matrix<Real> dL;
            leftInfo = RecursiveLLLWithQ( CL, ULNew, QRL, tL, dL, ctrl );

            auto UL = U( indL, indL );
            auto ULCopy( UL );
            Gemm( NORMAL, NORMAL, F(1), ULCopy, ULNew, F(0), UL );
        }
        else
        {
            Matrix<F> QRL, tL;
            Matrix<Real> dL;
            leftInfo = RecursiveLLLWithQ( CL, QRL, tL, dL, ctrl ); 
        }
        if( ctrl.time )
        {
            leftTime = timer.Stop(); 
            timer.Start();
        }

        double rightTime;
        LLLInfo<Real> rightInfo;
        if( maintainU )
        {
            Matrix<F> URNew, QRR, tR;
            Matrix<Real> dR;
            rightInfo = RecursiveLLLWithQ( CR, URNew, QRR, tR, dR, ctrl );

            auto UR = U( indR, indR );
            auto URCopy( UR );
            Gemm( NORMAL, NORMAL, F(1), URCopy, URNew, F(0), UR );
        }
        else
        {
            Matrix<F> QRR, tR;
            Matrix<Real> dR;
            rightInfo = RecursiveLLLWithQ( CR, QRR, tR, dR, ctrl );
        }
        if( ctrl.time )
            rightTime = timer.Stop();

        info.numSwaps += leftInfo.numSwaps + rightInfo.numSwaps;
        if( ctrl.progress || ctrl.time )
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
        // TODO: Incorporate norm of U if maintaining U
        if( ctrl.progress )
        {
            Output("  || C_L ||_1 = ",CLOneNorm);
            Output("  || C_R ||_1 = ",CROneNorm);
            Output("  || C_L ||_max = ",CLMaxNorm);
            Output("  || C_R ||_max = ",CRMaxNorm);
        }

        const Real COneNorm = Max(CLOneNorm,CROneNorm);
        const Real fudge = 1.5; // TODO: Make tunable
        const unsigned neededPrec = unsigned(Ceil(Log2(COneNorm)*fudge));
        if( ctrl.progress || ctrl.time )
        {
            Output("  || C ||_1 = ",COneNorm);
            Output("  Needed precision: ",neededPrec);
        }

        bool succeeded = false;
        Int numPrevSwaps = info.numSwaps;
        if( MantissaIsLonger<Real,float>::value &&
            MantissaBits<float>::value >= neededPrec )
        {
            try
            {
                info = LowerPrecisionMerge<F,float>
                  ( CL, CR, B, U, QR, t, d, maintainU, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            }
            catch( std::exception& e )
            { Output("e.what()=",e.what()); }
        }
        if( !succeeded && 
            MantissaIsLonger<Real,double>::value &&
            MantissaBits<double>::value >= neededPrec )
        {
            try
            {
                info = LowerPrecisionMerge<F,double>
                  ( CL, CR, B, U, QR, t, d, maintainU, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            }
            catch( std::exception& e )
            { Output("e.what()=",e.what()); }
        }
#ifdef EL_HAVE_QD
        if( !succeeded && 
            MantissaIsLonger<Real,DoubleDouble>::value &&
            MantissaBits<DoubleDouble>::value >= neededPrec )
        {
            try
            {
                info =
                  LowerPrecisionMerge<F,DoubleDouble>
                  ( CL, CR, B, U, QR, t, d, maintainU, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            }
            catch( std::exception& e )
            { Output("e.what()=",e.what()); }
        }
        if( !succeeded && 
            MantissaIsLonger<Real,QuadDouble>::value &&
            MantissaBits<QuadDouble>::value >= neededPrec )
        {
            try
            {
                info =
                  LowerPrecisionMerge<F,QuadDouble>
                  ( CL, CR, B, U, QR, t, d, maintainU, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            }
            catch( std::exception& e )
            { Output("e.what()=",e.what()); }
        }
#elif defined(EL_HAVE_QUAD)
        if( !succeeded && 
            MantissaIsLonger<Real,Quad>::value &&
            MantissaBits<Quad>::value >= neededPrec )
        {
            try
            {
                info =
                  LowerPrecisionMerge<F,Quad>
                  ( CL, CR, B, U, QR, t, d, maintainU, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            }
            catch( std::exception& e )
            { Output("e.what()=",e.what()); }
        }
#endif
#ifdef EL_HAVE_MPC
        // Only move down to a lower-precision MPFR type if the jump is
        // substantial. The current value has been naively chosen.
        const mpfr_prec_t minPrecDiff = 32;
        mpfr_prec_t inputPrec = mpc::Precision();
        if( !succeeded && neededPrec <= inputPrec-minPrecDiff )
        {
            mpc::SetPrecision( neededPrec );
            try {
                info = LowerPrecisionMerge<F,BigFloat>
                  ( CL, CR, B, U, QR, t, d, maintainU, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            }
            catch( std::exception& e )
            { Output("e.what()=",e.what()); }
            mpc::SetPrecision( inputPrec );
        }
#endif

        if( !succeeded )
        {
            // Interleave CL and CR to reform B before running LLL again
            for( Int jSub=0; jSub<n/2; ++jSub )
            {
                auto cL = CL( ALL, IR(jSub) );
                auto cR = CR( ALL, IR(jSub) ); 
                auto bL = B( ALL, IR(2*jSub) );
                auto bR = B( ALL, IR(2*jSub+1) );
                bL = cL;
                bR = cR;
            }
            if( firstHalf > n/2 )
            {
                auto cL = CL( ALL, IR(firstHalf-1) );
                auto bL = B( ALL, IR(n-1) ); 
                bL = cL;
            }

            if( maintainU )
            {
                auto UCopy( U );
                auto UL = U( ALL, indL );
                auto UR = U( ALL, indR );
                auto UCopyL = UCopy( ALL, indL );
                auto UCopyR = UCopy( ALL, indR );
                for( Int jSub=0; jSub<n/2; ++jSub )
                {
                    auto uCopyL = UCopyL( ALL, IR(jSub) );
                    auto uCopyR = UCopyR( ALL, IR(jSub) );
                    auto uL = UL( ALL, IR(2*jSub) );
                    auto uR = UR( ALL, IR(2*jSub+1) );
                    uL = uCopyL;
                    uR = uCopyR;
                }
                if( firstHalf > n/2 )
                {
                    auto uCopyL = UCopyL( ALL, IR(firstHalf-1) );
                    auto uL = U( ALL, IR(n-1) );
                    uL = uCopyL;
                }

                auto ctrlMod( ctrl );
                ctrlMod.jumpstart = true;
                ctrlMod.startCol = 0;
                info = LLLWithQ( B, U, QR, t, d, ctrlMod );
            }
            else
            {
                info = LLLWithQ( B, QR, t, d, ctrl );
            }
            info.numSwaps += numPrevSwaps;
        }
    }
    return info;
}

// Same as the above, but with the Complex<BigFloat> datatype avoided
template<typename Real>
LLLInfo<Real>
RecursiveHelper
( Matrix<Complex<Real>>& B,
  Matrix<Complex<Real>>& U,
  Matrix<Complex<Real>>& QR,
  Matrix<Complex<Real>>& t,
  Matrix<Real>& d,
  Int numShuffles,
  bool maintainU,
  const LLLCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("lll::RecursiveHelper"))

    typedef Complex<Real> F;
    const Int n = B.Width();
    if( n < ctrl.cutoff )
    {
        auto ctrlMod( ctrl );
        ctrlMod.recursive = false;
        if( maintainU )
            return LLLWithQ( B, U, QR, t, d, ctrlMod );
        else
            return LLLWithQ( B, QR, t, d, ctrlMod );
    }
    Timer timer;

    LLLInfo<Real> info;
    info.numSwaps = 0;
    for( Int shuffle=0; shuffle<=numShuffles; ++shuffle )
    {
        if( ctrl.progress || ctrl.time )
            Output("Shuffle=",shuffle);
        auto C( B ); 

        const Int firstHalf = n-(n/2);
        Range<Int> indL(0,firstHalf), indR(firstHalf,n);
        auto CL = C( ALL, indL );
        auto CR = C( ALL, indR );

        double leftTime;
        if( ctrl.time )
            timer.Start();
        LLLInfo<Real> leftInfo;
        if( maintainU )
        {
            Matrix<F> ULNew, QRL, tL;
            Matrix<Real> dL;
            leftInfo = RecursiveLLLWithQ( CL, ULNew, QRL, tL, dL, ctrl );

            auto UL = U( indL, indL );
            auto ULCopy( UL );
            Gemm( NORMAL, NORMAL, F(1), ULCopy, ULNew, F(0), UL );
        }
        else
        {
            Matrix<F> QRL, tL;
            Matrix<Real> dL;
            leftInfo = RecursiveLLLWithQ( CL, QRL, tL, dL, ctrl );
        }
        if( ctrl.time )
        {
            leftTime = timer.Stop();
            timer.Start();
        }

        double rightTime;
        LLLInfo<Real> rightInfo;
        if( maintainU )
        {
            Matrix<F> URNew, QRR, tR;
            Matrix<Real> dR;
            rightInfo = RecursiveLLLWithQ( CR, URNew, QRR, tR, dR, ctrl );

            auto UR = U( indR, indR );
            auto URCopy( UR );
            Gemm( NORMAL, NORMAL, F(1), URCopy, URNew, F(0), UR );
        }
        else
        {
            Matrix<F> QRR, tR;
            Matrix<Real> dR;
            rightInfo = RecursiveLLLWithQ( CR, QRR, tR, dR, ctrl );
        }
        if( ctrl.time )
            rightTime = timer.Stop();

        info.numSwaps += leftInfo.numSwaps + rightInfo.numSwaps;
        if( ctrl.progress || ctrl.time )
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
        // TODO: Incorporate U norm
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

        const Real COneNorm = Max(CLOneNorm,CROneNorm);
        const Real fudge = 1.5; // TODO: Make tunable
        const unsigned neededPrec = unsigned(Ceil(Log2(COneNorm)*fudge));
        if( ctrl.progress || ctrl.time )
        {
            Output("  || C ||_1 = ",COneNorm);
            Output("  Needed precision: ",neededPrec);
        }

        bool succeeded = false;
        Int numPrevSwaps = info.numSwaps;
        if( MantissaIsLonger<Real,float>::value &&
            MantissaBits<float>::value >= neededPrec )
        {
            try
            {
                info = LowerPrecisionMerge<F,float>
                  ( CL, CR, B, U, QR, t, d, maintainU, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            } catch( std::exception& e ) { Output("e.what()=",e.what()); }
        }
        if( !succeeded && 
            MantissaIsLonger<Real,double>::value &&
            MantissaBits<double>::value >= neededPrec )
        {
            try
            {
                info = LowerPrecisionMerge<F,double>
                  ( CL, CR, B, U, QR, t, d, maintainU, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            } catch( std::exception& e ) { Output("e.what()=",e.what()); }
        }
        // There is not yet support for Complex<{Quad,Double}Double>
#ifdef EL_HAVE_QUAD
        if( !succeeded &&
            MantissaIsLonger<Real,Quad>::value &&
            MantissaBits<Quad>::value >= neededPrec )
        {
            try
            {
                info = LowerPrecisionMerge<F,Quad>
                  ( CL, CR, B, U, QR, t, d, maintainU, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            } catch( std::exception& e ) { Output("e.what()=",e.what()); }
        }
#endif
        if( !succeeded )
        {
            // Interleave CL and CR to reform B before running LLL again
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

            if( maintainU )
            {
                auto UCopy( U );
                auto UL = U( ALL, indL );
                auto UR = U( ALL, indR );
                auto UCopyL = UCopy( ALL, indL );
                auto UCopyR = UCopy( ALL, indR );
                for( Int jSub=0; jSub<n/2; ++jSub )
                {
                    auto uCopyL = UCopyL( ALL, IR(jSub) );
                    auto uCopyR = UCopyR( ALL, IR(jSub) );
                    auto uL = UL( ALL, IR(2*jSub) );
                    auto uR = UR( ALL, IR(2*jSub+1) );
                    uL = uCopyL;
                    uR = uCopyR;
                }
                if( firstHalf > n/2 )
                {
                    auto uCopyL = UCopyL( ALL, IR(firstHalf-1) );
                    auto uL = U( ALL, IR(n-1) );
                    uL = uCopyL;
                }

                auto ctrlMod( ctrl );
                ctrlMod.jumpstart = true;
                ctrlMod.startCol = 0;
                info = LLLWithQ( B, U, QR, t, d, ctrlMod );
            }
            else
            {
                info = LLLWithQ( B, QR, t, d, ctrl );
            }
            info.numSwaps += numPrevSwaps;
        }
    }
    return info;
}

} // namespace lll

template<typename F>
LLLInfo<Base<F>>
RecursiveLLLWithQ
( Matrix<F>& B,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("RecursiveLLLWithQ"))
    if( ctrl.jumpstart && ctrl.startCol > 0 )
        LogicError("Cannot jumpstart LLL from this interface");

    // TODO: Make this runtime-tunable
    Int numShuffles = 1;
    Matrix<F> U;
    bool maintainU=false;
    return lll::RecursiveHelper( B, U, QR, t, d, numShuffles, maintainU, ctrl );
}

template<typename F>
LLLInfo<Base<F>>
RecursiveLLLWithQ
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("RecursiveLLL"))
    if( ctrl.jumpstart && ctrl.startCol > 0 )
        LogicError("Cannot jumpstart LLL from this interface");

    // TODO: Make this runtime-tunable
    Int numShuffles = 1;
    bool maintainU=true;
    return lll::RecursiveHelper
      ( B, U, QR, t, d, numShuffles, maintainU, ctrl );
}

template<typename F>
LLLInfo<Base<F>>
LLL
( Matrix<F>& B,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LLL"))
    if( ctrl.jumpstart && ctrl.startCol > 0 )
        LogicError("Cannot jumpstart LLL from this interface");
    Matrix<F> R;
    return LLL( B, R, ctrl );
}

template<typename F>
void DeepColSwap( Matrix<F>& B, Int i, Int k )
{
    const Int m = B.Height();
    auto bi = B( ALL, IR(i) );
    auto bk = B( ALL, IR(k) );
    auto bkCopy( bk );

    F* BBuf = B.Buffer();
    const Int BLDim = B.LDim();
    for( Int l=k-1; l>=i; --l )
        blas::Copy( m, &BBuf[l*BLDim], 1, &BBuf[(l+1)*BLDim], 1 );

    bi = bkCopy;
}

template<typename F>
void DeepRowSwap( Matrix<F>& B, Int i, Int k )
{
    const Int n = B.Width();
    auto bi = B( IR(i), ALL );
    auto bk = B( IR(k), ALL );
    auto bkCopy( bk );

    F* BBuf = B.Buffer();
    const Int BLDim = B.LDim();
    for( Int l=k-1; l>=i; --l )
        blas::Copy( n, &BBuf[l], BLDim, &BBuf[l+1], BLDim );

    bi = bkCopy;
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
    Matrix<F>& R, \
    const LLLCtrl<Base<F>>& ctrl ); \
  template LLLInfo<Base<F>> LLLWithQ \
  ( Matrix<F>& B, \
    Matrix<F>& QR, \
    Matrix<F>& t, \
    Matrix<Base<F>>& d, \
    const LLLCtrl<Base<F>>& ctrl ); \
  template LLLInfo<Base<F>> LLLWithQ \
  ( Matrix<F>& B, \
    Matrix<F>& U, \
    Matrix<F>& QR, \
    Matrix<F>& t, \
    Matrix<Base<F>>& d, \
    const LLLCtrl<Base<F>>& ctrl ); \
  template void DeepColSwap( Matrix<F>& B, Int i, Int k ); \
  template void DeepRowSwap( Matrix<F>& B, Int i, Int k ); \
  template std::pair<Base<F>,Base<F>> lll::Achieved \
  ( const Matrix<F>& R, const LLLCtrl<Base<F>>& ctrl ); \
  template Base<F> lll::LogVolume( const Matrix<F>& R );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
