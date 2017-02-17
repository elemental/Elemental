/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2016, Ron Estrin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LATTICE_LLL_HPP
#define EL_LATTICE_LLL_HPP

// Lenstra-Lenstra-Lovasz (LLL) lattice reduction
// ==============================================
// A reduced basis, say D, is an LLL(delta) reduction of an m x n matrix B if
//
//    B U = D = Q R,
//
// where U is unimodular (integer-valued with absolute determinant of 1)
// and Q R is a floating-point QR factorization of D that satisfies the three
//  properties:
//
//   1. R has non-negative diagonal
//
//   2. R is (eta) size-reduced:
//
//        | R(i,j) / R(i,i) | < phi(F) eta,  for all i < j, and
//
//      where phi(F) is 1 for a real field F or sqrt(2) for a complex
//      field F, and
//
//   3. R is (delta) Lovasz reduced:
//
//        delta R(i,i)^2 <= R(i+1,i+1)^2 + |R(i,i+1)|^2,  for all i.
//
// Please see
//
//   Henri Cohen, "A course in computational algebraic number theory"
//
// for more information on the "MLLL" variant of LLL used by Elemental to
// handle linearly dependent vectors (the algorithm was originally suggested by
// Mike Pohst).
//

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

template<typename Real>
struct LLLInfo
{
    Real delta;
    Real eta;
    Int rank;
    Int nullity;
    Int numSwaps=0;
    Int firstSwap;
    Real logVol;

    template<typename OtherReal>
    LLLInfo<Real>& operator=( const LLLInfo<OtherReal>& info )
    {
        delta = Real(info.delta);
        eta = Real(info.eta);
        rank = info.rank;
        nullity = info.nullity;
        numSwaps = info.numSwaps;
        firstSwap = info.firstSwap;
        logVol = Real(info.logVol);
        return *this;
    }

    LLLInfo() { }
    LLLInfo( const LLLInfo<Real>& ctrl ) { *this = ctrl; }
    template<typename OtherReal>
    LLLInfo( const LLLInfo<OtherReal>& ctrl ) { *this = ctrl; }
};

// Return the Gaussian estimate of the minimum-length vector
//
//   GH(L) = (1/sqrt(pi)) Gamma(n/2+1)^{1/n} |det(L)|^{1/n}.
//
// where n is the rank of the lattice L.
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real LatticeGaussianHeuristic( Int n, Real logVol )
{
    return Exp((LogGamma(Real(n)/Real(2)+Real(1))+logVol)/Real(n))/
           Sqrt(Pi<Real>());
}

enum LLLVariant {
  // A 'weak' LLL reduction only ensures that | R(i,i+1) / R(i,i) | is
  // bounded above by eta (or, for complex data, by sqrt(2) eta), but it often
  // produces much lower-quality basis vectors
  LLL_WEAK,
  LLL_NORMAL,
  // LLL with 'deep insertion' is no longer guaranteed to be polynomial time
  // but produces significantly higher quality bases than normal LLL.
  // See Schnorr and Euchner's "Lattice Basis Reduction: Improved Practical
  // Algorithms and Solving Subset Sum Problems".
  LLL_DEEP,
  // Going one step further, one can perform additional size reduction before
  // checking each deep insertion condition. See Schnorr's article
  // "Progress on LLL and Lattice Reduction" in the book "The LLL Algorithm",
  // edited by Nguyen and Vallee.
  LLL_DEEP_REDUCE
};

template<typename Real>
struct LLLCtrl
{
    Real delta=Real(3)/Real(4);
    Real eta=Real(1)/Real(2) + Pow(limits::Epsilon<Real>(),Real(0.9));

    LLLVariant variant=LLL_NORMAL;
    bool recursive=false;
    Int cutoff=10;

    // Fudge factor for determining whether to drop precision
    Real precisionFudge=Real(2);

    Int minColThresh = 0;

    // Ignore precision limits for QR factorization?
    bool unsafeSizeReduct=false;

    // Preprocessing with a "rank-obscuring" column-pivoted QR factorization
    // (in the manner suggested by Wubben et al.) can greatly decrease
    // the number of swaps within LLL in some circumstances
    bool presort=false;
    bool smallestFirst=true;

    // If the size-reduced column has a two-norm that is less than or
    // equal to `reorthogTol` times the  original two-norm, then reorthog.
    Real reorthogTol=0;

    // The number of times to execute the orthogonalization
    Int numOrthog=1;

    // If a size-reduced column has a two-norm less than or equal to 'zeroTol',
    // then it is interpreted as a zero vector (and forced to zero)
    Real zeroTol=Pow(limits::Epsilon<Real>(),Real(0.9));

    // Exploit the sparsity in the size reduction Axpy's unless at least the
    // following percentage of reductions were non-trivial
    float blockingThresh = 0.5f;

    bool progress=false;
    bool time=false;

    // If 'jumpstart' is true, start LLL under the assumption that the first
    // 'startCol' columns are already processed
    bool jumpstart=false;
    Int startCol=0;

    // In case of conversion from BigFloat to BigFloat with different precision
    LLLCtrl<Real>& operator=( const LLLCtrl<Real>& ctrl )
    {
        const Real eps = limits::Epsilon<Real>();
        const Real etaMin = Real(1)/Real(2)+Pow(eps,Real(0.9));
        const Real zeroTolMin = Pow(eps,Real(0.9));

        delta = Real(ctrl.delta);
        // NOTE: This does *not* seem to be equivalent to Max if the precisions
        //       are different
        eta = Real(ctrl.eta);
        if( eta < etaMin )
            eta = etaMin;
        precisionFudge = Real(ctrl.precisionFudge);
        minColThresh = ctrl.minColThresh;
        unsafeSizeReduct = ctrl.unsafeSizeReduct;
        variant = ctrl.variant;
        recursive = ctrl.recursive;
        cutoff = ctrl.cutoff;
        presort = ctrl.presort;
        smallestFirst = ctrl.smallestFirst;
        reorthogTol = Real(ctrl.reorthogTol);
        numOrthog = ctrl.numOrthog;
        // NOTE: This does *not* seem to be equivalent to Max if the precisions
        //       are different
        zeroTol = Real(ctrl.zeroTol);
        if( zeroTol < zeroTolMin )
            zeroTol = zeroTolMin;
        blockingThresh = ctrl.blockingThresh;
        progress = ctrl.progress;
        time = ctrl.time;
        jumpstart = ctrl.jumpstart;
        startCol = ctrl.startCol;
        return *this;
    }

    // We frequently need to convert datatypes, so make this easy
    template<typename OtherReal>
    LLLCtrl<Real>& operator=( const LLLCtrl<OtherReal>& ctrl )
    {
        const Real eps = limits::Epsilon<Real>();
        const Real etaMin = Real(1)/Real(2)+Pow(eps,Real(0.9));
        const Real zeroTolMin = Pow(eps,Real(0.9));

        delta = Real(ctrl.delta);
        eta = Max(etaMin,Real(ctrl.eta));
        precisionFudge = Real(ctrl.precisionFudge);
        minColThresh = ctrl.minColThresh;
        unsafeSizeReduct = ctrl.unsafeSizeReduct;
        variant = ctrl.variant;
        recursive = ctrl.recursive;
        cutoff = ctrl.cutoff;
        presort = ctrl.presort;
        smallestFirst = ctrl.smallestFirst;
        reorthogTol = Real(ctrl.reorthogTol);
        numOrthog = ctrl.numOrthog;
        zeroTol = Max(zeroTolMin,Real(ctrl.zeroTol));
        blockingThresh = ctrl.blockingThresh;
        progress = ctrl.progress;
        time = ctrl.time;
        jumpstart = ctrl.jumpstart;
        startCol = ctrl.startCol;
        return *this;
    }

    LLLCtrl() { }
    LLLCtrl( const LLLCtrl<Real>& ctrl ) { *this = ctrl; }
    template<typename OtherReal>
    LLLCtrl( const LLLCtrl<OtherReal>& ctrl ) { *this = ctrl; }
};

// TODO(poulson): Maintain B in BigInt form

template<typename Z,typename F=Z>
LLLInfo<Base<F>> LLL
( Matrix<Z>& B,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

template<typename Z,typename F=Z>
LLLInfo<Base<F>> LLL
( Matrix<Z>& B,
  Matrix<F>& R,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

template<typename Z,typename F=Z>
LLLInfo<Base<F>> LLL
( Matrix<Z>& B,
  Matrix<Z>& U,
  Matrix<F>& R,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

template<typename Z,typename F=Z>
LLLInfo<Base<F>> LLLWithQ
( Matrix<Z>& B,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

template<typename Z,typename F=Z>
LLLInfo<Base<F>> LLLWithQ
( Matrix<Z>& B,
  Matrix<Z>& U,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

namespace lll {

static Timer stepTimer, houseStepTimer,
       houseViewTimer, houseReflectTimer,
       applyHouseTimer, roundTimer,
       formSInvTimer;

// Return the achieved delta and eta reduction properties
template<typename F>
std::pair<Base<F>,Base<F>>
Achieved
( const Matrix<F>& R,
  const LLLCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
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
    // TODO(poulson): Decide if m < n requires checking the k=m-1 case.
    //
    Matrix<F> z;
    Real delta = limits::Max<Real>();
    for( Int i=0; i<minDim-1; ++i )
    {
        const Real rho_i_i = RealPart(R(i,i));
        if( Abs(rho_i_i) <= ctrl.zeroTol )
            continue;
        const Real rho_i_ip1 = Abs(R(i,i+1));
        const Real rho_ip1_ip1 = RealPart(R(i+1,i+1));

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
    // TODO(poulson): Decide if m < n requires checking the k=m-1 case.
    //
    // NOTE: phi(F) is 1 for real F and sqrt(2) for complex F.
    //
    Real eta = 0;
    if( ctrl.variant == LLL_WEAK )
    {
        for( Int i=0; i<minDim-1; ++i )
        {
            const F rho_ii = R(i,i);
            if( Abs(rho_ii) <= ctrl.zeroTol )
                continue;
            else
                eta = Max(eta,Abs(R(i,i+1)/rho_ii));
        }
    }
    else
    {
        for( Int i=0; i<minDim-1; ++i )
        {
            const F rho_ii = R(i,i);
            if( Abs(rho_ii) <= ctrl.zeroTol )
                continue;
            else
                for( Int j=i+1; j<n; ++j )
                    eta = Max(eta,Abs(R(i,j)/rho_ii));
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
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int m = R.Height();
    const Int n = R.Width();
    const Int minDim = Min(m,n);

    Real logVol = 0;
    for( Int j=0; j<minDim; ++j )
    {
        Real rho_j_j = RealPart(R(j,j));
        if( rho_j_j > Real(0) )
            logVol += Log(rho_j_j);
    }
    return logVol;
}

} // namespace lll

} // namespace El

#include <El/number_theory/lattice/LLL/Left.hpp>

namespace El {

template<typename Z, typename F>
LLLInfo<Base<F>> LLLWithQ
( Matrix<Z>& B,
  Matrix<Z>& U,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const LLLCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
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

    Int firstSwap = n;
    if( ctrl.presort )
    {
        if( ctrl.jumpstart )
            LogicError("Cannot combine jumpstarting with presorting");
        QRCtrl<Base<Z>> qrCtrl;
        qrCtrl.smallestFirst = ctrl.smallestFirst;

        auto BCopy = B;
        Matrix<Z> tPre;
        Matrix<Base<Z>> dPre;
        Permutation Omega;
        // TODO(poulson): Add support for qr::ProxyHouseholder as well
        El::QR( BCopy, tPre, dPre, Omega, qrCtrl );
        Omega.PermuteCols( B );
        Omega.PermuteCols( U );

        // TODO(poulson): Do not use such a pessimistic lower bound
        firstSwap = 0;
    }

    const bool formU = true;
    LLLInfo<Real> info;
    if( ctrl.variant == LLL_DEEP_REDUCE )
        info = lll::LeftDeepReduceAlg( B, U, QR, t, d, formU, ctrl );
    else if( ctrl.variant == LLL_DEEP )
        return lll::LeftDeepAlg( B, U, QR, t, d, formU, ctrl );
    else
        info = lll::LeftAlg( B, U, QR, t, d, formU, ctrl );

    info.firstSwap = Min(info.firstSwap,firstSwap);
    return info;
}

template<typename Z, typename F>
LLLInfo<Base<F>> LLL
( Matrix<Z>& B,
  Matrix<Z>& U,
  Matrix<F>& R,
  const LLLCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.jumpstart && ctrl.startCol > 0 )
        LogicError("Cannot jumpstart from this interface");
    typedef Base<F> Real;
    Matrix<F> t;
    Matrix<Real> d;
    auto info = LLLWithQ( B, U, R, t, d, ctrl );
    MakeTrapezoidal( UPPER, R );
    return info;
}

template<typename Z, typename F>
LLLInfo<Base<F>>
LLLWithQ
( Matrix<Z>& B,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const LLLCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
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

    Int firstSwap = n;
    if( ctrl.presort )
    {
        if( ctrl.jumpstart )
            LogicError("Cannot combine jumpstarting with presorting");

        QRCtrl<Base<Z>> qrCtrl;
        qrCtrl.smallestFirst = ctrl.smallestFirst;

        auto BCopy = B;
        Matrix<Z> tPre;
        Matrix<Base<Z>> dPre;
        Permutation Omega;
        // TODO(poulson): Add support for qr::ProxyHouseholder as well
        El::QR( BCopy, tPre, dPre, Omega, qrCtrl );
        Omega.PermuteCols( B );

        // TODO(poulson): Do not use such a pessimistic lower bound
        firstSwap = 0;
    }

    const bool formU = false;
    Matrix<Z> U;
    if( ctrl.variant == LLL_DEEP_REDUCE )
    {
        // Start with standard LLL
        auto infoReg = lll::LeftAlg( B, U, QR, t, d, formU, ctrl );
        firstSwap = Min(firstSwap,infoReg.firstSwap);

        // Move up from standard to deep
        auto ctrlMod( ctrl );
        ctrlMod.jumpstart = true;
        ctrlMod.startCol = 0;
        auto infoDeep = lll::LeftDeepAlg( B, U, QR, t, d, formU, ctrlMod );
        firstSwap = Min(firstSwap,infoDeep.firstSwap);

        // Move up from deep insertion to deep reduction
        auto infoDeepRed =
          lll::LeftDeepReduceAlg( B, U, QR, t, d, formU, ctrlMod );
        firstSwap = Min(firstSwap,infoDeepRed.firstSwap);

        infoDeepRed.numSwaps += infoReg.numSwaps + infoDeep.numSwaps;
        infoDeepRed.firstSwap = firstSwap;
        return infoDeepRed;
    }
    else if( ctrl.variant == LLL_DEEP )
    {
        // Start with standard LLL
        auto infoReg = lll::LeftAlg( B, U, QR, t, d, formU, ctrl );
        firstSwap = Min(firstSwap,infoReg.firstSwap);

        // Move from standard to deep insertion
        auto ctrlMod( ctrl );
        ctrlMod.jumpstart = true;
        ctrlMod.startCol = 0;
        auto infoDeep = lll::LeftDeepAlg( B, U, QR, t, d, formU, ctrlMod );
        firstSwap = Min(firstSwap,infoDeep.firstSwap);

        infoDeep.numSwaps += infoReg.numSwaps;
        infoDeep.firstSwap = firstSwap;
        return infoDeep;
    }
    else
        return lll::LeftAlg( B, U, QR, t, d, formU, ctrl );
}

template<typename Z, typename F>
LLLInfo<Base<F>>
LLL
( Matrix<Z>& B,
  Matrix<F>& R,
  const LLLCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
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

// TODO(poulson): Provide a way to display when changing precision without
// showing all of the backtracks and deep insertions.

template<typename F>
bool IsInteger( const Matrix<F>& A )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const F* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();

    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            if( ABuf[i+j*ALDim] != Round(ABuf[i+j*ALDim]) )
                return false;
    return true;
}

namespace lll {

template<typename Z,typename F,typename RealZLower,typename RealLower>
LLLInfo<RealLower>
LowerPrecisionMerge
( const Matrix<Z>& CL,
  const Matrix<Z>& CR,
        Matrix<Z>& B,
        Matrix<Z>& U,
        Matrix<F>& QR,
        Matrix<F>& t,
        Matrix<Base<F>>& d,
        bool maintainU,
  const LLLCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    typedef ConvertBase<F,RealLower> FLower;
    typedef ConvertBase<Z,RealZLower> ZLower;
    const string typeStringF = TypeName<RealLower>();
    const string typeStringZ = TypeName<RealZLower>();

    const Int n = B.Width();
    const Int firstHalf = n-(n/2);

    if( ctrl.progress || ctrl.time )
    {
        Output("  Dropping B  to " + typeStringZ);
        Output("  Dropping QR to " + typeStringF);
    }
    Matrix<ZLower> BLower;
    BLower.Resize( B.Height(), n );

    // Interleave CL and CR to reform B before running LLL again
    // NOTE: This does not seem to make a substantial difference
    for( Int jSub=0; jSub<n/2; ++jSub )
    {
        auto cL = CL( ALL, IR(jSub) );
        auto cR = CR( ALL, IR(jSub) );
        auto bL = BLower( ALL, IR(2*jSub) );
        auto bR = BLower( ALL, IR(2*jSub+1) );
        Copy( cL, bL );
        Copy( cR, bR );
    }
    if( firstHalf > n/2 )
    {
        auto cL = CL( ALL, IR(firstHalf-1) );
        auto bL = BLower( ALL, IR(n-1) );
        Copy( cL, bL );
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
    Matrix<ZLower> UNewLower;
    Matrix<FLower> tLower;
    Matrix<RealLower> dLower;
    if( maintainU )
        infoLower =
          LLLWithQ( BLower, UNewLower, QRLower, tLower, dLower, ctrlLower );
    else
        infoLower =
          LLLWithQ( BLower, QRLower, tLower, dLower, ctrlLower );
    if( ctrl.time )
        Output("  (" + typeStringZ + "," + typeStringF + ") LLL took ",timer.Stop()," seconds");
    Copy( BLower, B );
    Copy( QRLower, QR );
    Copy( tLower, t );
    Copy( dLower, d );

    if( maintainU )
    {
        Matrix<Z> UNew;
        Copy( UNewLower, UNew );
        auto UCopy( U );
        Gemm( NORMAL, NORMAL, Z(1), UCopy, UNew, Z(0), U );
    }

    return infoLower;
}

template<typename Z,typename F,typename RealLower>
bool TryLowerPrecisionMerge
( const Matrix<Z>& CL,
  const Matrix<Z>& CR,
        Matrix<Z>& B,
        Matrix<Z>& U,
        Matrix<F>& QR,
        Matrix<F>& t,
        Matrix<Base<F>>& d,
        bool maintainU,
  const LLLCtrl<Base<F>>& ctrl,
        unsigned neededPrec,
        LLLInfo<Base<F>>& info )
{
    bool succeeded = false;
    if( MantissaIsLonger<Base<Z>,RealLower>::value &&
        MantissaBits<RealLower>::value >= neededPrec )
    {
        try
        {
            if( MantissaIsLonger<Base<F>,RealLower>::value )
            {
                info = LowerPrecisionMerge<Z,F,RealLower,RealLower>
                  ( CL, CR, B, U, QR, t, d, maintainU, ctrl );
            }
            else
            {
                info = LowerPrecisionMerge<Z,F,RealLower,Base<F>>
                  ( CL, CR, B, U, QR, t, d, maintainU, ctrl );
            }
            succeeded = true;
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
    return succeeded;
}

#ifdef EL_HAVE_MPC
template<typename Z,typename F>
bool TryLowerPrecisionBigFloatMerge
( const Matrix<Z>& CL,
  const Matrix<Z>& CR,
        Matrix<Z>& B,
        Matrix<Z>& U,
        Matrix<F>& QR,
        Matrix<F>& t,
        Matrix<Base<F>>& d,
        bool maintainU,
  const LLLCtrl<Base<F>>& ctrl,
        unsigned neededPrec,
        LLLInfo<Base<F>>& info )
{
    bool succeeded = false;
    if( !IsFixedPrecision<Base<F>>::value )
    {
        // Only move down to a lower-precision MPFR type if the jump is
        // substantial. The current value has been naively chosen.
        const mpfr_prec_t minPrecDiff = 32;
        mpfr_prec_t inputPrec = mpfr::Precision();
        if( neededPrec <= inputPrec-minPrecDiff )
        {
            mpfr::SetPrecision( neededPrec );
            try
            {
                info = LowerPrecisionMerge<Z,F,BigFloat,Base<F>>
                  ( CL, CR, B, U, QR, t, d, maintainU, ctrl );
                succeeded = true;
            }
            catch( std::exception& e )
            {
                Output("e.what()=",e.what());
            }
            mpfr::SetPrecision( inputPrec );
        }
    }
    return succeeded;
}
#endif

template<typename Z,typename F>
LLLInfo<Base<F>>
RecursiveHelper
( Matrix<Z>& B,
  Matrix<Z>& U,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  Int numShuffles,
  bool maintainU,
  const LLLCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
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
    // NOTE: We can assume that the following lower bound should essentially
    //       always be tight.
    info.firstSwap = 0;
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
            Matrix<Z> ULNew;
            Matrix<F> QRL, tL;
            Matrix<Real> dL;
            leftInfo = RecursiveLLLWithQ( CL, ULNew, QRL, tL, dL, ctrl );

            auto UL = U( ALL, indL );
            auto ULCopy( UL );
            Gemm( NORMAL, NORMAL, Z(1), ULCopy, ULNew, Z(0), UL );
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
            Matrix<Z> URNew;
            Matrix<F> QRR, tR;
            Matrix<Real> dR;
            rightInfo = RecursiveLLLWithQ( CR, URNew, QRR, tR, dR, ctrl );

            auto UR = U( ALL, indR );
            auto URCopy( UR );
            Gemm( NORMAL, NORMAL, Z(1), URCopy, URNew, Z(0), UR );
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

        const bool isInteger = IsInteger( C );

        bool succeeded = false;
        Int numPrevSwaps = info.numSwaps;
        if( isInteger )
        {
            // Can we use QR for this without recomputing norms?
            Matrix<F> CLF;
            Copy(CL, CLF);
            Matrix<F> CRF;
            Copy(CR, CRF);
            const Real CLOneNorm = OneNorm( CLF );
            const Real CROneNorm = OneNorm( CRF );
            const Real CLMaxNorm = MaxNorm( CLF );
            const Real CRMaxNorm = MaxNorm( CRF );
            // TODO(poulson): Incorporate norm of U if maintaining U
            if( ctrl.progress )
            {
                Output("  || C_L ||_1 = ",CLOneNorm);
                Output("  || C_R ||_1 = ",CROneNorm);
                Output("  || C_L ||_max = ",CLMaxNorm);
                Output("  || C_R ||_max = ",CRMaxNorm);
            }
            const Real COneNorm = Max(CLOneNorm,CROneNorm);
            // TODO(poulson): Make tunable
            const Real fudge = ctrl.precisionFudge;
            const unsigned neededPrec = unsigned(Ceil(Log2(COneNorm)*fudge));
            if( ctrl.progress || ctrl.time )
            {
                Output("  || C ||_1 = ",COneNorm);
                Output("  Needed precision: ",neededPrec);
            }

            succeeded = TryLowerPrecisionMerge<Z,F,float>
              ( CL, CR, B, U, QR, t, d, maintainU, ctrl, neededPrec, info );
            if( !succeeded )
                succeeded = TryLowerPrecisionMerge<Z,F,double>
                  ( CL, CR, B, U, QR, t, d, maintainU, ctrl, neededPrec, info );
#ifdef EL_HAVE_QD
            if( !succeeded )
                succeeded = TryLowerPrecisionMerge<Z,F,DoubleDouble>
                  ( CL, CR, B, U, QR, t, d, maintainU, ctrl, neededPrec, info );
            if( !succeeded )
                succeeded = TryLowerPrecisionMerge<Z,F,QuadDouble>
                  ( CL, CR, B, U, QR, t, d, maintainU, ctrl, neededPrec, info );
#elif defined(EL_HAVE_QUAD)
            if( !succeeded )
                succeeded = TryLowerPrecisionMerge<Z,F,Quad>
                  ( CL, CR, B, U, QR, t, d, maintainU, ctrl, neededPrec, info );
#endif
#ifdef EL_HAVE_MPC
            if( !succeeded )
                succeeded = TryLowerPrecisionBigFloatMerge<Z,F>
                  ( CL, CR, B, U, QR, t, d, maintainU, ctrl, neededPrec, info );
#endif

            if( succeeded )
                info.numSwaps += numPrevSwaps;
        }
        // TODO(poulson): Allow for dropping with non-integer coefficients?

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

            auto ctrlMod( ctrl );
            ctrlMod.jumpstart = true;
            ctrlMod.startCol = 0;
            ctrlMod.recursive = false;
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

                info = LLLWithQ( B, U, QR, t, d, ctrlMod );
            }
            else
            {
                info = LLLWithQ( B, QR, t, d, ctrlMod );
            }
            info.numSwaps += numPrevSwaps;
        }
    }
    return info;
}

} // namespace lll

template<typename Z,typename F>
LLLInfo<Base<F>>
RecursiveLLLWithQ
( Matrix<Z>& B,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const LLLCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.jumpstart && ctrl.startCol > 0 )
        Output("Warning: Recursive LLL ignores jumpstarts");
    auto ctrlMod( ctrl );
    ctrlMod.jumpstart = false;

    // TODO(poulson): Make this runtime-tunable
    Int numShuffles = 1;
    Matrix<Z> U;
    bool maintainU=false;
    return
      lll::RecursiveHelper( B, U, QR, t, d, numShuffles, maintainU, ctrlMod );
}

template<typename Z,typename F>
LLLInfo<Base<F>>
RecursiveLLLWithQ
( Matrix<Z>& B,
  Matrix<Z>& U,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const LLLCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.jumpstart && ctrl.startCol > 0 )
        Output("Warning: Recursive LLL ignores jumpstarts");
    auto ctrlMod( ctrl );
    ctrlMod.jumpstart = false;

    // TODO(poulson): Make this runtime-tunable
    Int numShuffles = 1;
    bool maintainU=true;
    return lll::RecursiveHelper
      ( B, U, QR, t, d, numShuffles, maintainU, ctrlMod );
}

template<typename Z, typename F>
LLLInfo<Base<F>>
LLL
( Matrix<Z>& B,
  const LLLCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.jumpstart && ctrl.startCol > 0 )
        LogicError("Cannot jumpstart LLL from this interface");
    Matrix<F> R;
    return LLL( B, R, ctrl );
}

template<typename Z>
void DeepColSwap( Matrix<Z>& B, Int i, Int k )
{
    const Int m = B.Height();
    auto bi = B( ALL, IR(i) );
    auto bk = B( ALL, IR(k) );
    auto bkCopy( bk );

    Z* BBuf = B.Buffer();
    const Int BLDim = B.LDim();
    for( Int l=k-1; l>=i; --l )
        blas::Copy( m, &BBuf[l*BLDim], 1, &BBuf[(l+1)*BLDim], 1 );

    bi = bkCopy;
}

template<typename Z>
void DeepRowSwap( Matrix<Z>& B, Int i, Int k )
{
    const Int n = B.Width();
    auto bi = B( IR(i), ALL );
    auto bk = B( IR(k), ALL );
    auto bkCopy( bk );

    Z* BBuf = B.Buffer();
    const Int BLDim = B.LDim();
    for( Int l=k-1; l>=i; --l )
        blas::Copy( n, &BBuf[l], BLDim, &BBuf[l+1], BLDim );

    bi = bkCopy;
}

} // namespace El

#endif // ifndef EL_LATTICE_LLL_HPP
