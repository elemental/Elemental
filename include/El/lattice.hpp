/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LATTICE_HPP
#define EL_LATTICE_HPP

namespace El {

// Deep insertion
// ==============
template<typename F>
void DeepColSwap( Matrix<F>& B, Int i, Int k );
template<typename F>
void DeepRowSwap( Matrix<F>& B, Int i, Int k );

// Lattice Log Potential
// =====================
// NOTE: This operates on the Gaussian Normal Form (R from QR)
template<typename F>
Base<F> LatticeLogPotential( const Matrix<F>& R );

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

// TODO: Maintain B in BigInt form

template<typename F>
LLLInfo<Base<F>> LLL
( Matrix<F>& B,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

template<typename F>
LLLInfo<Base<F>> LLL
( Matrix<F>& B,
  Matrix<F>& R,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

template<typename F>
LLLInfo<Base<F>> LLL
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& R,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

template<typename F>
LLLInfo<Base<F>> LLLWithQ
( Matrix<F>& B,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

template<typename F>
LLLInfo<Base<F>> LLLWithQ
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

// Fill M with its (quasi-reduced) image of B, and fill K with the
// LLL-reduced basis for the kernel of B.
//
// This is essentially Algorithm 2.7.1 from Cohen's
// "A course in computational algebraic number theory". The main difference
// is that we avoid solving the normal equations and call a least squares
// solver.
// 
template<typename F>
void LatticeImageAndKernel
( const Matrix<F>& B,
        Matrix<F>& M,
        Matrix<F>& K,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

// Fill K with the LLL-reduced basis for the image of B.
// TODO: Make this faster than LatticeImageAndKernel
template<typename F>
void LatticeImage
( const Matrix<F>& B,
        Matrix<F>& M,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

// Fill K with the LLL-reduced basis for the kernel of B.
// This will eventually mirror Algorithm 2.7.2 from Cohen's
// "A course in computational algebraic number theory".
template<typename F>
void LatticeKernel
( const Matrix<F>& B,
        Matrix<F>& K,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

// Search for Z-dependence
// =======================
// Search for Z-dependence of a vector of real or complex numbers, z, via
// the quadratic form
//
//   Q(a) = || a ||_2^2 + N | z^T a |^2,
//
// which is generated by the basis matrix
//
//   
//   B = [I; sqrt(N) z^T],
//
// as Q(a) = a^T B^T B a = || B a ||_2^2. Cohen has advice for the choice of
// the (large) parameter N within subsection 2.7.2 within his book. 
//
// The number of (nearly) exact Z-dependences detected is returned.
//
template<typename F>
Int ZDependenceSearch
( const Matrix<F>& z,
        Base<F> NSqrt,
        Matrix<F>& B,
        Matrix<F>& U, 
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

// Search for an algebraic relation
// ================================
// Search for the (Gaussian) integer coefficients of a polynomial of alpha
// that is (nearly) zero.
template<typename F>
Int AlgebraicRelationSearch
( F alpha,
  Int n,
  Base<F> NSqrt,
  Matrix<F>& B,
  Matrix<F>& U, 
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

// Schnorr-Euchner enumeration
// ===========================

enum EnumType {
  FULL_ENUM,
  GNR_ENUM,
  YSPARSE_ENUM
};

template<typename Real>
struct EnumCtrl
{
    EnumType enumType=FULL_ENUM;

    bool disablePrecDrop=false;
    Real fudge=Real(1.5); // fudge factor for number of bits of precision

    bool time=false;
    bool progress=false;

    // For monitoring the core (bounded) enumeration procedure
    bool innerProgress=false;

    // Explicitly transpose 'N' to encourage unit-stride access
    bool explicitTranspose=true;

    // GNR_ENUM
    // --------
    // TODO: Add ability to further tune the bounding function
    bool linearBounding=false;
    Int numTrials=1000;

    // YSPARSE_ENUM
    // ------------
    Int phaseLength=10;
    double enqueueProb=1.;

    bool customStartIndex=false; 
    Int startIndex;

    bool customMinInfNorms=false;
    vector<Int> minInfNorms;

    bool customMaxInfNorms=false;
    vector<Int> maxInfNorms;

    bool customMinOneNorms=false;
    vector<Int> minOneNorms;

    bool customMaxOneNorms=false;
    vector<Int> maxOneNorms; 

    Int progressLevel=0;

    template<typename OtherReal>
    EnumCtrl<Real>& operator=( const EnumCtrl<OtherReal>& ctrl )
    {
        enumType = ctrl.enumType;
        disablePrecDrop = ctrl.disablePrecDrop;
        fudge = Real(ctrl.fudge);
        time = ctrl.time;
        progress = ctrl.progress;
        innerProgress = ctrl.innerProgress;
        explicitTranspose = ctrl.explicitTranspose;

        // GNR_ENUM
        // --------
        linearBounding = ctrl.linearBounding;
        numTrials = ctrl.numTrials;

        // YSPARSE_ENUM
        // ------------
        phaseLength = ctrl.phaseLength;
        enqueueProb = ctrl.enqueueProb;
        customStartIndex = ctrl.customStartIndex;
        startIndex = ctrl.startIndex;

        customMinInfNorms = ctrl.customMinInfNorms;
        minInfNorms = ctrl.minInfNorms;
        customMaxInfNorms = ctrl.customMaxInfNorms;
        maxInfNorms = ctrl.maxInfNorms;

        customMinOneNorms = ctrl.customMinOneNorms;
        minOneNorms = ctrl.minOneNorms;
        customMaxOneNorms = ctrl.customMaxOneNorms;
        maxOneNorms = ctrl.maxOneNorms;

        progressLevel = ctrl.progressLevel;

        return *this;
    }

    EnumCtrl() { }
    EnumCtrl( const EnumCtrl<Real>& ctrl ) { *this = ctrl; }
    template<typename OtherReal>
    EnumCtrl( const EnumCtrl<OtherReal>& ctrl ) { *this = ctrl; }
};

namespace svp {

// If successful, fills 'v' with the integer coordinates of the columns of 
// an m x n matrix B (represented by its n x n upper-triangular Gaussian Normal
// Form; the 'R' from the QR factorization) which had a norm profile
// underneath the vector 'u' of upper bounds (|| (B v)(0:j) ||_2 < u(j)).
// Notice that the inequalities are strict.
//
// If not successful, the return value is a value greater than u(n-1) and 
// the contents of 'v' should be ignored.
//
// NOTE: There is not currently a complex implementation, though algorithms
//       exist.
template<typename F>
Base<F> GNREnumeration
( const Matrix<Base<F>>& d,
  const Matrix<F>& N,
  const Matrix<Base<F>>& u,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );

// Convert to/from the so-called "y-sparse" representation of
//
//   Dan Ding, Guizhen Zhu, Yang Yu, and Zhongxiang Zheng,
//   "A fast phase-based enumeration algorithm for SVP challenge through
//    y-sparse representations of short lattice vectors"

template<typename F>
void CoordinatesToSparse
( const Matrix<F>& N, const Matrix<F>& v, Matrix<F>& y );
template<typename F>
void TransposedCoordinatesToSparse
( const Matrix<F>& NTrans, const Matrix<F>& v, Matrix<F>& y );

template<typename F>
void SparseToCoordinates
( const Matrix<F>& N, const Matrix<F>& y, Matrix<F>& v );
template<typename F>
void TransposedSparseToCoordinates
( const Matrix<F>& NTrans, const Matrix<F>& y, Matrix<F>& v );

template<typename F>
Base<F> CoordinatesToNorm
( const Matrix<Base<F>>& d, const Matrix<F>& N, const Matrix<F>& v );
template<typename F>
Base<F> TransposedCoordinatesToNorm
( const Matrix<Base<F>>& d, const Matrix<F>& NTrans, const Matrix<F>& v );

template<typename F>
Base<F> SparseToNorm
( const Matrix<Base<F>>& d, const Matrix<F>& N, const Matrix<F>& y );
template<typename F>
Base<F> TransposedSparseToNorm
( const Matrix<Base<F>>& d, const Matrix<F>& NTrans, const Matrix<F>& y );

template<typename Real>
Real PhaseEnumeration
( const Matrix<Real>& B,
  const Matrix<Real>& d,
  const Matrix<Real>& N,
        Real normUpperBound,
        Int startIndex,
        Int phaseLength,
        double enqueueProb,
  const vector<Int>& minInfNorms,
  const vector<Int>& maxInfNorms,
  const vector<Int>& minOneNorms,
  const vector<Int>& maxOneNorms,
        Matrix<Real>& v,
        Int progressLevel=0 );
template<typename Real>
std::pair<Real,Int>
PhaseEnumeration
( const Matrix<Real>& B,
  const Matrix<Real>& d,
  const Matrix<Real>& N,
  const Matrix<Real>& normUpperBounds,
        Int startIndex,
        Int phaseLength,
        double enqueueProb,
  const vector<Int>& minInfNorms,
  const vector<Int>& maxInfNorms,
  const vector<Int>& minOneNorms,
  const vector<Int>& maxOneNorms,
        Matrix<Real>& v,
        Int progressLevel=0 );

} // namespace svp

// Given a reduced lattice B and its Gaussian Normal Form, R, either find a
// member of the lattice (given by B v, with v the output) with norm less than
// the upper bound (and return its norm), or return a value greater than the
// upper bound.
template<typename F>
Base<F> ShortVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
        Base<F> normUpperBound,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );

template<typename F>
std::pair<Base<F>,Int>
MultiShortVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
  const Matrix<Base<F>>& normUpperBounds,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );

// Given a reduced lattice B and its Gaussian Normal Form, R, find the shortest
// member of the lattice (with the shortest vector given by B v).
//
// The return value is the norm of the (approximately) shortest vector.
template<typename F>
Base<F> ShortestVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );

// If an upper-bound on the shortest vector which is better than || b_0 ||_2 is
// available
template<typename F>
Base<F> ShortestVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
        Base<F> normUpperBound,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );

template<typename F>
std::pair<Base<F>,Int>
MultiShortestVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
  const Matrix<Base<F>>& normUpperBounds,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );

// If a shorter vector is found, insert it into the first position
// ---------------------------------------------------------------

// The return value is the norm of the (approximately) shortest vector.
template<typename F>
Base<F> ShortestVectorEnrichment
(       Matrix<F>& B,
  const Matrix<F>& R,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );
template<typename F>
Base<F> ShortestVectorEnrichment
(       Matrix<F>& B,
        Matrix<F>& U,
  const Matrix<F>& R,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );

// If an upper-bound on the shortest vector which is better than || b_0 ||_2 is
// available
template<typename F>
Base<F> ShortestVectorEnrichment
(       Matrix<F>& B,
  const Matrix<F>& R,
        Base<F> normUpperBound,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );
template<typename F>
Base<F> ShortestVectorEnrichment
(       Matrix<F>& B,
        Matrix<F>& U,
  const Matrix<F>& R,
        Base<F> normUpperBound,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );

template<typename F>
std::pair<Base<F>,Int>
MultiShortestVectorEnrichment
(       Matrix<F>& B,
  const Matrix<F>& R,
  const Matrix<Base<F>>& normUpperBounds,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );
template<typename F>
std::pair<Base<F>,Int>
MultiShortestVectorEnrichment
(       Matrix<F>& B,
        Matrix<F>& U,
  const Matrix<F>& R,
  const Matrix<Base<F>>& normUpperBounds,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );

// Block Korkin-Zolotarev (BKZ) reduction
// ======================================
// TODO: Tailor BKZInfo; it is currently a copy of LLLInfo
template<typename Real>
struct BKZInfo
{
    Real delta;
    Real eta; 
    Int rank;
    Int nullity; 
    Int numSwaps;
    Int numEnums;
    Int numEnumFailures;
    Real logVol;

    template<typename OtherReal>
    BKZInfo<Real>& operator=( const BKZInfo<OtherReal>& info )
    {
        delta = Real(info.delta);
        eta = Real(info.eta);
        rank = info.rank;
        nullity = info.nullity;
        numSwaps = info.numSwaps;
        numEnums = info.numEnums;
        numEnumFailures = info.numEnumFailures;
        logVol = Real(info.logVol);
        return *this;
    }

    BKZInfo() { }
    BKZInfo( const BKZInfo<Real>& ctrl ) { *this = ctrl; }
    template<typename OtherReal>
    BKZInfo( const BKZInfo<OtherReal>& ctrl ) { *this = ctrl; }
};

template<typename Real>
struct BKZCtrl
{
    Int blocksize=20;
    bool time=false;
    bool progress=false;

    bool earlyAbort=false;
    Int numEnumsBeforeAbort=1000; // only used if earlyAbort=true

    bool variableBlocksize=false;
    function<Int(Int)> blocksizeFunc;

    bool variableEnumType=false;
    function<EnumType(Int)> enumTypeFunc;

    // Y-sparse enumeration supports simultaneous searches for improving a
    // contiguous window of vectors
    Int multiEnumWindow=15;

    bool skipInitialLLL=false;
    bool jumpstart=false;
    Int startCol=0;

    EnumCtrl<Real> enumCtrl;

    // Rather than running LLL after a productive enumeration, one could run
    // BKZ with a smaller blocksize (perhaps with early abort)
    bool subBKZ=true;
    function<Int(Int)> subBlocksizeFunc =
      function<Int(Int)>( []( Int bsize )
      { return Max(Min(bsize/2,Int(20)),Int(2)); } );
    bool subEarlyAbort = true;
    Int subNumEnumsBeforeAbort = 100;

    // This seems to be *more* expensive but lead to higher quality (perhaps).
    // Note that this is different from GNR recursion and uses a tree method
    // with shuffling at each merge.
    bool recursive=false;

    bool logFailedEnums=false;
    std::string failedEnumFile="BKZFailedEnums.txt";

    bool logStreakSizes=false;
    std::string streakSizesFile="BKZStreakSizes.txt";

    bool logNontrivialCoords=false;
    std::string nontrivialCoordsFile="BKZNontrivialCoords.txt";

    bool logNorms=false;
    std::string normsFile="BKZNorms.txt";

    bool logProjNorms=false;
    std::string projNormsFile="BKZProjNorms.txt";

    bool checkpoint=false;
    FileFormat checkpointFormat=ASCII;
    std::string checkpointFileBase="BKZCheckpoint";
    std::string tourFileBase="BKZTour";

    LLLCtrl<Real> lllCtrl;

    // We frequently need to convert datatypes, so make this easy
    template<typename OtherReal>
    BKZCtrl<Real>& operator=( const BKZCtrl<OtherReal>& ctrl )
    {
        blocksize = ctrl.blocksize;
        time = ctrl.time;
        progress = ctrl.progress;

        earlyAbort = ctrl.earlyAbort;
        numEnumsBeforeAbort = ctrl.numEnumsBeforeAbort;

        variableBlocksize = ctrl.variableBlocksize;
        blocksizeFunc = ctrl.blocksizeFunc;

        variableEnumType = ctrl.variableEnumType;
        enumTypeFunc = ctrl.enumTypeFunc;

        multiEnumWindow = ctrl.multiEnumWindow;

        skipInitialLLL = ctrl.skipInitialLLL;
        jumpstart = ctrl.jumpstart;
        startCol = ctrl.startCol;

        enumCtrl = ctrl.enumCtrl;

        subBKZ = ctrl.subBKZ;
        subBlocksizeFunc = ctrl.subBlocksizeFunc;
        subEarlyAbort = ctrl.subEarlyAbort;
        subNumEnumsBeforeAbort = ctrl.subNumEnumsBeforeAbort;

        recursive = ctrl.recursive;

        logFailedEnums = ctrl.logFailedEnums;
        logStreakSizes = ctrl.logStreakSizes;
        logNontrivialCoords = ctrl.logNontrivialCoords;
        logNorms = ctrl.logNorms;
        logProjNorms = ctrl.logProjNorms;
        checkpoint = ctrl.checkpoint;
        failedEnumFile = ctrl.failedEnumFile;
        streakSizesFile = ctrl.streakSizesFile;
        nontrivialCoordsFile = ctrl.nontrivialCoordsFile;
        normsFile = ctrl.normsFile;
        projNormsFile = ctrl.projNormsFile;
        checkpointFileBase = ctrl.checkpointFileBase;
        tourFileBase = ctrl.tourFileBase;
        checkpointFormat = ctrl.checkpointFormat;

        lllCtrl = ctrl.lllCtrl;
        return *this;
    }
    
    BKZCtrl() { }
    BKZCtrl( const BKZCtrl<Real>& ctrl ) { *this = ctrl; }
    template<typename OtherReal>
    BKZCtrl( const BKZCtrl<OtherReal>& ctrl ) { *this = ctrl; }
};

template<typename F>
BKZInfo<Base<F>> BKZ
( Matrix<F>& B,
  const BKZCtrl<Base<F>>& ctrl=BKZCtrl<Base<F>>() );

template<typename F>
BKZInfo<Base<F>> BKZ
( Matrix<F>& B,
  Matrix<F>& R,
  const BKZCtrl<Base<F>>& ctrl=BKZCtrl<Base<F>>() );

template<typename F>
BKZInfo<Base<F>> BKZ
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& R,
  const BKZCtrl<Base<F>>& ctrl=BKZCtrl<Base<F>>() );

template<typename F>
BKZInfo<Base<F>> BKZWithQ
( Matrix<F>& B,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const BKZCtrl<Base<F>>& ctrl=BKZCtrl<Base<F>>() );

template<typename F>
BKZInfo<Base<F>> BKZWithQ
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const BKZCtrl<Base<F>>& ctrl=BKZCtrl<Base<F>>() );

// Lattice coordinates
// ===================
// Seek the coordinates x in Z^n of a vector y within a lattice B, i.e.,
//
//     B x = y.
//
// Return 'true' if such coordinates could be found and 'false' otherwise.
template<typename F>
bool LatticeCoordinates( const Matrix<F>& B, const Matrix<F>& y, Matrix<F>& x );

// Enrich a lattice with a particular vector
// =========================================
// Push B v into the first column of B via a unimodular transformation
template<typename F>
void EnrichLattice( Matrix<F>& B, const Matrix<F>& v );
template<typename F>
void EnrichLattice( Matrix<F>& B, Matrix<F>& U, const Matrix<F>& v );

// Closest vector problem
// ======================

template<typename F>
void NearestPlane
( const Matrix<F>& B,
  const Matrix<F>& T,
        Matrix<F>& X,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

template<typename F>
void NearestPlane
( const Matrix<F>& B,
  const Matrix<F>& QR,
  const Matrix<F>& t,
  const Matrix<Base<F>>& d,
  const Matrix<F>& T,
        Matrix<F>& X,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

} // namespace El

#include "El/lattice/LLL.hpp"
#include "El/lattice/NearestPlane.hpp"
#include "El/lattice/Enrich.hpp"
#include "El/lattice/BKZ.hpp"

#endif // ifndef EL_LATTICE_HPP
