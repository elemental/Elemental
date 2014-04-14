/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_PSEUDOSPECTRUM_UTIL_HPP
#define ELEM_PSEUDOSPECTRUM_UTIL_HPP

#include ELEM_MULTISHIFTQUASITRSM_INC
#include ELEM_MULTISHIFTTRSM_INC
#include ELEM_MULTISHIFTHESSSOLVE_INC
#include ELEM_ZERONORM_INC
#include ELEM_ONES_INC

#include "./Util/Rearrange.hpp"
#include "./Util/BasicMath.hpp"
#include "./Util/Snapshot.hpp"

namespace elem {

template<typename Real>
struct PseudospecCtrl
{
    // Preprocessing configuration
    bool schur; // begin with reduction to Schur form?
    bool forceComplexSchur;
    bool forceComplexPs;
    SdcCtrl<Real> sdcCtrl;

    // Convergence and deflation criteria
    Int maxIts;
    Real tol;
    bool deflate; 

    // (Implicitly Restarted) Arnoldi/Lanczos. If basisSize > 1, then
    // there is implicit restarting
    bool arnoldi;
    Int basisSize;
    bool reorthog; // only matters for IRL, which isn't currently used

    // Whether or not to print progress information at each iteration
    bool progress;

    SnapshotCtrl snapCtrl;

    PseudospecCtrl()
    : schur(true), forceComplexSchur(false), forceComplexPs(false), sdcCtrl(),
      maxIts(200), tol(1e-6), deflate(true),
      arnoldi(true), basisSize(10), reorthog(true),
      progress(false), snapCtrl()
    { }
};

} // namespace elem

#endif // ifndef ELEM_PSEUDOSPECTRUM_UTIL_HPP
