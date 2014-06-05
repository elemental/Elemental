/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_PSEUDOSPECTRUM_UTIL_HPP
#define EL_PSEUDOSPECTRUM_UTIL_HPP

#include EL_MULTISHIFTHESSSOLVE_INC
#include EL_ZERONORM_INC
#include EL_ONES_INC

#include "./Util/Rearrange.hpp"
#include "./Util/BasicMath.hpp"
#include "./Util/Snapshot.hpp"

namespace El {

enum PseudospecNorm {
  PS_TWO_NORM,
  PS_ONE_NORM
  /* For now, handle the infinity norm by using the adjoint matrix */
};

template<typename Real>
struct PseudospecCtrl
{
    PseudospecNorm norm;
    Int blockWidth; // block width for block 1-norm estimator

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
    : norm(PS_TWO_NORM), blockWidth(10),
      schur(true), forceComplexSchur(false), forceComplexPs(false), sdcCtrl(),
      maxIts(200), tol(1e-6), deflate(true),
      arnoldi(true), basisSize(10), reorthog(true),
      progress(false), snapCtrl()
    { }
};

} // namespace El

#endif // ifndef EL_PSEUDOSPECTRUM_UTIL_HPP
