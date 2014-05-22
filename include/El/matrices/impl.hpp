/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_MATRICES_IMPL_HPP
#define EL_MATRICES_IMPL_HPP

// Deterministic matrices
// ======================

#include "./BullsHead.hpp"
#include "./Cauchy.hpp"
#include "./CauchyLike.hpp"
#include "./Circulant.hpp"
#include "./Demmel.hpp"
#include "./Diagonal.hpp"
#include "./Egorov.hpp"
#include "./Ehrenfest.hpp"
#include "./ExtendedKahan.hpp"
#include "./Fiedler.hpp"
#include "./Forsythe.hpp"
#include "./FoxLi.hpp"
#include "./Fourier.hpp"
#include "./GCDMatrix.hpp"
#include "./Gear.hpp"
#include "./GKS.hpp"
#include "./Grcar.hpp"
#include "./Hankel.hpp"
#include "./Hanowa.hpp"
#include "./HatanoNelson.hpp"
#include "./Helmholtz.hpp"
#include "./HelmholtzPML.hpp"
#include "./HermitianFromEVD.hpp"
#include "./Hilbert.hpp"
#include "./Identity.hpp"
#include "./Jordan.hpp"
#include "./Kahan.hpp"
#include "./KMS.hpp"
#include "./Laplacian.hpp"
#include "./Lauchli.hpp"
#include "./Legendre.hpp"
#include "./Lehmer.hpp"
#include "./Lotkin.hpp"
#include "./MinIJ.hpp"
#include "./NormalFromEVD.hpp"
#include "./Ones.hpp"
#include "./OneTwoOne.hpp"
#include "./Parter.hpp"
#include "./Pei.hpp"
#include "./Redheffer.hpp"
#include "./Riemann.hpp"
#include "./Riffle.hpp"
#include "./Ris.hpp"
#include "./Toeplitz.hpp"
#include "./Trefethen.hpp"
#include "./Triangle.hpp"
#include "./TriW.hpp"
#include "./Walsh.hpp"
#include "./Whale.hpp"
#include "./Wilkinson.hpp"
#include "./Zeros.hpp"

// Random matrices
// ===============

// Uniform
#include "./HermitianUniformSpectrum.hpp"
#include "./NormalUniformSpectrum.hpp"
#include "./Uniform.hpp"
#include "./UniformHelmholtzGreens.hpp"

// Gaussian
#include "./Gaussian.hpp"
#include "./Wigner.hpp"
#include "./Haar.hpp"

#endif // ifndef EL_MATRICES_IMPL_HPP
