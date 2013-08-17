/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_IMPL_HPP
#define ELEM_MATRICES_IMPL_HPP

//
// Deterministic
//

#include "./matrices/Cauchy.hpp"
#include "./matrices/CauchyLike.hpp"
#include "./matrices/Circulant.hpp"
#include "./matrices/Diagonal.hpp"
#include "./matrices/Egorov.hpp"
#include "./matrices/ExtendedKahan.hpp"
#include "./matrices/Fiedler.hpp"
#include "./matrices/Forsythe.hpp"
#include "./matrices/Fourier.hpp"
#include "./matrices/GCDMatrix.hpp"
#include "./matrices/Gear.hpp"
#include "./matrices/GKS.hpp"
#include "./matrices/Grcar.hpp"
#include "./matrices/Hankel.hpp"
#include "./matrices/Hanowa.hpp"
#include "./matrices/Helmholtz.hpp"
#include "./matrices/HermitianFromEVD.hpp"
#include "./matrices/Hilbert.hpp"
#include "./matrices/Identity.hpp"
#include "./matrices/Jordan.hpp"
#include "./matrices/Kahan.hpp"
#include "./matrices/KMS.hpp"
#include "./matrices/Laplacian.hpp"
#include "./matrices/Lauchli.hpp"
#include "./matrices/Legendre.hpp"
#include "./matrices/Lehmer.hpp"
#include "./matrices/Lotkin.hpp"
#include "./matrices/MinIJ.hpp"
#include "./matrices/NormalFromEVD.hpp"
#include "./matrices/Ones.hpp"
#include "./matrices/OneTwoOne.hpp"
#include "./matrices/Parter.hpp"
#include "./matrices/Pei.hpp"
#include "./matrices/Redheffer.hpp"
#include "./matrices/Riemann.hpp"
#include "./matrices/Ris.hpp"
#include "./matrices/Toeplitz.hpp"
#include "./matrices/TriW.hpp"
#include "./matrices/Walsh.hpp"
#include "./matrices/Wilkinson.hpp"
#include "./matrices/Zeros.hpp"

//
// Random
//

// Uniform
#include "./matrices/HermitianUniformSpectrum.hpp"
#include "./matrices/NormalUniformSpectrum.hpp"
#include "./matrices/Uniform.hpp"

// Gaussian
#include "./matrices/Gaussian.hpp"
#include "./matrices/Wigner.hpp"
#include "./matrices/Haar.hpp"

#endif // ifndef ELEM_MATRICES_IMPL_HPP
