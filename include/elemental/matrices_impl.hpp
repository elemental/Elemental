/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_IMPL_HPP
#define MATRICES_IMPL_HPP

//
// Deterministic
//

#include "./matrices/Cauchy.hpp"
#include "./matrices/CauchyLike.hpp"
#include "./matrices/Circulant.hpp"
#include "./matrices/Diagonal.hpp"
#include "./matrices/DiscreteFourier.hpp"
#include "./matrices/Hankel.hpp"
#include "./matrices/Hilbert.hpp"
#include "./matrices/Identity.hpp"
#include "./matrices/Kahan.hpp"
#include "./matrices/Legendre.hpp"
#include "./matrices/Ones.hpp"
#include "./matrices/OneTwoOne.hpp"
#include "./matrices/Toeplitz.hpp"
#include "./matrices/Walsh.hpp"
#include "./matrices/Wilkinson.hpp"
#include "./matrices/Zeros.hpp"

//
// Random
//

#include "./matrices/Uniform.hpp"
#include "./matrices/HermitianUniformSpectrum.hpp"
#include "./matrices/NormalUniformSpectrum.hpp"

// TODO: Gaussian

#endif // ifndef MATRICES_IMPL_HPP
