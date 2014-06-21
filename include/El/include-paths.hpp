/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_INCLUDEPATHS_HPP
#define EL_INCLUDEPATHS_HPP

// Special matrices
// ================
#define EL_MATRICES_INC      "El/matrices.hpp"
#define EL_MATRICES_DECL_INC "El/matrices/decl.hpp"
#define EL_MATRICES_IMPL_INC "El/matrices/impl.hpp"

// Deterministic matrices
// ----------------------
#define EL_EGOROV_INC           "El/matrices/Egorov.hpp"
#define EL_ONETWOONE_INC        "El/matrices/OneTwoOne.hpp"
#define EL_PARTER_INC           "El/matrices/Parter.hpp"
#define EL_PEI_INC              "El/matrices/Pei.hpp"
#define EL_REDHEFFER_INC        "El/matrices/Redheffer.hpp"
#define EL_RIEMANN_INC          "El/matrices/Riemann.hpp"
#define EL_RIFFLE_INC           "El/matrices/Riffle.hpp"
#define EL_RIS_INC              "El/matrices/Ris.hpp"
#define EL_TOEPLITZ_INC         "El/matrices/Toeplitz.hpp"
#define EL_TREFETHEN_INC        "El/matrices/Trefethen.hpp"
#define EL_TRIANGLE_INC         "El/matrices/Triangle.hpp"
#define EL_TRIW_INC             "El/matrices/TriW.hpp"
#define EL_WALSH_INC            "El/matrices/Walsh.hpp"
#define EL_WHALE_INC            "El/matrices/Whale.hpp"
#define EL_WILKINSON_INC        "El/matrices/Wilkinson.hpp"

// Random matrices
// ---------------
#define EL_NORMALUNIFORMSPECTRUM_INC \
  "El/matrices/NormalUniformSpectrum.hpp"
#define EL_UNIFORMHELMHOLTZGREENS_INC \
  "El/matrices/UniformHelmholtzGreens.hpp"

#endif // ifndef EL_INCLUDEPATHS_HPP
