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

// LAPACK-like routines
// ====================
#define EL_LAPACK_INC      "El/lapack-like.hpp"
#define EL_LAPACK_DECL_INC "El/lapack-like/decl.hpp"
#define EL_LAPACK_IMPL_INC "El/lapack-like/impl.hpp"

// Matrix properties
// -----------------
#define EL_PROPS_INC      "El/lapack-like/props.hpp"
#define EL_PROPS_DECL_INC "El/lapack-like/props/decl.hpp"
#define EL_PROPS_IMPL_INC "El/lapack-like/props/impl.hpp"

#define EL_CONDITION_INC \
  "El/lapack-like/props/Condition.hpp"
#define EL_DETERMINANT_INC \
  "El/lapack-like/props/Determinant.hpp"
#define EL_INERTIA_INC \
  "El/lapack-like/props/Inertia.hpp"
#define EL_NORM_INC \
  "El/lapack-like/props/Norm.hpp"
#define EL_PSEUDOSPECTRUM_INC \
  "El/lapack-like/props/Pseudospectrum.hpp"
#define EL_TRACE_INC \
  "El/lapack-like/props/Trace.hpp"
// Specific condition-number routines
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#define EL_FROBENIUSCONDITION_INC \
  "El/lapack-like/props/Condition/Frobenius.hpp"
#define EL_INFINITYCONDITION_INC \
  "El/lapack-like/props/Condition/Infinity.hpp"
#define EL_MAXCONDITION_INC \
  "El/lapack-like/props/Condition/Max.hpp"
#define EL_ONECONDITION_INC \
  "El/lapack-like/props/Condition/One.hpp"
#define EL_TWOCONDITION_INC \
  "El/lapack-like/props/Condition/Two.hpp"
// Specific norm routines
// ^^^^^^^^^^^^^^^^^^^^^^
#define EL_ENTRYWISENORM_INC \
  "El/lapack-like/props/Norm/Entrywise.hpp"
#define EL_ENTRYWISEONENORM_INC \
  "El/lapack-like/props/Norm/EntrywiseOne.hpp"
#define EL_FROBENIUSNORM_INC \
  "El/lapack-like/props/Norm/Frobenius.hpp"
#define EL_INFINITYNORM_INC \
  "El/lapack-like/props/Norm/Infinity.hpp"
#define EL_KYFANNORM_INC \
  "El/lapack-like/props/Norm/KyFan.hpp"
#define EL_MAXNORM_INC \
  "El/lapack-like/props/Norm/Max.hpp"
#define EL_NUCLEARNORM_INC \
  "El/lapack-like/props/Norm/Nuclear.hpp"
#define EL_ONENORM_INC \
  "El/lapack-like/props/Norm/One.hpp"
#define EL_SCHATTENNORM_INC \
  "El/lapack-like/props/Norm/Schatten.hpp"
#define EL_TWONORMESTIMATE_INC \
  "El/lapack-like/props/Norm/TwoEstimate.hpp"
#define EL_TWONORM_INC \
  "El/lapack-like/props/Norm/Two.hpp"
#define EL_ZERONORM_INC \
  "El/lapack-like/props/Norm/Zero.hpp"

// Special matrices
// ================
#define EL_MATRICES_INC      "El/matrices.hpp"
#define EL_MATRICES_DECL_INC "El/matrices/decl.hpp"
#define EL_MATRICES_IMPL_INC "El/matrices/impl.hpp"

// Deterministic matrices
// ----------------------
#define EL_BULLSHEAD_INC        "El/matrices/BullsHead.hpp"
#define EL_CAUCHY_INC           "El/matrices/Cauchy.hpp"
#define EL_CAUCHYLIKE_INC       "El/matrices/CauchyLike.hpp"
#define EL_CIRCULANT_INC        "El/matrices/Circulant.hpp"
#define EL_DEMMEL_INC           "El/matrices/Demmel.hpp"
#define EL_DIAGONAL_INC         "El/matrices/Diagonal.hpp"
#define EL_EGOROV_INC           "El/matrices/Egorov.hpp"
#define EL_EHRENFEST_INC        "El/matrices/Ehrenfest.hpp"
#define EL_EXTENDEDKAHAN_INC    "El/matrices/ExtendedKahan.hpp"
#define EL_FIEDLER_INC          "El/matrices/Fiedler.hpp"
#define EL_FORSYTHE_INC         "El/matrices/Forsythe.hpp"
#define EL_FOURIER_INC          "El/matrices/Fourier.hpp"
#define EL_FOXLI_INC            "El/matrices/FoxLi.hpp"
#define EL_GCDMATRIX_INC        "El/matrices/GCDMatrix.hpp"
#define EL_GEAR_INC             "El/matrices/Gear.hpp"
#define EL_GKS_INC              "El/matrices/GKS.hpp"
#define EL_GRCAR_INC            "El/matrices/Grcar.hpp"
#define EL_HANKEL_INC           "El/matrices/Hankel.hpp"
#define EL_HANOWA_INC           "El/matrices/Hanowa.hpp"
#define EL_HELMHOLTZ_INC        "El/matrices/Helmholtz.hpp"
#define EL_HELMHOLTZPML_INC     "El/matrices/HelmholtzPML.hpp"
#define EL_HERMITIANFROMEVD_INC "El/matrices/HermitianFromEVD.hpp"
#define EL_HILBERT_INC          "El/matrices/Hilbert.hpp"
#define EL_IDENTITY_INC         "El/matrices/Identity.hpp"
#define EL_JORDAN_INC           "El/matrices/Jordan.hpp"
#define EL_KAHAN_INC            "El/matrices/Kahan.hpp"
#define EL_KMS_INC              "El/matrices/KMS.hpp"
#define EL_LAPLACIAN_INC        "El/matrices/Laplacian.hpp"
#define EL_LAUCHLI_INC          "El/matrices/Lauchli.hpp"
#define EL_LEGENDRE_INC         "El/matrices/Legendre.hpp"
#define EL_LEHMER_INC           "El/matrices/Lehmer.hpp"
#define EL_LOTKIN_INC           "El/matrices/Lotkin.hpp"
#define EL_MINIJ_INC            "El/matrices/MinIJ.hpp"
#define EL_NORMALFROMEVD_INC    "El/matrices/NormalFromEVD.hpp"
#define EL_ONES_INC             "El/matrices/Ones.hpp"
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
#define EL_ZEROS_INC            "El/matrices/Zeros.hpp"

// Random matrices
// ---------------
#define EL_GAUSSIAN_INC \
  "El/matrices/Gaussian.hpp"
#define EL_HAAR_INC \
  "El/matrices/Haar.hpp"
#define EL_HATANONELSON_INC \
  "El/matrices/HatanoNelson.hpp"
#define EL_HERMITIANUNIFORMSPECTRUM_INC \
  "El/matrices/HermitianUniformSpectrum.hpp"
#define EL_NORMALUNIFORMSPECTRUM_INC \
  "El/matrices/NormalUniformSpectrum.hpp"
#define EL_UNIFORM_INC \
  "El/matrices/Uniform.hpp"
#define EL_UNIFORMHELMHOLTZGREENS_INC \
  "El/matrices/UniformHelmholtzGreens.hpp"
#define EL_WIGNER_INC \
  "El/matrices/Wigner.hpp"

// Convex optimization
// ===================
#define EL_CONVEX_INC      "El/convex.hpp"
#define EL_CONVEX_DECL_INC "El/convex/decl.hpp"
#define EL_CONVEX_IMPL_INC "El/convex/impl.hpp"

// Utilities
// ---------
#define EL_CLIP_INC          "El/convex/Clip.hpp"
#define EL_COVARIANCE_INC    "El/convex/Covariance.hpp"
#define EL_LOGBARRIER_INC    "El/convex/LogBarrier.hpp"
#define EL_LOGDETDIV_INC     "El/convex/LogDetDiv.hpp"
#define EL_SOFTTHRESHOLD_INC "El/convex/SoftThreshold.hpp"
#define EL_SVT_INC           "El/convex/SVT.hpp"
// Specific SVT routines
// ^^^^^^^^^^^^^^^^^^^^^
#define EL_SVT_CROSS_INC     "El/convex/SVT/Cross.hpp"
#define EL_SVT_NORMAL_INC    "El/convex/SVT/Normal.hpp"
#define EL_SVT_PIVOTEDQR_INC "El/convex/SVT/PivotedQR.hpp"
#define EL_SVT_TSQR_INC      "El/convex/SVT/TSQR.hpp"

// ADMM
// ----
#define EL_BASISPURSUIT_INC     "El/convex/BasisPursuit.hpp"
#define EL_LINEARPROGRAM_INC    "El/convex/LinearProgram.hpp"
#define EL_QUADRATICPROGRAM_INC "El/convex/QuadraticProgram.hpp"
#define EL_SPARSEINVCOV_INC     "El/convex/SparseInvCov.hpp"

#endif // ifndef EL_INCLUDEPATHS_HPP
