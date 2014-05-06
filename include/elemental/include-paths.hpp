/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_INCLUDEPATHS_HPP
#define ELEM_INCLUDEPATHS_HPP

// BLAS-like routines
// ==================
#define ELEM_BLAS_INC      "elemental/blas-like.hpp"
#define ELEM_BLAS_DECL_INC "elemental/blas-like/decl.hpp"
#define ELEM_BLAS_IMPL_INC "elemental/blas-like/impl.hpp"

// Level 1 routines
// ----------------
#define ELEM_BLAS1_INC "elemental/blas-like/level1.hpp"

#define ELEM_ADJOINT_INC \
  "elemental/blas-like/level1/Adjoint.hpp"
#define ELEM_AXPY_INC \
  "elemental/blas-like/level1/Axpy.hpp"
#define ELEM_AXPYTRIANGLE_INC \
  "elemental/blas-like/level1/AxpyTriangle.hpp"
#define ELEM_CONJUGATE_INC \
  "elemental/blas-like/level1/Conjugate.hpp"
#define ELEM_COPY_INC \
  "elemental/blas-like/level1/Copy.hpp"
#define ELEM_DIAGONALSCALE_INC \
  "elemental/blas-like/level1/DiagonalScale.hpp"
#define ELEM_DIAGONALSCALETRAPEZOID_INC \
  "elemental/blas-like/level1/DiagonalScaleTrapezoid.hpp"
#define ELEM_DIAGONALSOLVE_INC \
  "elemental/blas-like/level1/DiagonalSolve.hpp"
#define ELEM_DOT_INC \
  "elemental/blas-like/level1/Dot.hpp"
#define ELEM_DOTU_INC \
  "elemental/blas-like/level1/Dotu.hpp"
#define ELEM_ENTRYWISEMAP_INC \
  "elemental/blas-like/level1/EntrywiseMap.hpp"
#define ELEM_HADAMARD_INC \
  "elemental/blas-like/level1/Hadamard.hpp"
#define ELEM_HILBERTSCHMIDT_INC \
  "elemental/blas-like/level1/HilbertSchmidt.hpp"
#define ELEM_MAKEHERMITIAN_INC \
  "elemental/blas-like/level1/MakeHermitian.hpp"
#define ELEM_MAKEREAL_INC \
  "elemental/blas-like/level1/MakeReal.hpp"
#define ELEM_MAKESYMMETRIC_INC \
  "elemental/blas-like/level1/MakeSymmetric.hpp"
#define ELEM_MAKETRAPEZOIDAL_INC \
  "elemental/blas-like/level1/MakeTrapezoidal.hpp"
#define ELEM_MAKETRIANGULAR_INC \
  "elemental/blas-like/level1/MakeTriangular.hpp"
#define ELEM_MAX_INC \
  "elemental/blas-like/level1/Max.hpp"
#define ELEM_MAXABS_INC \
  "elemental/blas-like/level1/MaxAbs.hpp"
#define ELEM_MIN_INC \
  "elemental/blas-like/level1/Min.hpp"
#define ELEM_MINABS_INC \
  "elemental/blas-like/level1/MinAbs.hpp"
#define ELEM_NRM2_INC \
  "elemental/blas-like/level1/Nrm2.hpp"
#define ELEM_QUASIDIAGONALSCALE_INC \
  "elemental/blas-like/level1/QuasiDiagonalScale.hpp"
#define ELEM_QUASIDIAGONALSOLVE_INC \
  "elemental/blas-like/level1/QuasiDiagonalSolve.hpp"
#define ELEM_SCALE_INC \
  "elemental/blas-like/level1/Scale.hpp"
#define ELEM_SCALETRAPEZOID_INC \
  "elemental/blas-like/level1/ScaleTrapezoid.hpp"
#define ELEM_SETDIAGONAL_INC \
  "elemental/blas-like/level1/SetDiagonal.hpp"
#define ELEM_SWAP_INC \
  "elemental/blas-like/level1/Swap.hpp"
#define ELEM_SYMMETRIC2X2INV_INC \
  "elemental/blas-like/level1/Symmetric2x2Inv.hpp"
#define ELEM_SYMMETRIC2X2SCALE_INC \
  "elemental/blas-like/level1/Symmetric2x2Scale.hpp"
#define ELEM_SYMMETRIC2X2SOLVE_INC \
  "elemental/blas-like/level1/Symmetric2x2Solve.hpp"
#define ELEM_TRANSPOSE_INC \
  "elemental/blas-like/level1/Transpose.hpp"
#define ELEM_UPDATEDIAGONAL_INC \
  "elemental/blas-like/level1/UpdateDiagonal.hpp"
#define ELEM_ZERO_INC \
  "elemental/blas-like/level1/Zero.hpp"

// Level 2 routines
// ----------------
#define ELEM_BLAS2_INC "elemental/blas-like/level2.hpp"

#define ELEM_GEMV_INC \
  "elemental/blas-like/level2/Gemv.hpp"
#define ELEM_GER_INC \
  "elemental/blas-like/level2/Ger.hpp"
#define ELEM_GERU_INC \
  "elemental/blas-like/level2/Geru.hpp"
#define ELEM_HEMV_INC \
  "elemental/blas-like/level2/Hemv.hpp"
#define ELEM_HER2_INC \
  "elemental/blas-like/level2/Her2.hpp"
#define ELEM_HER_INC \
  "elemental/blas-like/level2/Her.hpp"
#define ELEM_QUASITRSV_INC \
  "elemental/blas-like/level2/QuasiTrsv.hpp"
#define ELEM_SYMV_INC \
  "elemental/blas-like/level2/Symv.hpp"
#define ELEM_SYR2_INC \
  "elemental/blas-like/level2/Syr2.hpp"
#define ELEM_SYR_INC \
  "elemental/blas-like/level2/Syr.hpp"
#define ELEM_TRMV_INC \
  "elemental/blas-like/level2/Trmv.hpp"
#define ELEM_TRR2_INC \
  "elemental/blas-like/level2/Trr2.hpp"
#define ELEM_TRR_INC \
  "elemental/blas-like/level2/Trr.hpp"
#define ELEM_TRSV_INC \
  "elemental/blas-like/level2/Trsv.hpp"

// Level 3 routines
// ----------------
#define ELEM_BLAS3_INC "elemental/blas-like/level3.hpp"

#define ELEM_GEMM_INC         "elemental/blas-like/level3/Gemm.hpp"
#define ELEM_HEMM_INC         "elemental/blas-like/level3/Hemm.hpp"
#define ELEM_HER2K_INC        "elemental/blas-like/level3/Her2k.hpp"
#define ELEM_HERK_INC         "elemental/blas-like/level3/Herk.hpp"
#define ELEM_MULTISHIFTQUASITRSM_INC \
                            "elemental/blas-like/level3/MultiShiftQuasiTrsm.hpp"
#define ELEM_MULTISHIFTTRSM_INC \
                            "elemental/blas-like/level3/MultiShiftTrsm.hpp"
#define ELEM_QUASITRSM_INC    "elemental/blas-like/level3/QuasiTrsm.hpp"
#define ELEM_SYMM_INC         "elemental/blas-like/level3/Symm.hpp"
#define ELEM_SYR2K_INC        "elemental/blas-like/level3/Syr2k.hpp"
#define ELEM_SYRK_INC         "elemental/blas-like/level3/Syrk.hpp"
#define ELEM_TRDTRMM_INC      "elemental/blas-like/level3/Trdtrmm.hpp"
#define ELEM_TRMM_INC         "elemental/blas-like/level3/Trmm.hpp"
#define ELEM_TRSM_INC         "elemental/blas-like/level3/Trsm.hpp"
#define ELEM_TRSTRM_INC       "elemental/blas-like/level3/Trstrm.hpp"
#define ELEM_TRTRMM_INC       "elemental/blas-like/level3/Trtrmm.hpp"
#define ELEM_TWOSIDEDTRMM_INC "elemental/blas-like/level3/TwoSidedTrmm.hpp"
#define ELEM_TWOSIDEDTRSM_INC "elemental/blas-like/level3/TwoSidedTrsm.hpp"

// LAPACK-like routines
// ====================
#define ELEM_LAPACK_INC      "elemental/lapack-like.hpp"
#define ELEM_LAPACK_DECL_INC "elemental/lapack-like/decl.hpp"
#define ELEM_LAPACK_IMPL_INC "elemental/lapack-like/impl.hpp"

// Reduction to condensed form
// ---------------------------
#define ELEM_CONDENSE_INC      "elemental/lapack-like/condense.hpp"
#define ELEM_CONDENSE_DECL_INC "elemental/lapack-like/condense/decl.hpp"
#define ELEM_CONDENSE_IMPL_INC "elemental/lapack-like/condense/impl.hpp"

#define ELEM_BIDIAG_INC \
  "elemental/lapack-like/condense/Bidiag.hpp"
#define ELEM_HERMITIANTRIDIAG_INC \
  "elemental/lapack-like/condense/HermitianTridiag.hpp"
#define ELEM_HESSENBERG_INC \
  "elemental/lapack-like/condense/Hessenberg.hpp"

// Decompositions
// --------------
#define ELEM_DECOMP_INC      "elemental/lapack-like/decomp.hpp"
#define ELEM_DECOMP_DECL_INC "elemental/lapack-like/decomp/decl.hpp"
#define ELEM_DECOMP_IMPL_INC "elemental/lapack-like/decomp/impl.hpp"

#define ELEM_HERMITIANEIG_INC \
  "elemental/lapack-like/decomp/HermitianEig.hpp"
#define ELEM_HERMITIANGENDEFINITEEIG_INC \
  "elemental/lapack-like/decomp/HermitianGenDefiniteEig.hpp"
#define ELEM_HERMITIANTRIDIAGEIG_INC \
  "elemental/lapack-like/decomp/HermitianTridiagEig.hpp"
#define ELEM_POLAR_INC \
  "elemental/lapack-like/decomp/Polar.hpp"
#define ELEM_SCHUR_INC \
  "elemental/lapack-like/decomp/Schur.hpp"
#define ELEM_SKEWHERMITIANEIG_INC \
  "elemental/lapack-like/decomp/SkewHermitianEig.hpp"
#define ELEM_SVD_INC \
  "elemental/lapack-like/decomp/SVD.hpp"

// Specific HermitianEig routines
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#define ELEM_HERMITIANEIG_SDC_INC \
  "elemental/lapack-like/decomp/HermitianEig/SDC.hpp"

// Specific HermitianTridiagEig routines
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#define ELEM_HERMITIANTRIDIAGEIG_SORT_INC \
  "elemental/lapack-like/decomp/HermitianTridiagEig/Sort.hpp"

// Specific Polar routines
// ^^^^^^^^^^^^^^^^^^^^^^^
#define ELEM_POLAR_QDWH_INC \
  "elemental/lapack-like/decomp/Polar/QDWH.hpp"
#define ELEM_POLAR_SVD_INC \
  "elemental/lapack-like/decomp/Polar/SVD.hpp"

// Specific Schur routines
// ^^^^^^^^^^^^^^^^^^^^^^^
#define ELEM_SCHUR_INVERSEFREESDC_INC \
  "elemental/lapack-like/decomp/Schur/InverseFreeSDC.hpp"
#define ELEM_SCHUR_QR_INC \
  "elemental/lapack-like/decomp/Schur/QR.hpp"
#define ELEM_SCHUR_SDC_INC \
  "elemental/lapack-like/decomp/Schur/SDC.hpp"

// Factorizations
// --------------
#define ELEM_FACTOR_INC      "elemental/lapack-like/factor.hpp"
#define ELEM_FACTOR_DECL_INC "elemental/lapack-like/factor/decl.hpp"
#define ELEM_FACTOR_IMPL_INC "elemental/lapack-like/factor/impl.hpp"

#define ELEM_CHOLESKY_INC \
  "elemental/lapack-like/factor/Cholesky.hpp"
#define ELEM_LDL_INC \
  "elemental/lapack-like/factor/LDL.hpp"
#define ELEM_LU_INC \
  "elemental/lapack-like/factor/LU.hpp"

#define ELEM_LQ_INC \
  "elemental/lapack-like/factor/LQ.hpp"
#define ELEM_QR_INC \
  "elemental/lapack-like/factor/QR.hpp"
#define ELEM_RQ_INC \
  "elemental/lapack-like/factor/RQ.hpp"

#define ELEM_GQR_INC \
  "elemental/lapack-like/factor/GQR.hpp"
#define ELEM_GRQ_INC \
  "elemental/lapack-like/factor/GRQ.hpp"

#define ELEM_ID_INC \
  "elemental/lapack-like/factor/ID.hpp"
#define ELEM_SKELETON_INC \
  "elemental/lapack-like/factor/Skeleton.hpp"

// Specific LDL routines
// ^^^^^^^^^^^^^^^^^^^^^
#define ELEM_LDL_INERTIA_INC \
  "elemental/lapack-like/factor/LDL/Inertia.hpp"
#define ELEM_LDL_MULTIPLYAFTER_INC \
  "elemental/lapack-like/factor/LDL/MultiplyAfter.hpp"
#define ELEM_LDL_PIVOTED_INC \
  "elemental/lapack-like/factor/LDL/Pivoted.hpp"
#define ELEM_LDL_SOLVEAFTER_INC \
  "elemental/lapack-like/factor/LDL/SolveAfter.hpp"
#define ELEM_LDL_VAR3_INC \
  "elemental/lapack-like/factor/LDL/Var3.hpp"

// Specific LQ routines
// ^^^^^^^^^^^^^^^^^^^^
#define ELEM_LQ_APPLYQ_INC \
  "elemental/lapack-like/factor/LQ/ApplyQ.hpp"
#define ELEM_LQ_EXPLICIT_INC \
  "elemental/lapack-like/factor/LQ/Explicit.hpp"
#define ELEM_LQ_HOUSEHOLDER_INC \
  "elemental/lapack-like/factor/LQ/Householder.hpp"
#define ELEM_LQ_PANELHOUSEHOLDER_INC \
  "elemental/lapack-like/factor/LQ/PanelHouseholder.hpp"

// Specific LU routines
// ^^^^^^^^^^^^^^^^^^^^
#define ELEM_LU_FULL_INC \
  "elemental/lapack-like/factor/LU/Full.hpp"
#define ELEM_LU_LOCAL_INC \
  "elemental/lapack-like/factor/LU/Local.hpp"
#define ELEM_LU_PANEL_INC \
  "elemental/lapack-like/factor/LU/Panel.hpp"
#define ELEM_LU_SOLVEAFTER_INC \
  "elemental/lapack-like/factor/LU/SolveAfter.hpp"

// Specific QR routines
// ^^^^^^^^^^^^^^^^^^^^
#define ELEM_QR_APPLYQ_INC \
  "elemental/lapack-like/factor/QR/ApplyQ.hpp"
#define ELEM_QR_BUSINGERGOLUB_INC \
  "elemental/lapack-like/factor/QR/BusingerGolub.hpp"
#define ELEM_QR_CHOLESKY_INC \
  "elemental/lapack-like/factor/QR/Cholesky.hpp"
#define ELEM_QR_EXPLICIT_INC \
  "elemental/lapack-like/factor/QR/Explicit.hpp"
#define ELEM_QR_HOUSEHOLDER_INC \
  "elemental/lapack-like/factor/QR/Householder.hpp"
#define ELEM_QR_PANELHOUSEHOLDER_INC \
  "elemental/lapack-like/factor/QR/PanelHouseholder.hpp"
#define ELEM_QR_TS_INC \
  "elemental/lapack-like/factor/QR/TS.hpp"

// Specific RQ routines
// ^^^^^^^^^^^^^^^^^^^^
#define ELEM_RQ_APPLYQ_INC \
  "elemental/lapack-like/factor/RQ/ApplyQ.hpp"
#define ELEM_RQ_CHOLESKY_INC \
  "elemental/lapack-like/factor/RQ/Cholesky.hpp"
#define ELEM_RQ_HOUSEHOLDER_INC \
  "elemental/lapack-like/factor/RQ/Householder.hpp"
#define ELEM_RQ_PANELHOUSEHOLDER_INC \
  "elemental/lapack-like/factor/RQ/PanelHouseholder.hpp"

// Matrix functions
// ----------------
#define ELEM_FUNCS_INC      "elemental/lapack-like/funcs.hpp"
#define ELEM_FUNCS_DECL_INC "elemental/lapack-like/funcs/decl.hpp"
#define ELEM_FUNCS_IMPL_INC "elemental/lapack-like/funcs/impl.hpp"

#define ELEM_HERMITIANFUNCTION_INC \
  "elemental/lapack-like/funcs/HermitianFunction.hpp"
#define ELEM_INVERSE_INC \
  "elemental/lapack-like/funcs/Inverse.hpp"
#define ELEM_PSEUDOINVERSE_INC \
  "elemental/lapack-like/funcs/Pseudoinverse.hpp"
#define ELEM_SIGN_INC \
  "elemental/lapack-like/funcs/Sign.hpp"
#define ELEM_SQUAREROOT_INC \
  "elemental/lapack-like/funcs/SquareRoot.hpp"
// Specific inversion routines
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^
#define ELEM_TRIANGULARINVERSE_INC \
  "elemental/lapack-like/funcs/Inverse/Triangular.hpp"
#define ELEM_GENERALINVERSE_INC \
  "elemental/lapack-like/funcs/Inverse/General.hpp"
#define ELEM_HPDINVERSE_INC \
  "elemental/lapack-like/funcs/Inverse/HPD.hpp"
#define ELEM_SYMMETRICINVERSE_INC \
  "elemental/lapack-like/funcs/Inverse/Symmetric.hpp"
#define ELEM_HERMITIANINVERSE_INC \
  "elemental/lapack-like/funcs/Inverse/Hermitian.hpp"

// Permutations
// ------------
#define ELEM_APPLYCOLPIVOTS_INC \
  "elemental/lapack-like/perm/ApplyColPivots.hpp"
#define ELEM_APPLYROWPIVOTS_INC \
  "elemental/lapack-like/perm/ApplyRowPivots.hpp"
#define ELEM_APPLYSYMMETRICPIVOTS_INC \
  "elemental/lapack-like/perm/ApplySymmetricPivots.hpp"
#define ELEM_EXPLICITPERMUTATION_INC \
  "elemental/lapack-like/perm/ExplicitPermutation.hpp"
#define ELEM_INVERTPERMUTATION_INC \
  "elemental/lapack-like/perm/InvertPermutation.hpp"
#define ELEM_PERMUTATIONMETA_INC \
  "elemental/lapack-like/perm/PermutationMeta.hpp"
#define ELEM_PERMUTECOLS_INC \
  "elemental/lapack-like/perm/PermuteCols.hpp"
#define ELEM_PERMUTEROWS_INC \
  "elemental/lapack-like/perm/PermuteRows.hpp"
#define ELEM_PIVOTSTOPARTIALPERMUTATION_INC \
  "elemental/lapack-like/perm/PivotsToPartialPermutation.hpp"
#define ELEM_PIVOTSTOPERMUTATION_INC \
  "elemental/lapack-like/perm/PivotsToPermutation.hpp"

// Matrix properties
// -----------------
#define ELEM_PROPS_INC      "elemental/lapack-like/props.hpp"
#define ELEM_PROPS_DECL_INC "elemental/lapack-like/props/decl.hpp"
#define ELEM_PROPS_IMPL_INC "elemental/lapack-like/props/impl.hpp"

#define ELEM_CONDITION_INC \
  "elemental/lapack-like/props/Condition.hpp"
#define ELEM_DETERMINANT_INC \
  "elemental/lapack-like/props/Determinant.hpp"
#define ELEM_INERTIA_INC \
  "elemental/lapack-like/props/Inertia.hpp"
#define ELEM_NORM_INC \
  "elemental/lapack-like/props/Norm.hpp"
#define ELEM_PSEUDOSPECTRUM_INC \
  "elemental/lapack-like/props/Pseudospectrum.hpp"
#define ELEM_TRACE_INC \
  "elemental/lapack-like/props/Trace.hpp"
// Specific condition-number routines
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#define ELEM_FROBENIUSCONDITION_INC \
  "elemental/lapack-like/props/Condition/Frobenius.hpp"
#define ELEM_INFINITYCONDITION_INC \
  "elemental/lapack-like/props/Condition/Infinity.hpp"
#define ELEM_MAXCONDITION_INC \
  "elemental/lapack-like/props/Condition/Max.hpp"
#define ELEM_ONECONDITION_INC \
  "elemental/lapack-like/props/Condition/One.hpp"
#define ELEM_TWOCONDITION_INC \
  "elemental/lapack-like/props/Condition/Two.hpp"
// Specific norm routines
// ^^^^^^^^^^^^^^^^^^^^^^
#define ELEM_ENTRYWISENORM_INC \
  "elemental/lapack-like/props/Norm/Entrywise.hpp"
#define ELEM_ENTRYWISEONENORM_INC \
  "elemental/lapack-like/props/Norm/EntrywiseOne.hpp"
#define ELEM_FROBENIUSNORM_INC \
  "elemental/lapack-like/props/Norm/Frobenius.hpp"
#define ELEM_INFINITYNORM_INC \
  "elemental/lapack-like/props/Norm/Infinity.hpp"
#define ELEM_KYFANNORM_INC \
  "elemental/lapack-like/props/Norm/KyFan.hpp"
#define ELEM_MAXNORM_INC \
  "elemental/lapack-like/props/Norm/Max.hpp"
#define ELEM_NUCLEARNORM_INC \
  "elemental/lapack-like/props/Norm/Nuclear.hpp"
#define ELEM_ONENORM_INC \
  "elemental/lapack-like/props/Norm/One.hpp"
#define ELEM_SCHATTENNORM_INC \
  "elemental/lapack-like/props/Norm/Schatten.hpp"
#define ELEM_TWONORMESTIMATE_INC \
  "elemental/lapack-like/props/Norm/TwoEstimate.hpp"
#define ELEM_TWONORM_INC \
  "elemental/lapack-like/props/Norm/Two.hpp"
#define ELEM_ZERONORM_INC \
  "elemental/lapack-like/props/Norm/Zero.hpp"

// Solvers
// -------
#define ELEM_SOLVE_INC      "elemental/lapack-like/solve.hpp"
#define ELEM_SOLVE_DECL_INC "elemental/lapack-like/solve/decl.hpp"
#define ELEM_SOLVE_IMPL_INC "elemental/lapack-like/solve/impl.hpp"

#define ELEM_GAUSSIANELIMINATION_INC \
  "elemental/lapack-like/solve/GaussianElimination.hpp"
#define ELEM_HPDSOLVE_INC \
  "elemental/lapack-like/solve/HPDSolve.hpp"

#define ELEM_HERMITIANSOLVE_INC \
  "elemental/lapack-like/solve/HermitianSolve.hpp"
#define ELEM_SYMMETRICSOLVE_INC \
  "elemental/lapack-like/solve/SymmetricSolve.hpp"

#define ELEM_LEASTSQUARES_INC \
  "elemental/lapack-like/solve/LeastSquares.hpp"

#define ELEM_GLM_INC \
  "elemental/lapack-like/solve/GLM.hpp"
#define ELEM_LSE_INC \
  "elemental/lapack-like/solve/LSE.hpp"

#define ELEM_MULTISHIFTHESSSOLVE_INC \
  "elemental/lapack-like/solve/MultiShiftHessSolve.hpp"

// Utilities
// ---------
#define ELEM_LAPACKUTIL_INC      "elemental/lapack-like/util.hpp"
#define ELEM_LAPACKUTIL_DECL_INC "elemental/lapack-like/util/decl.hpp"
#define ELEM_LAPACKUTIL_IMPL_INC "elemental/lapack-like/util/impl.hpp"

#define ELEM_APPLYPACKEDREFLECTORS_INC \
  "elemental/lapack-like/util/ApplyPackedReflectors.hpp"
#define ELEM_EXPANDPACKEDREFLECTORS_INC \
  "elemental/lapack-like/util/ExpandPackedReflectors.hpp"
#define ELEM_HYPERBOLICREFLECTOR_INC \
  "elemental/lapack-like/util/HyperbolicReflector.hpp"
#define ELEM_MEDIAN_INC \
  "elemental/lapack-like/util/Median.hpp"
#define ELEM_PERMUTATIONPARITY_INC \
  "elemental/lapack-like/util/PermutationParity.hpp"
#define ELEM_PIVOTPARITY_INC \
  "elemental/lapack-like/util/PivotParity.hpp"
#define ELEM_REFLECTOR_INC \
  "elemental/lapack-like/util/Reflector.hpp"
#define ELEM_SORT_INC \
  "elemental/lapack-like/util/Sort.hpp"

// Special matrices
// ================
#define ELEM_MATRICES_INC      "elemental/matrices.hpp"
#define ELEM_MATRICES_DECL_INC "elemental/matrices/decl.hpp"
#define ELEM_MATRICES_IMPL_INC "elemental/matrices/impl.hpp"

// Deterministic matrices
// ----------------------
#define ELEM_BULLSHEAD_INC        "elemental/matrices/BullsHead.hpp"
#define ELEM_CAUCHY_INC           "elemental/matrices/Cauchy.hpp"
#define ELEM_CAUCHYLIKE_INC       "elemental/matrices/CauchyLike.hpp"
#define ELEM_CIRCULANT_INC        "elemental/matrices/Circulant.hpp"
#define ELEM_DEMMEL_INC           "elemental/matrices/Demmel.hpp"
#define ELEM_DIAGONAL_INC         "elemental/matrices/Diagonal.hpp"
#define ELEM_EGOROV_INC           "elemental/matrices/Egorov.hpp"
#define ELEM_EXTENDEDKAHAN_INC    "elemental/matrices/ExtendedKahan.hpp"
#define ELEM_FIEDLER_INC          "elemental/matrices/Fiedler.hpp"
#define ELEM_FORSYTHE_INC         "elemental/matrices/Forsythe.hpp"
#define ELEM_FOURIER_INC          "elemental/matrices/Fourier.hpp"
#define ELEM_FOXLI_INC            "elemental/matrices/FoxLi.hpp"
#define ELEM_GCDMATRIX_INC        "elemental/matrices/GCDMatrix.hpp"
#define ELEM_GEAR_INC             "elemental/matrices/Gear.hpp"
#define ELEM_GKS_INC              "elemental/matrices/GKS.hpp"
#define ELEM_GRCAR_INC            "elemental/matrices/Grcar.hpp"
#define ELEM_HANKEL_INC           "elemental/matrices/Hankel.hpp"
#define ELEM_HANOWA_INC           "elemental/matrices/Hanowa.hpp"
#define ELEM_HELMHOLTZ_INC        "elemental/matrices/Helmholtz.hpp"
#define ELEM_HELMHOLTZPML_INC     "elemental/matrices/HelmholtzPML.hpp"
#define ELEM_HERMITIANFROMEVD_INC "elemental/matrices/HermitianFromEVD.hpp"
#define ELEM_HILBERT_INC          "elemental/matrices/Hilbert.hpp"
#define ELEM_IDENTITY_INC         "elemental/matrices/Identity.hpp"
#define ELEM_JORDAN_INC           "elemental/matrices/Jordan.hpp"
#define ELEM_KAHAN_INC            "elemental/matrices/Kahan.hpp"
#define ELEM_KMS_INC              "elemental/matrices/KMS.hpp"
#define ELEM_LAPLACIAN_INC        "elemental/matrices/Laplacian.hpp"
#define ELEM_LAUCHLI_INC          "elemental/matrices/Lauchli.hpp"
#define ELEM_LEGENDRE_INC         "elemental/matrices/Legendre.hpp"
#define ELEM_LEHMER_INC           "elemental/matrices/Lehmer.hpp"
#define ELEM_LOTKIN_INC           "elemental/matrices/Lotkin.hpp"
#define ELEM_MINIJ_INC            "elemental/matrices/MinIJ.hpp"
#define ELEM_NORMALFROMEVD_INC    "elemental/matrices/NormalFromEVD.hpp"
#define ELEM_ONES_INC             "elemental/matrices/Ones.hpp"
#define ELEM_ONETWOONE_INC        "elemental/matrices/OneTwoOne.hpp"
#define ELEM_PARTER_INC           "elemental/matrices/Parter.hpp"
#define ELEM_PEI_INC              "elemental/matrices/Pei.hpp"
#define ELEM_REDHEFFER_INC        "elemental/matrices/Redheffer.hpp"
#define ELEM_RIEMANN_INC          "elemental/matrices/Riemann.hpp"
#define ELEM_RIS_INC              "elemental/matrices/Ris.hpp"
#define ELEM_TOEPLITZ_INC         "elemental/matrices/Toeplitz.hpp"
#define ELEM_TREFETHEN_INC        "elemental/matrices/Trefethen.hpp"
#define ELEM_TRIANGLE_INC         "elemental/matrices/Triangle.hpp"
#define ELEM_TRIW_INC             "elemental/matrices/TriW.hpp"
#define ELEM_WALSH_INC            "elemental/matrices/Walsh.hpp"
#define ELEM_WHALE_INC            "elemental/matrices/Whale.hpp"
#define ELEM_WILKINSON_INC        "elemental/matrices/Wilkinson.hpp"
#define ELEM_ZEROS_INC            "elemental/matrices/Zeros.hpp"

// Random matrices
// ---------------
#define ELEM_GAUSSIAN_INC \
  "elemental/matrices/Gaussian.hpp"
#define ELEM_HAAR_INC \
  "elemental/matrices/Haar.hpp"
#define ELEM_HATANONELSON_INC \
  "elemental/matrices/HatanoNelson.hpp"
#define ELEM_HERMITIANUNIFORMSPECTRUM_INC \
  "elemental/matrices/HermitianUniformSpectrum.hpp"
#define ELEM_NORMALUNIFORMSPECTRUM_INC \
  "elemental/matrices/NormalUniformSpectrum.hpp"
#define ELEM_UNIFORM_INC \
  "elemental/matrices/Uniform.hpp"
#define ELEM_UNIFORMHELMHOLTZGREENS_INC \
  "elemental/matrices/UniformHelmholtzGreens.hpp"
#define ELEM_WIGNER_INC \
  "elemental/matrices/Wigner.hpp"

// Input/output
// ============
#define ELEM_IO_INC      "elemental/io.hpp"
#define ELEM_IO_DECL_INC "elemental/io/decl.hpp"
#define ELEM_IO_IMPL_INC "elemental/io/impl.hpp"

#define ELEM_DISPLAY_INC "elemental/io/Display.hpp"
#define ELEM_PRINT_INC   "elemental/io/Print.hpp"
#define ELEM_READ_INC    "elemental/io/Read.hpp"
#define ELEM_SPY_INC     "elemental/io/Spy.hpp"
#define ELEM_WRITE_INC   "elemental/io/Write.hpp"

// Convex optimization
// ===================
#define ELEM_CONVEX_INC      "elemental/convex.hpp"
#define ELEM_CONVEX_DECL_INC "elemental/convex/decl.hpp"
#define ELEM_CONVEX_IMPL_INC "elemental/convex/impl.hpp"

// Utilities
// ---------
#define ELEM_CLIP_INC          "elemental/convex/Clip.hpp"
#define ELEM_COVARIANCE_INC    "elemental/convex/Covariance.hpp"
#define ELEM_LOGBARRIER_INC    "elemental/convex/LogBarrier.hpp"
#define ELEM_LOGDETDIV_INC     "elemental/convex/LogDetDiv.hpp"
#define ELEM_SOFTTHRESHOLD_INC "elemental/convex/SoftThreshold.hpp"
#define ELEM_SVT_INC           "elemental/convex/SVT.hpp"
// Specific SVT routines
// ^^^^^^^^^^^^^^^^^^^^^
#define ELEM_SVT_CROSS_INC     "elemental/convex/SVT/Cross.hpp"
#define ELEM_SVT_NORMAL_INC    "elemental/convex/SVT/Normal.hpp"
#define ELEM_SVT_PIVOTEDQR_INC "elemental/convex/SVT/PivotedQR.hpp"
#define ELEM_SVT_TSQR_INC      "elemental/convex/SVT/TSQR.hpp"

// ADMM
// ----
#define ELEM_BASISPURSUIT_INC     "elemental/convex/BasisPursuit.hpp"
#define ELEM_LINEARPROGRAM_INC    "elemental/convex/LinearProgram.hpp"
#define ELEM_QUADRATICPROGRAM_INC "elemental/convex/QuadraticProgram.hpp"
#define ELEM_SPARSEINVCOV_INC     "elemental/convex/SparseInvCov.hpp"


#endif // ifndef ELEM_INCLUDEPATHS_HPP
