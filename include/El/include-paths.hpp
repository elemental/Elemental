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

// BLAS-like routines
// ==================
#define EL_BLAS_INC      "El/blas-like.hpp"
#define EL_BLAS_DECL_INC "El/blas-like/decl.hpp"
#define EL_BLAS_IMPL_INC "El/blas-like/impl.hpp"

// Level 1 routines
// ----------------
#define EL_BLAS1_INC "El/blas-like/level1.hpp"

#define EL_ADJOINT_INC \
  "El/blas-like/level1/Adjoint.hpp"
#define EL_AXPY_INC \
  "El/blas-like/level1/Axpy.hpp"
#define EL_AXPYTRIANGLE_INC \
  "El/blas-like/level1/AxpyTriangle.hpp"
#define EL_CONJUGATE_INC \
  "El/blas-like/level1/Conjugate.hpp"
#define EL_COPY_INC \
  "El/blas-like/level1/Copy.hpp"
#define EL_DIAGONALSCALE_INC \
  "El/blas-like/level1/DiagonalScale.hpp"
#define EL_DIAGONALSCALETRAPEZOID_INC \
  "El/blas-like/level1/DiagonalScaleTrapezoid.hpp"
#define EL_DIAGONALSOLVE_INC \
  "El/blas-like/level1/DiagonalSolve.hpp"
#define EL_DOT_INC \
  "El/blas-like/level1/Dot.hpp"
#define EL_DOTU_INC \
  "El/blas-like/level1/Dotu.hpp"
#define EL_ENTRYWISEMAP_INC \
  "El/blas-like/level1/EntrywiseMap.hpp"
#define EL_HADAMARD_INC \
  "El/blas-like/level1/Hadamard.hpp"
#define EL_HILBERTSCHMIDT_INC \
  "El/blas-like/level1/HilbertSchmidt.hpp"
#define EL_MAKEHERMITIAN_INC \
  "El/blas-like/level1/MakeHermitian.hpp"
#define EL_MAKEREAL_INC \
  "El/blas-like/level1/MakeReal.hpp"
#define EL_MAKESYMMETRIC_INC \
  "El/blas-like/level1/MakeSymmetric.hpp"
#define EL_MAKETRAPEZOIDAL_INC \
  "El/blas-like/level1/MakeTrapezoidal.hpp"
#define EL_MAKETRIANGULAR_INC \
  "El/blas-like/level1/MakeTriangular.hpp"
#define EL_MAX_INC \
  "El/blas-like/level1/Max.hpp"
#define EL_MAXABS_INC \
  "El/blas-like/level1/MaxAbs.hpp"
#define EL_MIN_INC \
  "El/blas-like/level1/Min.hpp"
#define EL_MINABS_INC \
  "El/blas-like/level1/MinAbs.hpp"
#define EL_NRM2_INC \
  "El/blas-like/level1/Nrm2.hpp"
#define EL_QUASIDIAGONALSCALE_INC \
  "El/blas-like/level1/QuasiDiagonalScale.hpp"
#define EL_QUASIDIAGONALSOLVE_INC \
  "El/blas-like/level1/QuasiDiagonalSolve.hpp"
#define EL_SCALE_INC \
  "El/blas-like/level1/Scale.hpp"
#define EL_SCALETRAPEZOID_INC \
  "El/blas-like/level1/ScaleTrapezoid.hpp"
#define EL_SETDIAGONAL_INC \
  "El/blas-like/level1/SetDiagonal.hpp"
#define EL_SWAP_INC \
  "El/blas-like/level1/Swap.hpp"
#define EL_SYMMETRIC2X2INV_INC \
  "El/blas-like/level1/Symmetric2x2Inv.hpp"
#define EL_SYMMETRIC2X2SCALE_INC \
  "El/blas-like/level1/Symmetric2x2Scale.hpp"
#define EL_SYMMETRIC2X2SOLVE_INC \
  "El/blas-like/level1/Symmetric2x2Solve.hpp"
#define EL_TRANSPOSE_INC \
  "El/blas-like/level1/Transpose.hpp"
#define EL_UPDATEDIAGONAL_INC \
  "El/blas-like/level1/UpdateDiagonal.hpp"
#define EL_ZERO_INC \
  "El/blas-like/level1/Zero.hpp"

// Level 2 routines
// ----------------
#define EL_BLAS2_INC "El/blas-like/level2.hpp"

#define EL_GEMV_INC \
  "El/blas-like/level2/Gemv.hpp"
#define EL_GER_INC \
  "El/blas-like/level2/Ger.hpp"
#define EL_GERU_INC \
  "El/blas-like/level2/Geru.hpp"
#define EL_HEMV_INC \
  "El/blas-like/level2/Hemv.hpp"
#define EL_HER2_INC \
  "El/blas-like/level2/Her2.hpp"
#define EL_HER_INC \
  "El/blas-like/level2/Her.hpp"
#define EL_QUASITRSV_INC \
  "El/blas-like/level2/QuasiTrsv.hpp"
#define EL_SYMV_INC \
  "El/blas-like/level2/Symv.hpp"
#define EL_SYR2_INC \
  "El/blas-like/level2/Syr2.hpp"
#define EL_SYR_INC \
  "El/blas-like/level2/Syr.hpp"
#define EL_TRMV_INC \
  "El/blas-like/level2/Trmv.hpp"
#define EL_TRR2_INC \
  "El/blas-like/level2/Trr2.hpp"
#define EL_TRR_INC \
  "El/blas-like/level2/Trr.hpp"
#define EL_TRSV_INC \
  "El/blas-like/level2/Trsv.hpp"

// Level 3 routines
// ----------------
#define EL_BLAS3_INC "El/blas-like/level3.hpp"

#define EL_GEMM_INC         "El/blas-like/level3/Gemm.hpp"
#define EL_HEMM_INC         "El/blas-like/level3/Hemm.hpp"
#define EL_HER2K_INC        "El/blas-like/level3/Her2k.hpp"
#define EL_HERK_INC         "El/blas-like/level3/Herk.hpp"
#define EL_MULTISHIFTQUASITRSM_INC \
                            "El/blas-like/level3/MultiShiftQuasiTrsm.hpp"
#define EL_MULTISHIFTTRSM_INC \
                            "El/blas-like/level3/MultiShiftTrsm.hpp"
#define EL_QUASITRSM_INC    "El/blas-like/level3/QuasiTrsm.hpp"
#define EL_SYMM_INC         "El/blas-like/level3/Symm.hpp"
#define EL_SYR2K_INC        "El/blas-like/level3/Syr2k.hpp"
#define EL_SYRK_INC         "El/blas-like/level3/Syrk.hpp"
#define EL_TRDTRMM_INC      "El/blas-like/level3/Trdtrmm.hpp"
#define EL_TRMM_INC         "El/blas-like/level3/Trmm.hpp"
#define EL_TRSM_INC         "El/blas-like/level3/Trsm.hpp"
#define EL_TRSTRM_INC       "El/blas-like/level3/Trstrm.hpp"
#define EL_TRTRMM_INC       "El/blas-like/level3/Trtrmm.hpp"
#define EL_TWOSIDEDTRMM_INC "El/blas-like/level3/TwoSidedTrmm.hpp"
#define EL_TWOSIDEDTRSM_INC "El/blas-like/level3/TwoSidedTrsm.hpp"

// LAPACK-like routines
// ====================
#define EL_LAPACK_INC      "El/lapack-like.hpp"
#define EL_LAPACK_DECL_INC "El/lapack-like/decl.hpp"
#define EL_LAPACK_IMPL_INC "El/lapack-like/impl.hpp"

// Reduction to condensed form
// ---------------------------
#define EL_CONDENSE_INC      "El/lapack-like/condense.hpp"
#define EL_CONDENSE_DECL_INC "El/lapack-like/condense/decl.hpp"
#define EL_CONDENSE_IMPL_INC "El/lapack-like/condense/impl.hpp"

#define EL_BIDIAG_INC \
  "El/lapack-like/condense/Bidiag.hpp"
#define EL_HERMITIANTRIDIAG_INC \
  "El/lapack-like/condense/HermitianTridiag.hpp"
#define EL_HESSENBERG_INC \
  "El/lapack-like/condense/Hessenberg.hpp"

// Decompositions
// --------------
#define EL_DECOMP_INC      "El/lapack-like/decomp.hpp"
#define EL_DECOMP_DECL_INC "El/lapack-like/decomp/decl.hpp"
#define EL_DECOMP_IMPL_INC "El/lapack-like/decomp/impl.hpp"

#define EL_HERMITIANEIG_INC \
  "El/lapack-like/decomp/HermitianEig.hpp"
#define EL_HERMITIANGENDEFINITEEIG_INC \
  "El/lapack-like/decomp/HermitianGenDefiniteEig.hpp"
#define EL_HERMITIANTRIDIAGEIG_INC \
  "El/lapack-like/decomp/HermitianTridiagEig.hpp"
#define EL_POLAR_INC \
  "El/lapack-like/decomp/Polar.hpp"
#define EL_SCHUR_INC \
  "El/lapack-like/decomp/Schur.hpp"
#define EL_SKEWHERMITIANEIG_INC \
  "El/lapack-like/decomp/SkewHermitianEig.hpp"
#define EL_SVD_INC \
  "El/lapack-like/decomp/SVD.hpp"

// Specific HermitianEig routines
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#define EL_HERMITIANEIG_SDC_INC \
  "El/lapack-like/decomp/HermitianEig/SDC.hpp"

// Specific HermitianTridiagEig routines
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#define EL_HERMITIANTRIDIAGEIG_SORT_INC \
  "El/lapack-like/decomp/HermitianTridiagEig/Sort.hpp"

// Specific Polar routines
// ^^^^^^^^^^^^^^^^^^^^^^^
#define EL_POLAR_QDWH_INC \
  "El/lapack-like/decomp/Polar/QDWH.hpp"
#define EL_POLAR_SVD_INC \
  "El/lapack-like/decomp/Polar/SVD.hpp"

// Specific Schur routines
// ^^^^^^^^^^^^^^^^^^^^^^^
#define EL_SCHUR_INVERSEFREESDC_INC \
  "El/lapack-like/decomp/Schur/InverseFreeSDC.hpp"
#define EL_SCHUR_QR_INC \
  "El/lapack-like/decomp/Schur/QR.hpp"
#define EL_SCHUR_SDC_INC \
  "El/lapack-like/decomp/Schur/SDC.hpp"

// Factorizations
// --------------
#define EL_FACTOR_INC      "El/lapack-like/factor.hpp"
#define EL_FACTOR_DECL_INC "El/lapack-like/factor/decl.hpp"
#define EL_FACTOR_IMPL_INC "El/lapack-like/factor/impl.hpp"

#define EL_CHOLESKY_INC \
  "El/lapack-like/factor/Cholesky.hpp"
#define EL_LDL_INC \
  "El/lapack-like/factor/LDL.hpp"
#define EL_LU_INC \
  "El/lapack-like/factor/LU.hpp"

#define EL_LQ_INC \
  "El/lapack-like/factor/LQ.hpp"
#define EL_QR_INC \
  "El/lapack-like/factor/QR.hpp"
#define EL_RQ_INC \
  "El/lapack-like/factor/RQ.hpp"

#define EL_GQR_INC \
  "El/lapack-like/factor/GQR.hpp"
#define EL_GRQ_INC \
  "El/lapack-like/factor/GRQ.hpp"

#define EL_ID_INC \
  "El/lapack-like/factor/ID.hpp"
#define EL_SKELETON_INC \
  "El/lapack-like/factor/Skeleton.hpp"

// Specific LDL routines
// ^^^^^^^^^^^^^^^^^^^^^
#define EL_LDL_INERTIA_INC \
  "El/lapack-like/factor/LDL/Inertia.hpp"
#define EL_LDL_MULTIPLYAFTER_INC \
  "El/lapack-like/factor/LDL/MultiplyAfter.hpp"
#define EL_LDL_PIVOTED_INC \
  "El/lapack-like/factor/LDL/Pivoted.hpp"
#define EL_LDL_SOLVEAFTER_INC \
  "El/lapack-like/factor/LDL/SolveAfter.hpp"
#define EL_LDL_VAR3_INC \
  "El/lapack-like/factor/LDL/Var3.hpp"

// Specific LQ routines
// ^^^^^^^^^^^^^^^^^^^^
#define EL_LQ_APPLYQ_INC \
  "El/lapack-like/factor/LQ/ApplyQ.hpp"
#define EL_LQ_EXPLICIT_INC \
  "El/lapack-like/factor/LQ/Explicit.hpp"
#define EL_LQ_HOUSEHOLDER_INC \
  "El/lapack-like/factor/LQ/Householder.hpp"
#define EL_LQ_PANELHOUSEHOLDER_INC \
  "El/lapack-like/factor/LQ/PanelHouseholder.hpp"

// Specific LU routines
// ^^^^^^^^^^^^^^^^^^^^
#define EL_LU_FULL_INC \
  "El/lapack-like/factor/LU/Full.hpp"
#define EL_LU_LOCAL_INC \
  "El/lapack-like/factor/LU/Local.hpp"
#define EL_LU_PANEL_INC \
  "El/lapack-like/factor/LU/Panel.hpp"
#define EL_LU_SOLVEAFTER_INC \
  "El/lapack-like/factor/LU/SolveAfter.hpp"

// Specific QR routines
// ^^^^^^^^^^^^^^^^^^^^
#define EL_QR_APPLYQ_INC \
  "El/lapack-like/factor/QR/ApplyQ.hpp"
#define EL_QR_BUSINGERGOLUB_INC \
  "El/lapack-like/factor/QR/BusingerGolub.hpp"
#define EL_QR_CHOLESKY_INC \
  "El/lapack-like/factor/QR/Cholesky.hpp"
#define EL_QR_EXPLICIT_INC \
  "El/lapack-like/factor/QR/Explicit.hpp"
#define EL_QR_HOUSEHOLDER_INC \
  "El/lapack-like/factor/QR/Householder.hpp"
#define EL_QR_PANELHOUSEHOLDER_INC \
  "El/lapack-like/factor/QR/PanelHouseholder.hpp"
#define EL_QR_TS_INC \
  "El/lapack-like/factor/QR/TS.hpp"

// Specific RQ routines
// ^^^^^^^^^^^^^^^^^^^^
#define EL_RQ_APPLYQ_INC \
  "El/lapack-like/factor/RQ/ApplyQ.hpp"
#define EL_RQ_CHOLESKY_INC \
  "El/lapack-like/factor/RQ/Cholesky.hpp"
#define EL_RQ_HOUSEHOLDER_INC \
  "El/lapack-like/factor/RQ/Householder.hpp"
#define EL_RQ_PANELHOUSEHOLDER_INC \
  "El/lapack-like/factor/RQ/PanelHouseholder.hpp"

// Matrix functions
// ----------------
#define EL_FUNCS_INC      "El/lapack-like/funcs.hpp"
#define EL_FUNCS_DECL_INC "El/lapack-like/funcs/decl.hpp"
#define EL_FUNCS_IMPL_INC "El/lapack-like/funcs/impl.hpp"

#define EL_HERMITIANFUNCTION_INC \
  "El/lapack-like/funcs/HermitianFunction.hpp"
#define EL_INVERSE_INC \
  "El/lapack-like/funcs/Inverse.hpp"
#define EL_PSEUDOINVERSE_INC \
  "El/lapack-like/funcs/Pseudoinverse.hpp"
#define EL_SIGN_INC \
  "El/lapack-like/funcs/Sign.hpp"
#define EL_SQUAREROOT_INC \
  "El/lapack-like/funcs/SquareRoot.hpp"
// Specific inversion routines
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^
#define EL_TRIANGULARINVERSE_INC \
  "El/lapack-like/funcs/Inverse/Triangular.hpp"
#define EL_GENERALINVERSE_INC \
  "El/lapack-like/funcs/Inverse/General.hpp"
#define EL_HPDINVERSE_INC \
  "El/lapack-like/funcs/Inverse/HPD.hpp"
#define EL_SYMMETRICINVERSE_INC \
  "El/lapack-like/funcs/Inverse/Symmetric.hpp"
#define EL_HERMITIANINVERSE_INC \
  "El/lapack-like/funcs/Inverse/Hermitian.hpp"

// Permutations
// ------------
#define EL_APPLYCOLPIVOTS_INC \
  "El/lapack-like/perm/ApplyColPivots.hpp"
#define EL_APPLYROWPIVOTS_INC \
  "El/lapack-like/perm/ApplyRowPivots.hpp"
#define EL_APPLYSYMMETRICPIVOTS_INC \
  "El/lapack-like/perm/ApplySymmetricPivots.hpp"
#define EL_EXPLICITPERMUTATION_INC \
  "El/lapack-like/perm/ExplicitPermutation.hpp"
#define EL_INVERTPERMUTATION_INC \
  "El/lapack-like/perm/InvertPermutation.hpp"
#define EL_PERMUTATIONMETA_INC \
  "El/lapack-like/perm/PermutationMeta.hpp"
#define EL_PERMUTECOLS_INC \
  "El/lapack-like/perm/PermuteCols.hpp"
#define EL_PERMUTEROWS_INC \
  "El/lapack-like/perm/PermuteRows.hpp"
#define EL_PIVOTSTOPARTIALPERMUTATION_INC \
  "El/lapack-like/perm/PivotsToPartialPermutation.hpp"
#define EL_PIVOTSTOPERMUTATION_INC \
  "El/lapack-like/perm/PivotsToPermutation.hpp"

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

// Solvers
// -------
#define EL_SOLVE_INC      "El/lapack-like/solve.hpp"
#define EL_SOLVE_DECL_INC "El/lapack-like/solve/decl.hpp"
#define EL_SOLVE_IMPL_INC "El/lapack-like/solve/impl.hpp"

#define EL_GAUSSIANELIMINATION_INC \
  "El/lapack-like/solve/GaussianElimination.hpp"
#define EL_HPDSOLVE_INC \
  "El/lapack-like/solve/HPDSolve.hpp"

#define EL_HERMITIANSOLVE_INC \
  "El/lapack-like/solve/HermitianSolve.hpp"
#define EL_SYMMETRICSOLVE_INC \
  "El/lapack-like/solve/SymmetricSolve.hpp"

#define EL_LEASTSQUARES_INC \
  "El/lapack-like/solve/LeastSquares.hpp"

#define EL_GLM_INC \
  "El/lapack-like/solve/GLM.hpp"
#define EL_LSE_INC \
  "El/lapack-like/solve/LSE.hpp"

#define EL_MULTISHIFTHESSSOLVE_INC \
  "El/lapack-like/solve/MultiShiftHessSolve.hpp"

// Utilities
// ---------
#define EL_LAPACKUTIL_INC      "El/lapack-like/util.hpp"
#define EL_LAPACKUTIL_DECL_INC "El/lapack-like/util/decl.hpp"
#define EL_LAPACKUTIL_IMPL_INC "El/lapack-like/util/impl.hpp"

#define EL_APPLYPACKEDREFLECTORS_INC \
  "El/lapack-like/util/ApplyPackedReflectors.hpp"
#define EL_EXPANDPACKEDREFLECTORS_INC \
  "El/lapack-like/util/ExpandPackedReflectors.hpp"
#define EL_HYPERBOLICREFLECTOR_INC \
  "El/lapack-like/util/HyperbolicReflector.hpp"
#define EL_MEDIAN_INC \
  "El/lapack-like/util/Median.hpp"
#define EL_PERMUTATIONPARITY_INC \
  "El/lapack-like/util/PermutationParity.hpp"
#define EL_PIVOTPARITY_INC \
  "El/lapack-like/util/PivotParity.hpp"
#define EL_REFLECTOR_INC \
  "El/lapack-like/util/Reflector.hpp"
#define EL_SORT_INC \
  "El/lapack-like/util/Sort.hpp"

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

// Input/output
// ============
#define EL_IO_INC      "El/io.hpp"
#define EL_IO_DECL_INC "El/io/decl.hpp"
#define EL_IO_IMPL_INC "El/io/impl.hpp"

#define EL_DISPLAY_INC "El/io/Display.hpp"
#define EL_PRINT_INC   "El/io/Print.hpp"
#define EL_READ_INC    "El/io/Read.hpp"
#define EL_SPY_INC     "El/io/Spy.hpp"
#define EL_WRITE_INC   "El/io/Write.hpp"

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
