%module elem
%{
#include "elemental.hpp"
%}

%include "std_except.i"
%include "std_iostream.i"
%include "std_string.i"

/*
 * This code converts a Python list of strings to the input format requested for
 * Elemental's Initialize() statement. Note that it does not try to return the
 * modified argc/argv back to Python.
 */

%include "typemaps.i"
%typemap(in) ( int& argc, char**& argv ) {
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  $1 = (int*)malloc(sizeof(int));
  $1[0] = PyList_Size($input);
  $2 = (char***) malloc(sizeof(char**));
  $2[0] = (char **) malloc((*$1+1)*sizeof(char *));
  for (i = 0; i < $1[0]; i++) {
    PyObject *s = PyList_GetItem($input,i);
    if (!PyString_Check(s)) {
        free($2[0]);
        free($1);
        PyErr_SetString(PyExc_ValueError, "List items must be strings");
        return NULL;
    }
    $2[0][i] = PyString_AsString(s);
  }
  $2[0][i] = 0;
}
%typemap(freearg) ( int& argc, char**& argv ) {
  if ( $2 && $2[0] ) free($2[0]);
  if ( $2 ) free($2);
  if ( $1 ) free($1);
}

/*
 * Blanket exception handling.
 */

%exception {
  try {
    $action
  } catch (std::exception) {
    PyErr_SetString(PyExc_RuntimeError,"Exception caught from Elemental");
  }
}

%rename(copy) *::operator=;
%rename(MPI_Initialize) elem::mpi::Initialize;
%rename(MPI_Initialized) elem::mpi::Initialized;
%rename(MPI_Finalize) elem::mpi::Finalize;
%rename(MPI_Finalized) elem::mpi::Finalized;

%include "elemental/config.h"

// The equivalent of elemental/core.hpp
%include "elemental/core.hpp"
%include "elemental/core/timer_decl.hpp"
%include "elemental/core/memory_decl.hpp"
%include "elemental/core/complex_decl.hpp"
%include "elemental/core/types_decl.hpp"
%include "elemental/core/matrix_forward_decl.hpp"
%include "elemental/core/dist_matrix_forward_decl.hpp"
%include "elemental/core/view_decl.hpp"
%include "elemental/core/matrix.hpp"
%include "elemental/core/imports/mpi.hpp"
%include "elemental/core/grid_decl.hpp"
%include "elemental/core/dist_matrix.hpp"
%include "elemental/core/dist_matrix/abstract.hpp"
%include "elemental/core/dist_matrix/mc_mr.hpp"
%include "elemental/core/dist_matrix/mc_star.hpp"
%include "elemental/core/dist_matrix/md_star.hpp"
%include "elemental/core/dist_matrix/mr_mc.hpp"
%include "elemental/core/dist_matrix/mr_star.hpp"
%include "elemental/core/dist_matrix/star_mc.hpp"
%include "elemental/core/dist_matrix/star_md.hpp"
%include "elemental/core/dist_matrix/star_mr.hpp"
%include "elemental/core/dist_matrix/star_star.hpp"
%include "elemental/core/dist_matrix/star_vc.hpp"
%include "elemental/core/dist_matrix/star_vr.hpp"
%include "elemental/core/dist_matrix/vc_star.hpp"
%include "elemental/core/dist_matrix/vr_star.hpp"
%include "elemental/core/environment_decl.hpp"
%include "elemental/core/indexing_decl.hpp"

// The equivalent of elemental/matrices.hpp
%include "elemental/matrices/Cauchy.hpp"
%include "elemental/matrices/CauchyLike.hpp"
%include "elemental/matrices/Circulant.hpp"
%include "elemental/matrices/Diagonal.hpp"
%include "elemental/matrices/Egorov.hpp"
%include "elemental/matrices/ExtendedKahan.hpp"
%include "elemental/matrices/Fiedler.hpp"
%include "elemental/matrices/Forsythe.hpp"
%include "elemental/matrices/Fourier.hpp"
%include "elemental/matrices/GCDMatrix.hpp"
%include "elemental/matrices/Gear.hpp"
%include "elemental/matrices/GKS.hpp"
%include "elemental/matrices/Grcar.hpp"
%include "elemental/matrices/Hankel.hpp"
%include "elemental/matrices/Hanowa.hpp"
%include "elemental/matrices/Helmholtz.hpp"
%include "elemental/matrices/HermitianUniformSpectrum.hpp"
%include "elemental/matrices/Hilbert.hpp"
%include "elemental/matrices/Identity.hpp"
%include "elemental/matrices/Jordan.hpp"
%include "elemental/matrices/Kahan.hpp"
%include "elemental/matrices/KMS.hpp"
%include "elemental/matrices/Laplacian.hpp"
%include "elemental/matrices/Lauchli.hpp"
%include "elemental/matrices/Legendre.hpp"
%include "elemental/matrices/Lehmer.hpp"
%include "elemental/matrices/Lotkin.hpp"
%include "elemental/matrices/MinIJ.hpp"
%include "elemental/matrices/NormalUniformSpectrum.hpp"
%include "elemental/matrices/Ones.hpp"
%include "elemental/matrices/OneTwoOne.hpp"
%include "elemental/matrices/Parter.hpp"
%include "elemental/matrices/Pei.hpp"
%include "elemental/matrices/Redheffer.hpp"
%include "elemental/matrices/Riemann.hpp"
%include "elemental/matrices/Ris.hpp"
%include "elemental/matrices/Toeplitz.hpp"
%include "elemental/matrices/TriW.hpp"
%include "elemental/matrices/Uniform.hpp"
%include "elemental/matrices/Walsh.hpp"
%include "elemental/matrices/Wilkinson.hpp"
%include "elemental/matrices/Zeros.hpp"

namespace elem {
  %template(Complex_s) Complex<float>;
  %template(Complex_d) Complex<double>;

  %template(Base_i) Base<int>;
  %template(Base_s) Base<float>;
  %template(Base_d) Base<double>;
  %template(Base_c) Base<Complex<float> >;
  %template(Base_z) Base<Complex<double> >;

  %template(Matrix_i) Matrix<int,int>;
  %template(Matrix_s) Matrix<float,int>;
  %template(Matrix_d) Matrix<double,int>;
  %template(Matrix_c) Matrix<Complex<float>,int>;
  %template(Matrix_z) Matrix<Complex<double>,int>;

  %template(AbstractDistMatrix_i) AbstractDistMatrix<int,int>;
  %template(AbstractDistMatrix_s) AbstractDistMatrix<float,int>;
  %template(AbstractDistMatrix_d) AbstractDistMatrix<double,int>;
  %template(AbstractDistMatrix_c) AbstractDistMatrix<Complex<float>,int>;
  %template(AbstractDistMatrix_z) AbstractDistMatrix<Complex<double>,int>;

  %template(DistMatrix_i) DistMatrix<int,MC,MR,int>;
  %template(DistMatrix_s) DistMatrix<float,MC,MR,int>;
  %template(DistMatrix_d) DistMatrix<double,MC,MR,int>;
  %template(DistMatrix_c) DistMatrix<Complex<float>,MC,MR,int>;
  %template(DistMatrix_z) DistMatrix<Complex<double>,MC,MR,int>;

  //
  // BLAS-like routines (TODO)
  //

  //
  // LAPACK-like routines (TODO)
  //

  //
  // Special matrices
  //

  %template(Cauchy_s) Cauchy<float>;
  %template(Cauchy_s) Cauchy<float,MC,MR>;
  %template(Cauchy_d) Cauchy<double>;
  %template(Cauchy_d) Cauchy<double,MC,MR>;
  %template(Cauchy_c) Cauchy<Complex<float> >;
  %template(Cauchy_c) Cauchy<Complex<float>,MC,MR>;
  %template(Cauchy_z) Cauchy<Complex<double> >;
  %template(Cauchy_z) Cauchy<Complex<double>,MC,MR>;

  %template(CauchyLike_s) CauchyLike<float>;
  %template(CauchyLike_s) CauchyLike<float,MC,MR>;
  %template(CauchyLike_d) CauchyLike<double>;
  %template(CauchyLike_d) CauchyLike<double,MC,MR>;
  %template(CauchyLike_c) CauchyLike<Complex<float> >;
  %template(CauchyLike_c) CauchyLike<Complex<float>,MC,MR>;
  %template(CauchyLike_z) CauchyLike<Complex<double> >;
  %template(CauchyLike_z) CauchyLike<Complex<double>,MC,MR>;

  %template(Circulant_i) Circulant<int>;
  %template(Circulant_i) Circulant<int,MC,MR>;
  %template(Circulant_s) Circulant<float>;
  %template(Circulant_s) Circulant<float,MC,MR>;
  %template(Circulant_d) Circulant<double>;
  %template(Circulant_d) Circulant<double,MC,MR>;
  %template(Circulant_c) Circulant<Complex<float> >;
  %template(Circulant_c) Circulant<Complex<float>,MC,MR>;
  %template(Circulant_z) Circulant<Complex<double> >;
  %template(Circulant_z) Circulant<Complex<double>,MC,MR>;

  %template(Diagonal_i) Diagonal<int>;
  %template(Diagonal_i) Diagonal<int,MC,MR>;
  %template(Diagonal_s) Diagonal<float>;
  %template(Diagonal_s) Diagonal<float,MC,MR>;
  %template(Diagonal_d) Diagonal<double>;
  %template(Diagonal_d) Diagonal<double,MC,MR>;
  %template(Diagonal_c) Diagonal<Complex<float> >;
  %template(Diagonal_c) Diagonal<Complex<float>,MC,MR>;
  %template(Diagonal_z) Diagonal<Complex<double> >;
  %template(Diagonal_z) Diagonal<Complex<double>,MC,MR>;

  // Not sure how to do this yet...
  //%template(Egorov_c) Egorov<float>;
  //%template(Egorov_c) Egorov<float,MC,MR>;
  //%template(Egorov_z) Egorov<double>;
  //%template(Egorov_z) Egorov<double,MC,MR>;

  %template(ExtendedKahan_s) ExtendedKahan<float>;
  %template(ExtendedKahan_s) ExtendedKahan<float,MC,MR>;
  %template(ExtendedKahan_d) ExtendedKahan<double>;
  %template(ExtendedKahan_d) ExtendedKahan<double,MC,MR>;
  %template(ExtendedKahan_c) ExtendedKahan<Complex<float> >;
  %template(ExtendedKahan_c) ExtendedKahan<Complex<float>,MC,MR>;
  %template(ExtendedKahan_z) ExtendedKahan<Complex<double> >;
  %template(ExtendedKahan_z) ExtendedKahan<Complex<double>,MC,MR>;

  %template(Fiedler_i) Fiedler<int>;
  %template(Fiedler_i) Fiedler<int,MC,MR>;
  %template(Fiedler_s) Fiedler<float>;
  %template(Fiedler_s) Fiedler<float,MC,MR>;
  %template(Fiedler_d) Fiedler<double>;
  %template(Fiedler_d) Fiedler<double,MC,MR>;
  %template(Fiedler_c) Fiedler<Complex<float> >;
  %template(Fiedler_c) Fiedler<Complex<float>,MC,MR>;
  %template(Fiedler_z) Fiedler<Complex<double> >;
  %template(Fiedler_z) Fiedler<Complex<double>,MC,MR>;

  %template(Forsythe_i) Forsythe<int>;
  %template(Forsythe_i) Forsythe<int,MC,MR>;
  %template(Forsythe_s) Forsythe<float>;
  %template(Forsythe_s) Forsythe<float,MC,MR>;
  %template(Forsythe_d) Forsythe<double>;
  %template(Forsythe_d) Forsythe<double,MC,MR>;
  %template(Forsythe_c) Forsythe<Complex<float> >;
  %template(Forsythe_c) Forsythe<Complex<float>,MC,MR>;
  %template(Forsythe_z) Forsythe<Complex<double> >;
  %template(Forsythe_z) Forsythe<Complex<double>,MC,MR>;

  %template(Fourier_c) Fourier<float>;
  %template(Fourier_c) Fourier<float,MC,MR>;
  %template(Fourier_z) Fourier<double>;
  %template(Fourier_z) Fourier<double,MC,MR>;

  %template(GCDMatrix_i) GCDMatrix<int>;
  %template(GCDMatrix_i) GCDMatrix<int,MC,MR>;
  %template(GCDMatrix_s) GCDMatrix<float>;
  %template(GCDMatrix_s) GCDMatrix<float,MC,MR>;
  %template(GCDMatrix_d) GCDMatrix<double>;
  %template(GCDMatrix_d) GCDMatrix<double,MC,MR>;
  %template(GCDMatrix_c) GCDMatrix<Complex<float> >;
  %template(GCDMatrix_c) GCDMatrix<Complex<float>,MC,MR>;
  %template(GCDMatrix_z) GCDMatrix<Complex<double> >;
  %template(GCDMatrix_z) GCDMatrix<Complex<double>,MC,MR>;

  %template(Gear_i) Gear<int>;
  %template(Gear_i) Gear<int,MC,MR>;
  %template(Gear_s) Gear<float>;
  %template(Gear_s) Gear<float,MC,MR>;
  %template(Gear_d) Gear<double>;
  %template(Gear_d) Gear<double,MC,MR>;
  %template(Gear_c) Gear<Complex<float> >;
  %template(Gear_c) Gear<Complex<float>,MC,MR>;
  %template(Gear_z) Gear<Complex<double> >;
  %template(Gear_z) Gear<Complex<double>,MC,MR>;

  %template(GKS_s) GKS<float>;
  %template(GKS_s) GKS<float,MC,MR>;
  %template(GKS_d) GKS<double>;
  %template(GKS_d) GKS<double,MC,MR>;
  %template(GKS_c) GKS<Complex<float> >;
  %template(GKS_c) GKS<Complex<float>,MC,MR>;
  %template(GKS_z) GKS<Complex<double> >;
  %template(GKS_z) GKS<Complex<double>,MC,MR>;

  %template(Grcar_i) Grcar<int>;
  %template(Grcar_i) Grcar<int,MC,MR>;
  %template(Grcar_s) Grcar<float>;
  %template(Grcar_s) Grcar<float,MC,MR>;
  %template(Grcar_d) Grcar<double>;
  %template(Grcar_d) Grcar<double,MC,MR>;
  %template(Grcar_c) Grcar<Complex<float> >;
  %template(Grcar_c) Grcar<Complex<float>,MC,MR>;
  %template(Grcar_z) Grcar<Complex<double> >;
  %template(Grcar_z) Grcar<Complex<double>,MC,MR>;

  %template(Hankel_i) Hankel<int>;
  %template(Hankel_i) Hankel<int,MC,MR>;
  %template(Hankel_s) Hankel<float>;
  %template(Hankel_s) Hankel<float,MC,MR>;
  %template(Hankel_d) Hankel<double>;
  %template(Hankel_d) Hankel<double,MC,MR>;
  %template(Hankel_c) Hankel<Complex<float> >;
  %template(Hankel_c) Hankel<Complex<float>,MC,MR>;
  %template(Hankel_z) Hankel<Complex<double> >;
  %template(Hankel_z) Hankel<Complex<double>,MC,MR>;

  %template(Hanowa_i) Hanowa<int>;
  %template(Hanowa_i) Hanowa<int,MC,MR>;
  %template(Hanowa_s) Hanowa<float>;
  %template(Hanowa_s) Hanowa<float,MC,MR>;
  %template(Hanowa_d) Hanowa<double>;
  %template(Hanowa_d) Hanowa<double,MC,MR>;
  %template(Hanowa_c) Hanowa<Complex<float> >;
  %template(Hanowa_c) Hanowa<Complex<float>,MC,MR>;
  %template(Hanowa_z) Hanowa<Complex<double> >;
  %template(Hanowa_z) Hanowa<Complex<double>,MC,MR>;

  %template(Helmholtz_s) Helmholtz<float>;
  %template(Helmholtz_s) Helmholtz<float,MC,MR>;
  %template(Helmholtz_d) Helmholtz<double>;
  %template(Helmholtz_d) Helmholtz<double,MC,MR>;
  %template(Helmholtz_c) Helmholtz<Complex<float> >;
  %template(Helmholtz_c) Helmholtz<Complex<float>,MC,MR>;
  %template(Helmholtz_z) Helmholtz<Complex<double> >;
  %template(Helmholtz_z) Helmholtz<Complex<double>,MC,MR>;

  %template(HermitianUniformSpectrum_s) HermitianUniformSpectrum<float>;
  %template(HermitianUniformSpectrum_s) HermitianUniformSpectrum<float,MC,MR>;
  %template(HermitianUniformSpectrum_d) HermitianUniformSpectrum<double>;
  %template(HermitianUniformSpectrum_d) HermitianUniformSpectrum<double,MC,MR>;
  %template(HermitianUniformSpectrum_c) HermitianUniformSpectrum<Complex<float> >;
  %template(HermitianUniformSpectrum_c) HermitianUniformSpectrum<Complex<float>,MC,MR>;
  %template(HermitianUniformSpectrum_z) HermitianUniformSpectrum<Complex<double> >;
  %template(HermitianUniformSpectrum_z) HermitianUniformSpectrum<Complex<double>,MC,MR>;

  %template(Hilbert_s) Hilbert<float>;
  %template(Hilbert_s) Hilbert<float,MC,MR>;
  %template(Hilbert_d) Hilbert<double>;
  %template(Hilbert_d) Hilbert<double,MC,MR>;
  %template(Hilbert_c) Hilbert<Complex<float> >;
  %template(Hilbert_c) Hilbert<Complex<float>,MC,MR>;
  %template(Hilbert_z) Hilbert<Complex<double> >;
  %template(Hilbert_z) Hilbert<Complex<double>,MC,MR>;

  %template(Identity_i) Identity<int>;
  %template(Identity_i) Identity<int,MC,MR>;
  %template(Identity_s) Identity<float>;
  %template(Identity_s) Identity<float,MC,MR>;
  %template(Identity_d) Identity<double>;
  %template(Identity_d) Identity<double,MC,MR>;
  %template(Identity_c) Identity<Complex<float> >;
  %template(Identity_c) Identity<Complex<float>,MC,MR>;
  %template(Identity_z) Identity<Complex<double> >;
  %template(Identity_z) Identity<Complex<double>,MC,MR>;

  %template(Jordan_i) Jordan<int>;
  %template(Jordan_i) Jordan<int,MC,MR>;
  %template(Jordan_s) Jordan<float>;
  %template(Jordan_s) Jordan<float,MC,MR>;
  %template(Jordan_d) Jordan<double>;
  %template(Jordan_d) Jordan<double,MC,MR>;
  %template(Jordan_c) Jordan<Complex<float> >;
  %template(Jordan_c) Jordan<Complex<float>,MC,MR>;
  %template(Jordan_z) Jordan<Complex<double> >;
  %template(Jordan_z) Jordan<Complex<double>,MC,MR>;

  %template(Kahan_s) Kahan<float>;
  %template(Kahan_s) Kahan<float,MC,MR>;
  %template(Kahan_d) Kahan<double>;
  %template(Kahan_d) Kahan<double,MC,MR>;
  %template(Kahan_c) Kahan<Complex<float> >;
  %template(Kahan_c) Kahan<Complex<float>,MC,MR>;
  %template(Kahan_z) Kahan<Complex<double> >;
  %template(Kahan_z) Kahan<Complex<double>,MC,MR>;

  %template(KMS_i) KMS<int>;
  %template(KMS_i) KMS<int,MC,MR>;
  %template(KMS_s) KMS<float>;
  %template(KMS_s) KMS<float,MC,MR>;
  %template(KMS_d) KMS<double>;
  %template(KMS_d) KMS<double,MC,MR>;
  %template(KMS_c) KMS<Complex<float> >;
  %template(KMS_c) KMS<Complex<float>,MC,MR>;
  %template(KMS_z) KMS<Complex<double> >;
  %template(KMS_z) KMS<Complex<double>,MC,MR>;

  %template(Laplacian_s) Laplacian<float>;
  %template(Laplacian_s) Laplacian<float,MC,MR>;
  %template(Laplacian_d) Laplacian<double>;
  %template(Laplacian_d) Laplacian<double,MC,MR>;
  %template(Laplacian_c) Laplacian<Complex<float> >;
  %template(Laplacian_c) Laplacian<Complex<float>,MC,MR>;
  %template(Laplacian_z) Laplacian<Complex<double> >;
  %template(Laplacian_z) Laplacian<Complex<double>,MC,MR>;

  %template(Lauchli_i) Lauchli<int>;
  %template(Lauchli_i) Lauchli<int,MC,MR>;
  %template(Lauchli_s) Lauchli<float>;
  %template(Lauchli_s) Lauchli<float,MC,MR>;
  %template(Lauchli_d) Lauchli<double>;
  %template(Lauchli_d) Lauchli<double,MC,MR>;
  %template(Lauchli_c) Lauchli<Complex<float> >;
  %template(Lauchli_c) Lauchli<Complex<float>,MC,MR>;
  %template(Lauchli_z) Lauchli<Complex<double> >;
  %template(Lauchli_z) Lauchli<Complex<double>,MC,MR>;

  %template(Legendre_s) Legendre<float>;
  %template(Legendre_s) Legendre<float,MC,MR>;
  %template(Legendre_d) Legendre<double>;
  %template(Legendre_d) Legendre<double,MC,MR>;
  %template(Legendre_c) Legendre<Complex<float> >;
  %template(Legendre_c) Legendre<Complex<float>,MC,MR>;
  %template(Legendre_z) Legendre<Complex<double> >;
  %template(Legendre_z) Legendre<Complex<double>,MC,MR>;

  %template(Lehmer_s) Lehmer<float>;
  %template(Lehmer_s) Lehmer<float,MC,MR>;
  %template(Lehmer_d) Lehmer<double>;
  %template(Lehmer_d) Lehmer<double,MC,MR>;
  %template(Lehmer_c) Lehmer<Complex<float> >;
  %template(Lehmer_c) Lehmer<Complex<float>,MC,MR>;
  %template(Lehmer_z) Lehmer<Complex<double> >;
  %template(Lehmer_z) Lehmer<Complex<double>,MC,MR>;

  %template(Lotkin_s) Lotkin<float>;
  %template(Lotkin_s) Lotkin<float,MC,MR>;
  %template(Lotkin_d) Lotkin<double>;
  %template(Lotkin_d) Lotkin<double,MC,MR>;
  %template(Lotkin_c) Lotkin<Complex<float> >;
  %template(Lotkin_c) Lotkin<Complex<float>,MC,MR>;
  %template(Lotkin_z) Lotkin<Complex<double> >;
  %template(Lotkin_z) Lotkin<Complex<double>,MC,MR>;

  %template(MinIJ_i) MinIJ<int>;
  %template(MinIJ_i) MinIJ<int,MC,MR>;
  %template(MinIJ_s) MinIJ<float>;
  %template(MinIJ_s) MinIJ<float,MC,MR>;
  %template(MinIJ_d) MinIJ<double>;
  %template(MinIJ_d) MinIJ<double,MC,MR>;
  %template(MinIJ_c) MinIJ<Complex<float> >;
  %template(MinIJ_c) MinIJ<Complex<float>,MC,MR>;
  %template(MinIJ_z) MinIJ<Complex<double> >;
  %template(MinIJ_z) MinIJ<Complex<double>,MC,MR>;

  %template(NormalUniformSpectrum_c) NormalUniformSpectrum<float>;
  %template(NormalUniformSpectrum_c) NormalUniformSpectrum<float,MC,MR>;
  %template(NormalUniformSpectrum_z) NormalUniformSpectrum<double>;
  %template(NormalUniformSpectrum_z) NormalUniformSpectrum<double,MC,MR>;

  %template(Ones_i) Ones<int>;
  %template(Ones_i) Ones<int,MC,MR>;
  %template(Ones_s) Ones<float>;
  %template(Ones_s) Ones<float,MC,MR>;
  %template(Ones_d) Ones<double>;
  %template(Ones_d) Ones<double,MC,MR>;
  %template(Ones_c) Ones<Complex<float> >;
  %template(Ones_c) Ones<Complex<float>,MC,MR>;
  %template(Ones_z) Ones<Complex<double> >;
  %template(Ones_z) Ones<Complex<double>,MC,MR>;

  %template(OneTwoOne_i) OneTwoOne<int>;
  %template(OneTwoOne_i) OneTwoOne<int,MC,MR>;
  %template(OneTwoOne_s) OneTwoOne<float>;
  %template(OneTwoOne_s) OneTwoOne<float,MC,MR>;
  %template(OneTwoOne_d) OneTwoOne<double>;
  %template(OneTwoOne_d) OneTwoOne<double,MC,MR>;
  %template(OneTwoOne_c) OneTwoOne<Complex<float> >;
  %template(OneTwoOne_c) OneTwoOne<Complex<float>,MC,MR>;
  %template(OneTwoOne_z) OneTwoOne<Complex<double> >;
  %template(OneTwoOne_z) OneTwoOne<Complex<double>,MC,MR>;

  %template(Parter_s) Parter<float>;
  %template(Parter_s) Parter<float,MC,MR>;
  %template(Parter_d) Parter<double>;
  %template(Parter_d) Parter<double,MC,MR>;
  %template(Parter_c) Parter<Complex<float> >;
  %template(Parter_c) Parter<Complex<float>,MC,MR>;
  %template(Parter_z) Parter<Complex<double> >;
  %template(Parter_z) Parter<Complex<double>,MC,MR>;

  %template(Pei_i) Pei<int>;
  %template(Pei_i) Pei<int,MC,MR>;
  %template(Pei_s) Pei<float>;
  %template(Pei_s) Pei<float,MC,MR>;
  %template(Pei_d) Pei<double>;
  %template(Pei_d) Pei<double,MC,MR>;
  %template(Pei_c) Pei<Complex<float> >;
  %template(Pei_c) Pei<Complex<float>,MC,MR>;
  %template(Pei_z) Pei<Complex<double> >;
  %template(Pei_z) Pei<Complex<double>,MC,MR>;

  %template(Redheffer_i) Redheffer<int>;
  %template(Redheffer_i) Redheffer<int,MC,MR>;
  %template(Redheffer_s) Redheffer<float>;
  %template(Redheffer_s) Redheffer<float,MC,MR>;
  %template(Redheffer_d) Redheffer<double>;
  %template(Redheffer_d) Redheffer<double,MC,MR>;
  %template(Redheffer_c) Redheffer<Complex<float> >;
  %template(Redheffer_c) Redheffer<Complex<float>,MC,MR>;
  %template(Redheffer_z) Redheffer<Complex<double> >;
  %template(Redheffer_z) Redheffer<Complex<double>,MC,MR>;

  %template(Riemann_i) Riemann<int>;
  %template(Riemann_i) Riemann<int,MC,MR>;
  %template(Riemann_s) Riemann<float>;
  %template(Riemann_s) Riemann<float,MC,MR>;
  %template(Riemann_d) Riemann<double>;
  %template(Riemann_d) Riemann<double,MC,MR>;
  %template(Riemann_c) Riemann<Complex<float> >;
  %template(Riemann_c) Riemann<Complex<float>,MC,MR>;
  %template(Riemann_z) Riemann<Complex<double> >;
  %template(Riemann_z) Riemann<Complex<double>,MC,MR>;

  %template(Ris_s) Ris<float>;
  %template(Ris_s) Ris<float,MC,MR>;
  %template(Ris_d) Ris<double>;
  %template(Ris_d) Ris<double,MC,MR>;
  %template(Ris_c) Ris<Complex<float> >;
  %template(Ris_c) Ris<Complex<float>,MC,MR>;
  %template(Ris_z) Ris<Complex<double> >;
  %template(Ris_z) Ris<Complex<double>,MC,MR>;

  %template(Toeplitz_i) Toeplitz<int>;
  %template(Toeplitz_i) Toeplitz<int,MC,MR>;
  %template(Toeplitz_s) Toeplitz<float>;
  %template(Toeplitz_s) Toeplitz<float,MC,MR>;
  %template(Toeplitz_d) Toeplitz<double>;
  %template(Toeplitz_d) Toeplitz<double,MC,MR>;
  %template(Toeplitz_c) Toeplitz<Complex<float> >;
  %template(Toeplitz_c) Toeplitz<Complex<float>,MC,MR>;
  %template(Toeplitz_z) Toeplitz<Complex<double> >;
  %template(Toeplitz_z) Toeplitz<Complex<double>,MC,MR>;

  %template(TriW_i) TriW<int>;
  %template(TriW_i) TriW<int,MC,MR>;
  %template(TriW_s) TriW<float>;
  %template(TriW_s) TriW<float,MC,MR>;
  %template(TriW_d) TriW<double>;
  %template(TriW_d) TriW<double,MC,MR>;
  %template(TriW_c) TriW<Complex<float> >;
  %template(TriW_c) TriW<Complex<float>,MC,MR>;
  %template(TriW_z) TriW<Complex<double> >;
  %template(TriW_z) TriW<Complex<double>,MC,MR>;

  %template(Uniform_i) Uniform<int>;
  %template(Uniform_i) Uniform<int,MC,MR>;
  %template(Uniform_s) Uniform<float>;
  %template(Uniform_s) Uniform<float,MC,MR>;
  %template(Uniform_d) Uniform<double>;
  %template(Uniform_d) Uniform<double,MC,MR>;
  %template(Uniform_c) Uniform<Complex<float> >;
  %template(Uniform_c) Uniform<Complex<float>,MC,MR>;
  %template(Uniform_z) Uniform<Complex<double> >;
  %template(Uniform_z) Uniform<Complex<double>,MC,MR>;

  %template(Walsh_i) Walsh<int>;
  %template(Walsh_i) Walsh<int,MC,MR>;
  %template(Walsh_s) Walsh<float>;
  %template(Walsh_s) Walsh<float,MC,MR>;
  %template(Walsh_d) Walsh<double>;
  %template(Walsh_d) Walsh<double,MC,MR>;
  %template(Walsh_c) Walsh<Complex<float> >;
  %template(Walsh_c) Walsh<Complex<float>,MC,MR>;
  %template(Walsh_z) Walsh<Complex<double> >;
  %template(Walsh_z) Walsh<Complex<double>,MC,MR>;

  %template(Wilkinson_i) Wilkinson<int>;
  %template(Wilkinson_i) Wilkinson<int,MC,MR>;
  %template(Wilkinson_s) Wilkinson<float>;
  %template(Wilkinson_s) Wilkinson<float,MC,MR>;
  %template(Wilkinson_d) Wilkinson<double>;
  %template(Wilkinson_d) Wilkinson<double,MC,MR>;
  %template(Wilkinson_c) Wilkinson<Complex<float> >;
  %template(Wilkinson_c) Wilkinson<Complex<float>,MC,MR>;
  %template(Wilkinson_z) Wilkinson<Complex<double> >;
  %template(Wilkinson_z) Wilkinson<Complex<double>,MC,MR>;

  %template(Zeros_i) Zeros<int>;
  %template(Zeros_i) Zeros<int,MC,MR>;
  %template(Zeros_s) Zeros<float>;
  %template(Zeros_s) Zeros<float,MC,MR>;
  %template(Zeros_d) Zeros<double>;
  %template(Zeros_d) Zeros<double,MC,MR>;
  %template(Zeros_c) Zeros<Complex<float> >;
  %template(Zeros_c) Zeros<Complex<float>,MC,MR>;
  %template(Zeros_z) Zeros<Complex<double> >;
  %template(Zeros_z) Zeros<Complex<double>,MC,MR>;

  //
  // Convex routines (TODO)
  //
};
