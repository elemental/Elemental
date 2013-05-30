/*
   Copyright (c) 2009-2013, Jack Poulson
                      2013, Michael C. Grant
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
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

%include "elemental/blas-like_decl.hpp"

// The equivalent of elemental/blas-like/level1.hpp
%include "elemental/blas-like/level1/Adjoint.hpp"
%include "elemental/blas-like/level1/Axpy.hpp"
%include "elemental/blas-like/level1/AxpyTriangle.hpp"
%include "elemental/blas-like/level1/Conjugate.hpp"
%include "elemental/blas-like/level1/Copy.hpp"
%include "elemental/blas-like/level1/DiagonalScale.hpp"
%include "elemental/blas-like/level1/DiagonalSolve.hpp"
%include "elemental/blas-like/level1/Dot.hpp"
%include "elemental/blas-like/level1/Dotu.hpp"
%include "elemental/blas-like/level1/MakeHermitian.hpp"
%include "elemental/blas-like/level1/MakeReal.hpp"
%include "elemental/blas-like/level1/MakeSymmetric.hpp"
%include "elemental/blas-like/level1/MakeTrapezoidal.hpp"
%include "elemental/blas-like/level1/MakeTriangular.hpp"
%include "elemental/blas-like/level1/Nrm2.hpp"
%include "elemental/blas-like/level1/Scale.hpp"
%include "elemental/blas-like/level1/ScaleTrapezoid.hpp"
%include "elemental/blas-like/level1/SetDiagonal.hpp"
%include "elemental/blas-like/level1/Transpose.hpp"
%include "elemental/blas-like/level1/Zero.hpp"

// The equivalent of elemental/blas-like/level2.hpp
%include "elemental/blas-like/level2/Gemv.hpp"
%include "elemental/blas-like/level2/Ger.hpp"
%include "elemental/blas-like/level2/Geru.hpp"
%include "elemental/blas-like/level2/Hemv.hpp"
%include "elemental/blas-like/level2/Her.hpp"
%include "elemental/blas-like/level2/Her2.hpp"
%include "elemental/blas-like/level2/Symv.hpp"
%include "elemental/blas-like/level2/Syr.hpp"
%include "elemental/blas-like/level2/Syr2.hpp"
%include "elemental/blas-like/level2/Trmv.hpp"
%include "elemental/blas-like/level2/Trsv.hpp"

// The equivalent of elemental/blas-like/level3.hpp
%include "elemental/blas-like/level3/Gemm.hpp"
%include "elemental/blas-like/level3/Hemm.hpp"
%include "elemental/blas-like/level3/Her2k.hpp"
%include "elemental/blas-like/level3/Herk.hpp"
%include "elemental/blas-like/level3/Symm.hpp"
%include "elemental/blas-like/level3/Syr2k.hpp"
%include "elemental/blas-like/level3/Syrk.hpp"
%include "elemental/blas-like/level3/Trdtrmm.hpp"
%include "elemental/blas-like/level3/Trmm.hpp"
%include "elemental/blas-like/level3/Trsm.hpp"
%include "elemental/blas-like/level3/Trtrmm.hpp"
%include "elemental/blas-like/level3/Trtrsm.hpp"
%include "elemental/blas-like/level3/TwoSidedTrmm.hpp"
%include "elemental/blas-like/level3/TwoSidedTrsm.hpp"

// The equivalent of elemental/lapack-like/lapack-like.hpp
%include "elemental/lapack-like_decl.hpp"
%include "elemental/lapack-like/ConditionNumber.hpp"
%include "elemental/lapack-like/Determinant.hpp"

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

%typemap(in) Complex<float> {
  Py_complex ans1 = PyComplex_AsCComplex($input);
  $1 = elem::Complex<float>
    ( static_cast<float>(ans1.real), static_cast<float>(ans1.imag) );
}
%typecheck(SWIG_TYPECHECK_COMPLEX) Complex<float> {
  $1 = PyComplex_Check( $input ) || PyFloat_Check( $input ) || 
       PyInt_Check( $input ) || PyLong_Check( $input ) ? 1 : 0;
}
%typemap(out) Complex<float> {
  $result = PyComplex_FromDoubles( $1.real, $1.imag );
}
%typemap(in) Complex<double> {
  Py_complex ans1 = PyComplex_AsCComplex($input);
  $1 = elem::Complex<double>( ans1.real, ans1.imag );
}
%typecheck(SWIG_TYPECHECK_COMPLEX) Complex<double> {
  $1 = PyComplex_Check( $input ) || PyFloat_Check( $input ) || 
       PyInt_Check( $input ) || PyLong_Check( $input ) ? 1 : 0;
}
%typemap(out) Complex<double> {
  $result = PyComplex_FromDoubles( $1.real, $1.imag );
}
%typemap(out) SafeProduct<float> {
  $result = PyDict_New();
  PyDict_SetItemString( $result, "rho", PyFloat_FromDouble( $1.rho ) );
  PyDict_SetItemString( $result, "kappa", PyFloat_FromDouble( $1.kappa ) );
  PyDict_SetItemString( $result, "n", PyInt_FromLong( $1.n ) );
}
%typemap(out) SafeProduct<double> {
  $result = PyDict_New();
  PyDict_SetItemString( $result, "rho", PyFloat_FromDouble( $1.rho ) );
  PyDict_SetItemString( $result, "kappa", PyFloat_FromDouble( $1.kappa ) );
  PyDict_SetItemString( $result, "n", PyInt_FromLong( $1.n ) );
}
%typemap(out) SafeProduct<Complex<float> > {
  $result = PyDict_New();
  PyDict_SetItemString( $result, "rho", PyComplex_FromDoubles( $1.rho.real, $1.rho.imag ) );
  PyDict_SetItemString( $result, "kappa", PyFloat_FromDouble( $1.kappa ) );
  PyDict_SetItemString( $result, "n", PyInt_FromLong( $1.n ) );
}
%typemap(out) SafeProduct<Complex<double> > {
  $result = PyDict_New();
  PyDict_SetItemString( $result, "rho", PyComplex_FromDoubles( $1.rho.real, $1.rho.imag ) );
  PyDict_SetItemString( $result, "kappa", PyFloat_FromDouble( $1.kappa ) );
  PyDict_SetItemString( $result, "n", PyInt_FromLong( $1.n ) );
}
%apply int { elem::Base<int> };
%apply float { elem::Base<float>, elem::Base<elem::Complex<float> > };
%apply double { elem::Base<double>, elem::Base<elem::Complex<double> > };
  
%ignore Matrix<int,int>::Matrix(int,int,const int*,int);
%ignore Matrix<int,int>::Matrix(int,int,int*,int);
%ignore Matrix<float,int>::Matrix(int,int,const float*,int);
%ignore Matrix<float,int>::Matrix(int,int,float*,int);
%ignore Matrix<double,int>::Matrix(int,int,const double*,int);
%ignore Matrix<double,int>::Matrix(int,int,double*,int);
%ignore Matrix<elem::Complex<float>,int>::Matrix(int,int,const elem::Complex<float>*,int);
%ignore Matrix<elem::Complex<float>,int>::Matrix(int,int,elem::Complex<float>*,int);
%ignore Matrix<elem::Complex<double>,int>::Matrix(int,int,const elem::Complex<double>*,int);
%ignore Matrix<elem::Complex<double>,int>::Matrix(int,int,elem::Complex<double>*,int);
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

%ignore DistMatrix<int,MC,MR,int>::DistMatrix(int,int,int,int,const int*,int,const elem::Grid&);
%ignore DistMatrix<int,MC,MR,int>::DistMatrix(int,int,int,int,int*,int,const elem::Grid&);
%ignore DistMatrix<float,MC,MR,int>::DistMatrix(int,int,int,int,const float*,int,const elem::Grid&);
%ignore DistMatrix<float,MC,MR,int>::DistMatrix(int,int,int,int,float*,int,const elem::Grid&);
%ignore DistMatrix<double,MC,MR,int>::DistMatrix(int,int,int,int,const double*,int,const elem::Grid&);
%ignore DistMatrix<double,MC,MR,int>::DistMatrix(int,int,int,int,double*,int,const elem::Grid&);
%ignore DistMatrix<elem::Complex<float>,MC,MR,int>::DistMatrix(int,int,int,int,const elem::Complex<float>*,int,const elem::Grid&);
%ignore DistMatrix<elem::Complex<float>,MC,MR,int>::DistMatrix(int,int,int,int,elem::Complex<float>*,int,const elem::Grid&);
%ignore DistMatrix<elem::Complex<double>,MC,MR,int>::DistMatrix(int,int,int,int,const elem::Complex<double>*,int,const elem::Grid&);
%ignore DistMatrix<elem::Complex<double>,MC,MR,int>::DistMatrix(int,int,int,int,elem::Complex<double>*,int,const elem::Grid&);
%template(DistMatrix_i) DistMatrix<int,MC,MR,int>;
%template(DistMatrix_s) DistMatrix<float,MC,MR,int>;
%template(DistMatrix_d) DistMatrix<double,MC,MR,int>;
%template(DistMatrix_c) DistMatrix<Complex<float>,MC,MR,int>;
%template(DistMatrix_z) DistMatrix<Complex<double>,MC,MR,int>;

#define OVERLOAD0(X) \
%template(X) X<float>; \
%template(X) X<double>; \
%template(X) X<Complex<float> >; \
%template(X) X<Complex<double> >;

#define OVERLOAD0_cpx(X) \
%template(X) X<float>; \
%template(X) X<double>;
  
#define OVERLOAD0_int(X) \
%template(X) X<int>; \
%template(X) X<float>; \
%template(X) X<double>; \
%template(X) X<Complex<float> >; \
%template(X) X<Complex<double> >;
 
#define OVERLOAD1(X) \
%template(X) X<float,MC,MR>; \
%template(X) X<double,MC,MR>; \
%template(X) X<Complex<float>,MC,MR>; \
%template(X) X<Complex<double>,MC,MR>;

#define OVERLOAD1_cpx(X) \
%template(X) X<float,MC,MR>; \
%template(X) X<double,MC,MR>;

#define OVERLOAD1_int(X) \
%template(X) X<int,MC,MR>; \
%template(X) X<float,MC,MR>; \
%template(X) X<double,MC,MR>; \
%template(X) X<Complex<float>,MC,MR>; \
%template(X) X<Complex<double>,MC,MR>;
 
#define OVERLOAD2(X) \
%template(X) X<float,MC,MR,MC,MR>; \
%template(X) X<double,MC,MR,MC,MR>; \
%template(X) X<Complex<float>,MC,MR,MC,MR>; \
%template(X) X<Complex<double>,MC,MR,MC,MR>;

#define OVERLOAD2_cpx(X) \
%template(X) X<float,MC,MR,MC,MR>; \
%template(X) X<double,MC,MR,MC,MR>;

#define OVERLOAD2_int(X) \
%template(X) X<int,MC,MR,MC,MR>; \
%template(X) X<float,MC,MR,MC,MR>; \
%template(X) X<double,MC,MR,MC,MR>; \
%template(X) X<Complex<float>,MC,MR,MC,MR>; \
%template(X) X<Complex<double>,MC,MR,MC,MR>;

#define OVERLOAD01(X) \
OVERLOAD0(X) \
OVERLOAD1(X)

#define OVERLOAD01_cpx(X) \
OVERLOAD0_cpx(X) \
OVERLOAD1_cpx(X)

#define OVERLOAD01_int(X) \
OVERLOAD0_int(X) \
OVERLOAD1_int(X)

#define OVERLOAD02(X) \
OVERLOAD0(X) \
OVERLOAD2(X)

#define OVERLOAD02_cpx(X) \
OVERLOAD0_cpx(X) \
OVERLOAD2_cpx(X)

#define OVERLOAD02_int(X) \
OVERLOAD0_int(X) \
OVERLOAD2_int(X)

#define OVERLOAD012(X) \
OVERLOAD0(X) \
OVERLOAD1(X) \
OVERLOAD2(X)

#define OVERLOAD012_cpx(X) \
OVERLOAD0_cpx(X) \
OVERLOAD1_cpx(X) \
OVERLOAD2_cpx(X)

#define OVERLOAD012_int(X) \
OVERLOAD0_int(X) \
OVERLOAD1_int(X) \
OVERLOAD2_int(X)
 
//
// BLAS-like routines (TODO)
//

OVERLOAD02_int(Adjoint)
OVERLOAD01_int(Axpy)
OVERLOAD01_int(AxpyTriangle)
OVERLOAD012_int(Conjugate)
OVERLOAD02_int(Copy)
OVERLOAD02_int(DiagonalScale)
OVERLOAD02(DiagonalSolve)
OVERLOAD01_int(Dot)
OVERLOAD01_int(Dotu)
OVERLOAD0_int(MakeHermitian)
OVERLOAD01_int(MakeReal)
OVERLOAD0_int(MakeSymmetric)
OVERLOAD01_int(MakeTrapezoidal)
OVERLOAD01_int(MakeTriangular)
OVERLOAD0(Nrm2)
OVERLOAD01_int(Scale)
OVERLOAD01_int(ScaleTrapezoid)
OVERLOAD01_int(SetDiagonal)
OVERLOAD02_int(Transpose)
OVERLOAD01_int(Zero)

// Most of these could support ints (except for Trsv)
OVERLOAD0(Gemv)
OVERLOAD0(Ger)
OVERLOAD0(Geru)
OVERLOAD0(Hemv)
OVERLOAD0(Her)
OVERLOAD0(Her2)
OVERLOAD0(Symv)
OVERLOAD0(Syr)
OVERLOAD0(Syr2)
OVERLOAD0(Trsv)

// Most of these could support ints if necessary
OVERLOAD0(Gemm)
OVERLOAD0(Hemm)
OVERLOAD0(Her2k)
OVERLOAD0(Herk)
OVERLOAD0(Symm)
OVERLOAD0(Syr2k)
OVERLOAD0(Syrk)
OVERLOAD0(Trdtrmm)
OVERLOAD0(Trmm)
OVERLOAD0(Trsm)
OVERLOAD0(Trtrmm)
OVERLOAD0(Trtrsm)
OVERLOAD0(TwoSidedTrmm)
OVERLOAD0(TwoSidedTrsm)

// Why are these commented out?
// OVERLOAD0(SetLocalSymvBlocksize)
// OVERLOAD0(LocalSymvBlocksize)
// OVERLOAD0(SetLocalTrrkBlocksize)
// OVERLOAD0(LocalTrrkBlocksize)
// OVERLOAD0(SetLocalTrr2kBlocksize)
// OVERLOAD0(LocalTrr2kBlocksize)

//
// LAPACK-like routines (TODO)
//

OVERLOAD01(ConditionNumber)
OVERLOAD0(Determinant)
OVERLOAD0(HPDDeterminant)
OVERLOAD0(SafeDeterminant)
OVERLOAD0(SafeHPDDeterminant)

//
// Special matrices
//

OVERLOAD01(Cauchy)
OVERLOAD01(CauchyLike)
OVERLOAD01_int(Circulant)
OVERLOAD01_int(Diagonal)
// Not sure how to handle Egorov yet...
OVERLOAD01(ExtendedKahan)
OVERLOAD01_int(Fiedler)
OVERLOAD01_int(Forsythe)
OVERLOAD01_cpx(Fourier)
OVERLOAD01_int(GCDMatrix)
OVERLOAD01_int(Gear)
OVERLOAD01(GKS)
OVERLOAD01_int(Grcar)
OVERLOAD01_int(Hankel)
OVERLOAD01_int(Hanowa)
OVERLOAD01(Helmholtz)
OVERLOAD01(HermitianUniformSpectrum)
OVERLOAD01(Hilbert)
OVERLOAD01_int(Identity)
OVERLOAD01_int(Jordan)
OVERLOAD01(Kahan)
OVERLOAD01_int(KMS)
OVERLOAD01(Laplacian)
OVERLOAD01_int(Lauchli)
OVERLOAD01(Legendre)
OVERLOAD01(Lehmer)
OVERLOAD01(Lotkin)
OVERLOAD01_int(MinIJ)
OVERLOAD01_cpx(NormalUniformSpectrum)
OVERLOAD01_int(Ones)
OVERLOAD01_int(OneTwoOne)
OVERLOAD01(Parter)
OVERLOAD01_int(Pei)
OVERLOAD01_int(Redheffer)
OVERLOAD01_int(Riemann)
OVERLOAD01(Ris)
OVERLOAD01_int(Toeplitz)
OVERLOAD01_int(TriW)
OVERLOAD01_int(Uniform)
OVERLOAD01_int(Walsh)
OVERLOAD01_int(Wilkinson)
OVERLOAD01_int(Zeros)

//
// Convex routines (TODO)
//

};
