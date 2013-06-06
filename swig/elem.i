/*
   Copyright (c) 2009-2013, Jack Poulson
                      2013, Michael C. Grant
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
%module elem

#define RELEASE

%{
#include "elemental.hpp"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
static PyObject *numpy, *npmatrix;
static PyObject* create_npmatrix( int H, int W, int LDim, size_t DSize, int DType, bool writable, const void* data )
{
    PyObject *arr, *result;
    npy_intp dims[2], strides[2];
    dims[0] = H; dims[1] = W;
    strides[0] = DSize; strides[1] = DSize * LDim;
    arr = PyArray_New( &PyArray_Type, 2, &dims[0], DType, &strides[0], 
    const_cast<void*>(data), DSize, writable ? NPY_ARRAY_WRITEABLE : 0, NULL );
    result = PyObject_CallFunctionObjArgs( npmatrix, arr, NULL );
    PyArray_free( arr );
    return result;
}
%}

%import "std_except.i"
%import "std_iostream.i"
%import "std_string.i"

%init %{
import_array();
numpy = PyImport_ImportModule("numpy");
npmatrix = PyObject_GetAttrString(numpy,"matrix");
%}

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
%typemap(in) elem::Complex<float> {
  Py_complex ans1 = PyComplex_AsCComplex($input);
  $1 = elem::Complex<float>
    ( static_cast<float>(ans1.real), static_cast<float>(ans1.imag) );
}
%typecheck(SWIG_TYPECHECK_COMPLEX) elem::Complex<float> {
  $1 = PyComplex_Check( $input ) || PyFloat_Check( $input ) || 
       PyInt_Check( $input ) || PyLong_Check( $input ) ? 1 : 0;
}
%typemap(out) elem::Complex<float> {
  $result = PyComplex_FromDoubles( $1.real, $1.imag );
}
%typemap(in) elem::Complex<double> {
  Py_complex ans1 = PyComplex_AsCComplex($input);
  $1 = elem::Complex<double>( ans1.real, ans1.imag );
}
%typecheck(SWIG_TYPECHECK_COMPLEX) elem::Complex<double> {
  $1 = PyComplex_Check( $input ) || PyFloat_Check( $input ) || 
       PyInt_Check( $input ) || PyLong_Check( $input ) ? 1 : 0;
}
%typemap(out) elem::Complex<double> {
  $result = PyComplex_FromDoubles( $1.real, $1.imag );
}
%typemap(out) elem::SafeProduct<float> {
  $result = PyDict_New();
  PyDict_SetItemString( $result, "rho",   PyFloat_FromDouble( $1.rho ) );
  PyDict_SetItemString( $result, "kappa", PyFloat_FromDouble( $1.kappa ) );
  PyDict_SetItemString( $result, "n",     PyInt_FromLong( $1.n ) );
}
%typemap(out) elem::SafeProduct<double> {
  $result = PyDict_New();
  PyDict_SetItemString( $result, "rho",   PyFloat_FromDouble( $1.rho ) );
  PyDict_SetItemString( $result, "kappa", PyFloat_FromDouble( $1.kappa ) );
  PyDict_SetItemString( $result, "n",     PyInt_FromLong( $1.n ) );
}
%typemap(out) elem::SafeProduct<elem::Complex<float> > {
  $result = PyDict_New();
  PyDict_SetItemString( $result, "rho",   PyComplex_FromDoubles( $1.rho.real, $1.rho.imag ) );
  PyDict_SetItemString( $result, "kappa", PyFloat_FromDouble( $1.kappa ) );
  PyDict_SetItemString( $result, "n",     PyInt_FromLong( $1.n ) );
}
%typemap(out) elem::SafeProduct<elem::Complex<double> > {
  $result = PyDict_New();
  PyDict_SetItemString( $result, "rho",   PyComplex_FromDoubles( $1.rho.real, $1.rho.imag ) );
  PyDict_SetItemString( $result, "kappa", PyFloat_FromDouble( $1.kappa ) );
  PyDict_SetItemString( $result, "n",     PyInt_FromLong( $1.n ) );
}
%define TYPEMAP_MATRIX(T,U)
%typemap(out) elem::Matrix<T,int>&
{ $result = create_npmatrix( $1->Height(), $1->Width(), $1->LDim(), sizeof(T), U, true, $1->Buffer() ); }
%typemap(out) const elem::Matrix<T,int>&
{ $result = create_npmatrix( $1->Height(), $1->Width(), $1->LDim(), sizeof(T), U, false, $1->LockedBuffer() ); }
%enddef
TYPEMAP_MATRIX(float,NPY_FLOAT)
TYPEMAP_MATRIX(double,NPY_DOUBLE)
TYPEMAP_MATRIX(elem::Complex<float>,NPY_CFLOAT)
TYPEMAP_MATRIX(elem::Complex<double>,NPY_CDOUBLE)

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

%ignore *::operator=;
%rename(MPI_Initialize) elem::mpi::Initialize;
%rename(MPI_Initialized) elem::mpi::Initialized;
%rename(MPI_Finalize) elem::mpi::Finalize;
%rename(MPI_Finalized) elem::mpi::Finalized;

/*
 * TYPES, GRID, MPI
 */

// We do not need to %include complex_decl.hpp because we are using typemaps to convert
// Elemental complex values to native Python complex values. Using %import prevents SWIG
// from generating any wrappers.
%import  "elemental/core/complex_decl.hpp"
%include "elemental/core/types_decl.hpp"
%include "elemental/core/environment_decl.hpp"
%include "elemental/core/imports/mpi.hpp"
%include "elemental/core/grid_decl.hpp"

/*
 * MATRIX
 */

%ignore elem::Matrix::Matrix( Int, Int, const T*, Int );
%ignore elem::Matrix::Matrix( Int, Int, T*, Int );

%include "elemental/core/matrix.hpp"

namespace elem {
%template(Matrix_i) Matrix<int,int>;
%template(Matrix_s) Matrix<float,int>;
%template(Matrix_d) Matrix<double,int>;
%template(Matrix_c) Matrix<Complex<float>,int>;
%template(Matrix_z) Matrix<Complex<double>,int>;
};

/*
 * ABSTRACTDISTMATRIX
 */

%include "elemental/core/dist_matrix_forward_decl.hpp"
%include "elemental/core/dist_matrix.hpp"
%include "elemental/core/dist_matrix/abstract.hpp"

namespace elem {
%template(AbstractDistMatrix_i) AbstractDistMatrix<int,int>;
%template(AbstractDistMatrix_s) AbstractDistMatrix<float,int>;
%template(AbstractDistMatrix_d) AbstractDistMatrix<double,int>;
%template(AbstractDistMatrix_c) AbstractDistMatrix<Complex<float>,int>;
%template(AbstractDistMatrix_z) AbstractDistMatrix<Complex<double>,int>;
};

/*
 * DISTMATRIX
 */

%ignore elem::DistMatrix::DistMatrix( Int, Int, const T*, Int, const elem::Grid& );
%ignore elem::DistMatrix::DistMatrix( Int, Int, Int, const T*, Int, const elem::Grid& );
%ignore elem::DistMatrix::DistMatrix( Int, Int, Int, Int, const T*, Int, const elem::Grid& );
%ignore elem::DistMatrix::DistMatrix( Int, Int, T*, Int, elem::Grid& );
%ignore elem::DistMatrix::DistMatrix( Int, Int, Int, T*, Int, elem::Grid& );
%ignore elem::DistMatrix::DistMatrix( Int, Int, Int, Int, T*, Int, elem::Grid& );

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

namespace elem {
%template(DistMatrix_i) DistMatrix<int,MC,MR,int>;
%template(DistMatrix_s) DistMatrix<float,MC,MR,int>;
%template(DistMatrix_d) DistMatrix<double,MC,MR,int>;
%template(DistMatrix_c) DistMatrix<Complex<float>,MC,MR,int>;
%template(DistMatrix_z) DistMatrix<Complex<double>,MC,MR,int>;

%define DISTMATRIX(F,U,V,sfx) 
%template(DistMatrix_ ## sfx ## _ ## U ## _ ## V) DistMatrix<F,U,V,int>;
%enddef
%define DISTMATRIX_int(U,V) 
DISTMATRIX(int,U,V,i)
%enddef
%define DISTMATRIX_real(U,V)
DISTMATRIX(float,U,V,s)
DISTMATRIX(double,U,V,d)
%enddef
%define DISTMATRIX_cplx(U,V)
DISTMATRIX(Complex<float>,U,V,c)
DISTMATRIX(Complex<double>,U,V,z)
%enddef
%define DISTMATRIX_noint(U,V)
DISTMATRIX_real(U,V)
DISTMATRIX_cplx(U,V)
%enddef
%define DISTMATRIX_all(U,V)
DISTMATRIX_int(U,V)
DISTMATRIX_noint(U,V)
%enddef

DISTMATRIX_all(MC,STAR)
DISTMATRIX_all(MD,STAR)
DISTMATRIX_all(MR,MC)
DISTMATRIX_all(MR,STAR)
DISTMATRIX_all(STAR,MC)
DISTMATRIX_all(STAR,MD)
DISTMATRIX_all(STAR,MR)
DISTMATRIX_all(STAR,STAR)
DISTMATRIX_all(STAR,VC)
DISTMATRIX_all(STAR,VR)
DISTMATRIX_all(VC,STAR)
DISTMATRIX_all(VR,STAR)

};

/*
 * OVERLOADED FUNCTION MACROS 
 * These macros simplify the process of instantiating functions that are overloaded for
 * each Matrix and DistMatrix class. Since some functions are templated on all data
 * types, others only on the reals, others omitting integers, etc. we need a variety
 * of macros to handle all cases.
 */
 
%define OVERLOAD0_R(X,Y)
%template(Y) X<float>;
%template(Y) X<double>;
%template(Y) X<Complex<float> >;
%template(Y) X<Complex<double> >;
%enddef
%define OVERLOAD0_cpx_R(X,Y)
%template(Y) X<float>;
%template(Y) X<double>;
%enddef
%define OVERLOAD0_int_R(X,Y)
%template(Y) X<int>;
%template(Y) X<float>;
%template(Y) X<double>;
%template(Y) X<Complex<float> >;
%template(Y) X<Complex<double> >;
%enddef
%define OVERLOAD1_R(X,Y)
%template(Y) X<float,MC,MR>;
%template(Y) X<double,MC,MR>;
%template(Y) X<Complex<float>,MC,MR>;
%template(Y) X<Complex<double>,MC,MR>;
%enddef
%define OVERLOAD1_cpx_R(X,Y)
%template(Y) X<float,MC,MR>;
%template(Y) X<double,MC,MR>;
%enddef
// TODO: Extend this to other cases, needed for Display (and Norm)
%define OVERLOAD1_int_R_UV(X,Y,U,V)
%template(Y) X<int,U,V>;
%template(Y) X<float,U,V>;
%template(Y) X<double,U,V>;
%template(Y) X<Complex<float>,U,V>;
%template(Y) X<Complex<double>,U,V>;
%enddef
#define OVERLOAD1_int_UV(X,U,V) OVERLOAD1_int_R_UV(X,X,U,V)
#define OVERLOAD1_int_R(X,Y) OVERLOAD1_int_R_UV(X,Y,MC,MR)
#define OVERLOAD1_int(X) OVERLOAD1_int_R(X,X)
// TODO: Extend
%define OVERLOAD_VIEW_seq(X)
%template(X) X<int,int>;
%template(X) X<float,int>;
%template(X) X<double,int>;
%template(X) X<Complex<float>,int>;
%template(X) X<Complex<double>,int>;
%enddef
%define OVERLOAD_VIEW_UV(X,U,V)
%template(X) X<int,U,V,int>;
%template(X) X<float,U,V,int>;
%template(X) X<double,U,V,int>;
%template(X) X<Complex<float>,U,V,int>;
%template(X) X<Complex<double>,U,V,int>;
%enddef
#define OVERLOAD_VIEW(X) \
  OVERLOAD_VIEW_seq(X) \
  OVERLOAD_VIEW_UV(X,MC,MR) \

%define OVERLOAD2_R(X,Y)
%template(Y) X<float,MC,MR,MC,MR>;
%template(Y) X<double,MC,MR,MC,MR>;
%template(Y) X<Complex<float>,MC,MR,MC,MR>;
%template(Y) X<Complex<double>,MC,MR,MC,MR>;
%enddef
%define OVERLOAD2_cpx_R(X,Y)
%template(Y) X<float,MC,MR,MC,MR>;
%template(Y) X<double,MC,MR,MC,MR>;
%enddef
%define OVERLOAD2_int_R(X,Y)
%template(Y) X<int,MC,MR,MC,MR>;
%template(Y) X<float,MC,MR,MC,MR>;
%template(Y) X<double,MC,MR,MC,MR>;
%template(Y) X<Complex<float>,MC,MR,MC,MR>;
%template(Y) X<Complex<double>,MC,MR,MC,MR>;
%enddef

#define OVERLOAD01_R(X,Y)      OVERLOAD0_R(X,Y)     OVERLOAD1_R(X,Y)
#define OVERLOAD01_cpx_R(X,Y)  OVERLOAD0_cpx_R(X,Y) OVERLOAD1_cpx_R(X,Y)
#define OVERLOAD01_int_R(X,Y)  OVERLOAD0_int_R(X,Y) OVERLOAD1_int_R(X,Y)
#define OVERLOAD02_R(X,Y)      OVERLOAD0_R(X,Y)     OVERLOAD2_R(X,Y)
#define OVERLOAD02_cpx_R(X,Y)  OVERLOAD0_cpx_R(X,Y) OVERLOAD2_cpx_R(X,Y)
#define OVERLOAD02_int_R(X,Y)  OVERLOAD0_int_R(X,Y) OVERLOAD2_int_R(X,Y)
#define OVERLOAD012_R(X,Y)     OVERLOAD0_R(X,Y)     OVERLOAD1_R(X,Y)     OVERLOAD2_R(X,Y)
#define OVERLOAD012_cpx_R(X,Y) OVERLOAD0_cpx_R(X,Y) OVERLOAD1_cpx_R(X,Y) OVERLOAD2_cpx_R(X,Y)
#define OVERLOAD012_int_R(X,Y) OVERLOAD0_int_R(X,Y) OVERLOAD1_int_R(X,Y) OVERLOAD2_int_R(X,Y)

#define OVERLOAD0(X)       OVERLOAD0_R(X,X)
#define OVERLOAD0_cpx(X)   OVERLOAD0_cpx_R(X,X)
#define OVERLOAD0_int(X)   OVERLOAD0_int_R(X,X)
#define OVERLOAD1(X)       OVERLOAD1_R(X,X)
#define OVERLOAD1_cpx(X)   OVERLOAD1_cpx_R(X,X)
#define OVERLOAD1_int(X)   OVERLOAD1_int_R(X,X)
#define OVERLOAD2(X)       OVERLOAD2_R(X,X)
#define OVERLOAD2_cpx(X)   OVERLOAD2_cpx_R(X,X)
#define OVERLOAD2_int(X)   OVERLOAD2_int_R(X,X)
#define OVERLOAD01(X)      OVERLOAD01_R(X,X)
#define OVERLOAD01_cpx(X)  OVERLOAD01_cpx_R(X,X)
#define OVERLOAD01_int(X)  OVERLOAD01_int_R(X,X)
#define OVERLOAD02(X)      OVERLOAD02_R(X,X)
#define OVERLOAD02_cpx(X)  OVERLOAD02_cpx_R(X,X)
#define OVERLOAD02_int(X)  OVERLOAD02_int_R(X,X)
#define OVERLOAD012(X)     OVERLOAD012_R(X,X)
#define OVERLOAD012_cpx(X) OVERLOAD012_cpx_R(X,X)
#define OVERLOAD012_int(X) OVERLOAD012_int_R(X,X)

%define NO_OVERLOAD(X,...)
%rename(name ## _i) name<int>(__VA_ARGS__);
%rename(name ## _s) name<float>(__VA_ARGS__);
%rename(name ## _d) name<double>(__VA_ARGS__);
%rename(name ## _c) name<Complex<float> >(__VA_ARGS__);
%rename(name ## _z) name<Complex<double> >(__VA_ARGS__);
%enddef

/*
 * VIEWING
 */

%include "elemental/core/view_decl.hpp"
%include "elemental/core/partition_decl.hpp"

namespace elem {

OVERLOAD_VIEW(View)
OVERLOAD_VIEW(LockedView)
OVERLOAD_VIEW(View1x2)
OVERLOAD_VIEW(LockedView1x2)
OVERLOAD_VIEW(View2x1)
OVERLOAD_VIEW(LockedView2x1)
OVERLOAD_VIEW(View2x2)
OVERLOAD_VIEW(LockedView2x2)

OVERLOAD_VIEW(PartitionUp)
OVERLOAD_VIEW(LockedPartitionUp)
OVERLOAD_VIEW(PartitionDown)
OVERLOAD_VIEW(LockedPartitionDown)
OVERLOAD_VIEW(PartitionLeft)
OVERLOAD_VIEW(LockedPartitionLeft)
OVERLOAD_VIEW(PartitionRight)
OVERLOAD_VIEW(LockedPartitionRight)
OVERLOAD_VIEW(PartitionUpDiagonal)
OVERLOAD_VIEW(LockedPartitionUpDiagonal)
OVERLOAD_VIEW(PartitionUpLeftDiagonal)
OVERLOAD_VIEW(LockedPartitionUpLeftDiagonal)
OVERLOAD_VIEW(PartitionUpRightDiagonal)
OVERLOAD_VIEW(LockedPartitionUpRightDiagonal)
OVERLOAD_VIEW(PartitionDownDiagonal)
OVERLOAD_VIEW(LockedPartitionDownDiagonal)
OVERLOAD_VIEW(PartitionDownLeftDiagonal)
OVERLOAD_VIEW(LockedPartitionDownLeftDiagonal)
OVERLOAD_VIEW(PartitionDownRightDiagonal)
OVERLOAD_VIEW(LockedPartitionDownRightDiagonal)

};

/*
 * BLAS MISCELLANEOUS
 */
 
NO_OVERLOAD(SetLocalSymvBlocksize,int);
NO_OVERLOAD(LocalSymvBlocksize);
NO_OVERLOAD(SetLocalTrrkBlocksize,int);
NO_OVERLOAD(LocalTrrkBlocksize);
NO_OVERLOAD(SetLocalTrr2kBlocksize,int);
NO_OVERLOAD(LocalTrr2kBlocksize);
%include "elemental/blas-like_decl.hpp"

/*
 * BLAS LEVEL 1
 */

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

namespace elem {
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
};

/*
 * BLAS LEVEL 2
 */

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

// Most of these could support ints (except for Trsv)
namespace elem {
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
};

/* 
 * BLAS LEVEL 3
 */ 
 
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

namespace elem {
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
};

/*
 * INVARIANTS, INNER PRODUCTS, NORMS, ETC.
 */
 
%ignore PivotOp;
%ignore CreatePivotOp;
%ignore DestroyPivotOp;

// TODO: Instantiate non-symmetric/Hermitian Norms for all distributions
%define OVERLOAD_NORM(name)
OVERLOAD01(name ## Norm)
OVERLOAD0(Hermitian ## name ## Norm)
OVERLOAD0(Symmetric ## name ## Norm)
%enddef

%include "elemental/lapack-like/ConditionNumber.hpp"
%include "elemental/lapack-like/Determinant.hpp"
%include "elemental/lapack-like/Norm.hpp"
%include "elemental/lapack-like/Norm/Entrywise.hpp"
%include "elemental/lapack-like/Norm/EntrywiseOne.hpp"
%include "elemental/lapack-like/Norm/Frobenius.hpp"
%include "elemental/lapack-like/Norm/Infinity.hpp"
%include "elemental/lapack-like/Norm/KyFan.hpp"
%include "elemental/lapack-like/Norm/Max.hpp"
%include "elemental/lapack-like/Norm/One.hpp"
%include "elemental/lapack-like/Norm/Nuclear.hpp"
%include "elemental/lapack-like/Norm/Schatten.hpp"
%include "elemental/lapack-like/Norm/Two.hpp"
%include "elemental/lapack-like/Norm/TwoLowerBound.hpp"
%include "elemental/lapack-like/Norm/TwoUpperBound.hpp"
%include "elemental/lapack-like/Norm/Zero.hpp"
%include "elemental/lapack-like/Trace.hpp"
%include "elemental/lapack-like/Hadamard.hpp"
%include "elemental/lapack-like/HilbertSchmidt.hpp"

namespace elem {
OVERLOAD01(ConditionNumber)
OVERLOAD0(Determinant)
// SWIG doesn't like empty macro arguments. Who knows, maybe C++ doesn't either
// OVERLOAD_NORM()
OVERLOAD01(Norm)
OVERLOAD0(HermitianNorm)
OVERLOAD0(SymmetricNorm)
OVERLOAD_NORM(Entrywise)
OVERLOAD_NORM(EntrywiseOne)
OVERLOAD_NORM(Frobenius)
OVERLOAD_NORM(Infinity)
OVERLOAD_NORM(KyFan)
OVERLOAD_NORM(Max)
OVERLOAD_NORM(One)
OVERLOAD_NORM(Nuclear)
OVERLOAD_NORM(Schatten)
OVERLOAD_NORM(Two)
OVERLOAD0(TwoNormLowerBound)
OVERLOAD0(TwoNormUpperBound)
OVERLOAD01(ZeroNorm)
OVERLOAD0(Trace)
OVERLOAD01(Hadamard)
OVERLOAD01(HilbertSchmidt)
};

/*
 * FACTORIZATIONS
 */

%ignore elem::LocalCholesky;
%ignore elem::LocalLDL;
%ignore elem::LocalLU;

%include "elemental/lapack-like/Cholesky.hpp"
%include "elemental/lapack-like/LDL.hpp"
%include "elemental/lapack-like/LU.hpp"
%include "elemental/lapack-like/LQ.hpp"
%include "elemental/lapack-like/QR.hpp"
%include "elemental/lapack-like/QR/BusingerGolub.hpp"
%include "elemental/lapack-like/ID.hpp"
%include "elemental/lapack-like/Skeleton.hpp"

namespace elem {
OVERLOAD0(Cholesky)
OVERLOAD0(LDLH)
OVERLOAD0(LDLT)
OVERLOAD0(LU)
OVERLOAD0_cpx(LQ)
OVERLOAD0_cpx(QR)
namespace qr { 
OVERLOAD0_cpx(BusingerGolub) 
};
OVERLOAD0(ID)
OVERLOAD0(Skeleton)
};

/*
 * LINEAR SOLVERS
 */
 
%include "elemental/lapack-like/HPDSolve.hpp"
%include "elemental/lapack-like/GaussianElimination.hpp"
%include "elemental/lapack-like/LeastSquares.hpp"
%include "elemental/lapack-like/Cholesky/SolveAfter.hpp"
%include "elemental/lapack-like/LU/SolveAfter.hpp"

namespace elem {
OVERLOAD0(HPDSolve)
OVERLOAD0(RowEchelon)
OVERLOAD0(GaussianElimination)
OVERLOAD0_cpx(LeastSquares)
namespace cholesky { OVERLOAD0_R(SolveAfter,SolveAfter_Cholesky) };
namespace lu       { OVERLOAD0_R(SolveAfter,SolveAfter_LU) };
};

/*
 * INVERSION
 */
 
%ignore elem::LocalInverse;
%ignore elem::LocalHPDInverse;
%ignore elem::LocalTriangularInverse;
 
%include "elemental/lapack-like/Inverse.hpp"
%include "elemental/lapack-like/TriangularInverse.hpp"

namespace elem {
OVERLOAD0(Inverse)
OVERLOAD0(HPDInverse)
OVERLOAD0(TriangularInverse)
};

/*
 * REDUCTION TO CONDENSED FORM
 */

%include "elemental/lapack-like_decl.hpp"
%include "elemental/lapack-like/Bidiag.hpp"

namespace elem {
OVERLOAD0_cpx(HermitianTridiag)
OVERLOAD0_cpx(Bidiag)
};

/*
 * EIGENSOLVERS AND SVD
 */
 
%include "elemental/lapack-like/HermitianEig/Sort.hpp"
%include "elemental/lapack-like/SkewHermitianEig.hpp"
%include "elemental/lapack-like/HermitianGenDefiniteEig.hpp"
%include "elemental/lapack-like/SVD.hpp"
%include "elemental/lapack-like/Polar.hpp"
%include "elemental/lapack-like/Polar/Halley.hpp"
// %include "elemental/lapack-like/Polar/QDWH.hpp"
%include "elemental/lapack-like/HPSDCholesky.hpp"
%include "elemental/lapack-like/HPSDSquareRoot.hpp"

namespace elem {
OVERLOAD0(HermitianEig)
namespace hermitian_eig { OVERLOAD0_cpx_R(Sort,Sort_HermitianEig) };
OVERLOAD0_cpx(SkewHermitianEig)
OVERLOAD0(HermitianGenDefiniteEig)
OVERLOAD0(HermitianSVD)
OVERLOAD0(Polar)
namespace polar { 
OVERLOAD0(Halley)
//OVERLOAD0(QDWH)
};
OVERLOAD0(SVD)
};

/*
 * MATRIX FUNCTIONS
 */
 
%include "elemental/lapack-like/Pseudoinverse.hpp"
%include "elemental/lapack-like/HPSDSquareRoot.hpp"
%include "elemental/lapack-like/HPSDCholesky.hpp"
 
namespace elem {
OVERLOAD0(Pseudoinverse)
OVERLOAD0(HermitianPseudoinverse)
OVERLOAD0(HPSDSquareRoot)
OVERLOAD0_cpx(HPSDCholesky)
};

/*
 * CONVEX OPTIMIZATION
 */
 
%include "elemental/convex/LogBarrier.hpp"
%include "elemental/convex/LogDetDivergence.hpp"
%include "elemental/convex/SingularValueSoftThreshold.hpp"
%include "elemental/convex/SoftThreshold.hpp"
%include "elemental/convex/UnitaryCoherence.hpp" 

namespace elem {
OVERLOAD0(LogBarrier)
OVERLOAD0(LogDetDivergence)
OVERLOAD0_cpx(SingularValueSoftThreshold)
OVERLOAD0(SoftThreshold)
OVERLOAD0(UnitaryCoherence)
}

/*
 * SPECIAL MATRICES
 */

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

};

/*
 * I/O
 */

%include "elemental/io/Display.hpp"
#ifdef HAVE_QT5
%include "elemental/io/Spy.hpp"
#endif

namespace elem {

OVERLOAD0_cpx(Display)
OVERLOAD1_int_UV(Display,MC,MR)
OVERLOAD1_int_UV(Display,MC,STAR)
OVERLOAD1_int_UV(Display,MD,STAR)
OVERLOAD1_int_UV(Display,MR,MC)
OVERLOAD1_int_UV(Display,MR,STAR)
OVERLOAD1_int_UV(Display,STAR,MC)
OVERLOAD1_int_UV(Display,STAR,MD)
OVERLOAD1_int_UV(Display,STAR,MR)
OVERLOAD1_int_UV(Display,STAR,STAR)
OVERLOAD1_int_UV(Display,STAR,VC)
OVERLOAD1_int_UV(Display,STAR,VR)
OVERLOAD1_int_UV(Display,VC,STAR)
OVERLOAD1_int_UV(Display,VR,STAR)

#ifdef HAVE_QT5
OVERLOAD0_cpx(Spy)
OVERLOAD1_int_UV(Spy,MC,MR)
OVERLOAD1_int_UV(Spy,MC,STAR)
OVERLOAD1_int_UV(Spy,MD,STAR)
OVERLOAD1_int_UV(Spy,MR,MC)
OVERLOAD1_int_UV(Spy,MR,STAR)
OVERLOAD1_int_UV(Spy,STAR,MC)
OVERLOAD1_int_UV(Spy,STAR,MD)
OVERLOAD1_int_UV(Spy,STAR,MR)
OVERLOAD1_int_UV(Spy,STAR,STAR)
OVERLOAD1_int_UV(Spy,STAR,VC)
OVERLOAD1_int_UV(Spy,STAR,VR)
OVERLOAD1_int_UV(Spy,VC,STAR)
OVERLOAD1_int_UV(Spy,VR,STAR)
#endif

};
