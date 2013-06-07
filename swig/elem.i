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

%include <exception.i>

%{
#include "elemental.hpp"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#include "numpy/ndarraytypes.h"
#include <cstring>

template <typename T> struct NPY { static const int DType = NPY_VOID; };
template <> struct NPY<int> { static const int DType = sizeof(int) == NPY_SIZEOF_INT ? NPY_INT : NPY_LONG; };
template <> struct NPY<float> { static const int DType = NPY_FLOAT; };
template <> struct NPY<double> { static const int DType = NPY_DOUBLE; };
template <> struct NPY<elem::Complex<float> > { static const int DType = NPY_CFLOAT; };
template <> struct NPY<elem::Complex<double> > { static const int DType = NPY_CDOUBLE; };

static void Finalize_if()
{
    if ( elem::Initialized() )
        elem::Finalize();
}

class SwigException : public std::exception {
    int _type;
    const char* _msg;
public:
    SwigException( int type, const char* msg ) : _type(type), _msg(msg) {}
    const char* what() const throw() { return _msg; }
    int type() const throw() { return _type; }
};
    
template <typename T>
class ElemPyMatrix : public elem::Matrix<T>
{
    PyArrayObject *input, *darray;
    int HA, WA, LDimA;
public:
    ElemPyMatrix( PyObject*, bool );
    void UpdateArray();
};

template <typename T>
static PyObject* create_npmatrix( const elem::Matrix<T,int>& matrix, bool writable )
{
    // TO DO: replace PyArray_Type below with the type object for numpy.matrix
    npy_intp dims[2], strides[2];
    dims[0] = matrix.Height();
    dims[1] = matrix.Width();
    strides[0] = sizeof(T); 
    strides[1] = sizeof(T) * matrix.LDim();
    if ( matrix.Locked() ) writable = false;
    T* data = const_cast<T*>(matrix.LockedBuffer());
    return PyArray_NewFromDescr( 
        &PyArray_Type, PyArray_DescrFromType(NPY<T>::DType), 
        2, &dims[0], &strides[0], data,
        writable ? NPY_ARRAY_WRITEABLE : 0, NULL );
}

static bool check_elematrix( PyObject* obj, int DType, bool writable )
{
    if ( obj == 0 ) return false;
    if ( !PyArray_Check( obj ) ) return false;
    PyArrayObject* o = reinterpret_cast<PyArrayObject*>(obj);
    if ( PyArray_TYPE( o ) != DType ) return false;
    if ( writable && !PyArray_ISWRITEABLE( o ) ) return false;
    return true;
}

static bool get_HWL( PyArrayObject* obj, size_t DSize, int& HA, int& WA, int& LDimA )
{
    npy_intp ndim  = PyArray_NDIM( obj ), *dims, *strs;
    if ( ndim != 0 ) {
        dims = PyArray_DIMS( obj );
        strs = PyArray_STRIDES( obj );
    }
    bool recast = false;
    switch ( ndim ) {
    case 0:
        HA = WA = LDimA = 1;
        break;
    case 1:
        HA = 1;
        WA = dims[0];
        LDimA = strs[0] / DSize;
        recast = LDimA * DSize != strs[0] || LDimA < 1;
        break;
    case 2:
        HA = dims[0];
        WA = dims[1];
        LDimA = strs[1] / DSize;
        recast = LDimA * DSize != strs[1] || LDimA < HA || strs[0] != DSize;
        break;
    default:
        HA = 1;
        for ( int k = 0 ; k < ndim - 1 ; ++k ) {
            recast = recast || strs[k] != HA * DSize;
            HA *= dims[k];
        }
        WA = dims[ndim-1];
        LDimA = strs[ndim-1] / DSize;
        recast = recast || LDimA * DSize != strs[ndim-1] || LDimA < HA;
        break;
    }
    return recast;
}    

template <class T>
ElemPyMatrix<T>::ElemPyMatrix( PyObject* o, bool writable ) :
input(0), darray(0)
{
    if ( !PyArray_Check( o ) )
        throw SwigException( SWIG_TypeError, "NumPy array or matrix expected" );
    PyArrayObject* obj = reinterpret_cast<PyArrayObject*>(o);
    if ( PyArray_TYPE( obj ) != NPY<T>::DType )
        throw SwigException( SWIG_TypeError, "Incompatible NumPy data type encountered" );
    if ( !PyArray_ISWRITEABLE( obj ) ) {
        if ( writable )
            throw SwigException( SWIG_TypeError, "Incompatible NumPy data type encountered" );
        writable = false;
    }
    input = obj;
    bool recast = get_HWL( obj, sizeof(T), HA, WA, LDimA );
    // Why 3? Because the owner is 1, and the two layers of Python function calls that
    // sit between the owner and this line of code each add one more.
    bool owner = writable && PyArray_BASE( obj ) == NULL &&
        PyArray_CHKFLAGS( obj, NPY_ARRAY_OWNDATA ) != 0 &&
        PyArray_REFCOUNT( obj ) <= 3;
    if ( recast ) {
        darray = reinterpret_cast<PyArrayObject*>( PyArray_NewLikeArray( obj, NPY_FORTRANORDER, NULL, 0 ) );
        if ( darray == 0 )
            throw SwigException( SWIG_RuntimeError, "Unable to create the transposed matrix" );
        int result = PyArray_CopyInto( darray, input );
        if ( result < 0 ) {
            Py_DECREF( darray );
            throw SwigException( SWIG_RuntimeError, "Cannot copy data into the transposed matrix" );
        }
        obj = darray;
        LDimA = HA;
    }
    T* data = reinterpret_cast<T*>( PyArray_DATA( obj ) );
    if ( owner ) {
        this->Control( HA, WA, data, LDimA );
    } else if ( writable ) {
        this->Attach( HA, WA, data, LDimA );
    } else {
        this->LockedAttach( HA, WA, data, LDimA );
    }
}

template <class T>
void ElemPyMatrix<T>::UpdateArray()
{
    if ( !this->Viewing() ) {
        int HM = this->Height(), WM = this->Width(), LDimM = this->LDim();
        if ( HA != HM || WA != WM || LDimA != LDimM ) {
            T* ndata = reinterpret_cast<T*>( PyArray_DATA( input ) );
            T* odata = this->Buffer();
            bool need_copy = ndata != odata;
            PyArray_Dims ndims;
            npy_intp sdims[2];
            ndims.ptr = &sdims[0];
            ndims.len = 2;
            sdims[0] = LDimM;
            sdims[1] = WM;
            if ( PyArray_Resize( input, &ndims, 0, NPY_FORTRANORDER ) == NULL )
                throw SwigException( SWIG_RuntimeError, "Unable to modify the original NumPy matrix" );
            npy_intp* dims = PyArray_DIMS( input );
            npy_intp* strs = PyArray_STRIDES( input );
            dims[0] = HM;
            strs[0] = sizeof(T);
            strs[1] = sizeof(T) * LDimM;
            PyArray_UpdateFlags( input, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_F_CONTIGUOUS|NPY_ARRAY_ALIGNED );
            if ( need_copy ) {
                ndata = reinterpret_cast<T*>( PyArray_DATA( input ) );
                memcpy( ndata, odata, LDimM * WM * sizeof(T) );
            }
            return;
        }
    }
    if ( darray ) {
        int result = PyArray_CopyInto( input, darray );
        Py_DECREF( darray );
        if ( result < 0 )
            throw SwigException( SWIG_RuntimeError, "Unable to write results to the original NumPy matrix" );
        return;
    }
}

%}

%import "std_except.i"
%import "std_iostream.i"
%import "std_string.i"

%init %{
  import_array();
  PyObject *sys = PyImport_ImportModule("sys");
  PyObject *sysargv = PyObject_GetAttrString(sys,"argv");
  int argc = 0;
  char** argv = NULL;
  if ( sysargv )
    argc = PyList_Size( sysargv );
  if ( argc != 0 ) {
    argv = new char* [ argc + 1 ];
    for ( int i = 0 ; i != argc ; ++i ) {
        char *s = PyString_AsString( PyList_GetItem( sysargv, i ) );
        if ( s == NULL ) { argc = i; break; }
        argv[i] = s;
    }
    argv[argc] = 0;
  }
  elem::Initialize( argc, argv );
  Py_AtExit(Finalize_if);
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
  $1 = new int[1];
  $1[0] = PyList_Size($input);
  $2 = new char**[1];
  $2[0] = new char*[$1[0]+1];
  for (i = 0; i < $1[0]; i++) {
    PyObject *s = PyList_GetItem($input,i);
    if (!PyString_Check(s)) {
        delete[] $2[0];
        delete[] $1;
        PyErr_SetString(PyExc_ValueError, "List items must be strings");
        return NULL;
    }
    $2[0][i] = PyString_AsString(s);
  }
  $2[0][i] = 0;
}
%typemap(freearg) ( int& argc, char**& argv ) {
  if ( $2 && $2[0] ) delete[] $2[0];
  if ( $2 ) delete[] $2;
  if ( $1 ) delete[] $1;
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
%define TYPEMAPIN(T,ISCONST)
    try { 
        $1 = new ElemPyMatrix<T >( $input, ISCONST );
    } catch (SwigException exc) { 
        SWIG_exception( exc.type(), exc.what() ); 
    }
%enddef
%define TYPEMAPFREE(T)
    try {
        reinterpret_cast<ElemPyMatrix<T >*>($1)->UpdateArray();
        delete $1;
    } catch (SwigException exc) {
        delete $1;
        SWIG_exception( exc.type(), exc.what() );
    }
%enddef
%define TYPEMAP_MATRIX(T,V)
%typecheck(V)     const elem::Matrix<T    >& { $1 = check_elematrix( $input, NPY<T >::DType, false ); }
%typecheck(V)           elem::Matrix<T    >& { $1 = check_elematrix( $input, NPY<T >::DType, true  ); }
%typecheck(V)     const elem::Matrix<T,int>& { $1 = check_elematrix( $input, NPY<T >::DType, false ); }
%typecheck(V)           elem::Matrix<T,int>& { $1 = check_elematrix( $input, NPY<T >::DType, true  ); }
%typemap(in)      const elem::Matrix<T    >& { TYPEMAPIN(T,false) }
%typemap(in)            elem::Matrix<T    >& { TYPEMAPIN(T,true)  }
%typemap(in)      const elem::Matrix<T,int>& { TYPEMAPIN(T,false) }
%typemap(in)            elem::Matrix<T,int>& { TYPEMAPIN(T,true)  }
%typemap(freearg) const elem::Matrix<T    >& { TYPEMAPFREE(T) }
%typemap(freearg)       elem::Matrix<T    >& { TYPEMAPFREE(T) }
%typemap(freearg) const elem::Matrix<T,int>& { TYPEMAPFREE(T) }
%typemap(freearg)       elem::Matrix<T,int>& { TYPEMAPFREE(T) }
%typemap(out)     const elem::Matrix<T    >& { $result = create_npmatrix( *$1, false ); }
%typemap(out)           elem::Matrix<T    >& { $result = create_npmatrix( *$1, true );  }
%typemap(out)     const elem::Matrix<T,int>& { $result = create_npmatrix( *$1, false ); }
%typemap(out)           elem::Matrix<T,int>& { $result = create_npmatrix( *$1, true );  }
%enddef
TYPEMAP_MATRIX(int,SWIG_TYPECHECK_INT32_ARRAY)
TYPEMAP_MATRIX(float,SWIG_TYPECHECK_FLOAT_ARRAY)
TYPEMAP_MATRIX(double,SWIG_TYPECHECK_DOUBLE_ARRAY)
TYPEMAP_MATRIX(elem::Complex<float>,SWIG_TYPECHECK_CHAR_ARRAY)
TYPEMAP_MATRIX(elem::Complex<double>,SWIG_TYPECHECK_STRING_ARRAY)

/*
 * Blanket exception handling.
 */

%exception {
    try {
        $action
    } catch (std::exception exc) {
        const char* msg = exc.what();
        PyErr_SetString( PyExc_RuntimeError, msg ? msg : "Exception caught from Elemental" );
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

// We do not need to %include complex_decl.hpp or matrix.hpp, because we are using
// typemaps to convert the Elemental classes to equivalent Python and NumPy objects.
// Using %import prevents SWIG from generating any wrappers.
%import  "elemental/core/complex_decl.hpp"
%include "elemental/core/types_decl.hpp"
%include "elemental/core/environment_decl.hpp"
%include "elemental/core/imports/mpi.hpp"
%include "elemental/core/grid_decl.hpp"
%import "elemental/core/matrix.hpp"

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
  OVERLOAD_VIEW_UV(X,MC,MR) 
%define OVERLOAD_COPY(X,U1,V1,U2,V2)
%template(X) X<int,U1,V1,U2,V2>;
%template(X) X<float,U1,V1,U2,V2>;
%template(X) X<double,U1,V1,U2,V2>;
%template(X) X<Complex<float>,U1,V1,U2,V2>;
%template(X) X<Complex<double>,U1,V1,U2,V2>;
%enddef

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
OVERLOAD_VIEW_UV(View,MC,STAR)
OVERLOAD_VIEW_UV(View,MR,MC)
OVERLOAD_VIEW_UV(View,MR,STAR)
OVERLOAD_VIEW_UV(View,STAR,MC)
OVERLOAD_VIEW_UV(View,STAR,MR)
OVERLOAD_VIEW_UV(View,STAR,STAR)
OVERLOAD_VIEW_UV(View,STAR,VC)
OVERLOAD_VIEW_UV(View,STAR,VR)
OVERLOAD_VIEW_UV(View,VC,STAR)
OVERLOAD_VIEW_UV(View,VR,STAR)
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
 
%ignore elem::DiagonalSolve<T>(elem::left_or_right_wrapper::LeftOrRight,elem::orientation_wrapper::Orientation,elem::Matrix<Base<T>::type> const &,elem::Matrix<T> &,bool);
%ignore elem::DiagonalScale<T>(elem::left_or_right_wrapper::LeftOrRight,elem::orientation_wrapper::Orientation,elem::Matrix<Base<T>::type> const &,elem::Matrix<T> &,bool);

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
OVERLOAD0_int(Copy)
OVERLOAD_COPY(Copy,MC,MR,MC,MR)
OVERLOAD_COPY(Copy,MC,MR,MC,STAR)
OVERLOAD_COPY(Copy,MC,MR,MD,STAR)
OVERLOAD_COPY(Copy,MC,MR,MR,MC)
OVERLOAD_COPY(Copy,MC,MR,MR,STAR)
OVERLOAD_COPY(Copy,MC,MR,STAR,MC)
OVERLOAD_COPY(Copy,MC,MR,STAR,MD)
OVERLOAD_COPY(Copy,MC,MR,STAR,MR)
OVERLOAD_COPY(Copy,MC,MR,STAR,STAR)
OVERLOAD_COPY(Copy,MC,MR,STAR,VC)
OVERLOAD_COPY(Copy,MC,MR,STAR,VR)
OVERLOAD_COPY(Copy,MC,MR,VC,STAR)
OVERLOAD_COPY(Copy,MC,MR,VR,STAR)
OVERLOAD_COPY(Copy,MC,STAR,MC,MR)
OVERLOAD_COPY(Copy,MC,STAR,MC,STAR)
OVERLOAD_COPY(Copy,MC,STAR,MD,STAR)
OVERLOAD_COPY(Copy,MC,STAR,MR,MC)
OVERLOAD_COPY(Copy,MC,STAR,MR,STAR)
OVERLOAD_COPY(Copy,MC,STAR,STAR,MC)
OVERLOAD_COPY(Copy,MC,STAR,STAR,MD)
OVERLOAD_COPY(Copy,MC,STAR,STAR,MR)
OVERLOAD_COPY(Copy,MC,STAR,STAR,STAR)
OVERLOAD_COPY(Copy,MC,STAR,STAR,VC)
OVERLOAD_COPY(Copy,MC,STAR,STAR,VR)
OVERLOAD_COPY(Copy,MC,STAR,VC,STAR)
OVERLOAD_COPY(Copy,MC,STAR,VR,STAR)
OVERLOAD_COPY(Copy,MD,STAR,MC,MR)
OVERLOAD_COPY(Copy,MD,STAR,MC,STAR)
OVERLOAD_COPY(Copy,MD,STAR,MD,STAR)
OVERLOAD_COPY(Copy,MD,STAR,MR,MC)
OVERLOAD_COPY(Copy,MD,STAR,MR,STAR)
OVERLOAD_COPY(Copy,MD,STAR,STAR,MC)
OVERLOAD_COPY(Copy,MD,STAR,STAR,MD)
OVERLOAD_COPY(Copy,MD,STAR,STAR,MR)
OVERLOAD_COPY(Copy,MD,STAR,STAR,STAR)
OVERLOAD_COPY(Copy,MD,STAR,STAR,VC)
OVERLOAD_COPY(Copy,MD,STAR,STAR,VR)
OVERLOAD_COPY(Copy,MD,STAR,VC,STAR)
OVERLOAD_COPY(Copy,MD,STAR,VR,STAR)
OVERLOAD_COPY(Copy,MR,MC,MC,MR)
OVERLOAD_COPY(Copy,MR,MC,MC,STAR)
OVERLOAD_COPY(Copy,MR,MC,MD,STAR)
OVERLOAD_COPY(Copy,MR,MC,MR,MC)
OVERLOAD_COPY(Copy,MR,MC,MR,STAR)
OVERLOAD_COPY(Copy,MR,MC,STAR,MC)
OVERLOAD_COPY(Copy,MR,MC,STAR,MD)
OVERLOAD_COPY(Copy,MR,MC,STAR,MR)
OVERLOAD_COPY(Copy,MR,MC,STAR,STAR)
OVERLOAD_COPY(Copy,MR,MC,STAR,VC)
OVERLOAD_COPY(Copy,MR,MC,STAR,VR)
OVERLOAD_COPY(Copy,MR,MC,VC,STAR)
OVERLOAD_COPY(Copy,MR,MC,VR,STAR)
OVERLOAD_COPY(Copy,MR,STAR,MC,MR)
OVERLOAD_COPY(Copy,MR,STAR,MC,STAR)
OVERLOAD_COPY(Copy,MR,STAR,MD,STAR)
OVERLOAD_COPY(Copy,MR,STAR,MR,MC)
OVERLOAD_COPY(Copy,MR,STAR,MR,STAR)
OVERLOAD_COPY(Copy,MR,STAR,STAR,MC)
OVERLOAD_COPY(Copy,MR,STAR,STAR,MD)
OVERLOAD_COPY(Copy,MR,STAR,STAR,MR)
OVERLOAD_COPY(Copy,MR,STAR,STAR,STAR)
OVERLOAD_COPY(Copy,MR,STAR,STAR,VC)
OVERLOAD_COPY(Copy,MR,STAR,STAR,VR)
OVERLOAD_COPY(Copy,MR,STAR,VC,STAR)
OVERLOAD_COPY(Copy,MR,STAR,VR,STAR)
OVERLOAD_COPY(Copy,STAR,MC,MC,MR)
OVERLOAD_COPY(Copy,STAR,MC,MC,STAR)
OVERLOAD_COPY(Copy,STAR,MC,MD,STAR)
OVERLOAD_COPY(Copy,STAR,MC,MR,MC)
OVERLOAD_COPY(Copy,STAR,MC,MR,STAR)
OVERLOAD_COPY(Copy,STAR,MC,STAR,MC)
OVERLOAD_COPY(Copy,STAR,MC,STAR,MD)
OVERLOAD_COPY(Copy,STAR,MC,STAR,MR)
OVERLOAD_COPY(Copy,STAR,MC,STAR,STAR)
OVERLOAD_COPY(Copy,STAR,MC,STAR,VC)
OVERLOAD_COPY(Copy,STAR,MC,STAR,VR)
OVERLOAD_COPY(Copy,STAR,MC,VC,STAR)
OVERLOAD_COPY(Copy,STAR,MC,VR,STAR)
OVERLOAD_COPY(Copy,STAR,MD,MC,MR)
OVERLOAD_COPY(Copy,STAR,MD,MC,STAR)
OVERLOAD_COPY(Copy,STAR,MD,MD,STAR)
OVERLOAD_COPY(Copy,STAR,MD,MR,MC)
OVERLOAD_COPY(Copy,STAR,MD,MR,STAR)
OVERLOAD_COPY(Copy,STAR,MD,STAR,MC)
OVERLOAD_COPY(Copy,STAR,MD,STAR,MD)
OVERLOAD_COPY(Copy,STAR,MD,STAR,MR)
OVERLOAD_COPY(Copy,STAR,MD,STAR,STAR)
OVERLOAD_COPY(Copy,STAR,MD,STAR,VC)
OVERLOAD_COPY(Copy,STAR,MD,STAR,VR)
OVERLOAD_COPY(Copy,STAR,MD,VC,STAR)
OVERLOAD_COPY(Copy,STAR,MD,VR,STAR)
OVERLOAD_COPY(Copy,STAR,MR,MC,MR)
OVERLOAD_COPY(Copy,STAR,MR,MC,STAR)
OVERLOAD_COPY(Copy,STAR,MR,MD,STAR)
OVERLOAD_COPY(Copy,STAR,MR,MR,MC)
OVERLOAD_COPY(Copy,STAR,MR,MR,STAR)
OVERLOAD_COPY(Copy,STAR,MR,STAR,MC)
OVERLOAD_COPY(Copy,STAR,MR,STAR,MD)
OVERLOAD_COPY(Copy,STAR,MR,STAR,MR)
OVERLOAD_COPY(Copy,STAR,MR,STAR,STAR)
OVERLOAD_COPY(Copy,STAR,MR,STAR,VC)
OVERLOAD_COPY(Copy,STAR,MR,STAR,VR)
OVERLOAD_COPY(Copy,STAR,MR,VC,STAR)
OVERLOAD_COPY(Copy,STAR,MR,VR,STAR)
OVERLOAD_COPY(Copy,STAR,STAR,MC,MR)
OVERLOAD_COPY(Copy,STAR,STAR,MC,STAR)
OVERLOAD_COPY(Copy,STAR,STAR,MD,STAR)
OVERLOAD_COPY(Copy,STAR,STAR,MR,MC)
OVERLOAD_COPY(Copy,STAR,STAR,MR,STAR)
OVERLOAD_COPY(Copy,STAR,STAR,STAR,MC)
OVERLOAD_COPY(Copy,STAR,STAR,STAR,MD)
OVERLOAD_COPY(Copy,STAR,STAR,STAR,MR)
OVERLOAD_COPY(Copy,STAR,STAR,STAR,STAR)
OVERLOAD_COPY(Copy,STAR,STAR,STAR,VC)
OVERLOAD_COPY(Copy,STAR,STAR,STAR,VR)
OVERLOAD_COPY(Copy,STAR,STAR,VC,STAR)
OVERLOAD_COPY(Copy,STAR,STAR,VR,STAR)
OVERLOAD_COPY(Copy,STAR,VC,MC,MR)
OVERLOAD_COPY(Copy,STAR,VC,MC,STAR)
OVERLOAD_COPY(Copy,STAR,VC,MD,STAR)
OVERLOAD_COPY(Copy,STAR,VC,MR,MC)
OVERLOAD_COPY(Copy,STAR,VC,MR,STAR)
OVERLOAD_COPY(Copy,STAR,VC,STAR,MC)
OVERLOAD_COPY(Copy,STAR,VC,STAR,MD)
OVERLOAD_COPY(Copy,STAR,VC,STAR,MR)
OVERLOAD_COPY(Copy,STAR,VC,STAR,STAR)
OVERLOAD_COPY(Copy,STAR,VC,STAR,VC)
OVERLOAD_COPY(Copy,STAR,VC,STAR,VR)
OVERLOAD_COPY(Copy,STAR,VC,VC,STAR)
OVERLOAD_COPY(Copy,STAR,VC,VR,STAR)
OVERLOAD_COPY(Copy,STAR,VR,MC,MR)
OVERLOAD_COPY(Copy,STAR,VR,MC,STAR)
OVERLOAD_COPY(Copy,STAR,VR,MD,STAR)
OVERLOAD_COPY(Copy,STAR,VR,MR,MC)
OVERLOAD_COPY(Copy,STAR,VR,MR,STAR)
OVERLOAD_COPY(Copy,STAR,VR,STAR,MC)
OVERLOAD_COPY(Copy,STAR,VR,STAR,MD)
OVERLOAD_COPY(Copy,STAR,VR,STAR,MR)
OVERLOAD_COPY(Copy,STAR,VR,STAR,STAR)
OVERLOAD_COPY(Copy,STAR,VR,STAR,VC)
OVERLOAD_COPY(Copy,STAR,VR,STAR,VR)
OVERLOAD_COPY(Copy,STAR,VR,VC,STAR)
OVERLOAD_COPY(Copy,STAR,VR,VR,STAR)
OVERLOAD_COPY(Copy,VC,STAR,MC,MR)
OVERLOAD_COPY(Copy,VC,STAR,MC,STAR)
OVERLOAD_COPY(Copy,VC,STAR,MD,STAR)
OVERLOAD_COPY(Copy,VC,STAR,MR,MC)
OVERLOAD_COPY(Copy,VC,STAR,MR,STAR)
OVERLOAD_COPY(Copy,VC,STAR,STAR,MC)
OVERLOAD_COPY(Copy,VC,STAR,STAR,MD)
OVERLOAD_COPY(Copy,VC,STAR,STAR,MR)
OVERLOAD_COPY(Copy,VC,STAR,STAR,STAR)
OVERLOAD_COPY(Copy,VC,STAR,STAR,VC)
OVERLOAD_COPY(Copy,VC,STAR,STAR,VR)
OVERLOAD_COPY(Copy,VC,STAR,VC,STAR)
OVERLOAD_COPY(Copy,VC,STAR,VR,STAR)
OVERLOAD_COPY(Copy,VR,STAR,MC,MR)
OVERLOAD_COPY(Copy,VR,STAR,MC,STAR)
OVERLOAD_COPY(Copy,VR,STAR,MD,STAR)
OVERLOAD_COPY(Copy,VR,STAR,MR,MC)
OVERLOAD_COPY(Copy,VR,STAR,MR,STAR)
OVERLOAD_COPY(Copy,VR,STAR,STAR,MC)
OVERLOAD_COPY(Copy,VR,STAR,STAR,MD)
OVERLOAD_COPY(Copy,VR,STAR,STAR,MR)
OVERLOAD_COPY(Copy,VR,STAR,STAR,STAR)
OVERLOAD_COPY(Copy,VR,STAR,STAR,VC)
OVERLOAD_COPY(Copy,VR,STAR,STAR,VR)
OVERLOAD_COPY(Copy,VR,STAR,VC,STAR)
OVERLOAD_COPY(Copy,VR,STAR,VR,STAR)
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
%template(SafeProduct_i) SafeProduct<int>;
%template(SafeProduct_s) SafeProduct<float>;
%template(SafeProduct_d) SafeProduct<double>;
%template(SafeProduct_c) SafeProduct<Complex<float> >;
%template(SafeProduct_z) SafeProduct<Complex<double> >;
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
OVERLOAD1_int_UV(Uniform,MC,STAR)
OVERLOAD1_int_UV(Uniform,MR,MC)
OVERLOAD1_int_UV(Uniform,MR,STAR)
OVERLOAD1_int_UV(Uniform,STAR,MC)
OVERLOAD1_int_UV(Uniform,STAR,MR)
OVERLOAD1_int_UV(Uniform,STAR,STAR)
OVERLOAD1_int_UV(Uniform,STAR,VC)
OVERLOAD1_int_UV(Uniform,STAR,VR)
OVERLOAD1_int_UV(Uniform,VC,STAR)
OVERLOAD1_int_UV(Uniform,VR,STAR)
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
