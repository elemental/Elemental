/*
   Copyright (c) 2009-2013, Jack Poulson
                      2013, Michael C. Grant
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

%include <typemaps.i>
%include <exception.i>

%{
#include "elemental.hpp"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#include "numpy/ndarraytypes.h"
%}

%init %{
  import_array();
%}

%{

template <typename T> struct NPY { static const int DType = NPY_VOID; };
template <> struct NPY<int> { static const int DType = sizeof(int) == NPY_SIZEOF_INT ? NPY_INT : NPY_LONG; };
template <> struct NPY<float> { static const int DType = NPY_FLOAT; };
template <> struct NPY<double> { static const int DType = NPY_DOUBLE; };
template <> struct NPY<elem::Complex<float> > { static const int DType = NPY_CFLOAT; };
template <> struct NPY<elem::Complex<double> > { static const int DType = NPY_CDOUBLE; };

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

/*
 * TYPEMAPS 
 */
 
/*
 * elem::Complex<T>
 */
  
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

/*
 * elem::SafeProduct<T>
 */

%typemap(in) elem::SafeProduct<int> {
	$result = PyTuple_Pack( 3,
		PyFloat_FromDouble( $1.rho ), 
		PyFloat_FromDouble( $1.kappa ), 
		PyInteger_FromLong( $1.n ) );
}
%typemap(in) elem::SafeProduct<float> {
	$result = PyTuple_Pack( 3,
		PyFloat_FromDouble( $1.rho ), 
		PyFloat_FromDouble( $1.kappa ), 
		PyInteger_FromLong( $1.n ) );
}
%typemap(in) elem::SafeProduct<double> {
	$result = PyTuple_Pack( 3,
		PyFloat_FromDouble( $1.rho ), 
		PyFloat_FromDouble( $1.kappa ), 
		PyInteger_FromLong( $1.n ) );
}
%typemap(in) elem::SafeProduct<elem::Complex<float> > {
	$result = PyTuple_Pack( 3,
		PyComplex_FromDoubles( $1.rho.real, $1.rho.imag ), 
		PyFloat_FromDouble( $1.kappa ), 
		PyInteger_FromLong( $1.n ) );
}
%typemap(in) elem::SafeProduct<elem::Complex<double> > {
	$result = PyTuple_Pack( 3,
		PyComplex_FromDoubles( $1.rho.real, $1.rho.imag ), 
		PyFloat_FromDouble( $1.kappa ), 
		PyInteger_FromLong( $1.n ) );
}

/*
 * Elemental's Matrix class
 */
 
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
%define OVERLOAD_VIEW(X)
OVERLOAD_VIEW_seq(X)
OVERLOAD_VIEW_UV(X,MC,MR) 
%enddef
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
