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
static void Finalize_if()
{
    if ( elem::Initialized() )
        elem::Finalize();
}
%}

%pythoncode %{
from elem_blas import *
from elem_lapack import *
from elem_view import *
from elem_matrices import *
from elem_convex import *
from elem_io import *
import elem_mpi
%}

%init %{
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

%include "common.swg"

%ignore *::operator=;
%ignore SingularMatrixException;
%ignore NonHPDMatrixException;
%ignore NonHPSDMatrixException;

/*
 * TYPES, GRID, MPI
 */

// We do not need to %include complex_decl.hpp or matrix.hpp, because we are using
// typemaps to convert the Elemental classes to equivalent Python and NumPy objects.
// Using %import prevents SWIG from generating any wrappers.

%import  "std_except.i"
%import  "elemental/core/complex_decl.hpp"
%include "elemental/core/types_decl.hpp"
%include "elemental/core/environment_decl.hpp"
%import  "elemental/core/imports/mpi.hpp"
%include "elemental/core/grid_decl.hpp"
%import  "elemental/core/matrix.hpp"

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
%ignore elem::DistMatrix::DistMatrix( Int, Int, const T*, Int, const elem::Grid&, int );
%ignore elem::DistMatrix::DistMatrix( Int, Int, Int, const T*, Int, const elem::Grid& );
%ignore elem::DistMatrix::DistMatrix( Int, Int, Int, Int, const T*, Int, const elem::Grid& );
%ignore elem::DistMatrix::DistMatrix( Int, Int, T*, Int, const elem::Grid& ); 
%ignore elem::DistMatrix::DistMatrix( Int, Int, T*, Int, const elem::Grid&, int );
%ignore elem::DistMatrix::DistMatrix( Int, Int, Int, T*, Int, const elem::Grid& );
%ignore elem::DistMatrix::DistMatrix( Int, Int, Int, Int, T*, Int, const elem::Grid& );

%include "elemental/core/dist_matrix/circ_circ.hpp"
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

%define DISTMATRIX(F,U,V,sfx)
%template(DistMatrix_ ## sfx) DistMatrix<F,U,V,int>;
%extend DistMatrix<F,U,V,int> {
	const char *__str__() {
		std::string ans;
		std::ostringstream msg;
		elem::Print( *$self, ans, msg );
		ans = msg.str();
		std::size_t found = ans.find_last_not_of(" \t\f\v\n\r");
		if ( found != std::string::npos ) 
			ans.erase( found + 1 );
		return ans.c_str();
	}
}
%enddef
%define DISTMATRIX_all(U,V)
DISTMATRIX(int,U,V,i_ ## U ## _ ## V)
DISTMATRIX(float,U,V,s_ ## U ## _ ## V)
DISTMATRIX(double,U,V,d_ ## U ## _ ## V)
DISTMATRIX(Complex<float>,U,V,c_ ## U ## _ ## V)
DISTMATRIX(Complex<double>,U,V,z_ ## U ## _ ## V)
%enddef

namespace elem {
DISTMATRIX(int,MC,MR,i)
DISTMATRIX(float,MC,MR,s)
DISTMATRIX(double,MC,MR,d)
DISTMATRIX(Complex<float>,MC,MR,c)
DISTMATRIX(Complex<double>,MC,MR,z)
DISTMATRIX_all(CIRC,CIRC)
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
