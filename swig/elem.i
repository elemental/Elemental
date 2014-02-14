/*
   Copyright (c) 2009-2014, Jack Poulson
                      2013, Michael C. Grant
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
%module elem

%include "common.swg"

%{
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
from elem_control import *
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

%ignore *::operator=;
%ignore SingularMatrixException;
%ignore NonHPDMatrixException;
%ignore NonHPSDMatrixException;
%ignore *::Attach;
%ignore *::LockedAttach;
%ignore elem::Grid::FindFactor;
%ignore elem::PushCallStack;
%ignore elem::PopCallStack;
%ignore elem::DumpCallStack;
%ignore elem::CallStackEntry;
%ignore elem::ComplainIfDebug;

/*
 * TYPES, GRID, MPI
 */

// We do not need to %include complex/decl.hpp or matrix.hpp, because we are using
// typemaps to convert the Elemental classes to equivalent Python and NumPy objects.
// Using %import prevents SWIG from generating any wrappers.

#define GATTI(T) %attribute(elem::Grid,int,T,T)
#define GATTP(T) %attribute(elem::Grid,void*,T,T)
#define GATTB(T) %attribute(elem::Grid,bool,T,T)
GATTI(Row)
GATTI(Col)
GATTI(Rank)
GATTI(Height)
GATTI(Width)
GATTI(Size)
GATTP(ColComm)
GATTP(RowComm)
GATTP(Comm)
GATTI(MCRank)
GATTI(MRRank)
GATTI(VCRank)
GATTI(VRRank)
GATTI(MCSize)
GATTI(MRSize)
GATTI(VCSize)
GATTI(VRSize)
GATTP(MCComm)
GATTP(MRComm)
GATTP(VCComm)
GATTP(VRComm)
GATTI(GCD)
GATTI(LCM)
GATTB(InGrid)
GATTI(OwningRank)
GATTI(ViewingRank)
GATTP(OwningGroup)
GATTP(OwningComm)
GATTP(ViewingComm)
#undef GATTI
#undef GATTP
#undef GATTB

%include "elemental/core.hpp"
%import  "elemental/core/complex/decl.hpp"
%import  "elemental/core/imports/mpi.hpp"
%import  "elemental/core/imports/choice.hpp"
%include "elemental/core/imports/mpi_choice.hpp"
%include "elemental/core/types/decl.hpp"
%include "elemental/core/environment/decl.hpp"
%include "elemental/core/grid/decl.hpp"
%import  "elemental/core/matrix.hpp"

/*
 * ABSTRACTDISTMATRIX
 */

%ignore elem::AbstractDistMatrix::Buffer;
%ignore elem::AbstractDistMatrix::LockedBuffer;
%ignore elem::AbstractDistMatrix::DistData;

%define AbDMT(T,Type,N) 
%attribute_readonly(%arg(elem::AbstractDistMatrix<T >),Type,N,N,self_->N())
%enddef
%define AbDM(T)
AbDMT(T,Int,Height)
AbDMT(T,Int,Width)
AbDMT(T,Int,LocalHeight)
AbDMT(T,Int,LocalWidth)
AbDMT(T,Int,LDim)
AbDMT(T,size_t,AllocatedMemory)
%attribute_custom(%arg(elem::AbstractDistMatrix<T >),elem::Grid,Grid,Grid,SetGrid,&self_->Grid(),self_->SetGrid(*val_))
%attribute_readonly(%arg(elem::AbstractDistMatrix<T >),PyObject*,Matrix,Matrix,create_npmatrix(self_->Matrix(),true))
%attribute_readonly(%arg(elem::AbstractDistMatrix<T >),PyObject*,LockedMatrix,LockedMatrix,create_npmatrix(self_->LockedMatrix(),false))
AbDMT(T,bool,ColConstrained)
AbDMT(T,bool,RowConstrained)
AbDMT(T,Int,ColAlign)
AbDMT(T,Int,RowAlign)
AbDMT(T,Int,ColShift)
AbDMT(T,Int,RowShift)
AbDMT(T,Int,ColStride)
AbDMT(T,Int,RowStride)
AbDMT(T,Int,ColRank)
AbDMT(T,Int,RowRank)
AbDMT(T,bool,Locked)
AbDMT(T,bool,Viewing)
AbDMT(T,bool,Participating)
%enddef
AbDM(Int)
AbDM(float)
AbDM(double)
AbDM(elem::Complex<float>)
AbDM(elem::Complex<double>)
#undef AbDMT
#undef AbDM

%nodefaultctor DistData;

%include "elemental/core/dist_matrix/forward_decl.hpp"
%include "elemental/core/dist_matrix.hpp"
%include "elemental/core/dist_matrix/abstract.hpp"
%include "elemental/core/dist_matrix/general.hpp"

namespace elem {
%template(AbstractDistMatrix_i) AbstractDistMatrix<Int>;
%template(AbstractDistMatrix_s) AbstractDistMatrix<float>;
%template(AbstractDistMatrix_d) AbstractDistMatrix<double>;
%template(AbstractDistMatrix_c) AbstractDistMatrix<Complex<float> >;
%template(AbstractDistMatrix_z) AbstractDistMatrix<Complex<double> >;
};

/*
 * GeneralDistMatrix
 */

%define GENERALDISTMATRIX(F,U,V,sfx)
%template(GeneralDistMatrix_ ## sfx) GeneralDistMatrix<F,U,V>;
%enddef
%define GENERALDISTMATRIX_all(U,V)
GENERALDISTMATRIX(Int,U,V,i_ ## U ## _ ## V)
GENERALDISTMATRIX(float,U,V,s_ ## U ## _ ## V)
GENERALDISTMATRIX(double,U,V,d_ ## U ## _ ## V)
GENERALDISTMATRIX(Complex<float>,U,V,c_ ## U ## _ ## V)
GENERALDISTMATRIX(Complex<double>,U,V,z_ ## U ## _ ## V)
%enddef

namespace elem {
GENERALDISTMATRIX(Int,MC,MR,i)
GENERALDISTMATRIX(float,MC,MR,s)
GENERALDISTMATRIX(double,MC,MR,d)
GENERALDISTMATRIX(Complex<float>,MC,MR,c)
GENERALDISTMATRIX(Complex<double>,MC,MR,z)
GENERALDISTMATRIX_all(CIRC,CIRC)
GENERALDISTMATRIX_all(MC,STAR)
GENERALDISTMATRIX_all(MD,STAR)
GENERALDISTMATRIX_all(MR,MC)
GENERALDISTMATRIX_all(MR,STAR)
GENERALDISTMATRIX_all(STAR,MC)
GENERALDISTMATRIX_all(STAR,MD)
GENERALDISTMATRIX_all(STAR,MR)
GENERALDISTMATRIX_all(STAR,STAR)
GENERALDISTMATRIX_all(STAR,VC)
GENERALDISTMATRIX_all(STAR,VR)
GENERALDISTMATRIX_all(VC,STAR)
GENERALDISTMATRIX_all(VR,STAR)
};

/*
 * DISTMATRIX
 */

%ignore elem::DistMatrix::DistMatrix( Int, Int, const T*, Int, const elem::Grid& );
%ignore elem::DistMatrix::DistMatrix( Int, Int, const T*, Int, const elem::Grid&, Int );
%ignore elem::DistMatrix::DistMatrix( Int, Int, Int, const T*, Int, const elem::Grid& );
%ignore elem::DistMatrix::DistMatrix( Int, Int, Int, Int, const T*, Int, const elem::Grid& );
%ignore elem::DistMatrix::DistMatrix( Int, Int, T*, Int, const elem::Grid& ); 
%ignore elem::DistMatrix::DistMatrix( Int, Int, T*, Int, const elem::Grid&, Int );
%ignore elem::DistMatrix::DistMatrix( Int, Int, Int, T*, Int, const elem::Grid& );
%ignore elem::DistMatrix::DistMatrix( Int, Int, Int, Int, T*, Int, const elem::Grid& );

%define DCC(T)
%attribute_custom(%arg(elem::DistMatrix<T,elem::distribution_wrapper::CIRC,elem::distribution_wrapper::CIRC>),Int,Root,Root,SetRoot,self_->Root(),self_->SetRoot(val_))
%enddef
DCC(Int)
DCC(float)
DCC(double)
DCC(elem::Complex<float>)
DCC(elem::Complex<double>)

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
%template(DistMatrix_ ## sfx) DistMatrix<F,U,V>;
%extend DistMatrix<F,U,V> {
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
DISTMATRIX(Int,U,V,i_ ## U ## _ ## V)
DISTMATRIX(float,U,V,s_ ## U ## _ ## V)
DISTMATRIX(double,U,V,d_ ## U ## _ ## V)
DISTMATRIX(Complex<float>,U,V,c_ ## U ## _ ## V)
DISTMATRIX(Complex<double>,U,V,z_ ## U ## _ ## V)
%enddef

namespace elem {
DISTMATRIX(Int,MC,MR,i)
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
