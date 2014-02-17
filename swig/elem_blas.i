/*
   Copyright (c) 2009-2014, Jack Poulson
                      2013, Michael C. Grant
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

%module elem_blas

%include "common.swg"
%import "elem.i"

%define OVERLOAD_DIAGONAL_SEQ(X,Y)
%template(Y) X<Int,Int>;
%template(Y) X<float,float>;
%template(Y) X<double,double>;
%template(Y) X<Complex<float>,float>;
%template(Y) X<Complex<float>,Complex<float> >;
%template(Y) X<Complex<double>,double>;
%template(Y) X<Complex<double>,Complex<double> >;
%enddef
%define OVERLOAD_DIAGONAL_DIST(X,Y,U,V)
%template(Y) X<Int,Int,U,V>;
%template(Y) X<float,float,U,V>;
%template(Y) X<double,double,U,V>;
%template(Y) X<Complex<float>,float,U,V>;
%template(Y) X<Complex<float>,Complex<float>,U,V>;
%template(Y) X<Complex<double>,double,U,V>;
%template(Y) X<Complex<double>,Complex<double>,U,V>;
%enddef
%define OVERLOAD_DIAGONAL(X)
OVERLOAD_DIAGONAL_SEQ(X,X);
OVERLOAD_DIAGONAL_DIST(X,X,CIRC,CIRC);
OVERLOAD_DIAGONAL_DIST(X,X,MC,  MR  );
OVERLOAD_DIAGONAL_DIST(X,X,MC,  STAR);
OVERLOAD_DIAGONAL_DIST(X,X,MD,  STAR);
OVERLOAD_DIAGONAL_DIST(X,X,MR,  MC  );
OVERLOAD_DIAGONAL_DIST(X,X,MR,  STAR);
OVERLOAD_DIAGONAL_DIST(X,X,STAR,MC  );
OVERLOAD_DIAGONAL_DIST(X,X,STAR,MD  );
OVERLOAD_DIAGONAL_DIST(X,X,STAR,MR  );
OVERLOAD_DIAGONAL_DIST(X,X,STAR,STAR);
OVERLOAD_DIAGONAL_DIST(X,X,STAR,VC  );
OVERLOAD_DIAGONAL_DIST(X,X,STAR,VR  );
OVERLOAD_DIAGONAL_DIST(X,X,VC,  STAR);
OVERLOAD_DIAGONAL_DIST(X,X,VR,  STAR);
%enddef

// BLAS Miscellaneous
// ==================
 
NO_OVERLOAD(SetLocalSymvBlocksize,int);
NO_OVERLOAD(LocalSymvBlocksize);
NO_OVERLOAD(SetLocalTrrkBlocksize,int);
NO_OVERLOAD(LocalTrrkBlocksize);
NO_OVERLOAD(SetLocalTrr2kBlocksize,int);
NO_OVERLOAD(LocalTrr2kBlocksize);
%include "elemental/blas-like/decl.hpp"

// Level 1
// -------

%include "elemental/blas-like/level1/Adjoint.hpp"
%include "elemental/blas-like/level1/Axpy.hpp"
%include "elemental/blas-like/level1/AxpyTriangle.hpp"
%include "elemental/blas-like/level1/Conjugate.hpp"
%include "elemental/blas-like/level1/Copy.hpp"
%include "elemental/blas-like/level1/DiagonalScale.hpp"
%include "elemental/blas-like/level1/DiagonalSolve.hpp"
%include "elemental/blas-like/level1/Dot.hpp"
%include "elemental/blas-like/level1/Dotu.hpp"
%include "elemental/blas-like/level1/Hadamard.hpp"
%include "elemental/blas-like/level1/HilbertSchmidt.hpp"
%include "elemental/blas-like/level1/MakeHermitian.hpp"
%include "elemental/blas-like/level1/MakeReal.hpp"
%include "elemental/blas-like/level1/MakeSymmetric.hpp"
%include "elemental/blas-like/level1/MakeTrapezoidal.hpp"
%include "elemental/blas-like/level1/MakeTriangular.hpp"
// TODO
/*
%include "elemental/blas-like/level1/Max.hpp"
*/
%include "elemental/blas-like/level1/Nrm2.hpp"
// TODO
/*
%include "elemental/blas-like/level1/QuasiDiagonalScale.hpp"
%include "elemental/blas-like/level1/QuasiDiagonalSolve.hpp"
*/
%include "elemental/blas-like/level1/Scale.hpp"
%include "elemental/blas-like/level1/ScaleTrapezoid.hpp"
%include "elemental/blas-like/level1/SetDiagonal.hpp"
%include "elemental/blas-like/level1/Swap.hpp"
%include "elemental/blas-like/level1/Transpose.hpp"
%include "elemental/blas-like/level1/UpdateDiagonal.hpp"
%include "elemental/blas-like/level1/Zero.hpp"
 
namespace elem {
OVERLOAD02_int(Adjoint)
OVERLOAD_REDIST(Axpy)
OVERLOAD01_int(AxpyTriangle)
OVERLOAD012_int(Conjugate)
OVERLOAD0_int(Copy)
OVERLOAD_REDIST(Copy)
OVERLOAD02_Compat(DiagonalScale)
OVERLOAD02_Compat(DiagonalSolve)
OVERLOAD01_int(Dot)
OVERLOAD01_int(Dotu)
OVERLOAD01(Hadamard)
OVERLOAD01(HilbertSchmidt)
OVERLOAD0_int(MakeHermitian)
OVERLOAD01_int(MakeReal)
OVERLOAD0_int(MakeSymmetric)
OVERLOAD01_int(MakeTrapezoidal)
OVERLOAD01_int(MakeTriangular)
OVERLOAD0(Nrm2)
// TODO: QuasiDiagonal{Scale,Solve}
OVERLOAD01_int(Scale)
OVERLOAD01_int(ScaleTrapezoid)
OVERLOAD_REDIST(Swap)
OVERLOAD_DIAGONAL(SetDiagonal)
OVERLOAD02_int(Transpose)
OVERLOAD_DIAGONAL(UpdateDiagonal)
OVERLOAD01_int(Zero)
};

// Level 2
// -------

// The equivalent of elemental/blas-like/level2.hpp
// TODO
/*
%include "elemental/blas-like/level2/ApplyColumnPivots.hpp"
%include "elemental/blas-like/level2/ApplyRowPivots.hpp"
%include "elemental/blas-like/level2/ApplySymmetricPivots.hpp"
%include "elemental/blas-like/level2/ComposePivots.hpp"
*/
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
%include "elemental/blas-like/level2/Trr.hpp"
%include "elemental/blas-like/level2/Trr2.hpp"
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
OVERLOAD0(Trr)
OVERLOAD0(Trr2)
OVERLOAD0(Trsv)
};

// Level 3
// -------
 
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
%include "elemental/blas-like/level3/Trstrm.hpp"
%include "elemental/blas-like/level3/Trtrmm.hpp"
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
OVERLOAD0(Trstrm)
OVERLOAD0(Trtrmm)
OVERLOAD0(TwoSidedTrmm)
OVERLOAD0(TwoSidedTrsm)
};
