/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T,typename S>
void GetMappedDiagonal
( const Matrix<T>& A, Matrix<S>& d, std::function<S(T)> func, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GetMappedDiagonal"))
    const Int diagLength = A.DiagonalLength(offset);
    d.Resize( diagLength, 1 );

    const Int iStart = Max(-offset,0);
    const Int jStart = Max( offset,0);
    S* dBuf = d.Buffer();
    const T* buffer = A.LockedBuffer();
    const Int ldim = A.LDim();
    EL_PARALLEL_FOR
    for( Int k=0; k<diagLength; ++k )
    {
        const Int i = iStart + k;
        const Int j = jStart + k;
        dBuf[k] = func(buffer[i+j*ldim]);
    }
}

// TODO: SparseMatrix implementation

template<typename T,typename S,Dist U,Dist V>
void GetMappedDiagonal
( const DistMatrix<T,U,V>& A, AbstractDistMatrix<S>& dPre, 
  std::function<S(T)> func, Int offset )
{ 
    DEBUG_ONLY(
      CallStackEntry cse("GetMappedDiagonal");
      AssertSameGrids( A, dPre );
    )
    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = A.DiagonalAlign(offset);
    ctrl.rootConstrain = true;
    ctrl.root = A.DiagonalRoot(offset);
    auto dPtr = WriteProxy<S,DiagCol<U,V>(),DiagRow<U,V>()>(&dPre,ctrl);
    auto& d = *dPtr;

    d.Resize( A.DiagonalLength(offset), 1 );
    if( d.Participating() )
    {
        const Int diagShift = d.ColShift();
        const Int iStart = diagShift + Max(-offset,0);
        const Int jStart = diagShift + Max( offset,0);

        const Int colStride = A.ColStride();
        const Int rowStride = A.RowStride();
        const Int iLocStart = (iStart-A.ColShift()) / colStride;
        const Int jLocStart = (jStart-A.RowShift()) / rowStride;
        const Int iLocStride = d.ColStride() / colStride;
        const Int jLocStride = d.ColStride() / rowStride;

        const Int localDiagLength = d.LocalHeight();
        S* dBuf = d.Buffer();
        const T* buffer = A.LockedBuffer();
        const Int ldim = A.LDim();
        EL_PARALLEL_FOR
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLoc = iLocStart + k*iLocStride;
            const Int jLoc = jLocStart + k*jLocStride;
            dBuf[k] = func(buffer[iLoc+jLoc*ldim]);
        }
    }
}

// TODO: DistSparseMatrix implementation

#define PROTO_DIST_TYPES(S,T,U,V) \
  template void GetMappedDiagonal \
  ( const DistMatrix<T,U,V>& A, AbstractDistMatrix<S>& d, \
    std::function<S(T)> func, Int offset );

#define PROTO_TYPES(S,T) \
  template void GetMappedDiagonal \
  ( const Matrix<T>& A, Matrix<S>& B, std::function<S(T)> func, Int offset ); \
  PROTO_DIST_TYPES(S,T,CIRC,CIRC) \
  PROTO_DIST_TYPES(S,T,MC,  MR  ) \
  PROTO_DIST_TYPES(S,T,MC,  STAR) \
  PROTO_DIST_TYPES(S,T,MD,  STAR) \
  PROTO_DIST_TYPES(S,T,MR,  MC  ) \
  PROTO_DIST_TYPES(S,T,MR,  STAR) \
  PROTO_DIST_TYPES(S,T,STAR,MC  ) \
  PROTO_DIST_TYPES(S,T,STAR,MD  ) \
  PROTO_DIST_TYPES(S,T,STAR,MR  ) \
  PROTO_DIST_TYPES(S,T,STAR,STAR) \
  PROTO_DIST_TYPES(S,T,STAR,VC  ) \
  PROTO_DIST_TYPES(S,T,STAR,VR  ) \
  PROTO_DIST_TYPES(S,T,VC,  STAR) \
  PROTO_DIST_TYPES(S,T,VR,  STAR)

#define PROTO(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,float) \
  PROTO_TYPES(T,double) \
  PROTO_TYPES(T,Complex<float>) \
  PROTO_TYPES(T,Complex<double>)

#include "El/macros/Instantiate.h"

} // namespace El
