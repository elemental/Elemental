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
( const Matrix<T>& A, Matrix<S>& d, function<S(T)> func, Int offset )
{
    DEBUG_ONLY(CSE cse("GetMappedDiagonal"))
    const Int diagLength = A.DiagonalLength(offset);
    d.Resize( diagLength, 1 );

    const Int iStart = Max(-offset,0);
    const Int jStart = Max( offset,0);
    S* dBuf = d.Buffer();
    const T* ABuf = A.LockedBuffer();
    const Int ldim = A.LDim();
    EL_PARALLEL_FOR
    for( Int k=0; k<diagLength; ++k )
    {
        const Int i = iStart + k;
        const Int j = jStart + k;
        dBuf[k] = func(ABuf[i+j*ldim]);
    }
}

template<typename T,typename S,Dist U,Dist V>
void GetMappedDiagonal
( const DistMatrix<T,U,V>& A, ElementalMatrix<S>& dPre, 
  function<S(T)> func, Int offset )
{ 
    DEBUG_ONLY(
      CSE cse("GetMappedDiagonal");
      AssertSameGrids( A, dPre );
    )
    ElementalProxyCtrl ctrl;
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
        const T* ABuf = A.LockedBuffer();
        const Int ldim = A.LDim();
        EL_PARALLEL_FOR
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLoc = iLocStart + k*iLocStride;
            const Int jLoc = jLocStart + k*jLocStride;
            dBuf[k] = func(ABuf[iLoc+jLoc*ldim]);
        }
    }
}

template<typename T,typename S>
void GetMappedDiagonal
( const SparseMatrix<T>& A, Matrix<S>& d, function<S(T)> func, Int offset )
{
    DEBUG_ONLY(CSE cse("GetMappedDiagonal"))
    const Int m = A.Height();
    const Int n = A.Width();
    const T* valBuf = A.LockedValueBuffer();
    const Int* colBuf = A.LockedTargetBuffer();

    const Int iStart = Max(-offset,0);
    const Int jStart = Max( offset,0);

    const Int diagLength = El::DiagonalLength(m,n,offset);
    Zeros( d, diagLength, 1 );
    S* dBuf = d.Buffer();

    for( Int k=0; k<diagLength; ++k )
    {
        const Int i = iStart + k;
        const Int j = jStart + k;
        const Int thisOff = A.RowOffset(i);
        const Int nextOff = A.RowOffset(i+1);
        auto it = std::lower_bound( colBuf+thisOff, colBuf+nextOff, j );
        if( *it == j )
        {
            const Int e = it-colBuf;
            dBuf[Min(i,j)] = func(valBuf[e]);
        }
        else
            dBuf[Min(i,j)] = func(0);
    }
}

template<typename T,typename S>
void GetMappedDiagonal
( const DistSparseMatrix<T>& A, DistMultiVec<S>& d, 
  function<S(T)> func, Int offset )
{
    DEBUG_ONLY(CSE cse("GetMappedDiagonal"))
    const Int m = A.Height();
    const Int n = A.Width();
    const T* valBuf = A.LockedValueBuffer();
    const Int* colBuf = A.LockedTargetBuffer();

    if( m != n )
        LogicError("DistSparseMatrix GetMappedDiagonal assumes square matrix");
    if( offset != 0 )
        LogicError("DistSparseMatrix GetMappedDiagonal assumes offset=0");

    d.SetComm( A.Comm() );
    Ones( d, El::DiagonalLength(m,n,offset), 1 );
    S* dBuf = d.Matrix().Buffer();
    const Int dLocalHeight = d.LocalHeight();
    for( Int iLoc=0; iLoc<dLocalHeight; ++iLoc )
    {
        const Int i = d.GlobalRow(iLoc);
        const Int thisOff = A.RowOffset(iLoc);
        const Int nextOff = A.RowOffset(iLoc+1);
        auto it = std::lower_bound( colBuf+thisOff, colBuf+nextOff, i );
        if( *it == i )
        {
            const Int e = it-colBuf;
            dBuf[iLoc] = func(valBuf[e]);
        }
        else
            dBuf[iLoc] = func(0);
    }
}

#define PROTO_DIST_TYPES(S,T,U,V) \
  template void GetMappedDiagonal \
  ( const DistMatrix<T,U,V>& A, ElementalMatrix<S>& d, \
    function<S(T)> func, Int offset );

#define PROTO_TYPES(S,T) \
  template void GetMappedDiagonal \
  ( const Matrix<T>& A, Matrix<S>& B, function<S(T)> func, Int offset ); \
  template void GetMappedDiagonal \
  ( const SparseMatrix<T>& A, Matrix<S>& B, function<S(T)> func, Int offset ); \
  template void GetMappedDiagonal \
  ( const DistSparseMatrix<T>& A, DistMultiVec<S>& B, \
    function<S(T)> func, Int offset ); \
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

#ifdef EL_HAVE_QUAD

#define PROTO(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,float) \
  PROTO_TYPES(T,double) \
  PROTO_TYPES(T,Quad) \
  PROTO_TYPES(T,Complex<float>) \
  PROTO_TYPES(T,Complex<double>) \
  PROTO_TYPES(T,Complex<Quad>)

#else

#define PROTO(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,float) \
  PROTO_TYPES(T,double) \
  PROTO_TYPES(T,Complex<float>) \
  PROTO_TYPES(T,Complex<double>)

#endif

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
