/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>

namespace El {

// TODO(poulson): Think about using a more stable accumulation algorithm?

template<typename Ring>
Ring HilbertSchmidt( const Matrix<Ring>& A, const Matrix<Ring>& B )
{
    EL_DEBUG_CSE
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Matrices must be the same size");
    Ring innerProd(0);
    const Int width = A.Width();
    const Int height = A.Height();
    const Ring* ABuf = A.LockedBuffer();
    const Ring* BBuf = B.LockedBuffer();
    const Int ALDim = A.LDim();
    const Int BLDim = B.LDim();
    if( height == ALDim && height == BLDim )
    {
        innerProd += blas::Dot( height*width, ABuf, 1, BBuf, 1 );
    }
    else
    {
        for( Int j=0; j<width; ++j )
            for( Int i=0; i<height; ++i )
                innerProd += Conj(ABuf[i+j*ALDim])*BBuf[i+j*BLDim];
    }
    return innerProd;
}

template<typename Ring>
Ring HilbertSchmidt
( const ElementalMatrix<Ring>& A, const ElementalMatrix<Ring>& B )
{
    EL_DEBUG_CSE
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Matrices must be the same size");
    AssertSameGrids( A, B );
    // TODO(poulson): Add a general implementation using MatrixReadProxy
    if( A.DistData().colDist != B.DistData().colDist ||
        A.DistData().rowDist != B.DistData().rowDist )
        LogicError("A and B must have the same distribution");
    if( A.ColAlign() != B.ColAlign() || A.RowAlign() != B.RowAlign() )
        LogicError("Matrices must be aligned");

    Ring innerProd;
    if( A.Participating() )
    {
        Ring localInnerProd(0);
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        const Ring* ABuf = A.LockedBuffer();
        const Ring* BBuf = B.LockedBuffer();
        const Int ALDim = A.LDim();
        const Int BLDim = B.LDim();
        if( localHeight == ALDim && localHeight == BLDim )
        {
            localInnerProd +=
              blas::Dot( localHeight*localWidth, ABuf, 1, BBuf, 1 );
        }
        else
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    localInnerProd += Conj(ABuf[iLoc+jLoc*ALDim])*
                                           BBuf[iLoc+jLoc*BLDim];
        }
        innerProd = mpi::AllReduce( localInnerProd, A.DistComm() );
    }
    mpi::Broadcast( innerProd, A.Root(), A.CrossComm() );
    return innerProd;
}

template<typename Ring>
Ring HilbertSchmidt( const DistMultiVec<Ring>& A, const DistMultiVec<Ring>& B )
{
    EL_DEBUG_CSE
    if( !mpi::Congruent( A.Grid().Comm(), B.Grid().Comm() ) )
        LogicError("A and B must be congruent");
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("A and B must have the same dimensions");
    if( A.LocalHeight() != B.LocalHeight() )
        LogicError("A and B must have the same local heights");
    if( A.FirstLocalRow() != B.FirstLocalRow() )
        LogicError("A and B must own the same rows");

    Ring localInnerProd = 0;
    const Int localHeight = A.LocalHeight();
    const Int width = A.Width();
    const Ring* ABuf = A.LockedMatrix().LockedBuffer();
    const Ring* BBuf = B.LockedMatrix().LockedBuffer();
    const Int ALDim = A.LockedMatrix().LDim();
    const Int BLDim = B.LockedMatrix().LDim();
    for( Int j=0; j<width; ++j )
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            localInnerProd += Conj(ABuf[iLoc+j*ALDim])*BBuf[iLoc+j*BLDim];
    return mpi::AllReduce( localInnerProd, A.Grid().Comm() );
}

#define PROTO(Ring) \
  template Ring HilbertSchmidt \
  ( const Matrix<Ring>& A, const Matrix<Ring>& B ); \
  template Ring HilbertSchmidt \
  ( const ElementalMatrix<Ring>& A, const ElementalMatrix<Ring>& B ); \
  template Ring HilbertSchmidt \
  ( const DistMultiVec<Ring>& A, const DistMultiVec<Ring>& B );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
