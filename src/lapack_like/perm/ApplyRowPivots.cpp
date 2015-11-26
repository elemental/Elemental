/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void ApplyRowPivots
(       Matrix<T>& A,
  const Matrix<Int>& pivots,
        Int offset )
{
    DEBUG_ONLY(
      CSE cse("ApplyRowPivots");
      if( pivots.Width() != 1 )
          LogicError("p must be a column vector");
      if( pivots.Height() > A.Height() )
          LogicError("p cannot be larger than height of A");
    )
    const Int height = A.Height();
    const Int width = A.Width();
    if( height == 0 || width == 0 )
        return;

    const Int numPivots = pivots.Height();
    const Int* pivotsBuf = pivots.LockedBuffer();
    for( Int i=0; i<numPivots; ++i )
    {
        const Int k = pivotsBuf[i] - offset;
        RowSwap( A, i, k );
    }
}

template<typename T>
void ApplyInverseRowPivots
(       Matrix<T>& A,
  const Matrix<Int>& pivots,
        Int offset )
{
    DEBUG_ONLY(
      CSE cse("ApplyInverseRowPivots");
      if( pivots.Width() != 1 )
          LogicError("pivots must be a column vector");
      if( pivots.Height() > A.Height() )
          LogicError("pivots cannot be larger than height of A");
    )
    const Int height = A.Height();
    const Int width = A.Width();
    if( height == 0 || width == 0 )
        return;

    const Int numPivots = pivots.Height();
    const Int* pivotsBuf = pivots.LockedBuffer();
    for( Int i=numPivots-1; i>=0; --i )
    {
        const Int k = pivotsBuf[i] - offset;
        RowSwap( A, i, k );
    }
}

template<typename T>
void ApplyRowPivots
(       ElementalMatrix<T>& A,
  const ElementalMatrix<Int>& pivotsPre,
        Int offset )
{
    DEBUG_ONLY(CSE cse("ApplyRowPivots"))

    DistMatrixReadProxy<Int,Int,STAR,STAR> pivotsProx( pivotsPre );
    auto& pivots = pivotsProx.GetLocked();

    const Int numPivots = pivots.Height();
    const Int* pivotsBuf = pivots.LockedBuffer();
    for( Int i=0; i<numPivots; ++i )
    {
        const Int k = pivotsBuf[i] - offset;
        RowSwap( A, i, k );
    }
}

template<typename T>
void ApplyInverseRowPivots
(       ElementalMatrix<T>& A,
  const ElementalMatrix<Int>& pivotsPre,
        Int offset )
{
    DEBUG_ONLY(CSE cse("ApplyInverseRowPivots"))

    DistMatrixReadProxy<Int,Int,STAR,STAR> pivotsProx( pivotsPre );
    auto& pivots = pivotsProx.GetLocked();

    const Int numPivots = pivots.Height();
    const Int* pivotsBuf = pivots.LockedBuffer();
    for( Int i=numPivots-1; i>=0; --i )
    {
        const Int k = pivotsBuf[i] - offset;
        RowSwap( A, i, k );
    }
}

#define PROTO(T) \
  template void ApplyRowPivots \
  (       Matrix<T>& A, \
    const Matrix<Int>& pivots, \
          Int offset ); \
  template void ApplyRowPivots \
  (       ElementalMatrix<T>& A, \
    const ElementalMatrix<Int>& pivots, \
          Int offset ); \
  template void ApplyInverseRowPivots \
  (       Matrix<T>& A, \
    const Matrix<Int>& pivots, \
          Int offset ); \
  template void ApplyInverseRowPivots \
  (       ElementalMatrix<T>& A, \
    const ElementalMatrix<Int>& pivots, \
          Int offset );

#include "El/macros/Instantiate.h"

} // namespace El
