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
void ApplyColPivots
(       Matrix<T>& A,
  const Matrix<Int>& pivots,
        Int offset )
{
    DEBUG_ONLY(
      CSE cse("ApplyColPivots");
      if( pivots.Width() != 1 )
          LogicError("pivots must be a column vector");
      if( pivots.Height() > A.Width() )
          LogicError("pivots cannot be longer than width of A");
    )
    const Int height = A.Height();
    const Int width = A.Width();
    if( height == 0 || width == 0 )
        return;

    const Int numPivots = pivots.Height();
    const Int* pivotBuf = pivots.LockedBuffer();
    for( Int j=0; j<numPivots; ++j )
    {
        const Int k = pivotBuf[j] - offset;
        ColSwap( A, j, k );
    }
}

template<typename T>
void ApplyInverseColPivots
(       Matrix<T>& A,
  const Matrix<Int>& pivots,
        Int offset )
{
    DEBUG_ONLY(
      CSE cse("ApplyInverseColPivots");
      if( pivots.Width() != 1 )
          LogicError("pivots must be a column vector");
      if( pivots.Height() > A.Width() )
          LogicError("pivots cannot be larger than width of A");
    )
    const Int height = A.Height();
    const Int width = A.Width();
    if( height == 0 || width == 0 )
        return;

    const Int numPivots = pivots.Height();
    const Int* pivotBuf = pivots.LockedBuffer();
    for( Int j=numPivots-1; j>=0; --j )
    {
        const Int k = pivotBuf[j] - offset;
        ColSwap( A, j, k );
    }
}

template<typename T>
void ApplyColPivots
(       ElementalMatrix<T>& A,
  const ElementalMatrix<Int>& pivotsPre,
        Int offset )
{
    DEBUG_ONLY(CSE cse("ApplyColPivots"))
    auto pivotsPtr = ReadProxy<Int,STAR,STAR>( &pivotsPre );
    auto& pivots = *pivotsPtr;

    const Int numPivots = pivots.Height();
    const Int* pivotsBuf = pivots.LockedBuffer();
    for( Int j=0; j<numPivots; ++j )
    {
        const Int k = pivotsBuf[j] - offset;
        ColSwap( A, j, k );
    }
}

template<typename T>
void ApplyInverseColPivots
(       ElementalMatrix<T>& A,
  const ElementalMatrix<Int>& pivotsPre,
        Int offset )
{
    DEBUG_ONLY(CSE cse("ApplyInverseColPivots"))
    auto pivotsPtr = ReadProxy<Int,STAR,STAR>( &pivotsPre );
    auto& pivots = *pivotsPtr;

    const Int numPivots = pivots.Height();
    const Int* pivotsBuf = pivots.LockedBuffer();
    for( Int j=numPivots-1; j>=0; --j )
    {
        const Int k = pivotsBuf[j] - offset;
        ColSwap( A, j, k );
    }
}

#define PROTO(T) \
  template void ApplyColPivots \
  (       Matrix<T>& A, \
    const Matrix<Int>& pivots, \
    Int offset ); \
  template void ApplyColPivots \
  (       ElementalMatrix<T>& A, \
    const ElementalMatrix<Int>& pivots, \
          Int offset ); \
  template void ApplyInverseColPivots \
  (       Matrix<T>& A, \
    const Matrix<Int>& pivots, \
          Int offset ); \
  template void ApplyInverseColPivots \
  (       ElementalMatrix<T>& A, \
    const ElementalMatrix<Int>& pivots, \
          Int offset );

#include "El/macros/Instantiate.h"

} // namespace El
