/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace copy {

template<typename T,Dist U,Dist V>
void PartialRowSumScatter
( const DistMatrix<T,U,Partial<V>()>& A,
        DistMatrix<T,U,        V   >& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::PartialRowSumScatter"))
    AssertSameGrids( A, B );

    B.AlignAndResize
    ( A.ColAlign(), A.RowAlign(), A.Height(), A.Width(), false, false );
    Zeros( B.Matrix(), B.LocalHeight(), B.LocalWidth() );
    axpy::PartialRowSumScatter( T(1), A, B );
}

template<typename T,Dist U,Dist V>
void PartialRowSumScatter
( const BlockDistMatrix<T,U,Partial<V>()>& A,
        BlockDistMatrix<T,U,        V   >& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::PartialRowSumScatter"))
    AssertSameGrids( A, B );

    B.AlignAndResize
    ( A.BlockHeight(), A.BlockWidth(), 
      A.ColAlign(), A.RowAlign(), A.ColCut(), A.RowCut(), 
      A.Height(), A.Width(), false, false );
    Zeros( B.Matrix(), B.LocalHeight(), B.LocalWidth() );
    axpy::PartialRowSumScatter( T(1), A, B );
}

#define PROTO_DIST(T,U,V) \
  template void PartialRowSumScatter \
  ( const DistMatrix<T,U,Partial<V>()>& A, \
          DistMatrix<T,U,        V   >& B ); \
  template void PartialRowSumScatter \
  ( const BlockDistMatrix<T,U,Partial<V>()>& A, \
          BlockDistMatrix<T,U,        V   >& B );

#define PROTO(T) \
  PROTO_DIST(T,CIRC,CIRC) \
  PROTO_DIST(T,MC,  MR  ) \
  PROTO_DIST(T,MC,  STAR) \
  PROTO_DIST(T,MD,  STAR) \
  PROTO_DIST(T,MR,  MC  ) \
  PROTO_DIST(T,MR,  STAR) \
  PROTO_DIST(T,STAR,MC  ) \
  PROTO_DIST(T,STAR,MD  ) \
  PROTO_DIST(T,STAR,MR  ) \
  PROTO_DIST(T,STAR,STAR) \
  PROTO_DIST(T,STAR,VC  ) \
  PROTO_DIST(T,STAR,VR  ) \
  PROTO_DIST(T,VC,  STAR) \
  PROTO_DIST(T,VR,  STAR) 

#include "El/macros/Instantiate.h"

} // namespace copy
} // namespace El
