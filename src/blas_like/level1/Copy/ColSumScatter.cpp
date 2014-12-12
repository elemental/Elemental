/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace copy {

template<typename T,Dist U,Dist V>
void ColSumScatter
( const DistMatrix<T,Collect<U>(),V>& A,
        DistMatrix<T,        U,   V>& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::ColSumScatter"))
    AssertSameGrids( A, B );
    B.AlignRowsAndResize
    ( A.RowAlign(), A.Height(), A.Width(), false, false );
    // NOTE: This will be *slightly* slower than necessary due to the result
    //       of the MPI operations being added rather than just copied
    Zeros( B.Matrix(), B.LocalHeight(), B.LocalWidth() );
    axpy::ColSumScatter( T(1), A, B );
}

#define PROTO_DIST(T,U,V) \
  template void ColSumScatter \
  ( const DistMatrix<T,Collect<U>(),V>& A, \
          DistMatrix<T,        U,   V>& B );

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
