/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include EL_ZEROS_INC

namespace El {

void ExplicitPermutation( const Matrix<Int>& perm, Matrix<Int>& P )
{
    DEBUG_ONLY(
        CallStackEntry cse("ExplicitPermutation");
        if( perm.Width() != 1 )
            LogicError("perm must be a column vector");
    )
    const Int n = perm.Height();
    Zeros( P, n, n );
    for( Int i=0; i<n; ++i )
        P.Set( i, perm.Get(i,0), 1 );
}

template<Dist UPerm,Dist U,Dist V>
void ExplicitPermutation
( const DistMatrix<Int,UPerm,STAR>& perm, DistMatrix<Int,U,V>& P )
{
    DEBUG_ONLY(
        CallStackEntry cse("ExplicitPermutation");
        if( perm.Width() != 1 )
            LogicError("perm must be a column vector");
    )
    const Int n = perm.Height();
    Zeros( P, n, n );
    for( Int i=0; i<n; ++i )
        P.Set( i, perm.Get(i,0), 1 );
}

#define PROTO_DIST(U,V) \
  template void ExplicitPermutation \
  ( const DistMatrix<Int,MC,STAR>& perm, DistMatrix<Int,U,V>& P ); \
  template void ExplicitPermutation \
  ( const DistMatrix<Int,MD,STAR>& perm, DistMatrix<Int,U,V>& P ); \
  template void ExplicitPermutation \
  ( const DistMatrix<Int,MR,STAR>& perm, DistMatrix<Int,U,V>& P ); \
  template void ExplicitPermutation \
  ( const DistMatrix<Int,STAR,STAR>& perm, DistMatrix<Int,U,V>& P ); \
  template void ExplicitPermutation \
  ( const DistMatrix<Int,VC,STAR>& perm, DistMatrix<Int,U,V>& P ); \
  template void ExplicitPermutation \
  ( const DistMatrix<Int,VR,STAR>& perm, DistMatrix<Int,U,V>& P );

PROTO_DIST(CIRC,CIRC) 
PROTO_DIST(MC,  MR  ) 
PROTO_DIST(MC,  STAR) 
PROTO_DIST(MD,  STAR) 
PROTO_DIST(MR,  MC  ) 
PROTO_DIST(MR,  STAR) 
PROTO_DIST(STAR,MC  ) 
PROTO_DIST(STAR,MD  ) 
PROTO_DIST(STAR,MR  ) 
PROTO_DIST(STAR,STAR) 
PROTO_DIST(STAR,VC  ) 
PROTO_DIST(STAR,VR  ) 
PROTO_DIST(VR,  STAR) 
PROTO_DIST(VC,  STAR)

} // namespace El
