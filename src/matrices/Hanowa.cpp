/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void Hanowa( Matrix<T>& A, Int n, T mu )
{
    DEBUG_ONLY(CallStackEntry cse("Hanowa"))
    if( n % 2 != 0 )
        LogicError("n must be an even integer");
    A.Resize( n, n );
    const Int m = n/2;
    std::vector<T> d(m);

    for( Int j=0; j<m; ++j )
        d[j] = mu;
    auto ABlock = View( A, 0, 0, m, m );
    Diagonal( ABlock, d );
    ABlock = View( A, m, m, m, m );
    Diagonal( ABlock, d );

    for( Int j=0; j<m; ++j )
        d[j] = -(j+1);
    ABlock = View( A, 0, m, m, m );
    Diagonal( ABlock, d );

    for( Int j=0; j<m; ++j )
        d[j] = j+1;
    ABlock = View( A, m, 0, m, m );
    Diagonal( ABlock, d );
}

template<typename T,Dist U,Dist V>
void Hanowa( DistMatrix<T,U,V>& A, Int n, T mu )
{
    DEBUG_ONLY(CallStackEntry cse("Hanowa"))
    if( n % 2 != 0 )
        LogicError("n must be an even integer");
    A.Resize( n, n );
    const Int m = n/2;
    std::vector<T> d(m);

    for( Int j=0; j<m; ++j )
        d[j] = mu;
    auto ABlock = View( A, 0, 0, m, m );
    Diagonal( ABlock, d );
    ABlock = View( A, m, m, m, m );
    Diagonal( ABlock, d );

    for( Int j=0; j<m; ++j )
        d[j] = -(j+1);
    ABlock = View( A, 0, m, m, m );
    Diagonal( ABlock, d );

    for( Int j=0; j<m; ++j )
        d[j] = j+1;
    ABlock = View( A, m, 0, m, m );
    Diagonal( ABlock, d );
}

/*
template<typename T,Dist U,Dist V>
void Hanowa( BlockDistMatrix<T,U,V>& A, Int n, T mu )
{
    DEBUG_ONLY(CallStackEntry cse("Hanowa"))
    if( n % 2 != 0 )
        LogicError("n must be an even integer");
    A.Resize( n, n );
    const Int m = n/2;
    std::vector<T> d(m);

    for( Int j=0; j<m; ++j )
        d[j] = mu;
    auto ABlock = View( A, 0, 0, m, m );
    Diagonal( ABlock, d );
    ABlock = View( A, m, m, m, m );
    Diagonal( ABlock, d );

    for( Int j=0; j<m; ++j )
        d[j] = -(j+1);
    ABlock = View( A, 0, m, m, m );
    Diagonal( ABlock, d );

    for( Int j=0; j<m; ++j )
        d[j] = j+1;
    ABlock = View( A, m, 0, m, m );
    Diagonal( ABlock, d );
}
*/

#define PROTO_DIST(T,U,V) \
  template void Hanowa( DistMatrix<T,U,V>& A, Int n, T mu );
  //template void Hanowa( BlockDistMatrix<T,U,V>& A, Int n, T mu );

#define PROTO(T) \
  template void Hanowa( Matrix<T>& A, Int n, T mu ); \
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

PROTO(Int)
PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
