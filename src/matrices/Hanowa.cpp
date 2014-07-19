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

template<typename T>
void Hanowa( AbstractDistMatrix<T>& A, Int n, T mu )
{
    DEBUG_ONLY(CallStackEntry cse("Hanowa"))
    if( n % 2 != 0 )
        LogicError("n must be an even integer");
    A.Resize( n, n );
    const Int m = n/2;
    std::vector<T> d(m);

    for( Int j=0; j<m; ++j )
        d[j] = mu;
    std::unique_ptr<AbstractDistMatrix<T>> ABlock( A.Construct(A.Grid()) );
    View( *ABlock, A, 0, 0, m, m );
    Diagonal( *ABlock, d );
    View( *ABlock, A, m, m, m, m );
    Diagonal( *ABlock, d );

    for( Int j=0; j<m; ++j )
        d[j] = -(j+1);
    View( *ABlock, A, 0, m, m, m );
    Diagonal( *ABlock, d );

    for( Int j=0; j<m; ++j )
        d[j] = j+1;
    View( *ABlock, A, m, 0, m, m );
    Diagonal( *ABlock, d );
}

// TODO: AbstractBlockDistMatrix version

#define PROTO(T) \
  template void Hanowa( Matrix<T>& A, Int n, T mu ); \
  template void Hanowa( AbstractDistMatrix<T>& A, Int n, T mu );

#include "El/macros/Instantiate.h"

} // namespace El
