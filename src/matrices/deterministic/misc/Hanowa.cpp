/*
   Copyright (c) 2009-2016, Jack Poulson
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
    DEBUG_ONLY(CSE cse("Hanowa"))
    if( n % 2 != 0 )
        LogicError("n must be an even integer");
    A.Resize( n, n );
    const Int m = n/2;
    vector<T> d(m);

    for( Int j=0; j<m; ++j )
        d[j] = mu;
    auto ABlock = A( IR(0,m), IR(0,m) );
    Diagonal( ABlock, d );
    ABlock = A( IR(m,2*m), IR(m,2*m) );
    Diagonal( ABlock, d );

    for( Int j=0; j<m; ++j )
        d[j] = -(j+1);
    ABlock = A( IR(0,m), IR(m,2*m) );
    Diagonal( ABlock, d );

    for( Int j=0; j<m; ++j )
        d[j] = j+1;
    ABlock = A( IR(m,2*m), IR(0,m) );
    Diagonal( ABlock, d );
}

template<typename T>
void Hanowa( ElementalMatrix<T>& A, Int n, T mu )
{
    DEBUG_ONLY(CSE cse("Hanowa"))
    if( n % 2 != 0 )
        LogicError("n must be an even integer");
    A.Resize( n, n );
    const Int m = n/2;
    vector<T> d(m);

    for( Int j=0; j<m; ++j )
        d[j] = mu;
    unique_ptr<ElementalMatrix<T>> ABlock( A.Construct(A.Grid(),A.Root()) );
    View( *ABlock, A, IR(0,m), IR(0,m) );
    Diagonal( *ABlock, d );
    View( *ABlock, A, IR(m,2*m), IR(m,2*m) );
    Diagonal( *ABlock, d );

    for( Int j=0; j<m; ++j )
        d[j] = -(j+1);
    View( *ABlock, A, IR(0,m), IR(m,2*m) );
    Diagonal( *ABlock, d );

    for( Int j=0; j<m; ++j )
        d[j] = j+1;
    View( *ABlock, A, IR(m,2*m), IR(0,m) );
    Diagonal( *ABlock, d );
}

#define PROTO(T) \
  template void Hanowa( Matrix<T>& A, Int n, T mu ); \
  template void Hanowa( ElementalMatrix<T>& A, Int n, T mu );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
