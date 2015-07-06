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
void Fill( Matrix<T>& A, T alpha )
{
    DEBUG_ONLY(CSE cse("Fill"))
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            A.Set( i, j, alpha );
}

template<typename T>
void Fill( AbstractDistMatrix<T>& A, T alpha )
{
    DEBUG_ONLY(CSE cse("Fill"))
    Fill( A.Matrix(), alpha );
}

template<typename T>
void Fill( AbstractBlockDistMatrix<T>& A, T alpha )
{
    DEBUG_ONLY(CSE cse("Fill"))
    Fill( A.Matrix(), alpha );
}

template<typename T>
void Fill( DistMultiVec<T>& A, T alpha )
{
    DEBUG_ONLY(CSE cse("Fill"))
    Fill( A.Matrix(), alpha );
}

template<typename T>
void Fill( SparseMatrix<T>& A, T alpha )
{
    DEBUG_ONLY(CSE cse("Fill"))
    const Int m = A.Height();
    const Int n = A.Width();
    Zeros( A, m, n );
    if( alpha != T(0) )
    {
        A.Reserve( m*n ); 
        for( Int i=0; i<m; ++i )
            for( Int j=0; j<n; ++j )
                A.QueueUpdate( i, j, alpha );
        A.ProcessQueues();
    }
}

template<typename T>
void Fill( DistSparseMatrix<T>& A, T alpha )
{
    DEBUG_ONLY(CSE cse("Fill"))
    const Int m = A.Height();
    const Int n = A.Width();
    Zeros( A, m, n );
    if( alpha != T(0) ) 
    {
        const Int localHeight = A.LocalHeight();
        A.Reserve( localHeight*n );
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            for( Int j=0; j<n; ++j ) 
                A.QueueLocalUpdate( iLoc, j, alpha );
        A.ProcessLocalQueues();
    }
}

#define PROTO(T) \
  template void Fill( Matrix<T>& A, T alpha ); \
  template void Fill( AbstractDistMatrix<T>& A, T alpha ); \
  template void Fill( AbstractBlockDistMatrix<T>& A, T alpha ); \
  template void Fill( DistMultiVec<T>& A, T alpha ); \
  template void Fill( SparseMatrix<T>& A, T alpha ); \
  template void Fill( DistSparseMatrix<T>& A, T alpha );

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
