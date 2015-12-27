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
void SetSubmatrix
(       Matrix<T>& A, 
  const vector<Int>& I, 
  const vector<Int>& J, 
  const Matrix<T>& ASub )
{
    DEBUG_ONLY(CSE cse("SetSubmatrix"))
    const Int m = I.size();
    const Int n = J.size();

    // Fill in our locally-owned entries
    for( Int jSub=0; jSub<n; ++jSub )
    {
        const Int j = J[jSub];
        for( Int iSub=0; iSub<m; ++iSub )
        {
            const Int i = I[iSub];
            A.Set( i, j, ASub.Get(iSub,jSub) );
        }
    }
}

template<typename T>
void SetSubmatrix
(       AbstractDistMatrix<T>& A, 
  const vector<Int>& I, 
  const vector<Int>& J, 
  const AbstractDistMatrix<T>& ASub )
{
    DEBUG_ONLY(CSE cse("SetSubmatrix"))
    // Set the appropriate portion of A to zero before updating
    for( auto& i : I )
        if( A.IsLocalRow(i) )
            for( auto& j : J )
                if( A.IsLocalCol(j) )
                    A.Set( i, j, 0 );
    UpdateSubmatrix( A, I, J, T(1), ASub );
}

// TODO: DistMultiVec version similar to GetSubmatrix implementation

#define PROTO(T) \
  template void SetSubmatrix \
  (       Matrix<T>& A, \
    const vector<Int>& I, \
    const vector<Int>& J, \
    const Matrix<T>& ASub ); \
  template void SetSubmatrix \
  (       AbstractDistMatrix<T>& A, \
    const vector<Int>& I, \
    const vector<Int>& J, \
    const AbstractDistMatrix<T>& ASub );

#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
