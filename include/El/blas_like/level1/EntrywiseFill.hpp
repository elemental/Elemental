/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLAS_ENTRYWISEFILL_HPP
#define EL_BLAS_ENTRYWISEFILL_HPP

namespace El {

template<typename T>
void EntrywiseFill( Matrix<T>& A, function<T(void)> func )
{
    DEBUG_ONLY(CSE cse("EntrywiseFill"))
    const Int m = A.Height();
    const Int n = A.Width();
    T* ABuf = A.Buffer();
    const Int ALDim = A.LDim();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            ABuf[i+j*ALDim] = func();
}

template<typename T>
void EntrywiseFill( AbstractDistMatrix<T>& A, function<T(void)> func )
{ EntrywiseFill( A.Matrix(), func ); }

template<typename T>
void EntrywiseFill( DistMultiVec<T>& A, function<T(void)> func )
{ EntrywiseFill( A.Matrix(), func ); }

} // namespace El

#endif // ifndef EL_BLAS_ENTRYWISEFILL_HPP
