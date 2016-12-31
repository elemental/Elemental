/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_SHIFT_HPP
#define EL_BLAS_SHIFT_HPP

namespace El {

template<typename T,typename S>
void Shift( Matrix<T>& A, S alpha )
{
    EL_DEBUG_CSE
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            A(i,j) += alpha;
}

template<typename T,typename S>
void Shift( AbstractDistMatrix<T>& A, S alpha )
{
    EL_DEBUG_CSE
    Shift( A.Matrix(), alpha );
}

template<typename T,typename S>
void Shift( DistMultiVec<T>& A, S alpha )
{
    EL_DEBUG_CSE
    Shift( A.Matrix(), alpha );
}

} // namespace El

#endif // ifndef EL_BLAS_SHIFT_HPP
