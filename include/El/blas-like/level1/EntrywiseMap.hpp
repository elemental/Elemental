/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_ENTRYWISEMAP_HPP
#define EL_ENTRYWISEMAP_HPP

namespace El {

template<typename T,class Function>
inline void
EntrywiseMap( Matrix<T>& A, Function func )
{
    DEBUG_ONLY(CallStackEntry cse("EntrywiseMap"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, func(A.Get(i,j)) );
}

template<typename T,class Function>
inline void
EntrywiseMap( AbstractDistMatrix<T>& A, Function func )
{ EntrywiseMap( A.Matrix(), func ); }

template<typename T,class Function>
inline void
EntrywiseMap( AbstractBlockDistMatrix<T>& A, Function func )
{ EntrywiseMap( A.Matrix(), func ); }

} // namespace El

#endif // ifndef EL_ENTRYWISEMAP_HPP
