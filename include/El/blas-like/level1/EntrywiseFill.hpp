/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_ENTRYWISEFILL_HPP
#define EL_ENTRYWISEFILL_HPP

namespace El {

template<typename T,class Function>
inline void
EntrywiseFill( Matrix<T>& A, Function func )
{
    DEBUG_ONLY(CallStackEntry cse("EntrywiseFill"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, func() );
}

template<typename T,class Function>
inline void
EntrywiseFill( AbstractDistMatrix<T>& A, Function func )
{ EntrywiseFill( A.Matrix(), func ); }

template<typename T,class Function>
inline void
EntrywiseFill( AbstractBlockDistMatrix<T>& A, Function func )
{ EntrywiseFill( A.Matrix(), func ); }

} // namespace El

#endif // ifndef EL_ENTRYWISEFILL_HPP
