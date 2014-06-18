/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_PEI_HPP
#define EL_PEI_HPP

namespace El {

template<typename T> 
inline void
Pei( Matrix<T>& P, Int n, T alpha )
{
    DEBUG_ONLY(CallStackEntry cse("Pei"))
    Ones( P, n, n );
    UpdateDiagonal( P, alpha );
}

template<typename T>
inline void
Pei( AbstractDistMatrix<T>& P, Int n, T alpha )
{
    DEBUG_ONLY(CallStackEntry cse("Pei"))
    Ones( P, n, n );
    UpdateDiagonal( P, alpha );
}

template<typename T>
inline void
Pei( AbstractBlockDistMatrix<T>& P, Int n, T alpha )
{
    DEBUG_ONLY(CallStackEntry cse("Pei"))
    Ones( P, n, n );
    UpdateDiagonal( P, alpha );
}

} // namespace El

#endif // ifndef EL_PEI_HPP
