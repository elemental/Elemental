/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MAKEHERMITIAN_HPP
#define ELEM_MAKEHERMITIAN_HPP

#include ELEM_MAKESYMMETRIC_INC

namespace elem {

template<typename T>
inline void
MakeHermitian( UpperOrLower uplo, Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeHermitian"))
    MakeSymmetric( uplo, A, true );
}

template<typename T>
inline void
MakeHermitian( UpperOrLower uplo, DistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeHermitian"))
    MakeSymmetric( uplo, A, true );
}

} // namespace elem

#endif // ifndef ELEM_MAKEHERMITIAN_HPP
