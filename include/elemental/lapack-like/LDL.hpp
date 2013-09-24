/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_LDL_HPP
#define ELEM_LAPACK_LDL_HPP

namespace elem {
template<typename F>
void LocalLDL( Orientation orientation, DistMatrix<F,STAR,STAR>& A );
} // namespace elem

#include "./LDL/Pivoted.hpp"
#include "./LDL/Var3.hpp"

namespace elem {

template<typename F>
inline void
LocalLDL( Orientation orientation, DistMatrix<F,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("LocalLDL");
#endif
    ldl::Var3( orientation, A.Matrix() );
}

template<typename F>
inline void
LDLH( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("LDLH");
#endif
    ldl::Var3( ADJOINT, A );
}

template<typename F>
inline void
LDLH( Matrix<F>& A, Matrix<Int>& p )
{
#ifndef RELEASE
    CallStackEntry cse("LDLH");
#endif
    ldl::Pivoted( ADJOINT, A, p );
}

template<typename F>
inline void 
LDLH( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("LDLH");
#endif
    ldl::Var3( ADJOINT, A );
}

template<typename F>
inline void
LDLT( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("LDLT");
#endif
    ldl::Var3( TRANSPOSE, A );
}

template<typename F>
inline void
LDLT( Matrix<F>& A, Matrix<Int>& p )
{
#ifndef RELEASE
    CallStackEntry cse("LDLT");
#endif
    ldl::Pivoted( TRANSPOSE, A, p );
}

template<typename F>
inline void 
LDLT( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("LDLT");
#endif
    ldl::Var3( TRANSPOSE, A );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_LDL_HPP
