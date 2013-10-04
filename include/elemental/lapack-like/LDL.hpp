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
void LocalLDL( DistMatrix<F,STAR,STAR>& A, bool conjugate=false );
} // namespace elem

#include "./LDL/Var3.hpp"
#include "./LDL/Pivoted.hpp"

#include "./LDL/MultiplyAfter.hpp"
#include "./LDL/SolveAfter.hpp"

#include "./LDL/Inertia.hpp"

namespace elem {

template<typename F>
inline void
LocalLDL( DistMatrix<F,STAR,STAR>& A, bool conjugate )
{
#ifndef RELEASE
    CallStackEntry cse("LocalLDL");
#endif
    ldl::Var3( A.Matrix(), conjugate );
}

template<typename F>
inline void
LDLH( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("LDLH");
#endif
    ldl::Var3( A, true );
}

template<typename F>
inline void
LDLH( Matrix<F>& A, Matrix<F>& dSub, Matrix<Int>& p )
{
#ifndef RELEASE
    CallStackEntry cse("LDLH");
#endif
    ldl::Pivoted( A, dSub, p, true );
}

template<typename F>
inline void 
LDLH( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("LDLH");
#endif
    ldl::Var3( A, true );
}

template<typename F>
inline void
LDLH
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& dSub, DistMatrix<Int,VC,STAR>& p )
{
#ifndef RELEASE
    CallStackEntry cse("LDLH");
#endif
    ldl::Pivoted( A, dSub, p, true );
}

template<typename F>
inline void
LDLT( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("LDLT");
#endif
    ldl::Var3( A, false );
}

template<typename F>
inline void
LDLT( Matrix<F>& A, Matrix<F>& dSub, Matrix<Int>& p )
{
#ifndef RELEASE
    CallStackEntry cse("LDLT");
#endif
    ldl::Pivoted( A, dSub, p, false );
}

template<typename F>
inline void 
LDLT( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("LDLT");
#endif
    ldl::Var3( A, false );
}

template<typename F>
inline void
LDLT
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& dSub, DistMatrix<Int,VC,STAR>& p )
{
#ifndef RELEASE
    CallStackEntry cse("LDLT");
#endif
    ldl::Pivoted( A, dSub, p, false );
}

// The nonsingular qualifier is here because Bunch-Parlett (complete pivoting)
// should be used when the matrix might be singular. Not only are more 
// comparisons performed, but blocking is essentially impossible.
template<typename F>
inline elem::Inertia
NonsingularHermitianInertia( UpperOrLower uplo, Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("NonsingularHermitianInertia");
#endif
    if( uplo == UPPER )
        LogicError("This option not yet supported");
    Matrix<Int> p;
    Matrix<F> dSub;
    ldl::Pivoted( A, dSub, p, true ); 
    return ldl::Inertia( A.GetRealPartOfDiagonal(), dSub );
}

template<typename F>
inline elem::Inertia
NonsingularHermitianInertia( UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("NonsingularHermitianInertia");
#endif
    if( uplo == UPPER )
        LogicError("This option not yet supported");
    DistMatrix<Int,VC,STAR> p( A.Grid() );
    DistMatrix<F,MD,STAR> dSub( A.Grid() );
    ldl::Pivoted( A, dSub, p, true ); 
    return ldl::Inertia( A.GetRealPartOfDiagonal(), dSub );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_LDL_HPP
