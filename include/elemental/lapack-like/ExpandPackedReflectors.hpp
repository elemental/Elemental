/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_EXPANDPACKEDREFLECTORS_HPP
#define ELEM_LAPACK_EXPANDPACKEDREFLECTORS_HPP

#include "elemental/lapack-like/ApplyPackedReflectors/Util.hpp"

#include "./ExpandPackedReflectors/LV.hpp"

namespace elem {

template<typename F> 
inline void
ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation,
  int offset, Matrix<F>& H, const Matrix<F>& t )
{
#ifndef RELEASE
    CallStackEntry entry("ExpandPackedReflectors");
#endif
    if( uplo == LOWER && dir == VERTICAL )
        expand_packed_reflectors::LV( conjugation, offset, H, t );
    else
        throw std::logic_error("This option is not yet supported");
}

template<typename F> 
inline void
ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation,
  int offset, DistMatrix<F>& H, const DistMatrix<F,MD,STAR>& t )
{
#ifndef RELEASE
    CallStackEntry entry("ExpandPackedReflectors");
#endif
    if( uplo == LOWER && dir == VERTICAL )
        expand_packed_reflectors::LV( conjugation, offset, H, t );
    else
        throw std::logic_error("This option is not yet supported");
}

template<typename F> 
inline void
ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation,
  int offset, DistMatrix<F>& H, const DistMatrix<F,STAR,STAR>& t )
{
#ifndef RELEASE
    CallStackEntry entry("ExpandPackedReflectors");
#endif
    DistMatrix<F,MD,STAR> tDiag(H.Grid());
    tDiag.AlignWithDiagonal( H, offset );
    tDiag = t;
    ExpandPackedReflectors( uplo, dir, conjugation, offset, H, tDiag );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_EXPANDPACKEDREFLECTORS_HPP
