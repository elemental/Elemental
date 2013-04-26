/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_EXPANDPACKEDREFLECTORS_HPP
#define LAPACK_EXPANDPACKEDREFLECTORS_HPP

#include "elemental/lapack-like/ApplyPackedReflectors/Util.hpp"

#include "./ExpandPackedReflectors/LV.hpp"

namespace elem {

template<typename R> 
inline void
ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir,
  int offset, Matrix<R>& H )
{
#ifndef RELEASE
    CallStackEntry entry("ExpandPackedReflectors");
#endif
    // Since the complex version does not have the same argument list, there is
    // currently no good way to ensure that this version is not called with 
    // complex datatypes. Until C++11 compilers are commonplace, we cannot
    // use static_assert either.
    if( IsComplex<R>::val )
        throw std::logic_error("Called real routine with complex datatype");

    if( uplo == LOWER && dir == VERTICAL )
        expand_packed_reflectors::LV( offset, H );
    else
        throw std::logic_error("This option is not yet supported");
}

template<typename R> 
inline void
ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir, 
  int offset, DistMatrix<R>& H )
{
#ifndef RELEASE
    CallStackEntry entry("ExpandPackedReflectors");
#endif
    // Since the complex version does not have the same argument list, there is
    // currently no good way to ensure that this version is not called with 
    // complex datatypes. Until C++11 compilers are commonplace, we cannot
    // use static_assert either.
    if( IsComplex<R>::val )
        throw std::logic_error("Called real routine with complex datatype");

    if( uplo == LOWER && dir == VERTICAL )
        expand_packed_reflectors::LV( offset, H );
    else
        throw std::logic_error("This option is not yet supported");
}

template<typename R> 
inline void
ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation,
  int offset, Matrix<Complex<R> >& H, const Matrix<Complex<R> >& t )
{
#ifndef RELEASE
    CallStackEntry entry("ExpandPackedReflectors");
#endif
    if( uplo == LOWER && dir == VERTICAL )
        expand_packed_reflectors::LV( conjugation, offset, H, t );
    else
        throw std::logic_error("This option is not yet supported");
}

template<typename R> 
inline void
ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation,
  int offset, 
        DistMatrix<Complex<R> >& H, 
  const DistMatrix<Complex<R>,MD,STAR>& t )
{
#ifndef RELEASE
    CallStackEntry entry("ExpandPackedReflectors");
#endif
    if( uplo == LOWER && dir == VERTICAL )
        expand_packed_reflectors::LV( conjugation, offset, H, t );
    else
        throw std::logic_error("This option is not yet supported");
}

template<typename R> 
inline void
ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation,
  int offset,
        DistMatrix<Complex<R> >& H, 
  const DistMatrix<Complex<R>,STAR,STAR>& t )
{
#ifndef RELEASE
    CallStackEntry entry("ExpandPackedReflectors");
#endif
    DistMatrix<Complex<R>,MD,STAR> tDiag(H.Grid());
    tDiag.AlignWithDiagonal( H, offset );
    tDiag = t;
    ExpandPackedReflectors( uplo, dir, conjugation, offset, H, tDiag );
}

} // namespace elem

#endif // ifndef LAPACK_EXPANDPACKEDREFLECTORS_HPP
