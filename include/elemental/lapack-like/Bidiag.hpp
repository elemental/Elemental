/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_BIDIAG_HPP
#define LAPACK_BIDIAG_HPP

#include "./Bidiag/PanelL.hpp"
#include "./Bidiag/PanelU.hpp"
#include "./Bidiag/L.hpp"
#include "./Bidiag/LUnb.hpp"
#include "./Bidiag/U.hpp"
#include "./Bidiag/UUnb.hpp"

namespace elem {

template<typename R>
inline void
Bidiag( Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("Bidiag");
#endif
    if( IsComplex<R>::val )
        throw std::logic_error("Called real routine with complex datatype");
    if( A.Height() >= A.Width() )
        bidiag::U( A );
    else
        bidiag::L( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void Bidiag( DistMatrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("Bidiag");
#endif
    if( IsComplex<R>::val )
        throw std::logic_error("Called real routine with complex datatype");
    if( A.Height() >= A.Width() )
        bidiag::U( A );
    else
        bidiag::L( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void Bidiag
( Matrix<Complex<R> >& A,
  Matrix<Complex<R> >& tP,
  Matrix<Complex<R> >& tQ )
{
#ifndef RELEASE
    PushCallStack("Bidiag");
#endif
    if( A.Height() >= A.Width() )
        bidiag::U( A, tP, tQ );
    else
        bidiag::L( A, tP, tQ );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void Bidiag
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R>,STAR,STAR>& tP,
  DistMatrix<Complex<R>,STAR,STAR>& tQ )
{
#ifndef RELEASE
    PushCallStack("Bidiag");
#endif
    if( A.Height() >= A.Width() )
        bidiag::U( A, tP, tQ );
    else
        bidiag::L( A, tP, tQ );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef LAPACK_BIDIAG_HPP
