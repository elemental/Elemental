/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CONVEX_SVT_HPP
#define ELEM_CONVEX_SVT_HPP

#include "elemental/convex/SVT/Normal.hpp"

#include "elemental/convex/SVT/Cross.hpp"
#include "elemental/convex/SVT/PivotedQR.hpp"
#include "elemental/convex/SVT/TSQR.hpp"

namespace elem {

template<typename F>
inline Int
SVT( Matrix<F>& A, BASE(F) tau, bool relative=false )
{
#ifndef RELEASE
    CallStackEntry cse("SVT");
#endif
    return svt::Normal( A, tau, relative );
}

template<typename F>
inline Int
SVT( DistMatrix<F>& A, BASE(F) tau, bool relative=false )
{
#ifndef RELEASE
    CallStackEntry cse("SVT");
#endif
    // NOTE: This should be less accurate (but faster) than svt::Normal
    return svt::Cross( A, tau, relative );
}

template<typename F>
inline Int
SVT( Matrix<F>& A, BASE(F) tau, Int relaxedRank, bool relative=false )
{
#ifndef RELEASE
    CallStackEntry cse("SVT");
#endif
    // Preprocess with numSteps iterations of pivoted QR factorization
    return svt::PivotedQR( A, tau, relaxedRank, relative );
}

template<typename F>
inline Int
SVT( DistMatrix<F>& A, BASE(F) tau, Int relaxedRank, bool relative=false )
{
#ifndef RELEASE
    CallStackEntry cse("SVT");
#endif
    // Preprocess with numSteps iterations of pivoted QR factorization
    return svt::PivotedQR( A, tau, relaxedRank, relative );
}

// Singular-value soft-thresholding based on TSQR
template<typename F,Distribution U>
inline Int
SVT( DistMatrix<F,U,STAR>& A, BASE(F) tau, bool relative=false )
{
#ifndef RELEASE
    CallStackEntry cse("SVT");
#endif
    return svt::TSQR( A, tau, relative );
}

} // namespace elem

#endif // ifndef ELEM_CONVEX_SVT_HPP
