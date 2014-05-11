/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SVT_HPP
#define ELEM_SVT_HPP

#include "./SVT/Normal.hpp"
#include "./SVT/Cross.hpp"
#include "./SVT/PivotedQR.hpp"
#include "./SVT/TSQR.hpp"

namespace elem {

template<typename F>
inline Int
SVT( Matrix<F>& A, Base<F> tau, bool relative=false )
{
    DEBUG_ONLY(CallStackEntry cse("SVT"))
    return svt::Normal( A, tau, relative );
}

template<typename F>
inline Int
SVT( DistMatrix<F>& A, Base<F> tau, bool relative=false )
{
    DEBUG_ONLY(CallStackEntry cse("SVT"))
    // NOTE: This should be less accurate (but faster) than svt::Normal
    return svt::Cross( A, tau, relative );
}

template<typename F>
inline Int
SVT( Matrix<F>& A, Base<F> tau, Int relaxedRank, bool relative=false )
{
    DEBUG_ONLY(CallStackEntry cse("SVT"))
    // Preprocess with numSteps iterations of pivoted QR factorization
    return svt::PivotedQR( A, tau, relaxedRank, relative );
}

template<typename F>
inline Int
SVT( DistMatrix<F>& A, Base<F> tau, Int relaxedRank, bool relative=false )
{
    DEBUG_ONLY(CallStackEntry cse("SVT"))
    // Preprocess with numSteps iterations of pivoted QR factorization
    return svt::PivotedQR( A, tau, relaxedRank, relative );
}

// Singular-value soft-thresholding based on TSQR
template<typename F,Dist U>
inline Int
SVT( DistMatrix<F,U,STAR>& A, Base<F> tau, bool relative=false )
{
    DEBUG_ONLY(CallStackEntry cse("SVT"))
    return svt::TSQR( A, tau, relative );
}

} // namespace elem

#endif // ifndef ELEM_SVT_HPP
