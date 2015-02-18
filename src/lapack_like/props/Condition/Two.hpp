/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CONDITION_TWO_HPP
#define EL_CONDITION_TWO_HPP

namespace El {

template<typename F> 
Base<F> TwoCondition( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("TwoCondition"))
    typedef Base<F> Real;
    Matrix<F> B( A );
    Matrix<Real> s;
    SVD( B, s );

    Real cond = 1;
    const Int numVals = s.Height();
    if( numVals > 0 )
        cond = s.Get(0,0) / s.Get(numVals-1,0);
    return cond;
}

template<typename F> 
Base<F> TwoCondition( const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("TwoCondition"))
    typedef Base<F> Real;
    DistMatrix<F> B( A );
    DistMatrix<Real,VR,STAR> s( A.Grid() );
    SVD( B, s );

    Real cond = 1;
    const Int numVals = s.Height();
    if( numVals > 0 )
        cond = s.Get(0,0) / s.Get(numVals-1,0);
    return cond;
}

} // namespace El

#endif // ifndef EL_CONDITION_TWO_HPP
