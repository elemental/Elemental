/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CONDITION_MAX_HPP
#define EL_CONDITION_MAX_HPP

namespace El {

template<typename F> 
Base<F> MaxCondition( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MaxCondition"))
    typedef Base<F> Real;
    Matrix<F> B( A );
    const Real maxNorm = MaxNorm( B );
    try { Inverse( B ); }
    catch( SingularMatrixException& e ) 
    { return std::numeric_limits<Real>::infinity(); }
    const Real maxNormInv = MaxNorm( B );
    return maxNorm*maxNormInv;
}

template<typename F,Dist U,Dist V> 
Base<F> MaxCondition( const DistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MaxCondition"))
    typedef Base<F> Real;
    DistMatrix<F> B( A );
    const Real maxNorm = MaxNorm( B );
    try { Inverse( B ); }
    catch( SingularMatrixException& e ) 
    { return std::numeric_limits<Real>::infinity(); }
    const Real maxNormInv = MaxNorm( B );
    return maxNorm*maxNormInv;
}

} // namespace El

#endif // ifndef EL_CONDITION_MAX_HPP
