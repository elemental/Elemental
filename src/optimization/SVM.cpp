/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include EL_IDENTITY_INC
#include EL_ZEROS_INC

// NOTE: This is adapted from a MATLAB script written by AJ Friend.

namespace El {

template<typename F>
Int SVM
( const Matrix<F>& A, Base<F> gamma, Matrix<F>& z,
  Base<F> rho, Base<F> alpha, Int maxIter, Base<F> absTol, Base<F> relTol, 
  bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("SVM"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<F> P;
    if( m >= n )
    {
        Identity( P, n, n );        
        Herk( LOWER, ADJOINT, F(1), A, F(1), P );
    }
    else
    {
        Identity( P, m, m );
        Herk( LOWER, NORMAL, F(1), A, F(1), P );
    }
    if( inv )
        HPDInverse( LOWER, P );
    else
        Cholesky( LOWER, P ); 

    // Start the SVM
    Int numIter=0;
    Matrix<F> x, u, s, zOld, xHat;
    Zeros( x, n, 1 );
    Zeros( z, n, 1 );
    Zeros( u, n, 1 );
    while( numIter < maxIter )
    {
        // TODO
        LogicError("Not yet written");
    }
    if( maxIter == numIter )
        std::cout << "SVM failed to converge" << std::endl;
    return numIter;
}

template<typename F>
Int SVM
( const DistMatrix<F>& A, Base<F> gamma, DistMatrix<F>& z,
  Base<F> rho, Base<F> alpha, Int maxIter, Base<F> absTol, Base<F> relTol, 
  bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("SVM"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();

    DistMatrix<F> P(g);
    if( m >= n )
    {
        Identity( P, n, n );        
        Herk( LOWER, ADJOINT, F(1), A, F(1), P );
    }
    else
    {
        Identity( P, m, m );
        Herk( LOWER, NORMAL, F(1), A, F(1), P );
    }
    if( inv )
        HPDInverse( LOWER, P );
    else
        Cholesky( LOWER, P ); 

    // Start the SVM
    Int numIter=0;
    DistMatrix<F> x(g), u(g), s(g), zOld(g), xHat(g);
    Zeros( x, n, 1 );
    Zeros( z, n, 1 );
    Zeros( u, n, 1 );
    while( numIter < maxIter )
    {
        // TODO
        LogicError("Not yet written");
    }
    if( maxIter == numIter )
        std::cout << "SVM failed to converge" << std::endl;
    return numIter;
}


#define PROTO(F) \
  template Int SVM \
  ( const Matrix<F>& A, Base<F> gamma, \
    Matrix<F>& z, \
    Base<F> rho, Base<F> alpha, Int maxIter, Base<F> absTol, Base<F> relTol, \
    bool inv, bool progress ); \
  template Int SVM \
  ( const DistMatrix<F>& A, Base<F> gamma, \
    DistMatrix<F>& z, \
    Base<F> rho, Base<F> alpha, Int maxIter, Base<F> absTol, Base<F> relTol, \
    bool inv, bool progress ); \

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namepace elem
