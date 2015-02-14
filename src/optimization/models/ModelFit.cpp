/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// NOTE: 
// This abstract ADMM routine is adapted from a MATLAB script written by 
// AJ Friend which is available at:
// https://github.com/ajfriend/admm_model/blob/ab7cc6555f730f351434a53867de422143001d58/model_fit.m

namespace El {

template<typename Real>
Int ModelFit
( function<void(Matrix<Real>&,Real)> lossProx,
  function<void(Matrix<Real>&,Real)> regProx,
  const Matrix<Real>& A, const Matrix<Real>& b, Matrix<Real>& w, 
  Real rho, Int maxIter, bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("ModelFit"))
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<Real> P;
    if( m >= n )
    {
        Identity( P, n, n );        
        Herk( LOWER, ADJOINT, Real(1), A, Real(1), P );
    }
    else
    {
        Identity( P, m, m );
        Herk( LOWER, NORMAL, Real(1), A, Real(1), P );
    }
    if( inv )
        HPDInverse( LOWER, P );
    else
        Cholesky( LOWER, P ); 

    // Start the ADMM
    Int numIter=0;
    Matrix<Real> s, x0, x1, x2, ux,
                    y0, y1, y2, uy;

    Zeros( x2, n, 1 );
    Ones( y2, m, 1 );

    Zeros( ux, n, 1 );
    Zeros( uy, m, 1 );

    for( ; numIter<maxIter-1; ++numIter )
    //while( numIter < maxIter )
    {
        // Project onto A x1 + b = y1
        x0 = x2; 
        y0 = y2;
        Axpy( Real(-1), ux, x0 );
        Axpy( Real(-1), uy, y0 );
        x1 = x0;
        // Overwrite y0 to perform a single Gemv
        Axpy( Real(-1), b, y0 ); 
        Gemv( ADJOINT, Real(1), A, y0, Real(1), x1 );
        if( inv )
        {
            s = x1;
            Hemv( LOWER, Real(1), P, s, Real(0), x1 );
        }
        else
        {
            Trsv( LOWER, NORMAL, NON_UNIT, P, x1 );
            Trsv( LOWER, ADJOINT, NON_UNIT, P, x1 );
        }
        y1 = b;
        Gemv( NORMAL, Real(1), A, x1, Real(1), y1 );

        y0 = y1;
        Axpy( Real(1), uy, y0 );
        y2 = y0;
        lossProx( y2, rho );

        x0 = x1;
        Axpy( Real(1), ux, x0 );       
        x2 = x0;
        regProx( x2, rho );

        // Update dual variables
        Axpy( Real(1), x1, ux );
        Axpy( Real(1), y1, uy );
        Axpy( Real(-1), x2, ux );
        Axpy( Real(-1), y2, uy );
    }
    if( maxIter == numIter )
        cout << "Model fit failed to converge" << endl;
    w = x2;
    return numIter;
}

template<typename Real>
Int ModelFit
( function<void(DistMatrix<Real>&,Real)> lossProx,
  function<void(DistMatrix<Real>&,Real)> regProx,
  const AbstractDistMatrix<Real>& APre, const AbstractDistMatrix<Real>& bPre, 
        AbstractDistMatrix<Real>& wPre, 
  Real rho, Int maxIter, bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("ModelFit"))

    auto APtr = ReadProxy<Real,MC,MR>( &APre );  auto& A = *APtr;
    auto bPtr = ReadProxy<Real,MC,MR>( &bPre );  auto& b = *bPtr;
    auto wPtr = WriteProxy<Real,MC,MR>( &wPre ); auto& w = *wPtr;

    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();

    DistMatrix<Real> P(g);
    if( m >= n )
    {
        Identity( P, n, n );        
        Herk( LOWER, ADJOINT, Real(1), A, Real(1), P );
    }
    else
    {
        Identity( P, m, m );
        Herk( LOWER, NORMAL, Real(1), A, Real(1), P );
    }
    if( inv )
        HPDInverse( LOWER, P );
    else
        Cholesky( LOWER, P ); 

    // Start the ADMM
    Int numIter=0;
    DistMatrix<Real> s(g), x0(g), x1(g), x2(g), ux(g),
                           y0(g), y1(g), y2(g), uy(g);

    Zeros( x2, n, 1 );
    Ones( y2, m, 1 );

    Zeros( ux, n, 1 );
    Zeros( uy, m, 1 );

    for( ; numIter<maxIter-1; ++numIter )
    //while( numIter < maxIter )
    {
        // Project onto A x1 + b = y1
        x0 = x2; 
        y0 = y2;
        Axpy( Real(-1), ux, x0 );
        Axpy( Real(-1), uy, y0 );
        x1 = x0;
        // Overwrite y0 to perform a single Gemv
        Axpy( Real(-1), b, y0 ); 
        Gemv( ADJOINT, Real(1), A, y0, Real(1), x1 );
        if( inv )
        {
            s = x1;
            Hemv( LOWER, Real(1), P, s, Real(0), x1 );
        }
        else
        {
            Trsv( LOWER, NORMAL, NON_UNIT, P, x1 );
            Trsv( LOWER, ADJOINT, NON_UNIT, P, x1 );
        }
        y1 = b;
        Gemv( NORMAL, Real(1), A, x1, Real(1), y1 );

        y0 = y1;
        Axpy( Real(1), uy, y0 );
        y2 = y0;
        lossProx( y2, rho );

        x0 = x1;
        Axpy( Real(1), ux, x0 );       
        x2 = x0;
        regProx( x2, rho );

        // Update dual variables
        Axpy( Real(1), x1, ux );
        Axpy( Real(1), y1, uy );
        Axpy( Real(-1), x2, ux );
        Axpy( Real(-1), y2, uy );
    }
    if( maxIter == numIter )
        cout << "Model fit failed to converge" << endl;
    w = x2;
    return numIter;
}

#define PROTO(Real) \
  template Int ModelFit \
  ( function<void(Matrix<Real>&,Real)> lossProx, \
    function<void(Matrix<Real>&,Real)> regProx, \
    const Matrix<Real>& A, const Matrix<Real>& b, Matrix<Real>& w, \
    Real rho, Int maxIter, bool inv, bool progress ); \
  template Int ModelFit \
  ( function<void(DistMatrix<Real>&,Real)> lossProx, \
    function<void(DistMatrix<Real>&,Real)> regProx, \
    const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, \
          AbstractDistMatrix<Real>& w, \
    Real rho, Int maxIter, bool inv, bool progress );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namepace elem
