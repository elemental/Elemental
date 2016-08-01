/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// NOTE:
// This implementation is based upon a Newton strategy due to AJ Friend.
// The original implementation may be found here:
// https://github.com/ajfriend/admm_model/blob/ab7cc6555f730f351434a53867de422143001d58/prox_logistic.m

namespace El {

template<typename Real>
void LogisticProx( Matrix<Real>& A, Real tau, Int numIts )
{
    DEBUG_CSE
    auto logisticProx = 
      [=]( Real alpha ) -> Real
      {
        // Use an initial guess based upon the rough normal vector to
        // the logistic curve
        Real beta;
        if( alpha < -Real(5) ) 
            beta = alpha + 1/tau;
        else if( alpha < Real(5) )
            beta = (Real(1)/Real(2) + tau*alpha)/(tau+Real(1)/Real(10));
        else
            beta = alpha;

        // Run a fixed number of Newton steps based upon the entrywise 
        // derivative of the prox objective,
        //    f(x) = -1/(exp(x)+1) + tau(x-x0),
        // as well as its derivative,
        //    f'(x) = exp(x)/(exp(x)+1)^2 + tau
        for( Int j=0; j<numIts; ++j )
        {
            // TODO: Use a faster and/or more stable algorithm
            const Real gamma = Exp(beta);
            const Real gammaP1 = gamma+1;
            const Real gammaP1Sq = gammaP1*gammaP1;
            beta += (gammaP1 + tau*(alpha-beta)*gammaP1Sq) / 
                    (gamma + tau*gammaP1Sq);
        }
        return beta;
      };
    EntrywiseMap( A, function<Real(Real)>(logisticProx) );
}

template<typename Real>
void LogisticProx( AbstractDistMatrix<Real>& A, Real tau, Int numIts )
{
    DEBUG_CSE
    auto logisticProx = 
      [=]( Real alpha ) -> Real
      {
        // Use an initial guess based upon the rough normal vector to
        // the logistic curve
        Real beta;
        if( alpha < -Real(5) ) 
            beta = alpha + 1/tau;
        else if( alpha < Real(5) )
            beta = (Real(1)/Real(2) + tau*alpha)/(tau+Real(1)/Real(10));
        else
            beta = alpha;

        // Run a fixed number of Newton steps based upon the entrywise 
        // derivative of the prox objective,
        //    f(x) = -1/(exp(x)+1) + tau(x-x0),
        // as well as its derivative,
        //    f'(x) = exp(x)/(exp(x)+1)^2 + tau
        for( Int j=0; j<numIts; ++j )
        {
            // TODO: Use a faster and/or more stable algorithm
            const Real gamma = Exp(beta);
            const Real gammaP1 = gamma+1;
            const Real gammaP1Sq = gammaP1*gammaP1;
            beta += (gammaP1 + tau*(alpha-beta)*gammaP1Sq) / 
                    (gamma + tau*gammaP1Sq);
        }
        return beta;
      };
    EntrywiseMap( A, function<Real(Real)>(logisticProx) );
}

#define PROTO(Real) \
  template void LogisticProx \
  ( Matrix<Real>& A, Real tau, Int numIts ); \
  template void LogisticProx \
  ( AbstractDistMatrix<Real>& A, Real tau, Int numIts );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
