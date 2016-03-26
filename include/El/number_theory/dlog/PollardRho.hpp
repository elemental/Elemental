/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_NUMBER_THEORY_DLOG_POLLARD_RHO_HPP
#define EL_NUMBER_THEORY_DLOG_POLLARD_RHO_HPP

#ifdef EL_HAVE_MPC

namespace El {

namespace dlog {

// TODO: Add the ability to set a maximum number of iterations
inline BigInt PollardRho
( const BigInt& q,
  const BigInt& r,
  const BigInt& p,
  const PollardRhoCtrl& ctrl )
{
    BigInt pm1(p);
    pm1 -=1;

    BigInt pOneThird(p);
    pOneThird /= 3;

    BigInt pTwoThirds(p);
    pTwoThirds *= 2;
    pTwoThirds /= 3;

    if( ctrl.demandPrimitive )
    {
        bool rIsPrimitive =
          IsPrimitiveRoot( r, p, ctrl.progress, ctrl.factorCtrl );
        if( !rIsPrimitive )
            LogicError(r," was not a primitive root of ",p);
    }

    auto xAdvance =
      [&]( BigInt& x, BigInt& a, BigInt& b )
      {
          if( x <= pOneThird )
          {
              x *= q;
              x %= p;
              ++a;
              a %= pm1;
          }
          else if( x <= pTwoThirds )
          {
              x *= x;
              x %= p;
              a *= 2;
              a %= pm1;
              b *= 2;
              b %= pm1;
          }
          else
          {
              x *= r;
              x %= p;
              ++b;
              b %= pm1;
          }
      };

    // Initialize a_0, b_0, and x_0 = q^(a_0) * r^(b_0)
    BigInt ai=ctrl.a0, bi=ctrl.b0;
    BigInt xi;
    {
        PowMod( q, ai, p, xi );
        BigInt tmp;
        PowMod( r, bi, p, tmp );
        xi *= tmp;
    }
    
    BigInt a2i(ai), b2i(bi), x2i(xi);
    BigInt m, n, d, lambda, mu, k, exponent, theta, thetaExp, candidate;
    Int i=1; // it is okay for i to overflow since it is just for printing
    while( true )
    {
        // Advance xi once
        xAdvance( xi, ai, bi );

        // Advance x2i twice
        xAdvance( x2i, a2i, b2i );
        xAdvance( x2i, a2i, b2i );

        if( xi == x2i )
        {
            if( ctrl.progress )
                Output("Detected cycle at iteration ",i);

            m = ai;
            m -= a2i;
            m %= pm1;

            n = b2i;
            n -= bi;
            n %= pm1;

            ExtendedGCD( m, pm1, d, lambda, mu );
            if( ctrl.progress )
                Output("GCD(",m,",",pm1,")=",d);

            // Solve for k in lambda*n = d*k
            k = lambda;
            k *= n;
            k /= d;
            k %= pm1; // TODO: Decide if this mod is necessary

            // theta := pow( r, (p-1)/d ) 
            exponent = pm1;
            exponent /= d;
            PowMod( r, exponent, p, theta );

            // candidate := r^k
            PowMod( r, k, p, candidate );

            // Try each power of theta
            for( thetaExp=0; thetaExp<d; ++thetaExp )
            {
                if( candidate == q )
                {
                    BigInt discLog = k + thetaExp*exponent;
                    if( ctrl.progress )
                        Output("Returning ",discLog," at thetaExp=",thetaExp);
                    return discLog;
                }
                candidate *= theta;
                candidate %= p;
            }
            LogicError("No power of theta yielded q = r^k theta^i (mod p)");
        }
        ++i;
    }
    // This should never occur and is to prevent compiler warnings
    return BigInt(-1);
}

} // namespace dlog

} // namespace El

#endif // ifdef EL_HAVE_MPC

#endif // ifndef EL_NUMBER_THEORY_DLOG_POLLARD_RHO_HPP
