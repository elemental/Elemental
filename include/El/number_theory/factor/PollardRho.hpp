/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_NUMBER_THEORY_FACTOR_POLLARD_RHO_HPP
#define EL_NUMBER_THEORY_FACTOR_POLLARD_RHO_HPP

#ifdef EL_HAVE_MPC

namespace El {

namespace factor {

namespace pollard_rho {

// TODO: Add the ability to set a maximum number of iterations
inline BigInt FindFactor
( const BigInt& n,
  Int a,
  const PollardRhoCtrl& ctrl )
{
    if( a == 0 || a == -2 )
        Output("WARNING: Problematic choice of Pollard rho shift");
    BigInt tmp, gcd;
    BigInt one(1);

    auto xAdvance =
      [&]( BigInt& x )
      {
        if( ctrl.numSteps == 1 )
        {
            // TODO: Determine if there is a penalty to x *= x
            /*
            tmp = x;
            tmp *= x;
            tmp += a;
            x = tmp;
            x %= n;
            */
            x *= x;
            x += a;
            x %= n;
        }
        else
        {
            PowMod( x, 2*ctrl.numSteps, n, x );
            x += a;
            x %= n;
        }
      };

    auto QAdvance =
      [&]( const BigInt& x, const BigInt& x2, BigInt& Q )
      {
        tmp = x2;
        tmp -= x;
        Q *= tmp;
        Q %= n;
      };

    Int gcdDelay = ctrl.gcdDelay;
    BigInt xi=ctrl.x0;
    BigInt x2i(xi);
    BigInt xiSave=xi, x2iSave=x2i;
    BigInt Qi(1);
    Int k=1, i=1; // it is okay for i to overflow since it is just for printing
    while( true )
    {
        // Advance xi once
        xAdvance( xi );

        // Advance x2i twice
        xAdvance( x2i );
        xAdvance( x2i );

        // Advance Qi
        QAdvance( xi, x2i, Qi );

        if( k >= gcdDelay )
        {
            GCD( Qi, n, gcd );
            if( gcd > one )
            {
                // NOTE: This was not suggested by Pollard's original paper
                if( gcd == n )
                {
                    if( gcdDelay == 1 )
                    {
                        RuntimeError("(x) converged before (x mod p) at i=",i);
                    }
                    else
                    {
                        if( ctrl.progress )
                            Output("Backtracking at i=",i);
                        i = Max( i-(gcdDelay+1), 0 );
                        gcdDelay = 1;
                        xi = xiSave;
                        x2i = x2iSave;
                    }
                }
                else
                {
                    if( ctrl.progress )
                        Output("Found factor ",gcd," at i=",i); 
                    return gcd;
                }
            }

            // NOTE: This was not suggested by Pollard's original paper
            k = 0;
            xiSave = xi;
            x2iSave = x2i;
            Qi = 1;
        }
        ++k;
        ++i;
    }
}

} // namespace pollard_rho

inline vector<BigInt> PollardRho
( const BigInt& n,
  const PollardRhoCtrl& ctrl )
{
    vector<BigInt> factors;
    BigInt nRem = n;

    Timer timer;
    PushIndent();
    while( true )
    {
        // Try Miller-Rabin first
        if( ctrl.time )
            timer.Start();
        Primality primality = PrimalityTest( nRem, ctrl.numReps );
        if( primality == PRIME )
        {
            if( ctrl.time )
                Output(nRem," is prime (",timer.Stop()," seconds)");
            else if( ctrl.progress )
                Output(nRem," is prime");
            factors.push_back( nRem );       
            break;
        }
        else if( primality == PROBABLY_PRIME )
        {
            if( ctrl.time )
                Output(nRem," is probably prime (",timer.Stop()," seconds)");
            else if( ctrl.progress )
                Output(nRem," is probably prime");
            factors.push_back( nRem );
            break;
        }
        else
        {
            if( ctrl.time )
                Output(nRem," is composite (",timer.Stop()," seconds)");
            else if( ctrl.progress )
                Output(nRem," is composite");
        }

        if( ctrl.progress )
            Output("Attempting to factor ",nRem," with a=",ctrl.a0);
        if( ctrl.time )
            timer.Start();
        PushIndent();
        BigInt factor;
        try
        {
            factor = pollard_rho::FindFactor( nRem, ctrl.a0, ctrl );
        }
        catch( const exception& e ) // TODO: Introduce factor exception?
        {
            // Try again with a=ctrl.a1
            if( ctrl.progress )
                Output("Attempting to factor ",nRem," with a=",ctrl.a1);
            factor = pollard_rho::FindFactor( nRem, ctrl.a1, ctrl );
        }
        if( ctrl.time )
            Output("Pollard-rho: ",timer.Stop()," seconds");
        PopIndent();

        factors.push_back( factor );
        nRem /= factor;
    }
    PopIndent();
    sort( factors.begin(), factors.end() );
    return factors;
}

} // namespace factor

} // namespace El

#endif // ifdef EL_HAVE_MPC

#endif // ifndef EL_NUMBER_THEORY_FACTOR_POLLARD_RHO_HPP
