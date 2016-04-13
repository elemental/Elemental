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

namespace pollard_rho {

// For use within a Pohlig-Hellman decomposition
// NOTE: This implementation is meant to support subgroups of (Z/nZ)*, such
//       as the n=5 case with r=4 implies the subgroup {4,4^2=16=1} of order 2.
// TODO: Add the ability to set a maximum number of iterations
inline BigInt Subproblem
( const BigInt& q,
  const BigInt& r,
  const BigInt& n,
  const BigInt& subgroupOrder,
  const PollardRhoCtrl& ctrl )
{
    const BigInt& zero = BigIntZero();
    const BigInt& one = BigIntOne();

    // Ensure that q lives in (Z/nZ)*
    if( q < one || q >= n )
        LogicError(q," was not in [1,",n,")");
    if( GCD(q,n) != one )
        LogicError("GCD(",q,",",n,")=",GCD(q,n));

    // Ensure that r lives in (Z/nZ)*
    if( r < one || r >= n )
        LogicError(r," was not in [1,",n,")");
    if( GCD(r,n) != one )
        LogicError("GCD(",r,",",n,")=",GCD(r,n));

    // Check the (unlikely) case that r is one
    if( r == one )
    {
        if( q == one )
            return zero;
        else
            LogicError("One does not generate ",q);
    }

    BigInt nOneThird(n);
    nOneThird /= 3;

    BigInt nTwoThirds(n);
    nTwoThirds *= 2;
    nTwoThirds /= 3;

    auto xAdvance =
      [&]( BigInt& x, BigInt& a, BigInt& b )
      {
          if( x <= nOneThird )
          {
              x *= q;
              x %= n;
              ++a;
              a %= subgroupOrder;
          }
          else if( x <= nTwoThirds )
          {
              x *= x;
              x %= n;
              a *= 2;
              a %= subgroupOrder;
              b *= 2;
              b %= subgroupOrder;
          }
          else
          {
              x *= r;
              x %= n;
              ++b;
              b %= subgroupOrder;
          }
      };

    // Initialize a_0, b_0, and x_0 = q^(a_0) * r^(b_0)
    BigInt ai=ctrl.a0, bi=ctrl.b0;
    BigInt xi;
    {
        PowMod( q, ai, n, xi );
        BigInt tmp;
        PowMod( r, bi, n, tmp );
        xi *= tmp;
    }
    
    BigInt a2i(ai), b2i(bi), x2i(xi);
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

            BigInt aDiff = (ai - a2i) % subgroupOrder;
            BigInt bDiff = (b2i - bi) % subgroupOrder;
            // NOTE:
            // We should not necessarily throw an exception if bDiff=0;
            // consider the problem 1 = (n-1)^x (mod n), which will converge
            // at iteration 1 since (n-1)^2 = 1 (mod n) for any n. We will
            // instead attempt to detect degeneracy below.

            BigInt d, lambda, mu;
            ExtendedGCD( aDiff, subgroupOrder, d, lambda, mu );
            if( ctrl.progress )
                Output("GCD(",aDiff,",",subgroupOrder,")=",d);

            // Solve for k in lambda*bDiff = d*k.
            // Note that such a relationship of r^(lambda*bDiff) = r^(d*k)
            // need not exist if r does not generate q.
            BigInt k = (lambda*bDiff) / d;
            k %= subgroupOrder;

            // Q := q r^{-k}
            BigInt Q = PowMod( r, -k, n );
            Q *= q;
            Q %= n;

            // theta := pow( r, subgroupOrder/d ) 
            BigInt exponent(subgroupOrder);
            exponent /= d;
            BigInt theta = PowMod( r, exponent, n );

            // Test theta^i = Q for each i
            // (Also test theta^i = -Q, which implies theta^{i+d/2} = Q
            //  if r was a primitive root)
            BigInt thetaPow(theta);
            BigInt negQ(Q);
            negQ *= -1;
            negQ %= n;
            for( BigInt thetaExp=0; thetaExp<d; ++thetaExp )
            {
                if( thetaPow == Q )
                {
                    BigInt discLog = k + thetaExp*exponent;
                    if( ctrl.progress )
                        Output("Returning ",discLog," at thetaExp=",thetaExp);
                    return discLog;
                }
                else if( thetaPow == negQ )
                {
                    BigInt dHalf(d);
                    dHalf /= 2;
                    BigInt theta_dHalf = PowMod( theta, dHalf, n );
                    if( Mod(thetaPow*theta_dHalf,n) == Q )
                    {
                        BigInt discLog = k + (thetaExp+dHalf)*exponent;
                        if( ctrl.progress )
                            Output
                            ("Took -Q shortcut at thetaExp=",thetaExp,
                             " and found discLog=",discLog);
                        return discLog; 
                    }
                    else if( ctrl.progress )
                        Output("-Q shortcut failed at thetaExp=",thetaExp);
                } 
                if( thetaPow == one )
                {
                    LogicError
                    ("theta=r^(",subgroupOrder,"/",d,")=",theta,
                     " was a degenerate ",d,"'th root, as theta^",
                     thetaExp,"=1, and r does not generate q");
                }
                thetaPow *= theta;
                thetaPow %= n;
            }

            LogicError("This should not be possible");
        }
        ++i;
    }

    // This should never occur and is to prevent compiler warnings
    return BigInt(-1);
}

} // namespace pollard_rho

// TODO: Decide if the degeneracy detection within the core Pollard rho
//       algorithm is enough, or if there is a possibility of failure
//       within the Pohlig-Hellman decomposition logic as well
inline BigInt PollardRho
( const BigInt& q,
  const BigInt& r,
  const BigInt& n,
  const PollardRhoCtrl& ctrl )
{
    const BigInt& one = BigIntOne();

    // Ensure that q lives in (Z/nZ)*
    if( q < one || q >= n )
        LogicError(q," was not in [1,",n,")");
    if( GCD(q,n) != one )
        LogicError("GCD(",q,",",n,")=",GCD(q,n));

    // Ensure that r lives in (Z/nZ)*
    if( r < one || r >= n )
        LogicError(r," was not in [1,",n,")");
    if( GCD(r,n) != one )
        LogicError("GCD(",r,",",n,")=",GCD(r,n));

    // Check the (unlikely) case that r is one
    if( r == one )
    {
        if( q == one )
            return BigInt(0);
        else
            LogicError("One does not generate ",q);
    }

    if( ctrl.assumePrime && !ctrl.multistage )
    {
        // If we assume that n is prime, then \phi(n) = n-1 and we do not
        // need to attempt a prime factorization of n
        BigInt nTotient = n-1;
        if( ctrl.progress )
            Output("Assuming ",n," is prime");
        return pollard_rho::Subproblem( q, r, n, nTotient, ctrl );
    }

    // Since we are implicitly working within (Z/nZ)*, which has group order
    // 
    //     \phi(n) = \Prod_i p_i^(k-1) (p_i - 1),
    //
    // where n = \Prod_i p_i^k is the prime factorization of n, finding the
    // prime factorization of \phi(n) given a prime factorization of n is a  
    // matter of factoring each (p_i-1) and appending the factors to the 
    // list of factors from the p_i^(k-1). 
    auto nFactors = factor::PollardRho( n, ctrl.factorCtrl );
    if( ctrl.progress )
    {
        Output("Factored ",n," as:");
        for( const auto& factor : nFactors )
            Output("  ",factor);
    }
    vector<BigInt> nUniqueFactors;
    vector<Unsigned> nFactorPowers;
    nUniqueFactors.push_back( nFactors[0] );
    nFactorPowers.push_back( 1 );
    for( Unsigned index=1; index<nFactors.size(); ++index )
    {
        if( nFactors[index] == nFactors[index-1] )
        {
            ++nFactorPowers.back();
        }
        else
        {
            nUniqueFactors.push_back( nFactors[index] );
            nFactorPowers.push_back( 1 );
        }
    }

    // Compute Euler's totient function of n
    BigInt nTotient(1);
    for( Unsigned pIndex=0; pIndex<nUniqueFactors.size(); ++pIndex )
    {
        const BigInt& p = nUniqueFactors[pIndex];
        Unsigned power = nFactorPowers[pIndex];
        nTotient *= Pow(p,power-1);
        nTotient *= p-1;
    }
    if( ctrl.progress )
        Output("Computed phi(",n,")=",nTotient);

    if( !ctrl.multistage )
    {
        // We are now ready to call the one-shot Pollard rho algorithm
        if( ctrl.progress )
            Output("Running single-stage Pollard rho");
        return pollard_rho::Subproblem( q, r, n, nTotient, ctrl );
    }

    // Translate the list of unique prime powers factoring n into a list of
    // prime powers factoring \phi(n) by transforming one power of each p_i
    // into the prime factorization of (p_i-1).
    vector<BigInt> totientFactors;
    for( Unsigned pIndex=0; pIndex<nUniqueFactors.size(); ++pIndex )
    {
        const BigInt& p = nUniqueFactors[pIndex];    
        Unsigned power = nFactorPowers[pIndex];

        BigInt pm1 = p-1;
        auto pm1Factors = factor::PollardRho( pm1, ctrl.factorCtrl );
        for( const auto& factor : pm1Factors )
            totientFactors.push_back( factor );    

        for( Unsigned k=0; k<power-1; ++k )
            totientFactors.push_back( p );
    }
    sort( totientFactors.begin(), totientFactors.end() );
    if( ctrl.progress )
    {
        Output("Factored phi(",n,")=",nTotient," as:");
        for( const auto& factor : totientFactors )
            Output("  ",factor);
    }
    vector<BigInt> totientUniqueFactors;
    vector<Unsigned> totientFactorPowers;
    totientUniqueFactors.push_back( totientFactors[0] );
    totientFactorPowers.push_back( 1 );
    for( Unsigned index=1; index<totientFactors.size(); ++index )
    {
        if( totientFactors[index] == totientFactors[index-1] )
        {
            ++totientFactorPowers.back();
        }
        else
        {
            totientUniqueFactors.push_back( totientFactors[index] );
            totientFactorPowers.push_back( 1 );
        }
    }

    // We now apply the decomposition of Pohlig and Hellman to accelerate the
    // discrete logarithm using a prime factorization of Euler's totient
    // function, \phi(n). See the discussion in the proximity of Eq. (4.5) in
    // 
    //   Kevin S. McCurley, "The Discrete Logarithm Problem",
    //   Proceedings of Symposia in Applied Mathematics, Volume 42, 1990.
    //
    // for a concise overview of the algorithm.
    //
    BigInt x=0;
    for( Unsigned pIndex=0; pIndex<totientUniqueFactors.size(); ++pIndex )
    {
        const BigInt& p = totientUniqueFactors[pIndex];
        Unsigned power = totientFactorPowers[pIndex];
        if( ctrl.progress )
            Output("Prime power ",pIndex," is ",p,"^",power);

        BigInt rPrime = PowMod(r,nTotient/p,n);

        BigInt xp=0;
        vector<BigInt> radixDecomp(power);
        for( Unsigned i=0; i<power; ++i )
        {
            BigInt z(0);
            for( Unsigned j=0; j<i; ++j )
                z += radixDecomp[j]*Pow(p,j);
            
            BigInt qPrime = InvertMod(PowMod(r,z,n),n);
            qPrime *= q;
            qPrime %= n;
            qPrime = PowMod(qPrime,nTotient/Pow(p,i+1),n);

            if( ctrl.progress )
                Output
                ("  Solving subproblem ",qPrime," = ",rPrime,
                 "^x (mod ",n,") with subgroup order ",p);
            PushIndent();
            radixDecomp[i] =
              pollard_rho::Subproblem( qPrime, rPrime, n, p, ctrl );
            PopIndent();
            if( ctrl.progress )
                Output("  Subproblem index was ",radixDecomp[i]);
            xp += radixDecomp[i]*Pow(p,i);
        }
        if( ctrl.progress )
            Output("Index (mod ",p,"^",power,") was ",xp);

        // Add the contribution to the solution via Chinese Remainder Theorem
        BigInt pPower = Pow(p,power);
        BigInt pRem = nTotient/pPower;
        BigInt pRemInv = InvertMod(pRem,pPower);
        x += xp*pRem*pRemInv;
    }
    // TODO: Decide if this is necessary
    x %= nTotient;
    if( ctrl.progress )
        Output("Combined solution (via CRT) was ",x);
    //DEBUG_ONLY(
      if( PowMod(r,x,n) != q )
          LogicError
          ("Error in Pohlig-Hellman: ",r,"^",x," = ",PowMod(r,x,n)," != ",q);
    //)
    return x;
}

} // namespace dlog

} // namespace El

#endif // ifdef EL_HAVE_MPC

#endif // ifndef EL_NUMBER_THEORY_DLOG_POLLARD_RHO_HPP
