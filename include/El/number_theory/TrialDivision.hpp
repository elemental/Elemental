/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_NUMBER_THEORY_TRIAL_DIVISION_HPP
#define EL_NUMBER_THEORY_TRIAL_DIVISION_HPP

namespace El {

inline vector<unsigned long long>
TrialDivision( unsigned long long n, unsigned long long limit )
{
    limit = Min(limit,ISqrt(n));
    
    auto& sieve = TrialDivisionSieve();
    sieve.Generate( limit );
    
    vector<unsigned long long> factors;
    while( n % 2U == 0 )
    {
        factors.push_back( 2U );
        n /= 2U;
    }

    auto primeEnd =
      std::upper_bound
      ( sieve.oddPrimes.begin(),
        sieve.oddPrimes.end(), limit );
    for( auto iter=sieve.oddPrimes.begin(); iter<primeEnd; ++iter )
    {
        const auto& p = *iter; 
        while( n % p == 0 )
        {
            factors.push_back( p );
            n /= p;
        }
    }
    return factors;
}
#ifdef EL_HAVE_MPC
inline vector<unsigned long long>
TrialDivision( const BigInt& n, unsigned long long limit )
{
    // Implement Min carefully
    if( BigInt(limit) > ISqrt(n) )
        limit = static_cast<unsigned long long>(ISqrt(n));
    
    auto& sieve = TrialDivisionSieve();
    sieve.Generate( limit );
    
    BigInt nRem(n);
    vector<unsigned long long> factors;
    while( nRem % 2U == 0 )
    {
        factors.push_back( 2U );
        nRem /= 2U;
    }

    auto primeEnd =
      std::upper_bound
      ( sieve.oddPrimes.begin(),
        sieve.oddPrimes.end(), limit );
    for( auto iter=sieve.oddPrimes.begin(); iter<primeEnd; ++iter )
    {
        const auto& p = *iter; 
        // TODO: Combine modulus and division using mpz_fdiv_qr
        while( nRem % p == 0 )
        {
            factors.push_back( p );
            nRem /= p;
        }
    }
    return factors;
}
#endif

inline bool
HasTinyFactor( unsigned long long n, unsigned long long limit )
{
    limit = Min(limit,ISqrt(n));

    auto& sieve = TrialDivisionSieve();
    sieve.Generate( limit );
    
    if( n % 2U == 0 )
        return true;

    auto primeEnd =
      std::upper_bound
      ( sieve.oddPrimes.begin(),
        sieve.oddPrimes.end(), limit );
    for( auto iter=sieve.oddPrimes.begin(); iter<primeEnd; ++iter )
    {
        const auto& p = *iter; 
        if( n % p == 0 )
            return true;
    }

    return false;
}
#ifdef EL_HAVE_MPC
inline bool
HasTinyFactor( const BigInt& n, unsigned long long limit )
{
    // Implement Min carefully
    if( BigInt(limit) > ISqrt(n) )
        limit = static_cast<unsigned long long>(ISqrt(n));
    
    auto& sieve = TrialDivisionSieve();
    sieve.Generate( limit );
    
    vector<unsigned long long> factors;
    if( n % 2U == 0 )
        return true;

    auto primeEnd =
      std::upper_bound
      ( sieve.oddPrimes.begin(),
        sieve.oddPrimes.end(), limit );
    for( auto iter=sieve.oddPrimes.begin(); iter<primeEnd; ++iter )
    {
        const auto& p = *iter; 
        if( n % p == 0 )
            return true;
    }
    return false;
}
#endif

} // namespace El

#endif // ifndef EL_NUMBER_THEORY_TRIAL_DIVISION_HPP
