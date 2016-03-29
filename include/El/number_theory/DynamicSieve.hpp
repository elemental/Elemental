/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_NUMBER_THEORY_DYNAMIC_SIEVE_HPP
#define EL_NUMBER_THEORY_DYNAMIC_SIEVE_HPP

namespace El {

template<typename TUnsigned>
DynamicSieve<TUnsigned>::DynamicSieve( TUnsigned lowerBound )
{
    oddPrimes.resize( 9 ); 
    oddPrimes[0] = 3;
    oddPrimes[1] = 5;
    oddPrimes[2] = 7;
    oddPrimes[3] = 11;
    oddPrimes[4] = 13;
    oddPrimes[5] = 17;
    oddPrimes[6] = 19;
    oddPrimes[7] = 23;
    oddPrimes[8] = 29;

    oddOffset = lowerBound;
    halvedIndex = 0; // candidate for prime is oddOffset + 2*halvedIndex

    TUnsigned numPrimes = oddPrimes.size();
    starts.resize( numPrimes );
    for( Int i=0; i<numPrimes; ++i )
        starts[i] = ComputeStart( oddPrimes[i] );
        
    // Allocate a table which would use roughly the same amount of memory
    // as the list of primes (which we crudely approximate to have a cardinality
    // on the order of Sqrt(oddOffset))
    size_t minTableSize = 256;
    TUnsigned tableSize = Max( size_t(std::sqrt(oddOffset)), minTableSize );
    table.resize( tableSize );

    // Ensure that every prime factor of a number in the range covered by the
    // table, [oddOffset,oddOffset+2*newTableSize), is in our list of primes
    TUnsigned largestCandidate = oddOffset + 2*(tableSize-1);
    while( oddPrimes.back()*oddPrimes.back() < largestCandidate )
    {
        AugmentPrimes();
    }

    Sieve();
}

template<typename TUnsigned>
bool DynamicSieve<TUnsigned>::SeekPrecomputedPrime()
{
    if( freshTable )
        freshTable = false;
    else
        ++halvedIndex;

    // Set halvedIndex to the first location greater than or equal to the
    // current location that is marked prime
    const TUnsigned tableSize = table.size();
    for( TUnsigned i=halvedIndex; i<tableSize; ++i )
    {
        if( table[i] != 0 )
        {
            halvedIndex = i;
            return true;
        }
    }
    return false;
}

// Find the integer k such that oddOffset + 2*k is the first odd 
// multiple of p that is greater than or equal to Max(oddOffset,p^2)
template<typename TUnsigned>
TUnsigned DynamicSieve<TUnsigned>::ComputeStart( const TUnsigned& p ) const
{
    if( p*p >= oddOffset )
    {
        return (p*p - oddOffset) / 2;
    }
    else
    {
        TUnsigned complement = (p - (oddOffset % p)) % p;
        if( complement % 2 == 0 )
        {
            return complement / 2;
        }
        else
        {
            return (complement + p) / 2;
        }
    }
}

template<typename TUnsigned>
void DynamicSieve<TUnsigned>::AugmentPrimes()
{
    TUnsigned pIndex = oddPrimes.size();
    TUnsigned numPrimes = 2*pIndex;

    oddPrimes.resize( numPrimes );
    starts.resize( numPrimes );

    TUnsigned p = oddPrimes[pIndex-1];
    for( ; pIndex<numPrimes; ++pIndex )
    {
        // Find the next prime after p using trial division.
        auto begIter = oddPrimes.begin();
        auto boundIter = begIter;
        while( true )
        {
            p += 2;

            TUnsigned factorBound = TUnsigned(std::sqrt(p)) + 1;
            boundIter =
              std::lower_bound( boundIter, begIter+pIndex, factorBound );
            const TUnsigned indexEnd = TUnsigned(boundIter-begIter);

            TUnsigned remainder = 1; // this could be an arbitrary nonzero
            for( TUnsigned j=0; j<indexEnd; ++j )
            {
                remainder = p % oddPrimes[j];
                if( remainder == 0 )
                {
                    // p is composite with factor oddPrimes[j]
                    break;
                }
            }
            if( remainder )
            {
                // None of the stored oddPrimes factored p
                break;
            }
        }
        oddPrimes[pIndex] = p;
        starts[pIndex] = ComputeStart( p );
    }
}

template<typename TUnsigned>
void DynamicSieve<TUnsigned>::Sieve()
{
    halvedIndex = 0;
    std::fill( table.begin(), table.end(), byte(1) );

    TUnsigned tableSize = table.size();
    TUnsigned numPrimes = oddPrimes.size();
    for( TUnsigned j=0; j<numPrimes; ++j )
    {
        TUnsigned p = oddPrimes[j];

        // Iterate over the indices corresponding to odd multiples of p that
        // are at least as large as Max(oddOffset,p^2)
        TUnsigned k;
        for( k=starts[j]; k<tableSize; k+=p )
            table[k] = 0;

        // The next time we use this value, it will be for a new table starting
        // at oddOffsetNew = oddOffset + 2*tableSize. Then the adjusted value
        // should be as follows, where k is the minimum member of
        //     starts[j] + p N
        // that is greater than or equal to tableSize.
        starts[j] = k - tableSize;
    }

    freshTable = true;
    halvedIndex = 0;
}

template<typename TUnsigned>
void DynamicSieve<TUnsigned>::FormNewTable()
{
    // The current table handles the odd integers in the inclusive interval
    // [oddOffset,oddOffset+2*(tableSize-1)], so the next table should begin at
    // 2*tableSize.
    oddOffset += 2*table.size();

    // Allocate a new table which would use roughly the same amount of memory
    // as the list of primes (which we crudely approximate to have a cardinality
    // on the order of Sqrt(oddOffset))
    const TUnsigned newTableSize =
      Max( size_t(std::sqrt(oddOffset)), table.size() );
    table.resize( newTableSize );

    // Ensure that every prime factor of a number in the range covered by the
    // table, [oddOffset,oddOffset+2*newTableSize), is in our list of primes
    const TUnsigned largestCandidate = oddOffset + 2*(newTableSize-1);
    while( oddPrimes.back()*oddPrimes.back() < largestCandidate )
    {
        AugmentPrimes();
    }

    Sieve();
}

template<typename TUnsigned>
TUnsigned DynamicSieve<TUnsigned>::CurrentPrime() const
{ return oddOffset + 2*halvedIndex; }

template<typename TUnsigned>
TUnsigned DynamicSieve<TUnsigned>::NextPrime()
{
    if( !SeekPrecomputedPrime() )
    {
        FormNewTable();
        bool foundPrecomputed = SeekPrecomputedPrime();
        DEBUG_ONLY(
          if( !foundPrecomputed )
              LogicError("Bug in implementation");
        )
    }
    return CurrentPrime();
}

} // namespace El

#endif // ifndef EL_NUMBER_THEORY_DYNAMIC_SIEVE_HPP
