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
DynamicSieve<TUnsigned>::DynamicSieve( TUnsigned lowerBound, bool keepAll )
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

    lowerBound_ = segmentOffset_ = ( keepAll ? 3 : lowerBound );

    TUnsigned numPrimes = oddPrimes.size();
    
    // Allocate a table which would use roughly the same amount of memory
    // as the list of primes (which we crudely approximate to have a
    // cardinality on the order of Sqrt(segmentOffset_))
    size_t minTableSize = 256;
    TUnsigned tableSize =
      Max( size_t(std::sqrt(segmentOffset_)), minTableSize );
    segmentTable_.resize( tableSize );

    // Ensure that every prime factor of a number in the range covered by
    // the table, [segmentOffset_,segmentOffset_+2*newTableSize), 
    // is in our list of primes
    TUnsigned largestCandidate = segmentOffset_ + 2*(tableSize-1);
    while( oddPrimes.back()*oddPrimes.back() < largestCandidate )
    {
        AugmentPrimes( 2*oddPrimes.size() );
    }

    SieveSegment();

    if( keepAll )
    {
        Generate( lowerBound );
        lowerBound_ = lowerBound;
    }
}

template<typename TUnsigned>
void DynamicSieve<TUnsigned>::SetStorage( bool keepAll )
{
    if( keepAll && !keepAll_ )
    {
        // NOTE: This could be made substantially faster but is provided as
        //       a convenience. Advanced users should set this mode upon 
        //       construction.
        DynamicSieve<TUnsigned> tmpSieve( oddPrimes.back()+1 );
        tmpSieve.oddPrimes = oddPrimes;

        TUnsigned pntEstimate =
          segmentOffset_ / TUnsigned(Log2(double(segmentOffset_)));
        TUnsigned fudge = 10;
        oddPrimes.reserve( pntEstimate + fudge );

        while( tmpSieve.oddPrimes.back() < segmentOffset_ )
            oddPrimes.push_back( tmpSieve.NextPrime() );
    }
    keepAll_ = keepAll;
}

template<typename TUnsigned>
void DynamicSieve<TUnsigned>::Generate( TUnsigned upperBound )
{
    while( oddPrimes.back() < upperBound )
        NextPrime();
}

template<typename TUnsigned>
bool DynamicSieve<TUnsigned>::SeekSegmentPrime()
{
    // Set segmentIndex_ to the first location greater than or equal to the
    // current location that is marked prime
    const TUnsigned tableSize = segmentTable_.size();
    const TUnsigned start = 
      ( segmentOffset_ >= lowerBound_ ?
        0 :
        (lowerBound_-segmentOffset_+1) / 2 );
    for( TUnsigned i=start; i<tableSize; ++i )
    {
        if( segmentTable_[i] != 0 )
        {
            segmentIndex_ = i;
            return true;
        }
    }
    return false;
}

// Find the integer k such that segmentOffset_ + 2*k is the first odd 
// multiple of p that is greater than or equal to Max(segmentOffset_,p^2)
template<typename TUnsigned>
TUnsigned DynamicSieve<TUnsigned>::ComputeSegmentStart
( const TUnsigned& p ) const
{
    if( p*p >= segmentOffset_ )
    {
        return (p*p - segmentOffset_) / 2;
    }
    else
    {
        TUnsigned complement = (p - (segmentOffset_ % p)) % p;
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
void DynamicSieve<TUnsigned>::AugmentPrimes( TUnsigned numPrimes )
{
    TUnsigned pIndex = oddPrimes.size();
    oddPrimes.resize( numPrimes );

    // TODO: Consider switching to another scheme if a sufficient number of
    //       primes are to be kept
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
    }
}

template<typename TUnsigned>
void DynamicSieve<TUnsigned>::SieveSegment()
{
    std::fill( segmentTable_.begin(), segmentTable_.end(), byte(1) );

    TUnsigned tableSize = segmentTable_.size();
    TUnsigned numPrimes = oddPrimes.size();

    TUnsigned largestCandidate = segmentOffset_ + 2*(tableSize-1);
    TUnsigned factorBound = TUnsigned(std::sqrt(largestCandidate)) + 1;
    auto boundIter =
      std::lower_bound( oddPrimes.begin(), oddPrimes.end(), factorBound );
    const TUnsigned indexEnd = TUnsigned(boundIter-oddPrimes.begin());

    for( TUnsigned j=0; j<indexEnd; ++j )
    {
        TUnsigned p = oddPrimes[j];

        // Iterate over the indices corresponding to odd multiples of p that
        // are at least as large as Max(segmentOffset_,p^2)
        TUnsigned start = ComputeSegmentStart( p );
        for( TUnsigned k=start; k<tableSize; k+=p )
            segmentTable_[k] = 0;
    }
}

template<typename TUnsigned>
void DynamicSieve<TUnsigned>::FormNewSegment()
{
    // The current table handles the odd integers in the inclusive interval
    // [segmentOffset,segmentOffset+2*(tableSize-1)], so the next table should
    // begin at 2*tableSize (unless lowerBound_ is larger).
    segmentOffset_ += 2*segmentTable_.size();
    segmentOffset_ = Max(segmentOffset_,lowerBound_);

    // Allocate a new table which would use roughly the same amount of memory
    // as the list of primes (which we crudely approximate to have a cardinality
    // on the order of Sqrt(segmentOffset_))
    const TUnsigned newTableSize =
      Max( size_t(std::sqrt(segmentOffset_)), segmentTable_.size() );
    segmentTable_.resize( newTableSize );

    // Ensure that every prime factor of a number in the range covered by the
    // table, [segmentOffset_,segmentOffset_+2*newTableSize), is in our list of
    //  primes
    const TUnsigned largestCandidate = segmentOffset_ + 2*(newTableSize-1);
    if( keepAll_ )
    {
        // We should have all of the primes below segmentOffset_ precomputed
        // TODO: Decide a sharp way to assert this
        TUnsigned pntEstimate =
          largestCandidate / TUnsigned(Log2(double(largestCandidate)));
        TUnsigned fudge = 10;
        oddPrimes.reserve( pntEstimate + fudge );
    }
    else
    {
        while( oddPrimes.back()*oddPrimes.back() < largestCandidate )
        {
            AugmentPrimes( 2*oddPrimes.size() );
        }
    }

    SieveSegment();
}

template<typename TUnsigned>
TUnsigned DynamicSieve<TUnsigned>::NextPrime()
{
    // Use a stored prime if a sufficiently large one exists
    if( oddPrimes.back() >= lowerBound_ )
    {
        auto iter =
          std::lower_bound( oddPrimes.begin(), oddPrimes.end(), lowerBound_ );
        TUnsigned currentPrime = *iter;
        lowerBound_ = currentPrime + 1;
        return currentPrime;
    }

    // Fall back to the segment table
    if( !SeekSegmentPrime() )
    {
        FormNewSegment();
        bool foundPrecomputed = SeekSegmentPrime();
        DEBUG_ONLY(
          if( !foundPrecomputed )
              LogicError("Error in DynamicSieve's table computation");
        )
    }

    // We have guaranteed that we are currently pointing to a prime
    TUnsigned currentPrime = segmentOffset_ + 2*segmentIndex_;
    if( keepAll_ && currentPrime > oddPrimes.back() )
        oddPrimes.push_back( currentPrime );
    DEBUG_ONLY(
        if( currentPrime < lowerBound_ )
            LogicError("Error in DynamicSieve's lower bound mechanism");
    )
    lowerBound_ = currentPrime + 1;
    return currentPrime;
}

} // namespace El

#endif // ifndef EL_NUMBER_THEORY_DYNAMIC_SIEVE_HPP
