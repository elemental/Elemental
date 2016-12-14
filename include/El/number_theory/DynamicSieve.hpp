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

template<typename T,typename TSmall>
DynamicSieve<T,TSmall>::DynamicSieve
( T lowerBound,
  bool keepAll,
  TSmall segmentSize )
{
    keepAll_ = false;
    oddPrimeBound_ = 0;

    // Ensure that the lower bound is odd
    if( lowerBound % 2 == 0 )
        ++lowerBound;

    // Ensure that the segment size is even
    if( segmentSize % 2 == 1 )
        ++segmentSize;

    // This could be an arbitrary number of the first several odd primes,
    // but testing up to 53 being a good default for trial division is common
    // folklore
    oddPrimes.resize( 15 );
    oddPrimes[0] = 3;
    oddPrimes[1] = 5;
    oddPrimes[2] = 7;
    oddPrimes[3] = 11;
    oddPrimes[4] = 13;
    oddPrimes[5] = 17;
    oddPrimes[6] = 19;
    oddPrimes[7] = 23;
    oddPrimes[8] = 29;
    oddPrimes[9] = 31;
    oddPrimes[10] = 37;
    oddPrimes[11] = 41;
    oddPrimes[12] = 43;
    oddPrimes[13] = 47;
    oddPrimes[14] = 53;

    segmentOffset_ = lowerBound_ = ( keepAll ? 3 : lowerBound );
    segmentOffset_ = std::max( segmentOffset_, oddPrimes.back()+2 );
    segmentSize_ = segmentSize;
    segmentTable_.resize( segmentSize );

    // Initialize the sieve table
    T largestCandidate = segmentOffset_ + 2*(segmentSize-1);
    while( oddPrimes.back()*oddPrimes.back() < largestCandidate )
    {
        AugmentPrimes( 2*oddPrimes.size() );
    }
    SieveSegment();

    if( keepAll )
    {
        // Since it has been requested that we start generating primes at
        // 'lowerBound', to keep all primes, we should now generate all
        // primes before 'lowerBound'
        Generate( lowerBound );
        lowerBound_ = lowerBound;
    }
}

template<typename T,typename TSmall>
T DynamicSieve<T,TSmall>::PrimeCountingEstimate( T n ) const
{ return 2*(n/T(Log(double(n)))); }

template<typename T,typename TSmall>
void DynamicSieve<T,TSmall>::SetLowerBound( T lowerBound )
{
    // Ensure that the lower bound is odd
    if( lowerBound % 2 == 0 )
        ++lowerBound;
    if( lowerBound != lowerBound_ )
    {
        lowerBound_ = lowerBound;
        // TODO: Decide if we always need to pay this price
        MoveSegmentOffset( lowerBound_ );
    }
}

template<typename T,typename TSmall>
void DynamicSieve<T,TSmall>::MoveSegmentOffset( T segmentOffset )
{
    // TODO: Decide if this should only be enabled in debug mode
    if( segmentOffset % 2 == 0 )
        LogicError("Cannot choose an even segment offset");

    if( segmentOffset_ != segmentOffset )
    {
        segmentOffset_ = segmentOffset;
        for( TSmall i=0; i<oddPrimeBound_; ++i )
            oddPrimeStarts_[i] = ComputeSegmentStart( oddPrimes[i] );
    }
}

template<typename T,typename TSmall>
void DynamicSieve<T,TSmall>::SetStorage( bool keepAll )
{
    if( keepAll && !keepAll_ )
    {
        // Ensure that we have all of the primes below lowerBound_ stored
        if( lowerBound_ > oddPrimes.back()+2 )
        {
            T lowerBoundSave = lowerBound_;
            lowerBound_ = oddPrimes.back()+2;
            MoveSegmentOffset( lowerBound_ );

            keepAll_ = true;
            while( lowerBound_ < lowerBoundSave )
                NextPrime(); 
        }
    }
    keepAll_ = keepAll;
}

template<typename T,typename TSmall>
void DynamicSieve<T,TSmall>::Generate( T upperBound )
{
    SetStorage( true );

    // Compute all of the needed "small" odd primes for testing candidates
    // up to a (loose) upper bound of upperBound + 2*(segmentSize-1)
    T largestCandidate = upperBound + 2*(segmentSize_-1);
    oddPrimes.reserve( PrimeCountingEstimate(largestCandidate) );
    while( oddPrimes.back()*oddPrimes.back() < largestCandidate )
    {
        AugmentPrimes( 2*oddPrimes.size() );
    }

    if( segmentOffset_ < oddPrimes.back()+2 )
    {
        MoveSegmentOffset( oddPrimes.back()+2 );
    }

    for( ; segmentOffset_<upperBound; segmentOffset_+=2*segmentSize_ )
    {
        SieveSegment();
        for( TSmall k=0; k<segmentSize_; ++k )
        {
            if( segmentTable_[k] )
            {
                T p = segmentOffset_ + 2*k;
                oddPrimes.push_back( p );
            }
        }
    }

    SetStorage( false );
}

template<typename T,typename TSmall>
bool DynamicSieve<T,TSmall>::SeekSegmentPrime()
{
    // Set segmentIndex_ to the first location greater than or equal to the
    // current location that is marked prime
    const TSmall start =
      ( segmentOffset_ >= lowerBound_ ?
        0 :
        (lowerBound_-segmentOffset_+1) / 2 );
    for( TSmall i=start; i<segmentSize_; ++i )
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
template<typename T,typename TSmall>
TSmall DynamicSieve<T,TSmall>::ComputeSegmentStart
( const TSmall& p ) const
{
    if( p*p >= segmentOffset_ )
    {
        return (p*p - segmentOffset_) / 2;
    }
    else
    {
        auto complement = (p - (segmentOffset_ % p)) % p;
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

template<typename T,typename TSmall>
void DynamicSieve<T,TSmall>::AugmentPrimes
( T numPrimes )
{
    T pIndex = oddPrimes.size();
    oddPrimes.resize( numPrimes );

    // TODO:
    // Consider switching to another scheme if a sufficient number of primes
    // are to be kept
    T p = oddPrimes[pIndex-1];
    for( ; pIndex<numPrimes; ++pIndex )
    {
        // Find the next prime after p using trial division.
        auto begIter = oddPrimes.begin();
        auto boundIter = begIter;
        while( true )
        {
            p += 2;

            T factorBound = T(std::sqrt(p)) + 1;
            auto endIter = begIter;
            endIter += pIndex;
            boundIter = std::lower_bound( boundIter, endIter, factorBound );
            const T indexEnd = T(boundIter-begIter);

            T remainder = 1; // this could be an arbitrary nonzero
            for( T j=0; j<indexEnd; ++j )
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

template<typename T,typename TSmall>
void DynamicSieve<T,TSmall>::SieveSegment()
{
    std::fill( segmentTable_.begin(), segmentTable_.end(), char(1) );

    // Since sqrt(n) is an upper-bound on the smallest prime factor of 
    // a number n, we can a priori compute which of our stored primes is
    // the last that we need to sieve with in this segment
    T largestCandidate = segmentOffset_ + 2*(segmentSize_-1);
    auto factorBound = TSmall(std::sqrt(largestCandidate)) + 1;
    auto boundIter =
      std::lower_bound( oddPrimes.begin(), oddPrimes.end(), factorBound );
    auto indexEnd = TSmall(boundIter-oddPrimes.begin());

    // Compute the starting index in this segment of all of the newly
    // incorporated small primes 
    while( oddPrimeBound_ < indexEnd )
    {
        oddPrimeStarts_.push_back
        ( ComputeSegmentStart(oddPrimes[oddPrimeBound_]) );
        ++oddPrimeBound_;
    }

    for( TSmall j=0; j<indexEnd; ++j )
    {
        TSmall p = oddPrimes[j];

        // Iterate over the indices corresponding to odd multiples of p that
        // are at least as large as Max(segmentOffset_,p^2)
        auto k = oddPrimeStarts_[j];
        for( ; k<segmentSize_; k+=p )
        {
            segmentTable_[k] = 0;
        }

        // Store the next starting index of the j'th small prime
        // (under the assumption that the next segment will begin where the
        // current one ends)
        oddPrimeStarts_[j] = k - segmentSize_;
    }
}

template<typename T,typename TSmall>
void DynamicSieve<T,TSmall>::FormNewSegment()
{
    // The current table handles the odd integers in the inclusive interval
    // [segmentOffset,segmentOffset+2*(segmentSize-1)], so the next table
    // should begin at 2*segmentSize (unless lowerBound_ is larger).
    segmentOffset_ += 2*segmentSize_;

    // Ensure that every prime factor of a number in the range covered by the
    // table, [segmentOffset_,segmentOffset_+2*segmentSize), is in our list of
    //  primes
    const T largestCandidate = segmentOffset_ + 2*(segmentSize_-1);
    if( keepAll_ )
    {
        oddPrimes.reserve( PrimeCountingEstimate(largestCandidate) );
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

template<typename T,typename TSmall>
T DynamicSieve<T,TSmall>::NextPrime()
{
    // Use a stored prime if a sufficiently large one exists
    if( oddPrimes.back() >= lowerBound_ )
    {
        auto iter =
          std::lower_bound( oddPrimes.begin(), oddPrimes.end(), lowerBound_ );
        T currentPrime = *iter;
        lowerBound_ = currentPrime + 2;
        return currentPrime;
    }

    // Fall back to the segment table
    while( !SeekSegmentPrime() )
    {
        FormNewSegment();
    }

    // We have guaranteed that we are currently pointing to a prime
    T currentPrime = segmentOffset_ + 2*segmentIndex_;
    if( keepAll_ && currentPrime > oddPrimes.back() )
    {
        oddPrimes.push_back( currentPrime );
    }
    EL_DEBUG_ONLY(
      if( currentPrime < lowerBound_ )
          Output
          ("Current prime, ",currentPrime,
           ", smaller than lower bound of ",lowerBound_);
    )
    lowerBound_ = currentPrime + 2;
    return currentPrime;
}

} // namespace El

#endif // ifndef EL_NUMBER_THEORY_DYNAMIC_SIEVE_HPP
