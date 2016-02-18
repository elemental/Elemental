/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

namespace svp {

// NOTE: This accepts the N from the Q(DN), not R from QR, where R=DN
template<typename F>
void CoordinatesToSparse( const Matrix<F>& N, const Matrix<F>& v, Matrix<F>& y )
{
    DEBUG_ONLY(CSE cse("svp::CoordinatesToSparse"))
    y = v;
    Trmv( UPPER, NORMAL, UNIT, N, y );
    Round( y );
}
template<typename F>
void TransposedCoordinatesToSparse
( const Matrix<F>& NTrans, const Matrix<F>& v, Matrix<F>& y )
{
    DEBUG_ONLY(CSE cse("svp::TransposedCoordinatesToSparse"))
    y = v;
    Trmv( LOWER, TRANSPOSE, UNIT, NTrans, y );
    Round( y );
}

template<typename F>
void BatchCoordinatesToSparse
( const Matrix<F>& N, const Matrix<F>& V, Matrix<F>& Y )
{
    DEBUG_ONLY(CSE cse("svp::BatchCoordinatesToSparse"))
    Y = V;
    Trmm( LEFT, UPPER, NORMAL, UNIT, F(1), N, Y );
    Round( Y );
}
template<typename F>
void BatchTransposedCoordinatesToSparse
( const Matrix<F>& NTrans, const Matrix<F>& V, Matrix<F>& Y )
{
    DEBUG_ONLY(CSE cse("svp::BatchTransposedCoordinatesToSparse"))
    Y = V;
    Trmm( LEFT, LOWER, TRANSPOSE, UNIT, F(1), NTrans, Y );
    Round( Y );
}

template<typename F>
void SparseToCoordinates( const Matrix<F>& N, const Matrix<F>& y, Matrix<F>& v )
{
    DEBUG_ONLY(CSE cse("svp::SparseToCoordinates"))
    const Int n = N.Height();

    v = y;

    // A custom rounded analogue of an upper-triangular solve
    const F* NBuf = N.LockedBuffer();
          F* vBuf = v.Buffer();
    const Int NLDim = N.LDim();
    for( Int j=n-1; j>=0; --j )
    {
        F tau = 0;
        for( Int k=j+1; k<n; ++k ) 
            tau += NBuf[j+k*NLDim]*vBuf[k];
        vBuf[j] -= Round(tau);
    }
}
template<typename F>
void TransposedSparseToCoordinates
( const Matrix<F>& NTrans, const Matrix<F>& y, Matrix<F>& v )
{
    DEBUG_ONLY(CSE cse("svp::TransposedSparseToCoordinates"))
    const Int n = NTrans.Height();

    v = y;

    // A custom rounded analogue of an upper-triangular solve
    const F* NTransBuf = NTrans.LockedBuffer();
          F* vBuf = v.Buffer();
    const Int NTransLDim = NTrans.LDim();
    for( Int j=n-1; j>=0; --j )
    {
        const F* nBuf = &NTransBuf[j*NTransLDim];

        F tau = 0;
        for( Int k=j+1; k<n; ++k ) 
            tau += nBuf[k]*vBuf[k];
        vBuf[j] -= Round(tau);
    }
}

// TODO: Optimize this routine by changing the loop order?
template<typename F>
void BatchSparseToCoordinatesUnblocked
( const Matrix<F>& N, const Matrix<F>& Y, Matrix<F>& V )
{
    DEBUG_ONLY(CSE cse("svp::BatchSparseToCoordinatesUnblocked"))
    const Int n = N.Height();
    const Int numRHS = Y.Width();

    // A custom rounded analogue of an upper-triangular solve
    const F* NBuf = N.LockedBuffer();
    const F* YBuf = Y.LockedBuffer();
          F* VBuf = V.Buffer();
    const Int NLDim = N.LDim();
    const Int YLDim = Y.LDim();
    const Int VLDim = V.LDim();
    for( Int l=0; l<numRHS; ++l )
    {
              F* vBuf = &VBuf[l*VLDim];
        const F* yBuf = &YBuf[l*YLDim];
        for( Int j=n-1; j>=0; --j )
        {
            F tau = 0;
            for( Int k=j+1; k<n; ++k ) 
                tau += NBuf[j+k*NLDim]*vBuf[k];
            vBuf[j] = yBuf[j] + Round(vBuf[j]-tau); 
        }
    }
}
template<typename F>
void BatchTransposedSparseToCoordinatesUnblocked
( const Matrix<F>& NTrans, const Matrix<F>& Y, Matrix<F>& V )
{
    DEBUG_ONLY(CSE cse("svp::BatchTransposedSparseToCoordinatesUnblocked"))
    const Int n = NTrans.Height();
    const Int numRHS = Y.Width();

    // A custom rounded analogue of an upper-triangular solve
    const F* NTransBuf = NTrans.LockedBuffer();
    const F* YBuf = Y.LockedBuffer();
          F* VBuf = V.Buffer();
    const Int NTransLDim = NTrans.LDim();
    const Int YLDim = Y.LDim();
    const Int VLDim = V.LDim();
    for( Int l=0; l<numRHS; ++l )
    {
              F* vBuf = &VBuf[l*VLDim];
        const F* yBuf = &YBuf[l*YLDim];
        for( Int j=n-1; j>=0; --j )
        {
            const F* nBuf = &NTransBuf[j*NTransLDim];

            F tau = 0;
            for( Int k=j+1; k<n; ++k ) 
                tau += nBuf[k]*vBuf[k];
            vBuf[j] = yBuf[j] + Round(vBuf[j]-tau); 
        }
    }
}

// TODO: Optimize this routine
template<typename F>
void BatchSparseToCoordinates
( const Matrix<F>& N, const Matrix<F>& Y, Matrix<F>& V, Int blocksize )
{
    DEBUG_ONLY(CSE cse("svp::BatchSparseToCoordinates"))
    const Int n = N.Height();
    const Int numRHS = Y.Width();

    // TODO: Test the relative performance of this branch
    if( numRHS == 1 )
    {
        SparseToCoordinates( N, Y, V );
        return;
    }

    Zeros( V, n, numRHS );
    for( Int i=0; i<n; i+=blocksize )
    { 
        const Int nb = Min(n-i,blocksize);
        const Range<Int> ind0(0,i), ind1(i,i+nb);
        
        auto N01 = N( ind0, ind1 );
        auto N11 = N( ind1, ind1 );
        auto Y1 = Y( ind1, ALL );
        auto V0 = V( ind0, ALL );
        auto V1 = V( ind1, ALL );

        BatchSparseToCoordinatesUnblocked( N11, Y1, V1 );
        if( i+nb < n )
            Gemm( NORMAL, NORMAL, F(-1), N01, V1, F(1), V0 );
    }
}
template<typename F>
void BatchTransposedSparseToCoordinates
( const Matrix<F>& NTrans, const Matrix<F>& Y, Matrix<F>& V, Int blocksize )
{
    DEBUG_ONLY(CSE cse("svp::BatchTransposedSparseToCoordinates"))
    const Int n = NTrans.Height();
    const Int numRHS = Y.Width();

    // TODO: Test the relative performance of this branch
    if( numRHS == 1 )
    {
        TransposedSparseToCoordinates( NTrans, Y, V );
        return;
    }

    Zeros( V, n, numRHS );
    for( Int i=0; i<n; i+=blocksize )
    { 
        const Int nb = Min(n-i,blocksize);
        const Range<Int> ind0(0,i), ind1(i,i+nb);
        
        auto NTrans10 = NTrans( ind1, ind0 );
        auto NTrans11 = NTrans( ind1, ind1 );
        auto Y1 = Y( ind1, ALL );
        auto V0 = V( ind0, ALL );
        auto V1 = V( ind1, ALL );

        BatchTransposedSparseToCoordinatesUnblocked( NTrans11, Y1, V1 );
        if( i+nb < n )
            Gemm( TRANSPOSE, NORMAL, F(-1), NTrans10, V1, F(1), V0 );
    }
}

template<typename F>
Base<F> CoordinatesToNorm
( const Matrix<Base<F>>& d, const Matrix<F>& N, const Matrix<F>& v )
{
    DEBUG_ONLY(CSE cse("svp::CoordinatesToNorm"))
    Matrix<F> z( v );
    Trmv( UPPER, NORMAL, UNIT, N, z );
    DiagonalScale( LEFT, NORMAL, d, z );
    return FrobeniusNorm( z );
}
template<typename F>
Base<F> TransposedCoordinatesToNorm
( const Matrix<Base<F>>& d, const Matrix<F>& NTrans, const Matrix<F>& v )
{
    DEBUG_ONLY(CSE cse("svp::TransposedCoordinatesToNorm"))
    Matrix<F> z( v );
    Trmv( LOWER, TRANSPOSE, UNIT, NTrans, z );
    DiagonalScale( LEFT, NORMAL, d, z );
    return FrobeniusNorm( z );
}

template<typename F>
Matrix<Base<F>> BatchCoordinatesToNorms
( const Matrix<Base<F>>& d, const Matrix<F>& N, const Matrix<F>& V )
{
    DEBUG_ONLY(CSE cse("svp::BatchCoordinatesToNorms"))
    typedef Base<F> Real;
    Matrix<F> Z( V );
    // TODO: Decide whether this branch is necessary or not...
    if( V.Width() == 1 )
        Trmv( UPPER, NORMAL, UNIT, N, Z );
    else
        Trmm( LEFT, UPPER, NORMAL, UNIT, F(1), N, Z );
    DiagonalScale( LEFT, NORMAL, d, Z );

    Matrix<Real> colNorms;
    ColumnTwoNorms( Z, colNorms );
    return colNorms;
}
template<typename F>
Matrix<Base<F>> BatchTransposedCoordinatesToNorms
( const Matrix<Base<F>>& d, const Matrix<F>& NTrans, const Matrix<F>& V )
{
    DEBUG_ONLY(CSE cse("svp::BatchTransposedCoordinatesToNorms"))
    typedef Base<F> Real;
    Matrix<F> Z( V );
    // TODO: Decide whether this branch is necessary or not...
    if( V.Width() == 1 )
        Trmv( LOWER, TRANSPOSE, UNIT, NTrans, Z );
    else
        Trmm( LEFT, LOWER, TRANSPOSE, UNIT, F(1), NTrans, Z );
    DiagonalScale( LEFT, NORMAL, d, Z );

    Matrix<Real> colNorms;
    ColumnTwoNorms( Z, colNorms );
    return colNorms;
}

template<typename F>
Base<F> SparseToNormLowerBound
( const Matrix<Base<F>>& d, const Matrix<F>& y )
{
    DEBUG_ONLY(CSE cse("svp::SparseToNormLowerBound"))
    typedef Base<F> Real;
    const Int n = d.Height();
    const Real* dBuf = d.LockedBuffer();
    const F* yBuf = y.LockedBuffer();

    const Real oneHalf = Real(1)/Real(2);

    Real lowerBoundSquared = 0;
    for( Int j=0; j<n; ++j )
    {
        if( yBuf[j] != F(0) )
        {
            const Real arg = (Abs(yBuf[j])-oneHalf)*dBuf[j];
            lowerBoundSquared += Pow(arg,Real(2));
        }
    }
    return Sqrt(lowerBoundSquared);
}

template<typename F>
Matrix<Base<F>> BatchSparseToNormLowerBound
( const Matrix<Base<F>>& d, const Matrix<F>& Y )
{
    DEBUG_ONLY(CSE cse("svp::BatchSparseToNormLowerBound"))
    typedef Base<F> Real;
    const Int numRHS = Y.Width();
    Matrix<Real> normBounds;
    Zeros( normBounds, numRHS, 1 );
    for( Int j=0; j<numRHS; ++j )
        normBounds.Set( j, 0, SparseToNormLowerBound(d,Y(ALL,IR(j))) );
    return normBounds;
}

template<typename F>
Base<F> SparseToNorm
( const Matrix<Base<F>>& d, const Matrix<F>& N, const Matrix<F>& y )
{
    DEBUG_ONLY(CSE cse("svp::SparseToNorm"))
    Matrix<F> v;
    SparseToCoordinates( N, y, v );
    return CoordinatesToNorm( d, N, v );
}
template<typename F>
Base<F> TransposedSparseToNorm
( const Matrix<Base<F>>& d, const Matrix<F>& NTrans, const Matrix<F>& y )
{
    DEBUG_ONLY(CSE cse("svp::TransposedSparseToNorm"))
    Matrix<F> v;
    TransposedSparseToCoordinates( NTrans, y, v );
    return TransposedCoordinatesToNorm( d, NTrans, v );
}

template<typename F>
Matrix<Base<F>> BatchSparseToNorm
( const Matrix<Base<F>>& d,
  const Matrix<F>& N,
  const Matrix<F>& Y,
        Int blocksize )
{
    DEBUG_ONLY(CSE cse("svp::BatchSparseToNorm"))
    Matrix<F> V;
    BatchSparseToCoordinates( N, Y, V, blocksize );
    return BatchCoordinatesToNorms( d, N, V );
}
template<typename F>
Matrix<Base<F>> BatchTransposedSparseToNorm
( const Matrix<Base<F>>& d,
  const Matrix<F>& NTrans,
  const Matrix<F>& Y,
        Int blocksize )
{
    DEBUG_ONLY(CSE cse("svp::BatchTransposedSparseToNorm"))
    Matrix<F> V;
    BatchTransposedSparseToCoordinates( NTrans, Y, V, blocksize );
    return BatchTransposedCoordinatesToNorms( d, NTrans, V );
}

template<typename Real>
class PhaseEnumerationCache
{
private:
    const Matrix<Real>& B_;
    const Matrix<Real>& d_;
    const Matrix<Real>& N_;
          Matrix<Real> NTrans_;
    Real normUpperBound_;
    bool foundVector_=false;
    bool useTranspose_=true;

    Int numQueued_=0;
    Matrix<Real> Y_;
    Matrix<Real> VCand_;
    Matrix<Real> v_;
    Int blocksize_=32;

public:
    PhaseEnumerationCache
    ( const Matrix<Real>& B,
      const Matrix<Real>& d,
      const Matrix<Real>& N,
      Real normUpperBound, 
      Int batchSize=1000,
      Int blocksize=16,
      bool useTranspose=true )
    : B_(B),
      d_(d),
      N_(N),
      normUpperBound_(normUpperBound),
      foundVector_(false),
      numQueued_(0),
      blocksize_(blocksize),
      useTranspose_(useTranspose)
    { 
        Zeros( Y_, N.Height(), batchSize );   
        if( useTranspose )
            Transpose( N, NTrans_ );
    }

    bool FoundVector() const { return foundVector_; }
    const Matrix<Real>& BestVector() const { return v_; }
    Real NormUpperBound() const { return normUpperBound_; }

    void Flush()
    {
        if( numQueued_ == 0 )
            return;

        auto YActive = Y_( ALL, IR(0,numQueued_) );

        Matrix<Real> colNorms;
        if( useTranspose_ )
        {
            BatchTransposedSparseToCoordinates
            ( NTrans_, YActive, VCand_, blocksize_ );
            colNorms = BatchTransposedCoordinatesToNorms( d_, NTrans_, VCand_ );
        }
        else
        {
            BatchSparseToCoordinates( N_, YActive, VCand_, blocksize_ );
            colNorms = BatchCoordinatesToNorms( d_, N_, VCand_ );
        }

        for( Int j=0; j<numQueued_; ++j )
        {
            const Real bNorm = colNorms.Get(j,0);
            if( bNorm < normUpperBound_ && bNorm != Real(0) )
            {
                auto y = YActive(ALL,IR(j));
                auto vCand = VCand_(ALL,IR(j));

                Output("normUpperBound=",normUpperBound_,", bNorm=",bNorm);
                Print( y, "y" );
                Print( d_, "d" );
                Print( N_, "N" );

                // Check that the reverse transformation holds
                Matrix<Real> yCheck;
                CoordinatesToSparse( N_, vCand, yCheck );
                yCheck -= y;
                if( FrobeniusNorm(yCheck) != Real(0) )
                {
                    Print( vCand, "vCand" );
                    Print( yCheck, "eCheck" );
                    LogicError("Invalid sparse transformation");
                }

                Copy( vCand, v_ );

                Print( v_, "v" );
                Print( B_, "B" );

                Matrix<Real> b;
                Zeros( b, B_.Height(), 1 );
                Gemv( NORMAL, Real(1), B_, v_, Real(0), b );
                Print( b, "b" );

                normUpperBound_ = bNorm;
                foundVector_ = true;
            }
        }
        numQueued_ = 0;
    }

    void Enqueue( const Matrix<Real>& y )
    {
        MemCopy( Y_.Buffer(0,numQueued_), y.LockedBuffer(), y.Height() ); 

        ++numQueued_;
        if( numQueued_ == blocksize_ )
            Flush();
    }

    ~PhaseEnumerationCache() { }
};

template<typename Real>
void PhaseEnumerationBottom
( PhaseEnumerationCache<Real>& cache,
  Matrix<Real>& y,
  const Int& beg,
  const Int& minInf,
  const Int& maxInf,
  const Int& minOne,
  const Int& maxOne,
  const Int& baseOneNorm,
  const Int& baseInfNorm,
  const bool& zeroSoFar )
{
    DEBUG_ONLY(CSE cse("svp::PhaseEnumerationBottom"))
    const Int n = y.Height();
          Real* yBuf = y.Buffer();

    const Int bound = Min(maxInf,maxOne-baseOneNorm);
    for( Int absBeta=1; absBeta<=bound; ++absBeta )
    {
        const Int newOneNorm = baseOneNorm + absBeta;
        const Int newInfNorm = Max(baseInfNorm,absBeta);

        if( newOneNorm >= minOne && newInfNorm >= minInf )
        {
            // We can insert +-|beta| into any position and be admissible,
            // so enqueue each possibility
            for( Int i=beg; i<n; ++i )
            {
                yBuf[i] = absBeta;
                cache.Enqueue( y );

                if( !zeroSoFar )
                {
                    yBuf[i] = -absBeta;
                    cache.Enqueue( y );
                }

                yBuf[i] = 0;
            }
        }

        if( newOneNorm < maxOne )
        {
            // We can insert +-|beta| into any position and still have room
            // left for the one norm bound, so do so and then recurse
            for( Int i=beg; i<n-1; ++i )
            {
                yBuf[i] = absBeta;
                PhaseEnumerationBottom
                ( cache, y,
                  i+1, minInf, maxInf, minOne, maxOne,
                  newOneNorm, newInfNorm, false );

                if( !zeroSoFar )
                {
                    yBuf[i] = -absBeta;
                    PhaseEnumerationBottom
                    ( cache, y,
                      i+1, minInf, maxInf, minOne, maxOne,
                      newOneNorm, newInfNorm, false );
                }

                yBuf[i] = 0;
            }
        }
    }
}

template<typename Real>
void PhaseEnumerationInner
(       PhaseEnumerationCache<Real>& cache,
        Matrix<Real>& y,
  const Int& phaseLength,
  const Int& phase,
  const Int& beg,
  const Int& end,
  const vector<Int>& minInfNorms,
  const vector<Int>& maxInfNorms,
  const vector<Int>& minOneNorms,
  const vector<Int>& maxOneNorms,
  const Int& baseOneNorm,
  const Int& baseInfNorm,
  const bool& zeroSoFar,
  const Int& progressLevel )
{
    DEBUG_ONLY(CSE cse("svp::PhaseEnumerationInner"))
    const Int n = y.Height();
    Real* yBuf = y.Buffer();

    const Int phaseBeg = Max(end-phaseLength,0);
    const Int nextPhaseEnd = Min(end+phaseLength,n);

    if( end >= n )
    {
        const Int minInf = minInfNorms.back(); 
        const Int maxInf = maxInfNorms.back();
        const Int minOne = minOneNorms.back();
        const Int maxOne = maxOneNorms.back();
        // Enqueue the zero phase if it is admissible
        if( minInf == Int(0) && minOne == Int(0) )
            cache.Enqueue( y );

        PhaseEnumerationBottom
        ( cache, y, beg,
          minInf, maxInf, minOne, maxOne,
          Int(0), Int(0), zeroSoFar );
        return;
    }

    if( beg == phaseBeg && minInfNorms[phase] == 0 && minOneNorms[phase] == 0 )
    {
        // This phase can be all zeroes, so move to the next phase
        PhaseEnumerationInner
        ( cache, y, phaseLength,
          phase+1, end, nextPhaseEnd,
          minInfNorms, maxInfNorms,
          minOneNorms, maxOneNorms,
          Int(0), Int(0), zeroSoFar, progressLevel );
    }

    const Int bound = Min(maxInfNorms[phase],maxOneNorms[phase]-baseOneNorm);
    for( Int absBeta=1; absBeta<=bound; ++absBeta ) 
    {
        const Int phaseOneNorm = baseOneNorm + absBeta;
        const Int phaseInfNorm = Max( baseInfNorm, absBeta );

        if( phaseOneNorm >= minOneNorms[phase] &&
            phaseInfNorm >= minInfNorms[phase] )
        {
            // We could insert +-|beta| into any position and still have
            // an admissible choice for the current phase, so procede to
            // the next phase with each such choice

            for( Int i=beg; i<end; ++i )
            {
                // Fix y[i] = +|beta| and move to the next phase
                yBuf[i] = absBeta;
                if( phase < progressLevel )
                    Output("phase ",phase,": y[",i,"]=",absBeta);
                PhaseEnumerationInner
                ( cache, y, phaseLength,
                  phase+1, end, nextPhaseEnd,
                  minInfNorms, maxInfNorms, minOneNorms, maxOneNorms,
                  Int(0), Int(0), false, progressLevel );

                if( !zeroSoFar )
                {
                    // Fix y[i] = -|beta| and move to the next phase
                    yBuf[i] = -absBeta;
                    if( phase < progressLevel )
                        Output("phase ",phase,": y[",i,"]=",-absBeta);
                    PhaseEnumerationInner
                    ( cache, y, phaseLength,
                      phase+1, end, nextPhaseEnd,
                      minInfNorms, maxInfNorms, minOneNorms, maxOneNorms,
                      Int(0), Int(0), false, progressLevel );
                }
                
                yBuf[i] = 0;
            }
        }

        if( phaseOneNorm < maxOneNorms[phase] )
        {
            // Inserting +-|beta| into any position still leaves us with
            // room in the current phase

            for( Int i=beg; i<end; ++i )
            {
                // Fix y[i] = +|beta| and move to y[i+1]
                yBuf[i] = absBeta;
                PhaseEnumerationInner
                ( cache, y, phaseLength,
                  phase, i+1, end,
                  minInfNorms, maxInfNorms, minOneNorms, maxOneNorms,
                  phaseOneNorm, phaseInfNorm, false, progressLevel );

                if( !zeroSoFar )
                {
                    // Fix y[i] = -|beta| and move to y[i+1]
                    yBuf[i] = -absBeta;
                    PhaseEnumerationInner
                    ( cache, y, phaseLength,
                      phase, i+1, end,
                      minInfNorms, maxInfNorms, minOneNorms, maxOneNorms,
                      phaseOneNorm, phaseInfNorm, false, progressLevel );
                }
               
                yBuf[i] = 0;
            }
        }
    }
}

template<typename Real>
Real PhaseEnumeration
( const Matrix<Real>& B,
  const Matrix<Real>& d,
  const Matrix<Real>& N,
        Real normUpperBound,
        Int startIndex,
        Int phaseLength,
  const vector<Int>& maxInfNorms,
  const vector<Int>& maxOneNorms,
        Matrix<Real>& v,
        Int progressLevel )
{
    DEBUG_ONLY(CSE cse("svp::PhaseEnumeration"))
    const Int n = N.Height();
    if( n <= 1 )
        return 2*normUpperBound + 1;

    // TODO: Make starting index modifiable
    const Int numPhases = ((n-startIndex)+phaseLength-1)/phaseLength;
    if( numPhases != Int(maxInfNorms.size()) )
        LogicError("Invalid length of maxInfNorms");
    if( numPhases != Int(maxOneNorms.size()) )
        LogicError("Invalid length of maxOneNorms");

    // TODO: Make these values modifiable
    vector<Int> minInfNorms(numPhases,0), minOneNorms(numPhases,0);

    // TODO: Loop and increase bands for min and max one and inf norms?

    // NOTE: This does not seem to have a substantial impact on performance...
    const Int batchSize = 2000;
    //const Int batchSize = 1;
    const Int blocksize = 16;
    const bool useTranspose = true;
    PhaseEnumerationCache<Real>
      cache( B, d, N, normUpperBound, batchSize, blocksize, useTranspose );

    Int phase=0;
    Int beg=startIndex;
    Int end=Min(beg+phaseLength,n);
    Matrix<Real> y;
    Zeros( y, n, 1 );
    bool zeroSoFar = true;
    PhaseEnumerationInner
    ( cache, y, phaseLength, phase, beg, end,
      minInfNorms, maxInfNorms, minOneNorms, maxOneNorms,
      Int(0), Int(0), zeroSoFar, progressLevel );

    cache.Flush();

    if( cache.FoundVector() )
    {
        y = cache.BestVector();
        return cache.NormUpperBound();
    }
    else
    {
        return 2*normUpperBound+1;
    }
}

} // namespace svp

// TODO: Instantiate batched variants?
#define PROTO(F) \
  template void svp::CoordinatesToSparse \
  ( const Matrix<F>& N, const Matrix<F>& v, Matrix<F>& y ); \
  template void svp::TransposedCoordinatesToSparse \
  ( const Matrix<F>& NTrans, const Matrix<F>& v, Matrix<F>& y ); \
  template void svp::SparseToCoordinates \
  ( const Matrix<F>& N, const Matrix<F>& y, Matrix<F>& v ); \
  template void svp::TransposedSparseToCoordinates \
  ( const Matrix<F>& NTrans, const Matrix<F>& y, Matrix<F>& v ); \
  template Base<F> svp::CoordinatesToNorm \
  ( const Matrix<Base<F>>& d, const Matrix<F>& N, const Matrix<F>& v ); \
  template Base<F> svp::TransposedCoordinatesToNorm \
  ( const Matrix<Base<F>>& d, const Matrix<F>& NTrans, const Matrix<F>& v ); \
  template Base<F> svp::SparseToNorm \
  ( const Matrix<Base<F>>& d, const Matrix<F>& N, const Matrix<F>& y ); \
  template Base<F> svp::TransposedSparseToNorm \
  ( const Matrix<Base<F>>& d, const Matrix<F>& NTrans, const Matrix<F>& y ); \
  template Base<F> svp::PhaseEnumeration \
  ( const Matrix<F>& B, \
    const Matrix<Base<F>>& d, \
    const Matrix<F>& N, \
          Base<F> normUpperBound, \
          Int startIndex, \
          Int phaseLength, \
    const vector<Int>& maxInfNorms, \
    const vector<Int>& maxOneNorms, \
          Matrix<F>& v, \
          Int progressLevel );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
