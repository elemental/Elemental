/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

namespace svp {

// NOTE: This accepts the N from the Q(DN), not R from QR, where R=DN
template<typename Field>
void CoordinatesToSparse
( const Matrix<Field>& N,
  const Matrix<Field>& v,
        Matrix<Field>& y )
{
    EL_DEBUG_CSE
    y = v;
    Trmv( UPPER, NORMAL, UNIT, N, y );
    Round( y );
}

template<typename Field>
void TransposedCoordinatesToSparse
( const Matrix<Field>& NTrans,
  const Matrix<Field>& v,
        Matrix<Field>& y )
{
    EL_DEBUG_CSE
    y = v;
    Trmv( LOWER, TRANSPOSE, UNIT, NTrans, y );
    Round( y );
}

template<typename Field>
void BatchCoordinatesToSparse
( const Matrix<Field>& N,
  const Matrix<Field>& V,
        Matrix<Field>& Y )
{
    EL_DEBUG_CSE
    Y = V;
    Trmm( LEFT, UPPER, NORMAL, UNIT, Field(1), N, Y );
    Round( Y );
}

template<typename Field>
void BatchTransposedCoordinatesToSparse
( const Matrix<Field>& NTrans,
  const Matrix<Field>& V,
        Matrix<Field>& Y )
{
    EL_DEBUG_CSE
    Y = V;
    Trmm( LEFT, LOWER, TRANSPOSE, UNIT, Field(1), NTrans, Y );
    Round( Y );
}

template<typename Field>
void SparseToCoordinates
( const Matrix<Field>& N,
  const Matrix<Field>& y,
        Matrix<Field>& v )
{
    EL_DEBUG_CSE
    const Int n = N.Height();

    v = y;

    // A custom rounded analogue of an upper-triangular solve
    for( Int j=n-1; j>=0; --j )
    {
        Field tau = 0;
        for( Int k=j+1; k<n; ++k ) 
            tau += N(j,k)*v(k);
        v(j) -= Round(tau);
    }
}

template<typename Field>
void TransposedSparseToCoordinates
( const Matrix<Field>& NTrans,
  const Matrix<Field>& y,
        Matrix<Field>& v )
{
    EL_DEBUG_CSE
    const Int n = NTrans.Height();

    v = y;

    // A custom rounded analogue of an upper-triangular solve
    for( Int j=n-1; j>=0; --j )
    {
        const Field* nBuf = &NTrans(0,j);

        Field tau = 0;
        for( Int k=j+1; k<n; ++k ) 
            tau += nBuf[k]*v(k);
        v(j) -= Round(tau);
    }
}

// TODO(poulson): Optimize this routine by changing the loop order?
template<typename Field>
void BatchSparseToCoordinatesUnblocked
( const Matrix<Field>& N,
  const Matrix<Field>& Y,
        Matrix<Field>& V )
{
    EL_DEBUG_CSE
    const Int n = N.Height();
    const Int numRHS = Y.Width();

    // A custom rounded analogue of an upper-triangular solve
    for( Int l=0; l<numRHS; ++l )
    {
              Field* vBuf = &V(0,l);
        const Field* yBuf = &Y(0,l);
        for( Int j=n-1; j>=0; --j )
        {
            Field tau = 0;
            for( Int k=j+1; k<n; ++k ) 
                tau += N(j,k)*vBuf[k];
            vBuf[j] = yBuf[j] + Round(vBuf[j]-tau); 
        }
    }
}

template<typename Field>
void BatchTransposedSparseToCoordinatesUnblocked
( const Matrix<Field>& NTrans,
  const Matrix<Field>& Y,
        Matrix<Field>& V )
{
    EL_DEBUG_CSE
    const Int n = NTrans.Height();
    const Int numRHS = Y.Width();

    // A custom rounded analogue of an upper-triangular solve
    for( Int l=0; l<numRHS; ++l )
    {
              Field* vBuf = &V(0,l);
        const Field* yBuf = &Y(0,l);
        for( Int j=n-1; j>=0; --j )
        {
            const Field* nBuf = &NTrans(0,j);

            Field tau = 0;
            for( Int k=j+1; k<n; ++k ) 
                tau += nBuf[k]*vBuf[k];
            vBuf[j] = yBuf[j] + Round(vBuf[j]-tau); 
        }
    }
}

// TODO(poulson): Optimize this routine
template<typename Field>
void BatchSparseToCoordinates
( const Matrix<Field>& N,
  const Matrix<Field>& Y,
        Matrix<Field>& V,
        Int blocksize )
{
    EL_DEBUG_CSE
    const Int n = N.Height();
    const Int numRHS = Y.Width();

    // TODO(poulson): Test the relative performance of this branch
    if( numRHS == 1 )
    {
        SparseToCoordinates( N, Y, V );
        return;
    }

    Zeros( V, n, numRHS );
    const Int iLast = LastOffset( n, blocksize );
    for( Int i=iLast; i>=0; i-=blocksize )
    { 
        const Int nb = Min(n-i,blocksize);
        const Range<Int> ind0(0,i), ind1(i,i+nb);
        
        auto N01 = N( ind0, ind1 );
        auto N11 = N( ind1, ind1 );
        auto Y1 = Y( ind1, ALL );
        auto V0 = V( ind0, ALL );
        auto V1 = V( ind1, ALL );

        BatchSparseToCoordinatesUnblocked( N11, Y1, V1 );
        Gemm( NORMAL, NORMAL, Field(-1), N01, V1, Field(1), V0 );
    }
}

template<typename Field>
void BatchTransposedSparseToCoordinates
( const Matrix<Field>& NTrans,
  const Matrix<Field>& Y,
        Matrix<Field>& V,
        Int blocksize )
{
    EL_DEBUG_CSE
    const Int n = NTrans.Height();
    const Int numRHS = Y.Width();

    // TODO(poulson): Test the relative performance of this branch
    if( numRHS == 1 )
    {
        TransposedSparseToCoordinates( NTrans, Y, V );
        return;
    }

    Zeros( V, n, numRHS );
    const Int iLast = LastOffset( n, blocksize );
    for( Int i=iLast; i>=0; i-=blocksize )
    { 
        const Int nb = Min(n-i,blocksize);
        const Range<Int> ind0(0,i), ind1(i,i+nb);
        
        auto NTrans10 = NTrans( ind1, ind0 );
        auto NTrans11 = NTrans( ind1, ind1 );
        auto Y1 = Y( ind1, ALL );
        auto V0 = V( ind0, ALL );
        auto V1 = V( ind1, ALL );

        BatchTransposedSparseToCoordinatesUnblocked( NTrans11, Y1, V1 );
        Gemm( TRANSPOSE, NORMAL, Field(-1), NTrans10, V1, Field(1), V0 );
    }
}

template<typename Field>
Base<Field> CoordinatesToNorm
( const Matrix<Base<Field>>& d,
  const Matrix<Field>& N,
  const Matrix<Field>& v )
{
    EL_DEBUG_CSE
    Matrix<Field> z( v );
    Trmv( UPPER, NORMAL, UNIT, N, z );
    DiagonalScale( LEFT, NORMAL, d, z );
    return FrobeniusNorm( z );
}

template<typename Field>
Base<Field> TransposedCoordinatesToNorm
( const Matrix<Base<Field>>& d,
  const Matrix<Field>& NTrans,
  const Matrix<Field>& v )
{
    EL_DEBUG_CSE
    Matrix<Field> z( v );
    Trmv( LOWER, TRANSPOSE, UNIT, NTrans, z );
    DiagonalScale( LEFT, NORMAL, d, z );
    return FrobeniusNorm( z );
}

template<typename Field>
Matrix<Base<Field>>
NestedColumnTwoNorms( const Matrix<Field>& Z, Int numNested=1 )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = Z.Height();
    const Int numRHS = Z.Width();

    Matrix<Real> colNorms(numRHS,numNested);
    // Compute nested norms in linear time
    for( Int j=0; j<numRHS; ++j )
    {
        const Field* zBuf = Z.LockedBuffer(0,j);
        Real scale=0, scaledSquare=1;          
        for( Int i=n-1; i>=numNested; --i )
            UpdateScaledSquare( zBuf[i], scale, scaledSquare );
        for( Int i=numNested-1; i>=0; --i )
        {
            UpdateScaledSquare( zBuf[i], scale, scaledSquare );
            colNorms(j,i) = scale*Sqrt(scaledSquare);
        }
    }
    return colNorms;
}

template<typename Field>
Matrix<Base<Field>> BatchCoordinatesToNorms
( const Matrix<Base<Field>>& d,
  const Matrix<Field>& N,
  const Matrix<Field>& V,
        Int numNested=1 )
{
    EL_DEBUG_CSE
    Matrix<Field> Z( V );
    // TODO(poulson): Decide whether this branch is necessary or not...
    if( V.Width() == 1 )
        Trmv( UPPER, NORMAL, UNIT, N, Z );
    else
        Trmm( LEFT, UPPER, NORMAL, UNIT, Field(1), N, Z );
    DiagonalScale( LEFT, NORMAL, d, Z );

    return NestedColumnTwoNorms( Z, numNested );
}

template<typename Field>
Matrix<Base<Field>> BatchTransposedCoordinatesToNorms
( const Matrix<Base<Field>>& d,
  const Matrix<Field>& NTrans,
  const Matrix<Field>& V,
        Int numNested=1 )
{
    EL_DEBUG_CSE
    Matrix<Field> Z( V );
    // TODO(poulson): Decide whether this branch is necessary or not...
    if( V.Width() == 1 )
        Trmv( LOWER, TRANSPOSE, UNIT, NTrans, Z );
    else
        Trmm( LEFT, LOWER, TRANSPOSE, UNIT, Field(1), NTrans, Z );
    DiagonalScale( LEFT, NORMAL, d, Z );

    return NestedColumnTwoNorms( Z, numNested );
}

template<typename Field>
Base<Field> SparseToNormLowerBound
( const Matrix<Base<Field>>& d,
  const Matrix<Field>& y )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = d.Height();
    const Real* dBuf = d.LockedBuffer();
    const Field* yBuf = y.LockedBuffer();

    const Real oneHalf = Real(1)/Real(2);

    Real lowerBoundSquared = 0;
    for( Int j=0; j<n; ++j )
    {
        if( yBuf[j] != Field(0) )
        {
            const Real arg = (Abs(yBuf[j])-oneHalf)*dBuf[j];
            lowerBoundSquared += Pow(arg,Real(2));
        }
    }
    return Sqrt(lowerBoundSquared);
}

template<typename Field>
Matrix<Base<Field>> BatchSparseToNormLowerBound
( const Matrix<Base<Field>>& d,
  const Matrix<Field>& Y )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int numRHS = Y.Width();
    Matrix<Real> normBounds;
    Zeros( normBounds, numRHS, 1 );
    for( Int j=0; j<numRHS; ++j )
        normBounds(j) = SparseToNormLowerBound(d,Y(ALL,IR(j)));
    return normBounds;
}

template<typename Field>
Base<Field> SparseToNorm
( const Matrix<Base<Field>>& d,
  const Matrix<Field>& N,
  const Matrix<Field>& y )
{
    EL_DEBUG_CSE
    Matrix<Field> v;
    SparseToCoordinates( N, y, v );
    return CoordinatesToNorm( d, N, v );
}

template<typename Field>
Base<Field> TransposedSparseToNorm
( const Matrix<Base<Field>>& d,
  const Matrix<Field>& NTrans,
  const Matrix<Field>& y )
{
    EL_DEBUG_CSE
    Matrix<Field> v;
    TransposedSparseToCoordinates( NTrans, y, v );
    return TransposedCoordinatesToNorm( d, NTrans, v );
}

template<typename Field>
Matrix<Base<Field>> BatchSparseToNorm
( const Matrix<Base<Field>>& d,
  const Matrix<Field>& N,
  const Matrix<Field>& Y,
        Int blocksize,
        Int numNested=1 )
{
    EL_DEBUG_CSE
    Matrix<Field> V;
    BatchSparseToCoordinates( N, Y, V, blocksize );
    return BatchCoordinatesToNorms( d, N, V, numNested );
}

template<typename Field>
Matrix<Base<Field>> BatchTransposedSparseToNorm
( const Matrix<Base<Field>>& d,
  const Matrix<Field>& NTrans,
  const Matrix<Field>& Y,
        Int blocksize,
        Int numNested=1 )
{
    EL_DEBUG_CSE
    Matrix<Field> V;
    BatchTransposedSparseToCoordinates( NTrans, Y, V, blocksize );
    return BatchTransposedCoordinatesToNorms( d, NTrans, V, numNested );
}

template<typename Field>
class PhaseEnumerationCache
{
private:
    const Matrix<Field>& B_;
    const Matrix<Base<Field>>& d_;
    const Matrix<Field>& N_;
          Matrix<Field> NTrans_;
          Matrix<Base<Field>> normUpperBounds_;
    bool foundVector_=false;

    Int numQueued_=0;
    Matrix<Field> Y_;
    Matrix<Field> VCand_;

    Int insertionBound_;
    Matrix<Field> v_;

    Int blocksize_=32;
    bool useTranspose_=true;

public:
    PhaseEnumerationCache
    ( const Matrix<Field>& B,
      const Matrix<Base<Field>>& d,
      const Matrix<Field>& N,
      const Matrix<Base<Field>>& normUpperBounds,
            Int batchSize=256,
            Int blocksize=32,
            bool useTranspose=true )
    : B_(B),
      d_(d),
      N_(N),
      normUpperBounds_(normUpperBounds),
      foundVector_(false),
      numQueued_(0),
      insertionBound_(normUpperBounds.Height()),
      blocksize_(blocksize),
      useTranspose_(useTranspose)
    { 
        Zeros( Y_, N.Height(), batchSize );   
        if( useTranspose )
            Transpose( N, NTrans_ );
    }

    Int Height() const { return N_.Height(); }

    bool FoundVector() const { return foundVector_; }
    const Matrix<Field>& BestVector() const { return v_; }

    Int InsertionIndex() const
    {
        if( foundVector_ )
        {
            return insertionBound_-1;
        }
        else
        {
            LogicError("Did not find a shorter vector");
            return -1;
        }
    }

    Base<Field> InsertionNorm() const
    {
        const Int insertionIndex = InsertionIndex();
        return normUpperBounds_(insertionIndex);
    }

    void Flush()
    {
        if( numQueued_ == 0 )
            return;

        auto YActive = Y_( ALL, IR(0,numQueued_) );

        Matrix<Base<Field>> colNorms;
        if( useTranspose_ )
        {
            // TODO(poulson): Add this as an option
            /*
            Timer timer;
            timer.Start();
            BatchTransposedSparseToCoordinates
            ( NTrans_, YActive, VCand_, blocksize_ );
            const double transformTime = timer.Stop();
            const double n = YActive.Height();
            const double transformGflops =
              double(numQueued_)*n*n/(1.e9*transformTime);
            Output
            (numQueued_," transforms: ",timer.Stop()," seconds (",
             transformGflops," GFlop/s");
            timer.Start();
            colNorms =
              BatchTransposedCoordinatesToNorms
              ( d_, NTrans_, VCand_, insertionBound_ );
            const double normTime = timer.Stop();
            const double normGflops = double(numQueued_)*n*n/(1.e9*normTime);
            Output
            (numQueued_," norms: ",timer.Stop()," seconds (",
             normGflops," GFlop/s");
            */

            BatchTransposedSparseToCoordinates
            ( NTrans_, YActive, VCand_, blocksize_ );
            colNorms =
              BatchTransposedCoordinatesToNorms
              ( d_, NTrans_, VCand_, insertionBound_ );
        }
        else
        {
            BatchSparseToCoordinates( N_, YActive, VCand_, blocksize_ );
            colNorms =
              BatchCoordinatesToNorms( d_, N_, VCand_, insertionBound_ );
        }

        for( Int j=0; j<numQueued_; ++j )
        {
            for( Int k=0; k<insertionBound_; ++k )
            {
                const Base<Field> bNorm = colNorms(j,k);
                if( bNorm < normUpperBounds_(k) && bNorm != Base<Field>(0) )
                {
                    const Range<Int> subInd(k,END);

                    auto y = YActive(subInd,IR(j));
                    auto vCand = VCand_(subInd,IR(j));

                    Output
                    ("normUpperBound=",normUpperBounds_(k),
                     ", bNorm=",bNorm,", k=",k);
                    Print( y, "y" );

                    // Check that the reverse transformation holds
                    Matrix<Field> yCheck;
                    CoordinatesToSparse( N_(subInd,subInd), vCand, yCheck );
                    yCheck -= y;
                    if( FrobeniusNorm(yCheck) != Base<Field>(0) )
                    {
                        Print( B_(ALL,subInd), "B" );
                        Print( d_(subInd,ALL), "d" );
                        Print( N_(subInd,subInd), "N" );
                        Print( vCand, "vCand" );
                        Print( yCheck, "eCheck" );
                        LogicError("Invalid sparse transformation");
                    }

                    Copy( vCand, v_ );
                    Print( v_, "v" );

                    Matrix<Field> b;
                    Zeros( b, B_.Height(), 1 );
                    Gemv( NORMAL, Field(1), B_(ALL,subInd), v_, Field(0), b );
                    Print( b, "b" );

                    normUpperBounds_(k) = bNorm;
                    foundVector_ = true;
                    insertionBound_ = k+1;
                }
                // TODO(poulson): Keep track of 'stock' vectors?
            }
        }
        numQueued_ = 0;
        Zero( Y_ );
    }

    void Enqueue( const Matrix<Field>& y )
    {
        MemCopy( Y_.Buffer(0,numQueued_), y.LockedBuffer(), y.Height() ); 

        ++numQueued_;
        if( numQueued_ == Y_.Width() )
            Flush();
    }

    void Enqueue( const vector<pair<Int,Field>>& y )
    {
        Field* yBuf = Y_.Buffer(0,numQueued_);

        const Int numEntries = y.size();
        for( Int e=0; e<numEntries; ++e )
            yBuf[y[e].first] = y[e].second;

        ++numQueued_;
        if( numQueued_ == Y_.Width() )
            Flush();
    }

    void MaybeEnqueue( const vector<pair<Int,Field>>& y, double enqueueProb=1. )
    {
        if( enqueueProb >= 1. || SampleUniform<double>(0,1) <= enqueueProb )
            Enqueue( y );
    }

    ~PhaseEnumerationCache() { }
};

struct PhaseEnumerationCtrl
{
  const vector<Int>& phaseOffsets;
  const vector<Int>& minInfNorms;
  const vector<Int>& maxInfNorms;
  const vector<Int>& minOneNorms;
  const vector<Int>& maxOneNorms;
  const double enqueueProb=1.;
  const bool earlyExit=false;
  const Int progressLevel=4;

  PhaseEnumerationCtrl
  ( const vector<Int>& offsets,
    const vector<Int>& minInf,
    const vector<Int>& maxInf,
    const vector<Int>& minOne,
    const vector<Int>& maxOne,
    double accept=1.,
    bool early=false,
    Int level=4 )
  : phaseOffsets(offsets),
    minInfNorms(minInf), maxInfNorms(maxInf),
    minOneNorms(minOne), maxOneNorms(maxOne),
    enqueueProb(accept),
    earlyExit(early),
    progressLevel(level)
  { }
};

template<typename Field>
void PhaseEnumerationLeafInner
(       PhaseEnumerationCache<Field>& cache,
  const PhaseEnumerationCtrl& ctrl,
        vector<pair<Int,Field>>& y,
  const Int beg,
  const Int baseInf,
  const Int baseOne )
{
    EL_DEBUG_CSE
    const Int n = ctrl.phaseOffsets.back();
    const Int minInf = ctrl.minInfNorms.back();
    const Int maxInf = ctrl.maxInfNorms.back();
    const Int minOne = ctrl.minOneNorms.back();
    const Int maxOne = ctrl.maxOneNorms.back();
    const bool constrained = (y.size() == 0);

    SpiralState<Field> spiral;
    spiral.Initialize( constrained );
    while( true )
    {
        const Field beta = spiral.Step();
        const Int betaInf = Int(MaxAbs(beta));
        const Int betaOne = Int(OneAbs(beta));
        const Int newInf = Max(betaInf,baseInf);
        const Int newOne = betaOne + baseOne;
        if( newInf > maxInf || newOne > maxOne )
            break;

        if( newOne >= minOne && newInf >= minInf )
        {
            for( Int i=beg; i<n; ++i )
            {
                if( ctrl.enqueueProb >= 1. ||
                    SampleUniform<double>(0,1) <= ctrl.enqueueProb )
                {
                    y.emplace_back( i, beta );
                    cache.Enqueue( y );
                    y.pop_back();
                }
            }
        }

        if( newOne < maxOne )
        {
            // We can insert beta into any position and still have room
            // left for the one norm bound, so do so and then recurse
            for( Int i=beg; i<n-1; ++i )
            {
                y.emplace_back( i, beta );
                PhaseEnumerationLeafInner
                ( cache, ctrl, y, i+1, newInf, newOne );
                y.pop_back();
            }
        }
    }
}

template<typename Field>
void PhaseEnumerationLeaf
(       PhaseEnumerationCache<Field>& cache,
  const PhaseEnumerationCtrl& ctrl,
        vector<pair<Int,Field>>& y )
{
    EL_DEBUG_CSE

    // Enqueue the zero phase if it is admissible
    if( ctrl.minInfNorms.back() == Int(0) &&
        ctrl.minOneNorms.back() == Int(0) )
        cache.MaybeEnqueue( y, ctrl.enqueueProb );

    const Int beg = ctrl.phaseOffsets[ctrl.phaseOffsets.size()-2];
    const Int baseInf = 0;
    const Int baseOne = 0;
    PhaseEnumerationLeafInner( cache, ctrl, y, beg, baseInf, baseOne );
}

template<typename Field>
void PhaseEnumerationNodeInner
(       PhaseEnumerationCache<Field>& cache,
  const PhaseEnumerationCtrl& ctrl,
        vector<pair<Int,Field>>& y,
  const Int phase,
  const Int beg,
  const Int baseInf,
  const Int baseOne )
{
    EL_DEBUG_CSE
    const Int n = cache.Height();
    if( ctrl.earlyExit && cache.FoundVector() )
        return;

    const Int phaseBeg = ctrl.phaseOffsets[phase];
    const Int nextPhaseBeg = ctrl.phaseOffsets[phase+1];
    const Int minInf = ctrl.minInfNorms[phase];
    const Int maxInf = ctrl.maxInfNorms[phase];
    const Int minOne = ctrl.minOneNorms[phase];
    const Int maxOne = ctrl.maxOneNorms[phase];
    const bool constrained = (y.size() == 0);

    if( nextPhaseBeg >= n )
    {
        PhaseEnumerationLeaf( cache, ctrl, y );
        return;
    }

    if( beg == phaseBeg && minInf == 0 && minOne == 0 )
    {
        // This phase can be all zeroes, so move to the next phase
        PhaseEnumerationNode( cache, ctrl, y, phase+1 );
    }

    SpiralState<Field> spiral;
    spiral.Initialize( constrained );
    while( true )
    {
        const Field beta = spiral.Step();
        const Int betaInf = Int(MaxAbs(beta));
        const Int betaOne = Int(OneAbs(beta));
        const Int phaseInf = Max( betaInf, baseInf );
        const Int phaseOne = betaOne + baseOne;
        if( phaseOne > maxOne || phaseInf > maxInf )
            break;

        if( phaseOne >= minOne && phaseInf >= minInf )
        {
            // We could insert beta into any position and still have
            // an admissible choice for the current phase, so procede to
            // the next phase with each such choice

            for( Int i=beg; i<nextPhaseBeg; ++i )
            {
                // Fix y[i] = beta and move to the next phase
                y.emplace_back( i, beta );
                if( phase < ctrl.progressLevel )
                    Output("phase ",phase,": y[",i,"]=",beta);
                PhaseEnumerationNode( cache, ctrl, y, phase+1 );
                y.pop_back();
            }
        }

        if( phaseOne < maxOne )
        {
            // Inserting beta into any position still leaves us with
            // room in the current phase

            for( Int i=beg; i<nextPhaseBeg; ++i )
            {
                // Fix y[i] = beta and move to y[i+1]
                y.emplace_back( i, beta ); 
                PhaseEnumerationNodeInner
                ( cache, ctrl, y, 
                  phase, i+1,
                  phaseInf, phaseOne );
                y.pop_back();
            }
        }
    }
}

template<typename Field>
void PhaseEnumerationNode
(       PhaseEnumerationCache<Field>& cache,
  const PhaseEnumerationCtrl& ctrl,
        vector<pair<Int,Field>>& y,
  const Int phase )
{
    EL_DEBUG_CSE
    const Int baseInfNorm = 0;
    const Int baseOneNorm = 0;
    PhaseEnumerationNodeInner
    ( cache, ctrl, y,
      phase, ctrl.phaseOffsets[phase],
      baseInfNorm, baseOneNorm ); 
}

template<typename Field>
pair<Base<Field>,Int>
PhaseEnumeration
( const Matrix<Field>& B,
  const Matrix<Base<Field>>& d,
  const Matrix<Field>& N,
  const Matrix<Base<Field>>& normUpperBounds,
        Int startIndex,
        Int phaseLength,
        double enqueueProb,
  const vector<Int>& minInfNorms,
  const vector<Int>& maxInfNorms,
  const vector<Int>& minOneNorms,
  const vector<Int>& maxOneNorms,
        Matrix<Field>& v,
        Int progressLevel )
{
    EL_DEBUG_CSE
    const Int n = N.Height();
    if( n <= 1 )
        return pair<Base<Field>,Int>(2*normUpperBounds(0)+1,0);

    // TODO(poulson): Make starting index modifiable
    const Int numPhases = ((n-startIndex)+phaseLength-1)/phaseLength;
    if( numPhases != Int(maxInfNorms.size()) )
        LogicError("Invalid length of maxInfNorms");
    if( numPhases != Int(maxOneNorms.size()) )
        LogicError("Invalid length of maxOneNorms");

    // TODO(poulson): Loop and increase bands for min and max one and inf norms?

    // NOTE: The blocking doesn't seem to help the performance (yet)
    const Int batchSize = 512;
    const Int blocksize = 32;
    const bool useTranspose = true;
    PhaseEnumerationCache<Field>
      cache( B, d, N, normUpperBounds, batchSize, blocksize, useTranspose );

    const bool earlyExit = false;

    // TODO(poulson): Allow phaseOffsets to be an input
    vector<Int> phaseOffsets(numPhases+1);
    for( Int phase=0; phase<=numPhases; ++phase )
        phaseOffsets[phase] = Min(startIndex+phase*phaseLength,n);

    PhaseEnumerationCtrl
      ctrl
      (phaseOffsets,
       minInfNorms,maxInfNorms,
       minOneNorms,maxOneNorms,
       enqueueProb,
       earlyExit,progressLevel);

    vector<pair<Int,Field>> y;
    Int phase=0;
    PhaseEnumerationNode( cache, ctrl, y, phase );

    cache.Flush();

    if( cache.FoundVector() )
    {
        v = cache.BestVector();
        const Base<Field> insertionNorm = cache.InsertionNorm();
        const Int insertionIndex = cache.InsertionIndex();
        return pair<Base<Field>,Int>(insertionNorm,insertionIndex);
    }
    else
    {
        return pair<Base<Field>,Int>(2*normUpperBounds(0)+1,0);
    }
}

template<typename Field>
Base<Field> PhaseEnumeration
( const Matrix<Field>& B,
  const Matrix<Base<Field>>& d,
  const Matrix<Field>& N,
        Base<Field> normUpperBound,
        Int startIndex,
        Int phaseLength,
        double enqueueProb,
  const vector<Int>& minInfNorms,
  const vector<Int>& maxInfNorms,
  const vector<Int>& minOneNorms,
  const vector<Int>& maxOneNorms,
        Matrix<Field>& v,
        Int progressLevel )
{
    EL_DEBUG_CSE
    Matrix<Base<Field>> normUpperBounds(1,1);
    normUpperBounds(0) = normUpperBound;
    auto pair = 
      PhaseEnumeration
      ( B, d, N, normUpperBounds,
        startIndex, phaseLength, enqueueProb,
        minInfNorms, maxInfNorms,
        minOneNorms, maxOneNorms,
        v, progressLevel );
    return pair.first;
}

} // namespace svp

// TODO(poulson): Instantiate batched variants?
#define PROTO(Field) \
  template void svp::CoordinatesToSparse \
  ( const Matrix<Field>& N, \
    const Matrix<Field>& v, \
          Matrix<Field>& y ); \
  template void svp::TransposedCoordinatesToSparse \
  ( const Matrix<Field>& NTrans, \
    const Matrix<Field>& v, \
          Matrix<Field>& y ); \
  template void svp::SparseToCoordinates \
  ( const Matrix<Field>& N, \
    const Matrix<Field>& y, \
          Matrix<Field>& v ); \
  template void svp::TransposedSparseToCoordinates \
  ( const Matrix<Field>& NTrans, \
    const Matrix<Field>& y, \
          Matrix<Field>& v ); \
  template Base<Field> svp::CoordinatesToNorm \
  ( const Matrix<Base<Field>>& d, \
    const Matrix<Field>& N, \
    const Matrix<Field>& v ); \
  template Base<Field> svp::TransposedCoordinatesToNorm \
  ( const Matrix<Base<Field>>& d, \
    const Matrix<Field>& NTrans, \
    const Matrix<Field>& v ); \
  template Base<Field> svp::SparseToNorm \
  ( const Matrix<Base<Field>>& d, \
    const Matrix<Field>& N, \
    const Matrix<Field>& y ); \
  template Base<Field> svp::TransposedSparseToNorm \
  ( const Matrix<Base<Field>>& d, \
    const Matrix<Field>& NTrans, \
    const Matrix<Field>& y ); \
  template Base<Field> svp::PhaseEnumeration \
  ( const Matrix<Field>& B, \
    const Matrix<Base<Field>>& d, \
    const Matrix<Field>& N, \
          Base<Field> normUpperBound, \
          Int startIndex, \
          Int phaseLength, \
          double enqueueProb, \
    const vector<Int>& minInfNorms, \
    const vector<Int>& maxInfNorms, \
    const vector<Int>& minOneNorms, \
    const vector<Int>& maxOneNorms, \
          Matrix<Field>& v, \
          Int progressLevel ); \
  template pair<Base<Field>,Int> svp::PhaseEnumeration \
  ( const Matrix<Field>& B, \
    const Matrix<Base<Field>>& d, \
    const Matrix<Field>& N, \
    const Matrix<Base<Field>>& normUpperBounds, \
          Int startIndex, \
          Int phaseLength, \
          double enqueueProb, \
    const vector<Int>& minInfNorms, \
    const vector<Int>& maxInfNorms, \
    const vector<Int>& minOneNorms, \
    const vector<Int>& maxOneNorms, \
          Matrix<Field>& v, \
          Int progressLevel );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
