/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

// TODO: Move these into BLAS1?
template<typename F>
Base<F> MinAbsNonzero( const Matrix<F>& A, Base<F> upperBound )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Real minAbs = upperBound;
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            const Real alphaAbs = Abs(A(i,j));
            if( alphaAbs > Real(0) )
                minAbs = Min(minAbs,alphaAbs);
        }
    }
    return minAbs;
}

template<typename F>
Base<F> MinAbsNonzero( const SparseMatrix<F>& A, Base<F> upperBound )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int numEntries = A.NumEntries();
    Real minAbs = upperBound;
    for( Int e=0; e<numEntries; ++e )
    {
        const Real absVal = Abs(A.Value(e));
        if( absVal > Real(0) )
            minAbs = Min(minAbs,absVal);
    }
    return minAbs;
}

template<typename F>
Base<F> MinAbsNonzero
( const ElementalMatrix<F>& A, Base<F> upperBound )
{
    DEBUG_CSE
    typedef Base<F> Real;
    Real minAbs;
    if( A.Participating() )
    {
        const Real minLocAbs = MinAbsNonzero( A.LockedMatrix(), upperBound );
        minAbs = mpi::AllReduce( minLocAbs, mpi::MIN, A.DistComm() ); 
    }
    mpi::Broadcast( minAbs, A.Root(), A.CrossComm() );
    return minAbs;
}

template<typename F>
Base<F> MinAbsNonzero( const DistSparseMatrix<F>& A, Base<F> upperBound )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int numEntries = A.NumLocalEntries();
    Real minLocAbs = upperBound;
    for( Int e=0; e<numEntries; ++e )
    {
        const Real absVal = Abs(A.Value(e));
        if( absVal > Real(0) )
            minLocAbs = Min(minLocAbs,absVal);
    }
    return mpi::AllReduce( minLocAbs, mpi::MIN, A.Comm() );
}

template<typename F,Dist U,Dist V>
void GeometricColumnScaling
( const DistMatrix<F,      U,V   >& A, 
        DistMatrix<Base<F>,V,STAR>& geomScaling )
{
    DEBUG_CSE
    typedef Base<F> Real;
    DistMatrix<Real,V,STAR> maxScaling(A.Grid());
    ColumnMaxNorms( A, maxScaling );
    ColumnMinAbsNonzero( A, maxScaling, geomScaling );
    const Int nLocal = A.LocalWidth();
    auto& maxScalingLoc = maxScaling.Matrix();
    auto& geomScalingLoc = geomScaling.Matrix();
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Real maxAbs = maxScalingLoc(jLoc);
        const Real minAbs = geomScalingLoc(jLoc);
        geomScalingLoc(jLoc) = Sqrt(minAbs*maxAbs);
    }
}

template<typename F,Dist U,Dist V>
void StackedGeometricColumnScaling
( const DistMatrix<F,      U,V   >& A, 
  const DistMatrix<F,      U,V   >& B,
        DistMatrix<Base<F>,V,STAR>& geomScaling )
{
    DEBUG_CSE
    // NOTE: Assuming A.ColComm() == B.ColComm() and that the row alignments
    //       are equal
    typedef Base<F> Real;

    DistMatrix<Real,V,STAR> maxScalingA(A.Grid()),
                            maxScalingB(A.Grid());
    ColumnMaxNorms( A, maxScalingA );
    ColumnMaxNorms( B, maxScalingB );

    const Int mLocalA = A.LocalHeight();
    const Int mLocalB = B.LocalHeight();
    const Int nLocal = A.LocalWidth();
    geomScaling.AlignWith( maxScalingA );
    geomScaling.Resize( A.Width(), 1 );
    auto& ALoc = A.LockedMatrix();
    auto& BLoc = B.LockedMatrix();
    auto& geomScalingLoc = geomScaling.Matrix();
    auto& maxScalingALoc = maxScalingA.Matrix();
    auto& maxScalingBLoc = maxScalingB.Matrix();
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        Real minAbs = Max(maxScalingALoc(jLoc),maxScalingBLoc(jLoc));
        for( Int iLoc=0; iLoc<mLocalA; ++iLoc )
        {
            const Real absVal = Abs(ALoc(iLoc,jLoc));  
            if( absVal > 0 && absVal < minAbs )
                minAbs = Min(minAbs,absVal);
        }
        for( Int iLoc=0; iLoc<mLocalB; ++iLoc )
        {
            const Real absVal = Abs(BLoc(iLoc,jLoc));  
            if( absVal > 0 && absVal < minAbs )
                minAbs = Min(minAbs,absVal);
        }
        geomScalingLoc(jLoc) = minAbs;
    }
    mpi::AllReduce( geomScaling.Buffer(), nLocal, mpi::MIN, A.ColComm() );

    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Real maxAbsA = maxScalingALoc(jLoc);
        const Real maxAbsB = maxScalingBLoc(jLoc);
        const Real maxAbs = Max(maxAbsA,maxAbsB);
        const Real minAbs = geomScalingLoc(jLoc);
        geomScalingLoc(jLoc) = Sqrt(minAbs*maxAbs);
    }
}

template<typename F,Dist U,Dist V>
void GeometricRowScaling
( const DistMatrix<F,      U,V   >& A, 
        DistMatrix<Base<F>,U,STAR>& geomScaling )
{
    DEBUG_CSE
    typedef Base<F> Real;
    DistMatrix<Real,U,STAR> maxScaling(A.Grid());
    RowMaxNorms( A, maxScaling );
    RowMinAbsNonzero( A, maxScaling, geomScaling );
    const Int mLocal = A.LocalHeight();
    auto& maxScalingLoc = maxScaling.Matrix();
    auto& geomScalingLoc = geomScaling.Matrix();
    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
    {
        const Real maxAbs = maxScalingLoc(iLoc);
        const Real minAbs = geomScalingLoc(iLoc);
        geomScalingLoc(iLoc) = Sqrt(minAbs*maxAbs);
    }
}

} // namespace El
