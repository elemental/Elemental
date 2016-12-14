/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./HermitianEig/SDC.hpp"

// The targeted number of pieces to break the eigenvectors into during the
// redistribution from the [* ,VR] distribution after PMRRR to the [MC,MR]
// distribution needed for backtransformation.
#define TARGET_CHUNKS 20

namespace El {

// TODO(poulson): Decide if this should be lifted
template<typename T>
bool IsDistributed( const AbstractDistMatrix<T>& A )
{
    return A.ColStride() != 1 || A.RowStride() != 1 || A.CrossSize() != 1;
}

// Forward declarations:

namespace herm_eig {

template<typename Real>
void SortAndFilter
( Matrix<Real>& w,
  const HermitianTridiagEigCtrl<Real>& ctrl );

template<typename F>
void SortAndFilter
( Matrix<Base<F>>& w,
  Matrix<F>& Q,
  const HermitianTridiagEigCtrl<Base<F>>& ctrl );

template<typename Real>
void SortAndFilter
( AbstractDistMatrix<Real>& w,
  const HermitianTridiagEigCtrl<Real>& ctrl );

template<typename F>
void SortAndFilter
( AbstractDistMatrix<Base<F>>& w,
  AbstractDistMatrix<F>& Q,
  const HermitianTridiagEigCtrl<Base<F>>& ctrl );

} // namespace herm_eig

namespace herm_tridiag_eig {

template<typename Real>
Int MRRREstimate
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        mpi::Comm wColComm,
        Real vl,
        Real vu );

// Q is assumed to be sufficiently large and properly aligned
template<typename Real>
HermitianTridiagEigInfo
MRRRPostEstimate
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        AbstractDistMatrix<Real>& w,
        AbstractDistMatrix<Real>& Q,
        SortType sort,
        Real vl,
        Real vu );

} // namespace herm_tridiag_eig

namespace herm_eig {

// We create specialized redistribution routines for redistributing the 
// real eigenvectors of the symmetric tridiagonal matrix at the core of our 
// eigensolver in order to minimize the temporary memory usage.
template<typename F>
void InPlaceRedist( DistMatrix<F>& Q, Int rowAlign, const Base<F>* readBuffer )
{
    typedef Base<F> Real;
    const Grid& g = Q.Grid();
    const Int height = Q.Height();
    const Int width = Q.Width();

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = r * c;
    const Int row = g.Row();
    const Int col = g.Col();
    const Int rowShift = Q.RowShift();
    const Int colAlign = Q.ColAlign();
    const Int localWidth = Length(width,g.VRRank(),rowAlign,p);

    const Int maxHeight = MaxLength(height,r);
    const Int maxWidth = MaxLength(width,p);
    const Int portionSize = mpi::Pad( maxHeight*maxWidth );
    
    // Allocate our send/recv buffers
    vector<Real> buffer(2*r*portionSize);
    Real* sendBuffer = &buffer[0];
    Real* recvBuffer = &buffer[r*portionSize];

    // Pack
    EL_OUTER_PARALLEL_FOR
    for( Int k=0; k<r; ++k )
    {
        Real* data = &sendBuffer[k*portionSize];

        const Int thisColShift = Shift(k,colAlign,r);
        const Int thisLocalHeight = Length(height,thisColShift,r);

        EL_INNER_PARALLEL_FOR_COLLAPSE2
        for( Int j=0; j<localWidth; ++j )
            for( Int i=0; i<thisLocalHeight; ++i )
                data[i+j*thisLocalHeight] = 
                    readBuffer[thisColShift+i*r+j*height];
    }

    // Communicate
    mpi::AllToAll
    ( sendBuffer, portionSize,
      recvBuffer, portionSize, g.ColComm() );

    // Unpack
    const Int localHeight = Length(height,row,colAlign,r);
    EL_OUTER_PARALLEL_FOR
    for( Int k=0; k<r; ++k )
    {
        const Real* data = &recvBuffer[k*portionSize];

        const Int thisRank = col+k*c;
        const Int thisRowShift = Shift(thisRank,rowAlign,p);
        const Int thisRowOffset = (thisRowShift-rowShift) / c;
        const Int thisLocalWidth = Length(width,thisRowShift,p);

        EL_INNER_PARALLEL_FOR
        for( Int j=0; j<thisLocalWidth; ++j )
        {
            const Real* dataCol = &(data[j*localHeight]);
            Real* thisCol = (Real*)Q.Buffer(0,thisRowOffset+j*r);
            if( IsComplex<F>::value )
            {
                for( Int i=0; i<localHeight; ++i )
                {
                    thisCol[2*i] = dataCol[i];
                    thisCol[2*i+1] = 0;
                }
            }
            else
            {
                MemCopy( thisCol, dataCol, localHeight );
            }
        }
    }
}

} // namespace herm_eig

// Compute eigenvalues
// ===================

namespace herm_eig {

template<typename F>
HermitianEigInfo
BlackBox
( UpperOrLower uplo,
  Matrix<F>& A, 
  Matrix<Base<F>>& w,
  const HermitianEigCtrl<F>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    HermitianEigInfo info;
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");
    auto subset = ctrl.tridiagEigCtrl.subset;

    if( subset.indexSubset && subset.rangeSubset )
        LogicError("Cannot mix index and range subsets");
    if( (subset.rangeSubset && (subset.lowerBound >= subset.upperBound)) ||
        (subset.indexSubset && (subset.lowerIndex > subset.upperIndex)) )
    {
        w.Resize(0,1);
        return info;
    }

    // Check if we need to rescale the matrix, and do so if necessary
    const Real maxNormA = HermitianMaxNorm( uplo, A );
    const Real underflowThreshold = limits::Min<Real>();
    const Real overflowThreshold = limits::Max<Real>();
    // TODO(poulson): Better norm bounds?
    const Real normMax = overflowThreshold;
    const Real normMin = underflowThreshold;
    bool scaledDown=false, scaledUp=false;
    if( maxNormA == Real(0) )
    {
        const Int n = A.Height();
        Int numValid = n;
        if( subset.indexSubset )
        {   
            numValid = subset.upperIndex-subset.lowerIndex+1;
        }
        else if( subset.rangeSubset )
        {
            if( subset.lowerBound >= Real(0) || subset.upperBound < Real(0) )
            {   
                numValid = 0;
            }
        }
        Zeros( w, numValid, 1 );
        return info;
    }
    if( maxNormA > normMax )
    {
        scaledDown = true;
        SafeScaleTrapezoid( maxNormA, normMax, uplo, A );
    }
    else if( maxNormA < normMin )
    {
        scaledUp = true;
        SafeScaleTrapezoid( maxNormA, normMin, uplo, A );
    }

    // TODO(poulson): Extend interface to support accepting ctrl.tridiagCtrl
    herm_tridiag::ExplicitCondensed( uplo, A );

    auto d = GetRealPartOfDiagonal(A);
    auto dSub = GetDiagonal( A, (uplo==LOWER?-1:1) );
    info.tridiagEigInfo =
      HermitianTridiagEig( d, dSub, w, ctrl.tridiagEigCtrl );

    // Rescale the eigenvalues if necessary
    if( scaledDown )
    {
        SafeScale( normMax, maxNormA, w );
    }
    else if( scaledUp )
    {
        SafeScale( normMin, maxNormA, w );
    }

    return info;
}

} // namespace herm_eig

template<typename F>
HermitianEigInfo
HermitianEig
( UpperOrLower uplo,
  Matrix<F>& A, 
  Matrix<Base<F>>& w,
  const HermitianEigCtrl<F>& ctrl )
{
    EL_DEBUG_CSE
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");
    if( ctrl.useSDC )
    {
        HermitianEigInfo info;
        herm_eig::SDC( uplo, A, w, ctrl.sdcCtrl );
        herm_eig::SortAndFilter( w, ctrl.tridiagEigCtrl );
        return info;
    }
    return herm_eig::BlackBox( uplo, A, w, ctrl );
}

namespace herm_eig {

template<typename F>
HermitianEigInfo
SequentialHelper
( UpperOrLower uplo,
  AbstractDistMatrix<F>& A, 
  AbstractDistMatrix<Base<F>>& w,
  const HermitianEigCtrl<F>& ctrl )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( IsDistributed(A) )
          LogicError("A should not`have been distributed");
      if( IsDistributed(w) )
          LogicError("w should not`have been distributed");
    )
    const Int n = A.Height();
    w.SetGrid( A.Grid() );
    HermitianEigInfo info;

    auto subset = ctrl.tridiagEigCtrl.subset;
    if( subset.indexSubset || subset.rangeSubset )
    {
        Matrix<Base<F>> wProx;
        wProx.Resize( n, 1 );
        info = HermitianEig( uplo, A.Matrix(), wProx, ctrl );
        w.Resize( wProx.Height(), 1 );
        w.Matrix() = wProx;
    }
    else
    {
        w.Resize( n, 1 );
        info = HermitianEig( uplo, A.Matrix(), w.Matrix(), ctrl );
    }

    return info;
}


#ifdef EL_HAVE_SCALAPACK
template<typename F,typename=EnableIf<IsBlasScalar<Base<F>>>>
HermitianEigInfo
ScaLAPACKHelper
( UpperOrLower uplo,
  AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<Base<F>>& wPre,
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() )
{
    EL_DEBUG_CSE

    DistMatrixReadProxy<F,F,MC,MR,BLOCK> AProx( APre );
    auto& A = AProx.Get();

    DistMatrixWriteProxy<Base<F>,Base<F>,STAR,STAR> wProx( wPre );
    auto& w = wProx.Get();

    HermitianEigInfo info;
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");
    w.SetGrid( A.Grid() );

    auto subset = ctrl.tridiagEigCtrl.subset;
    if( subset.indexSubset && subset.rangeSubset )
        LogicError("Cannot mix index and range subsets");
    if( (subset.rangeSubset && (subset.lowerBound >= subset.upperBound)) ||
        (subset.indexSubset && (subset.lowerIndex > subset.upperIndex)) )
    {
        w.Resize(0,1);
        return info;
    }

    const Int n = A.Height();
    w.Resize( n, 1 );

    if( subset.rangeSubset )
    {
        LogicError("This option is not yet supported");
    }
    else if( subset.indexSubset )
    {
        LogicError("This option is not yet supported");
    }
    else
    {
        const char uploChar = UpperOrLowerToChar( uplo );
        auto descA = FillDesc( A );
        scalapack::HermitianEig
        ( uploChar, n, A.Buffer(), descA.data(), w.Buffer() );
    }

    return info;
}

template<typename F,typename=DisableIf<IsBlasScalar<Base<F>>>,typename=void>
HermitianEigInfo
ScaLAPACKHelper
( UpperOrLower uplo,
  AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<Base<F>>& wPre,
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() )
{
    EL_DEBUG_CSE
    LogicError("This routine should not be called");
    HermitianEigInfo info;
    return info;
}
#endif // ifdef EL_HAVE_SCALAPACK

template<typename F>
HermitianEigInfo
BlackBox
( UpperOrLower uplo,
  AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<Base<F>>& w,
  const HermitianEigCtrl<F>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    HermitianEigInfo info;

    auto subset = ctrl.tridiagEigCtrl.subset;
    if( subset.indexSubset && subset.rangeSubset )
        LogicError("Cannot mix index and range subsets");
    if( (subset.rangeSubset && (subset.lowerBound >= subset.upperBound)) ||
        (subset.indexSubset && (subset.lowerIndex > subset.upperIndex)) )
    {
        w.Resize(0,1);
        return info;
    }

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    // Check if we need to rescale the matrix, and do so if necessary
    const Real maxNormA = HermitianMaxNorm( uplo, A );
    const Real underflowThreshold = limits::Min<Real>();
    const Real overflowThreshold = limits::Max<Real>();
    // TODO(poulson): Better norm bounds?
    const Real normMax = overflowThreshold;
    const Real normMin = underflowThreshold;
    bool scaledDown=false, scaledUp=false;
    if( maxNormA == Real(0) )
    {
        const Int n = A.Height();
        Int numValid = n;
        if( subset.indexSubset )
        {
            numValid = subset.upperIndex-subset.lowerIndex+1;
        }
        else if( subset.rangeSubset )
        {
            if( subset.lowerBound >= Real(0) || subset.upperBound < Real(0) )
            {   
                numValid = 0;
            }
        }
        Zeros( w, numValid, 1 );
        return info;
    }
    if( maxNormA > normMax )
    {
        scaledDown = true;
        SafeScaleTrapezoid( maxNormA, normMax, uplo, A );
    }
    else if( maxNormA < normMin )
    {
        scaledUp = true;
        SafeScaleTrapezoid( maxNormA, normMin, uplo, A );
    }

    Timer timer;
    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
            timer.Start();
    }
   
    // Tridiagonalize A
    herm_tridiag::ExplicitCondensed( uplo, A, ctrl.tridiagCtrl );

    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
        {
            Output("  Condense:      ",timer.Stop()," secs");
            timer.Start();
        }
    }

    // Solve the symmetric tridiagonal EVP
    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
    auto d = GetRealPartOfDiagonal(A);
    auto e = GetRealPartOfDiagonal(A,subdiagonal);
    info.tridiagEigInfo = HermitianTridiagEig( d, e, w, ctrl.tridiagEigCtrl );

    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
            Output("  TridiagEig:    ",timer.Stop()," secs");
    }

    // Rescale the eigenvalues if necessary
    if( scaledDown )
    {
        SafeScale( normMax, maxNormA, w );
    }
    else if( scaledUp )
    {
        SafeScale( normMin, maxNormA, w );
    }

    return info;
}

} // namespace herm_eig

template<typename F>
HermitianEigInfo
HermitianEig
( UpperOrLower uplo,
  AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<Base<F>>& w,
  const HermitianEigCtrl<F>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    if( APre.Height() != APre.Width() )
        LogicError("Hermitian matrices must be square");

    if( ctrl.useScaLAPACK && IsBlasScalar<Real>::value )
    {
#ifdef EL_HAVE_SCALAPACK
        return herm_eig::ScaLAPACKHelper( uplo, APre, w, ctrl );
#endif
    }

    if( !IsDistributed(APre) && !IsDistributed(w) )
    {
        return herm_eig::SequentialHelper( uplo, APre, w, ctrl );
    }

    if( ctrl.useSDC )
    {
        HermitianEigInfo info;
        herm_eig::SDC( uplo, APre, w, ctrl.sdcCtrl );
        herm_eig::SortAndFilter( w, ctrl.tridiagEigCtrl );
        return info;
    }

    return herm_eig::BlackBox( uplo, APre, w, ctrl );
}

// Compute eigenpairs
// ==================

namespace herm_eig {

template<typename F>
HermitianEigInfo
BlackBox
( UpperOrLower uplo,
  Matrix<F>& A,
  Matrix<Base<F>>& w,
  Matrix<F>& Q, 
  const HermitianEigCtrl<F>& ctrl )
{
    EL_DEBUG_CSE
    HermitianEigInfo info;

    // TODO(poulson): Extend interface to support ctrl.tridiagCtrl
    Matrix<F> householderScalars;
    HermitianTridiag( uplo, A, householderScalars );

    auto d = GetRealPartOfDiagonal(A);
    auto dSub = GetDiagonal( A, (uplo==LOWER?-1:1) );
    info.tridiagEigInfo =
      HermitianTridiagEig( d, dSub, w, Q, ctrl.tridiagEigCtrl );

    herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, householderScalars, Q );

    return info;
}

template<typename F>
HermitianEigInfo
BlackBox
( UpperOrLower uplo,
  AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<Base<F>>& w,
  AbstractDistMatrix<F>& QPre, 
  const HermitianEigCtrl<F>& ctrl )
{
    EL_DEBUG_CSE
    const Grid& g = APre.Grid();
    HermitianEigInfo info;

    DistMatrixReadProxy<F,F,MC,MR> AProx( APre ); 
    auto& A = AProx.Get();

    // TODO(poulson): Extend interface to support ctrl.tridiagCtrl
    DistMatrix<F,VC,STAR> householderScalars(g);
    HermitianTridiag( uplo, A, householderScalars );

    auto d = GetRealPartOfDiagonal(A);
    auto dSub = GetDiagonal( A, (uplo==LOWER?-1:1) );

    if( ctrl.tridiagEigCtrl.accumulateEigVecs )
    {
        DistMatrixReadWriteProxy<F,F,MC,MR> QProx( QPre );
        auto& Q = QProx.Get();

        info.tridiagEigInfo =
          HermitianTridiagEig( d, dSub, w, Q, ctrl.tridiagEigCtrl );
        herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, householderScalars, Q );
    }
    else
    {
        DistMatrixWriteProxy<F,F,MC,MR> QProx( QPre );
        auto& Q = QProx.Get();

        info.tridiagEigInfo =
          HermitianTridiagEig( d, dSub, w, Q, ctrl.tridiagEigCtrl );
        herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, householderScalars, Q );
    }

    return info;
}

} // namespace herm_eig

template<typename F>
HermitianEigInfo
HermitianEig
( UpperOrLower uplo,
  Matrix<F>& A,
  Matrix<Base<F>>& w,
  Matrix<F>& Q, 
  const HermitianEigCtrl<F>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int n = A.Height();
    auto subset = ctrl.tridiagEigCtrl.subset;
    HermitianEigInfo info;

    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    if( subset.indexSubset && subset.rangeSubset )
        LogicError("Cannot mix index and range subsets");
    if( (subset.rangeSubset && (subset.lowerBound >= subset.upperBound)) ||
        (subset.indexSubset && (subset.lowerIndex > subset.upperIndex)) )
    {
        w.Resize(0,1);
        Q.Resize(n,0);
        return info;
    }

    // Check if we need to rescale the matrix, and do so if necessary
    const Real maxNormA = HermitianMaxNorm( uplo, A );
    const Real underflowThreshold = limits::Min<Real>();
    const Real overflowThreshold = limits::Max<Real>();
    // TODO(poulson): Better norm bounds?
    const Real normMax = overflowThreshold;
    const Real normMin = underflowThreshold;
    bool scaledDown=false, scaledUp=false;
    if( maxNormA == Real(0) )
    {
        const Int n = A.Height();
        Int numValid = n;
        if( subset.indexSubset )
        {   
            numValid = subset.upperIndex-subset.lowerIndex+1;
        }
        else if( subset.rangeSubset )
        {
            if( subset.lowerBound >= Real(0) || subset.upperBound < Real(0) )
            {   
                numValid = 0;
            }
        }
        Zeros( w, numValid, 1 );
        return info;
    }
    if( maxNormA > normMax )
    {
        scaledDown = true;
        SafeScaleTrapezoid( maxNormA, normMax, uplo, A );
    }
    else if( maxNormA < normMin )
    {
        scaledUp = true;
        SafeScaleTrapezoid( maxNormA, normMin, uplo, A );
    }

    if( ctrl.useSDC )
    {
        herm_eig::SDC( uplo, A, w, Q, ctrl.sdcCtrl );
        herm_eig::SortAndFilter( w, Q, ctrl.tridiagEigCtrl );
    }
    else
    {
        info = herm_eig::BlackBox( uplo, A, w, Q, ctrl );
    }

    // Rescale the eigenvalues if necessary
    if( scaledDown )
    {
        SafeScale( normMax, maxNormA, w );
    }
    else if( scaledUp )
    {
        SafeScale( normMin, maxNormA, w );
    }

    return info;
}

namespace herm_eig {

template<typename F>
HermitianEigInfo
SequentialHelper
( UpperOrLower uplo,
  AbstractDistMatrix<F>& A, 
  AbstractDistMatrix<Base<F>>& w,
  AbstractDistMatrix<F>& Q, 
  const HermitianEigCtrl<F>& ctrl )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( IsDistributed(A) )
          LogicError("A should not`have been distributed");
      if( IsDistributed(w) )
          LogicError("w should not`have been distributed");
      if( IsDistributed(Q) )
          LogicError("Q should not`have been distributed");
    )

    const Int n = A.Height();
    w.SetGrid( A.Grid() );
    Q.SetGrid( A.Grid() );
    HermitianEigInfo info;

    auto subset = ctrl.tridiagEigCtrl.subset;
    if( subset.indexSubset || subset.rangeSubset )
    {
        // TODO(poulson): Avoid the need for an extra matrix here
        Matrix<Base<F>> wProx;
        Matrix<F> QProx;
        wProx.Resize( n, 1 );
        QProx.Resize( n, n ); 

        info = HermitianEig( uplo, A.Matrix(), wProx, QProx, ctrl );

        w.Resize( wProx.Height(), 1 );
        w.Matrix() = wProx;

        Q.Resize( QProx.Height(), QProx.Width() );
        Q.Matrix() = QProx;
    }
    else
    {
        w.Resize( n, 1 );
        Q.Resize( n, n );
        info = HermitianEig( uplo, A.Matrix(), w.Matrix(), Q.Matrix(), ctrl );
    }

    return info;
}

#ifdef EL_HAVE_SCALAPACK
template<typename F,typename=EnableIf<IsBlasScalar<Base<F>>>>
HermitianEigInfo
ScaLAPACKHelper
( UpperOrLower uplo,
  AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<Base<F>>& wPre,
  AbstractDistMatrix<F>& QPre,
  const HermitianEigCtrl<F>& ctrl )
{
    EL_DEBUG_CSE

    DistMatrixReadProxy<F,F,MC,MR,BLOCK> AProx( APre );
    auto& A = AProx.Get();

    DistMatrixWriteProxy<Base<F>,Base<F>,STAR,STAR> wProx( wPre );
    auto& w = wProx.Get();

    DistMatrixWriteProxy<F,F,MC,MR,BLOCK> QProx( QPre );
    auto& Q = QProx.Get();

    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");
    HermitianEigInfo info;

    const Int n = A.Height();
    auto subset = ctrl.tridiagEigCtrl.subset;
    if( subset.indexSubset && subset.rangeSubset )
        LogicError("Cannot mix index and range subsets");
    if( (subset.rangeSubset && (subset.lowerBound >= subset.upperBound)) ||
        (subset.indexSubset && (subset.lowerIndex > subset.upperIndex)) )
    {
        w.SetGrid( A.Grid() );
        Q.SetGrid( A.Grid() );
        w.Resize(0,1);
        Q.Resize(n,0);
        return info;
    }

    w.Resize( n, 1 );
    Q.SetGrid( A.Grid() );
    Q.AlignWith( A );
    Q.Resize( n, n );


    if( subset.rangeSubset )
    {
        LogicError("This option is not yet supported");
    }
    else if( subset.indexSubset )
    {
        LogicError("This option is not yet supported");
    }
    else
    {
        const char uploChar = UpperOrLowerToChar( uplo );
        auto descA = FillDesc( A );
        auto descQ = FillDesc( Q );
        scalapack::HermitianEig
        ( uploChar, n,
          A.Buffer(), descA.data(),
          w.Buffer(),
          Q.Buffer(), descQ.data() );
    }

    return info;
}

template<typename F,typename=DisableIf<IsBlasScalar<Base<F>>>,typename=void>
HermitianEigInfo
ScaLAPACKHelper
( UpperOrLower uplo,
  AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<Base<F>>& w,
  AbstractDistMatrix<F>& Q,
  const HermitianEigCtrl<F>& ctrl )
{
    EL_DEBUG_CSE
    LogicError("This routine should not be called");
    HermitianEigInfo info;
    return info;
}
#endif // ifdef EL_HAVE_SCALAPACK

template<typename F>
HermitianEigInfo
MRRR
( UpperOrLower uplo,
  AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<Base<F>>& w,
  AbstractDistMatrix<F>& QPre,
  const HermitianEigCtrl<F>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int n = APre.Height();
    const Grid& g = APre.Grid();
    auto subset = ctrl.tridiagEigCtrl.subset;
    HermitianEigInfo info;
    Timer timer;

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    // Tridiagonalize A
    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
            timer.Start();
    }
    DistMatrix<F,STAR,STAR> householderScalars(g);
    HermitianTridiag( uplo, A, householderScalars, ctrl.tridiagCtrl );
    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
        {
            Output("  Condense time:      ",timer.Stop()," secs");
            timer.Start();
        }
    }

    Int kEst;
    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
    auto d = GetRealPartOfDiagonal(A);
    auto e = GetRealPartOfDiagonal(A,subdiagonal);
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<Real,STAR,STAR> e_STAR_STAR( g );
    e_STAR_STAR.Resize( n-1, 1, n );
    e_STAR_STAR = e;
    if( subset.rangeSubset )
    {
        // Get an upper-bound on the number of local eigenvalues in the range
        kEst = herm_tridiag_eig::MRRREstimate
          ( d, e, g.VRComm(), subset.lowerBound, subset.upperBound );
    }
    else if( subset.indexSubset )
        kEst = subset.upperIndex-subset.lowerIndex+1;
    else
        kEst = n;

    // We will use the same buffer for Q in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Q's dimensions slightly.
    const Int N = MaxLength(n,g.Height())*g.Height();
    const Int K = MaxLength(kEst,g.Size())*g.Size(); 

    ElementalProxyCtrl proxCtrl;
    proxCtrl.colConstrain = true;
    proxCtrl.rowConstrain = true;
    proxCtrl.colAlign = 0;
    proxCtrl.rowAlign = 0;
    
    DistMatrixWriteProxy<F,F,MC,MR> QProx( QPre, proxCtrl );
    auto& Q = QProx.Get();

    Q.Resize( N, K );
    DistMatrix<Real,STAR,VR> Q_STAR_VR(g);
    {
        // Grab a slice of size Q_STAR_VR_BufferSize from the very end
        // of QBuf so that we can later redistribute in place
        Real* QBuf = (Real*)Q.Buffer();
        const Int QBufSize =
            ( IsComplex<F>::value ? 2*Q.LDim()*Q.LocalWidth()
                                  :   Q.LDim()*Q.LocalWidth() );
        const Int Q_STAR_VR_LocalWidth = Length(kEst,g.VRRank(),g.Size());
        const Int Q_STAR_VR_BufSize = n*Q_STAR_VR_LocalWidth;
        Real* Q_STAR_VR_Buf = &QBuf[QBufSize-Q_STAR_VR_BufSize];
        Q_STAR_VR.Attach( n, kEst, g, 0, 0, Q_STAR_VR_Buf, n );
    }
    // NOTE: We should be guaranteeing that Q_STAR_VR does not need to
    //       reallocate a buffer
    if( subset.rangeSubset )
        info.tridiagEigInfo = herm_tridiag_eig::MRRRPostEstimate
        ( d_STAR_STAR, e_STAR_STAR, w, Q_STAR_VR, UNSORTED,
          subset.lowerBound, subset.upperBound );
    else
        info.tridiagEigInfo = HermitianTridiagEig
        ( d_STAR_STAR, e_STAR_STAR, w, Q_STAR_VR, ctrl.tridiagEigCtrl );

    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
        {
            Output("  TridiagEig:    ",timer.Stop()," secs");
            timer.Start();
        }
    }

    const Int k = w.Height();
    {
        // Redistribute Q piece-by-piece in place. This is to keep the 
        // send/recv buffer memory usage low.
        const Int p = g.Size();
        const Int numEqualPanels = K/p;
        const Int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const Int nbProp = numPanelsPerComm*p;

        // Manually maintain information about the implicit Q[* ,VR] stored 
        // at the end of the Q[MC,MR] buffers.
        Int alignment = 0;
        const Real* readBuffer = Q_STAR_VR.LockedBuffer();
        for( Int j=0; j<k; j+=nbProp )
        {
            const Int nb = Min(nbProp,k-j);
            auto Q1 = Q( IR(0,n), IR(j,j+nb) );

            // Redistribute Q1[MC,MR] <- Q1[* ,VR] in place.
            // NOTE: This assumes that Q_STAR_VR did not reallocate within
            //       herm_tridiag_eig::MRRR[PostEstimate] above
            herm_eig::InPlaceRedist( Q1, alignment, readBuffer );

            // Update the Q1[* ,VR] information
            const Int localWidth = nb/p;
            readBuffer = &readBuffer[localWidth*n];
            alignment = (alignment+nb) % p;
        }
    }
    Q.Resize( n, k ); // We can simply shrink matrices

    // Backtransform the tridiagonal eigenvectors, Q
    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
        {
            Output("  Redist:        ",timer.Stop()," secs");
            timer.Start();
        }
    }
    herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, householderScalars, Q );
    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
            Output("  Backtransform: ",timer.Stop()," secs");
    }

    return info;
}

} // namespace herm_eig

template<typename F>
HermitianEigInfo
HermitianEig
( UpperOrLower uplo,
  AbstractDistMatrix<F>& A,
  AbstractDistMatrix<Base<F>>& w,
  AbstractDistMatrix<F>& Q,
  const HermitianEigCtrl<F>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int n = A.Height();
    auto subset = ctrl.tridiagEigCtrl.subset;
    HermitianEigInfo info;
    Timer timer;

    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    if( subset.indexSubset && subset.rangeSubset )
        LogicError("Cannot mix index and range subsets");
    if( (subset.rangeSubset && (subset.lowerBound >= subset.upperBound)) ||
        (subset.indexSubset && (subset.lowerIndex > subset.upperIndex)) )
    {
        w.SetGrid( A.Grid() );
        w.Resize(0,1);
        Q.SetGrid( A.Grid() );
        Q.Resize(n,0);
        return info;
    }

    if( ctrl.useScaLAPACK && IsBlasScalar<Real>::value )
    {
#ifdef EL_HAVE_SCALAPACK
        return herm_eig::ScaLAPACKHelper( uplo, A, w, Q, ctrl );
#endif
    }
    if( !IsDistributed(A) && !IsDistributed(w) && !IsDistributed(Q) )
    {
        return herm_eig::SequentialHelper( uplo, A, w, Q, ctrl );
    }

    // Check if we need to rescale the matrix, and do so if necessary
    const Real maxNormA = HermitianMaxNorm( uplo, A );
    const Real underflowThreshold = limits::Min<Real>();
    const Real overflowThreshold = limits::Max<Real>();
    // TODO(poulson): Better norm bounds?
    const Real normMax = overflowThreshold;
    const Real normMin = underflowThreshold;
    bool scaledDown=false, scaledUp=false;
    if( maxNormA == Real(0) )
    {
        const Int n = A.Height();
        Int numValid = n;
        if( subset.indexSubset )
        {
            numValid = subset.upperIndex-subset.lowerIndex+1;
        }
        else if( subset.rangeSubset )
        {
            if( subset.lowerBound >= Real(0) || subset.upperBound < Real(0) )
            {
                numValid = 0;
            }
        }
        Zeros( w, numValid, 1 );
        Zeros( Q, n, numValid );
        return info;
    }
    if( maxNormA > normMax )
    {
        scaledDown = true;
        SafeScaleTrapezoid( maxNormA, normMax, uplo, A );
    }
    else if( maxNormA < normMin )
    {
        scaledUp = true;
        SafeScaleTrapezoid( maxNormA, normMin, uplo, A );
    }

    if( ctrl.useSDC )
    {
        herm_eig::SDC( uplo, A, w, Q, ctrl.sdcCtrl );
        herm_eig::SortAndFilter( w, Q, ctrl.tridiagEigCtrl );
    }
    else if( ctrl.tridiagEigCtrl.alg == HERM_TRIDIAG_EIG_MRRR )
    {
        info = herm_eig::MRRR( uplo, A, w, Q, ctrl );
    }
    else
    {
        info = herm_eig::BlackBox( uplo, A, w, Q, ctrl );
    }

    // Rescale the eigenvalues if necessary
    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
            timer.Start();
    }
    if( scaledDown )
    {
        SafeScale( normMax, maxNormA, w );
    }
    else if( scaledUp )
    {
        SafeScale( normMin, maxNormA, w );
    }

    auto sortPairs = TaggedSort( w, ctrl.tridiagEigCtrl.sort );
    for( Int j=0; j<n; ++j )
        w.Set( j, 0, sortPairs[j].value );
    ApplyTaggedSortToEachRow( sortPairs, Q );

    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
            Output("  Scale+sort:    ",timer.Stop()," secs");
    }

    return info;
}

#define EIGVAL_PROTO(F) \
  template HermitianEigInfo HermitianEig\
  ( UpperOrLower uplo, \
    Matrix<F>& A, \
    Matrix<Base<F>>& w, \
    const HermitianEigCtrl<F>& ctrl ); \
  template HermitianEigInfo HermitianEig\
  ( UpperOrLower uplo, \
    AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<Base<F>>& w, \
    const HermitianEigCtrl<F>& ctrl );

#define EIGPAIR_PROTO(F) \
  template HermitianEigInfo HermitianEig\
  ( UpperOrLower uplo, \
    Matrix<F>& A, \
    Matrix<Base<F>>& w, \
    Matrix<F>& Q,\
    const HermitianEigCtrl<F>& ctrl ); \
  template HermitianEigInfo HermitianEig\
  ( UpperOrLower uplo, \
    AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<Base<F>>& w, \
    AbstractDistMatrix<F>& Q, \
    const HermitianEigCtrl<F>& ctrl );

#define PROTO(F) \
  EIGVAL_PROTO(F) \
  EIGPAIR_PROTO(F)

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
