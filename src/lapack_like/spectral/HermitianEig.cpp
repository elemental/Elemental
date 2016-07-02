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

template<typename F>
bool CheckScale( UpperOrLower uplo, DistMatrix<F>& A, Base<F>& scale )
{
    typedef Base<F> Real;

    scale = 1;
    const Real maxNormOfA = HermitianMaxNorm( uplo, A );
    const Real underflowThreshold = limits::Min<Real>();
    const Real overflowThreshold = limits::Max<Real>();
    if( maxNormOfA > 0 && maxNormOfA < underflowThreshold )
    {
        scale = underflowThreshold / maxNormOfA;
        return true;
    }
    else if( maxNormOfA > overflowThreshold )
    {
        // This branch should never be activated with most datatypes since
        // overflowing implies an infinite result
        // (and thus rescaling would be meaningless)
        scale = overflowThreshold / maxNormOfA;
        return true;
    }
    else
        return false;
}

} // namespace herm_eig

// Compute eigenvalues
// ===================

template<typename F>
HermitianEigInfo
HermitianEig
( UpperOrLower uplo,
  Matrix<F>& A, 
  Matrix<Base<F>>& w,
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_CSE
    HermitianEigInfo info;
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");
    auto subset = ctrl.tridiagEigCtrl.subset;
    if( ctrl.useSDC )
    {
        herm_eig::SDC( uplo, A, w, ctrl.sdcCtrl );
        herm_eig::SortAndFilter( w, ctrl.tridiagEigCtrl );
        return info;
    }

    if( subset.indexSubset && subset.rangeSubset )
        LogicError("Cannot mix index and range subsets");
    if( (subset.rangeSubset && (subset.lowerBound >= subset.upperBound)) ||
        (subset.indexSubset && (subset.lowerIndex > subset.upperIndex)) )
    {
        w.Resize(0,1);
        return info;
    }

    // TODO(poulson): Extend interface to support accepting ctrl.tridiagCtrl
    herm_tridiag::ExplicitCondensed( uplo, A );

    auto d = GetRealPartOfDiagonal(A);
    auto dSub = GetDiagonal( A, (uplo==LOWER?-1:1) );
    info.tridiagEigInfo =
      HermitianTridiagEig( d, dSub, w, ctrl.tridiagEigCtrl );

    return info;
}

template<typename F>
HermitianEigInfo
HermitianEig
( UpperOrLower uplo,
  DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<Base<F>,STAR,STAR>& w,
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_CSE
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
        info = HermitianEig( uplo, A.Matrix(), w.Matrix(), ctrl );
    }

    return info;
}

template<typename F>
HermitianEigInfo
HermitianEig
( UpperOrLower uplo,
  ElementalMatrix<F>& APre,
  ElementalMatrix<Base<F>>& w,
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_CSE
    HermitianEigInfo info;
    if( APre.Height() != APre.Width() )
        LogicError("Hermitian matrices must be square");

    if( ctrl.useSDC )
    {
        herm_eig::SDC( uplo, APre, w, ctrl.sdcCtrl );
        // TODO(poulson): Add support for this function
        //herm_eig::SortAndFilter( w, ctrl.tridiagEigCtrl );
        return info;
    }

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
    Base<F> scale;
    const bool needRescaling = herm_eig::CheckScale( uplo, A, scale );
    if( needRescaling )
        ScaleTrapezoid( F(scale), uplo, A );

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
    if( needRescaling ) 
        w *= 1/scale;

    return info;
}

namespace herm_eig {

template<typename F,typename=EnableIf<IsBlasScalar<Base<F>>>>
HermitianEigInfo
Helper
( UpperOrLower uplo,
  DistMatrix<F,MC,MR,BLOCK>& A,
  DistMatrix<Base<F>,STAR,STAR>& w,
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() )
{
    DEBUG_CSE
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");
    HermitianEigInfo info;
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

#ifdef EL_HAVE_SCALAPACK
    const Int n = A.Height();
    w.Resize( n, 1 );

    const int bHandle = blacs::Handle( A );
    const int context = blacs::GridInit( bHandle, A );
    auto descA = FillDesc( A, context );
    
    const char uploChar = UpperOrLowerToChar( uplo );
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
        scalapack::HermitianEig
        ( uploChar, n, A.Buffer(), descA.data(), w.Buffer() );
    }
    return info;
#else
    DistMatrix<F,MC,MR> AElem( A );
    return HermitianEig( uplo, AElem, w, ctrl );
#endif
}

template<typename F,typename=DisableIf<IsBlasScalar<Base<F>>>,typename=void>
HermitianEigInfo
Helper
( UpperOrLower uplo,
  DistMatrix<F,MC,MR,BLOCK>& A,
  DistMatrix<Base<F>,STAR,STAR>& w,
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() )
{
    DEBUG_CSE
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");
    HermitianEigInfo info;
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

    DistMatrix<F,MC,MR> AElem( A );
    return HermitianEig( uplo, AElem, w, ctrl );
}

} // namespace herm_eig

template<typename F>
HermitianEigInfo
HermitianEig
( UpperOrLower uplo,
  DistMatrix<F,MC,MR,BLOCK>& A,
  DistMatrix<Base<F>,STAR,STAR>& w,
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_CSE
    return herm_eig::Helper( uplo, A, w, ctrl );
}

// Compute eigenpairs
// ==================

template<typename F>
HermitianEigInfo
HermitianEig
( UpperOrLower uplo,
  Matrix<F>& A,
  Matrix<Base<F>>& w,
  Matrix<F>& Q, 
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_CSE
    const Int n = A.Height();
    HermitianEigInfo info;
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");
    auto subset = ctrl.tridiagEigCtrl.subset;
    if( ctrl.useSDC )
    {
        herm_eig::SDC( uplo, A, w, Q, ctrl.sdcCtrl );
        herm_eig::SortAndFilter( w, Q, ctrl.tridiagEigCtrl );
        return info;
    }

    if( subset.indexSubset && subset.rangeSubset )
        LogicError("Cannot mix index and range subsets");
    if( (subset.rangeSubset && (subset.lowerBound >= subset.upperBound)) ||
        (subset.indexSubset && (subset.lowerIndex > subset.upperIndex)) )
    {
        w.Resize(0,1);
        Q.Resize(n,0);
        return info;
    }

    Matrix<F> phase;
    // TODO(poulson): Extend interface to support ctrl.tridiagCtrl
    HermitianTridiag( uplo, A, phase );

    auto d = GetRealPartOfDiagonal(A);
    auto dSub = GetDiagonal( A, (uplo==LOWER?-1:1) );
    info.tridiagEigInfo =
      HermitianTridiagEig( d, dSub, w, Q, ctrl.tridiagEigCtrl );

    herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, phase, Q );

    return info;
}

template<typename F>
HermitianEigInfo
HermitianEig
( UpperOrLower uplo,
  DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<Base<F>,STAR,STAR>& w,
  DistMatrix<F,STAR,STAR>& Q, 
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_CSE
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
        info = HermitianEig( uplo, A.Matrix(), w.Matrix(), Q.Matrix(), ctrl );
    }

    return info;
}

template<typename F>
HermitianEigInfo
HermitianEig
( UpperOrLower uplo,
  ElementalMatrix<F>& APre,
  ElementalMatrix<Base<F>>& w,
  ElementalMatrix<F>& QPre,
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int n = APre.Height();
    auto subset = ctrl.tridiagEigCtrl.subset;
    HermitianEigInfo info;
    if( APre.Height() != APre.Width() )
        LogicError("Hermitian matrices must be square");

    if( ctrl.useSDC )
    {
        herm_eig::SDC( uplo, APre, w, QPre, ctrl.sdcCtrl );
        // TODO(poulson): Add support for this function
        //herm_eig::SortAndFilter( w, QPre, ctrl.tridiagEigCtrl );
        return info;
    }

    if( subset.indexSubset && subset.rangeSubset )
        LogicError("Cannot mix index and range subsets");
    if( (subset.rangeSubset && (subset.lowerBound >= subset.upperBound)) ||
        (subset.indexSubset && (subset.lowerIndex > subset.upperIndex)) )
    {
        w.SetGrid( APre.Grid() );
        w.Resize(0,1);
        QPre.SetGrid( APre.Grid() );
        QPre.Resize(n,0);
        return info;
    }

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    // Check if we need to rescale the matrix, and do so if necessary
    Real scale;
    const bool needRescaling = herm_eig::CheckScale( uplo, A, scale );
    if( needRescaling )
        ScaleTrapezoid( F(scale), uplo, A );

    Timer timer;
    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
            timer.Start();
    }

    // Tridiagonalize A
    const Grid& g = A.Grid();
    DistMatrix<F,STAR,STAR> phase(g);
    HermitianTridiag( uplo, A, phase, ctrl.tridiagCtrl );

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
        kEst = HermitianTridiagEigEstimate
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
        info.tridiagEigInfo = HermitianTridiagEigPostEstimate
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
            //       HermitianTridiagEig[PostEstimate] above
            herm_eig::InPlaceRedist( Q1, alignment, readBuffer );

            // Update the Q1[* ,VR] information
            const Int localWidth = nb/p;
            readBuffer = &readBuffer[localWidth*n];
            alignment = (alignment+nb) % p;
        }
    }
    Q.Resize( n, k ); // We can simply shrink matrices

    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
        {
            Output("  Redist:        ",timer.Stop()," secs");
            timer.Start();
        }
    }

    // Backtransform the tridiagonal eigenvectors, Q
    herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, phase, Q );

    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
        {
            Output("  Backtransform: ",timer.Stop()," secs");
            timer.Start();
        }
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        w *= 1/scale;

    herm_eig::Sort( w, Q, ctrl.tridiagEigCtrl.sort );

    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
            Output("  Scale+sort:    ",timer.Stop()," secs");
    }

    return info;
}

namespace herm_eig {

template<typename F,typename=EnableIf<IsBlasScalar<Base<F>>>>
HermitianEigInfo
Helper
( UpperOrLower uplo,
  DistMatrix<F,MC,MR,BLOCK>& A,
  DistMatrix<Base<F>,STAR,STAR>& w,
  DistMatrix<F,MC,MR,BLOCK>& Q,
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_CSE
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

#ifdef EL_HAVE_SCALAPACK
    w.Resize( n, 1 );
    Q.SetGrid( A.Grid() );
    Q.AlignWith( A );
    Q.Resize( n, n );

    const int bHandle = blacs::Handle( A );
    const int context = blacs::GridInit( bHandle, A );
    auto descA = FillDesc( A, context );
    auto descQ = FillDesc( Q, context );
    
    const char uploChar = UpperOrLowerToChar( uplo );
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
        scalapack::HermitianEig
        ( uploChar, n,
          A.Buffer(), descA.data(),
          w.Buffer(),
          Q.Buffer(), descQ.data() );
    }
    return info;
#else
    DistMatrix<F,MC,MR> AElem( A ), QElem( A.Grid() );
    info = HermitianEig( uplo, AElem, w, QElem, ctrl );
    Q = QElem;
    return info;
#endif
}

template<typename F,typename=DisableIf<IsBlasScalar<Base<F>>>,typename=void>
HermitianEigInfo
Helper
( UpperOrLower uplo,
  DistMatrix<F,MC,MR,BLOCK>& A,
  DistMatrix<Base<F>,STAR,STAR>& w,
  DistMatrix<F,MC,MR,BLOCK>& Q,
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_CSE
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

    DistMatrix<F,MC,MR> AElem( A ), QElem( A.Grid() );
    info = HermitianEig( uplo, AElem, w, QElem, ctrl );
    Q = QElem;
    return info;
}

} // namespace herm_eig

template<typename F>
HermitianEigInfo
HermitianEig
( UpperOrLower uplo,
  DistMatrix<F,MC,MR,BLOCK>& A,
  DistMatrix<Base<F>,STAR,STAR>& w,
  DistMatrix<F,MC,MR,BLOCK>& Q,
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_CSE
    return herm_eig::Helper( uplo, A, w, Q, ctrl );
}

#define EIGVAL_PROTO(F) \
  template HermitianEigInfo HermitianEig\
  ( UpperOrLower uplo, \
    Matrix<F>& A, \
    Matrix<Base<F>>& w, \
    const HermitianEigCtrl<F>& ctrl ); \
  template HermitianEigInfo HermitianEig\
  ( UpperOrLower uplo, \
    DistMatrix<F,STAR,STAR>& A,\
    DistMatrix<Base<F>,STAR,STAR>& w, \
    const HermitianEigCtrl<F>& ctrl ); \
  template HermitianEigInfo HermitianEig\
  ( UpperOrLower uplo, \
    ElementalMatrix<F>& A, \
    ElementalMatrix<Base<F>>& w, \
    const HermitianEigCtrl<F>& ctrl ); \
  template HermitianEigInfo HermitianEig\
  ( UpperOrLower uplo, \
    DistMatrix<F,MC,MR,BLOCK>& A, \
    DistMatrix<Base<F>,STAR,STAR>& w, \
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
    DistMatrix<F,STAR,STAR>& A,\
    DistMatrix<Base<F>,STAR,STAR>& w, \
    DistMatrix<F,STAR,STAR>& Q,\
    const HermitianEigCtrl<F>& ctrl ); \
  template HermitianEigInfo HermitianEig\
  ( UpperOrLower uplo, \
    ElementalMatrix<F>& A, \
    ElementalMatrix<Base<F>>& w, \
    ElementalMatrix<F>& Q, \
    const HermitianEigCtrl<F>& ctrl ); \
  template HermitianEigInfo HermitianEig\
  ( UpperOrLower uplo, \
    DistMatrix<F,MC,MR,BLOCK>& A, \
    DistMatrix<Base<F>,STAR,STAR>& w, \
    DistMatrix<F,MC,MR,BLOCK>& Q, \
    const HermitianEigCtrl<F>& ctrl );

// Spectral Divide and Conquer
#define SDC_PROTO(F) \
  template void herm_eig::SDC \
  ( UpperOrLower uplo, \
    Matrix<F>& A, \
    Matrix<Base<F>>& w, \
    const HermitianSDCCtrl<Base<F>> ctrl ); \
  template void herm_eig::SDC \
  ( UpperOrLower uplo, \
    ElementalMatrix<F>& A, \
    ElementalMatrix<Base<F>>& w, \
    const HermitianSDCCtrl<Base<F>> ctrl ); \
  template void herm_eig::SDC \
  ( UpperOrLower uplo, \
    Matrix<F>& A, \
    Matrix<Base<F>>& w, \
    Matrix<F>& Q, \
    const HermitianSDCCtrl<Base<F>> ctrl ); \
  template void herm_eig::SDC \
  ( UpperOrLower uplo, \
    ElementalMatrix<F>& A, \
    ElementalMatrix<Base<F>>& w, \
    ElementalMatrix<F>& Q, \
    const HermitianSDCCtrl<Base<F>> ctrl );

#define PROTO(F) \
  EIGVAL_PROTO(F) \
  EIGPAIR_PROTO(F) \
  SDC_PROTO(F)

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
