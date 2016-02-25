/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

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
void InPlaceRedist( DistMatrix<F>& Z, Int rowAlign, const Base<F>* readBuffer )
{
    typedef Base<F> Real;
    const Grid& g = Z.Grid();
    const Int height = Z.Height();
    const Int width = Z.Width();

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = r * c;
    const Int row = g.Row();
    const Int col = g.Col();
    const Int rowShift = Z.RowShift();
    const Int colAlign = Z.ColAlign();
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
            Real* thisCol = (Real*)Z.Buffer(0,thisRowOffset+j*r);
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
    const Real underflowThreshold = lapack::MachineUnderflowThreshold<Real>();
    const Real overflowThreshold = lapack::MachineOverflowThreshold<Real>();
    if( maxNormOfA > 0 && maxNormOfA < underflowThreshold )
    {
        scale = underflowThreshold / maxNormOfA;
        return true;
    }
    else if( maxNormOfA > overflowThreshold )
    {
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
void HermitianEig
( UpperOrLower uplo,
  Matrix<F>& A, 
  Matrix<Base<F>>& w,
  SortType sort,
  const HermitianEigSubset<Base<F>> subset,
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_ONLY(CSE cse("HermitianEig"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");
    if( ctrl.useSDC )
    {
        herm_eig::SDC( uplo, A, w, ctrl.sdcCtrl );
        return;
    }

    if( subset.indexSubset && subset.rangeSubset )
        LogicError("Cannot mix index and range subsets");
    if( (subset.rangeSubset && (subset.lowerBound >= subset.upperBound)) ||
        (subset.indexSubset && (subset.lowerIndex > subset.upperIndex)) )
    {
        w.Resize(0,1);
        return;
    }

    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    w.Resize( n, 1 );
    if( subset.rangeSubset )
    {
        const Int numEigs = lapack::HermitianEig
          ( uploChar, BlasInt(n), A.Buffer(), BlasInt(A.LDim()), w.Buffer(), 
            subset.lowerBound, subset.upperBound );
        w.Resize( numEigs, 1 );
    }
    else if( subset.indexSubset )
    {
        lapack::HermitianEig
        ( uploChar, BlasInt(n), A.Buffer(), BlasInt(A.LDim()), w.Buffer(), 
          BlasInt(subset.lowerIndex), BlasInt(subset.upperIndex) );
        w.Resize( subset.upperIndex-subset.lowerIndex+1, 1 );
    }
    else
        lapack::HermitianEig
        ( uploChar, BlasInt(n), A.Buffer(), BlasInt(A.LDim()), w.Buffer() );
    Sort( w, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<Base<F>,STAR,STAR>& w,
  SortType sort,
  const HermitianEigSubset<Base<F>> subset,
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_ONLY(CSE cse("HermitianEig"))
    const Int n = A.Height();
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");
    if( ctrl.useSDC )
    {
        w.SetGrid( A.Grid() );
        w.Resize( n, 1 );
        herm_eig::SDC( uplo, A.Matrix(), w.Matrix(), ctrl.sdcCtrl );
        return;
    }

    if( subset.indexSubset && subset.rangeSubset )
        LogicError("Cannot mix index and range subsets");
    if( (subset.rangeSubset && (subset.lowerBound >= subset.upperBound)) ||
        (subset.indexSubset && (subset.lowerIndex > subset.upperIndex)) )
    {
        w.Resize(0,1);
        return;
    }

    const char uploChar = UpperOrLowerToChar( uplo );
    w.Resize( n, 1 );
    if( subset.rangeSubset )
    {
        const Int numEigs = lapack::HermitianEig
          ( uploChar, BlasInt(n), A.Buffer(), BlasInt(A.LDim()), w.Buffer(), 
            subset.lowerBound, subset.upperBound );
        w.Resize( numEigs, 1 );
    }
    else if( subset.indexSubset )
    {
        lapack::HermitianEig
        ( uploChar, BlasInt(n), A.Buffer(), BlasInt(A.LDim()), w.Buffer(), 
          BlasInt(subset.lowerIndex), BlasInt(subset.upperIndex) );
        w.Resize( subset.upperIndex-subset.lowerIndex+1, 1 );
    }
    else
        lapack::HermitianEig
        ( uploChar, BlasInt(n), A.Buffer(), BlasInt(A.LDim()), w.Buffer() );
    Sort( w, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo,
  ElementalMatrix<F>& APre,
  ElementalMatrix<Base<F>>& w,
  SortType sort, 
  const HermitianEigSubset<Base<F>> subset,
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_ONLY(CSE cse("HermitianEig"))
    if( APre.Height() != APre.Width() )
        LogicError("Hermitian matrices must be square");

    if( ctrl.useSDC )
    {
        herm_eig::SDC( uplo, APre, w, ctrl.sdcCtrl );
        return;
    }

    if( subset.indexSubset && subset.rangeSubset )
        LogicError("Cannot mix index and range subsets");
    if( (subset.rangeSubset && (subset.lowerBound >= subset.upperBound)) ||
        (subset.indexSubset && (subset.lowerIndex > subset.upperIndex)) )
    {
        w.Resize(0,1);
        return;
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
            cout << "  Condense time:      " << timer.Stop() << " secs" << endl;
            timer.Start();
        }
    }

    // Solve the symmetric tridiagonal EVP
    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
    auto d = GetRealPartOfDiagonal(A);
    auto e = GetRealPartOfDiagonal(A,subdiagonal);
    HermitianTridiagEig( d, e, w, sort, subset );

    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
            cout << "  TridiagEig time:    " << timer.Stop() << " secs" << endl;
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling ) 
        w *= 1/scale;
}

template<typename F>
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<F,MC,MR,BLOCK>& A,
  Matrix<Base<F>>& w,
  const HermitianEigSubset<Base<F>> subset )
{
    DEBUG_ONLY(CSE cse("HermitianEig"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    if( subset.indexSubset && subset.rangeSubset )
        LogicError("Cannot mix index and range subsets");
    if( (subset.rangeSubset && (subset.lowerBound >= subset.upperBound)) ||
        (subset.indexSubset && (subset.lowerIndex > subset.upperIndex)) )
    {
        w.Resize(0,1);
        return;
    }

    AssertScaLAPACKSupport();
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
#endif
}

// Compute eigenpairs
// ==================

template<typename F>
void HermitianEig
( UpperOrLower uplo,
  Matrix<F>& A,
  Matrix<Base<F>>& w,
  Matrix<F>& Z, 
  SortType sort,
  const HermitianEigSubset<Base<F>> subset,
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_ONLY(CSE cse("HermitianEig"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");
    if( ctrl.useSDC )
    {
        herm_eig::SDC( uplo, A, w, Z, ctrl.sdcCtrl );
        herm_eig::Sort( w, Z, sort );
        return;
    }

    if( subset.indexSubset && subset.rangeSubset )
        LogicError("Cannot mix index and range subsets");
    const Int n = A.Height();
    if( (subset.rangeSubset && (subset.lowerBound >= subset.upperBound)) ||
        (subset.indexSubset && (subset.lowerIndex > subset.upperIndex)) )
    {
        w.Resize(0,1);
        Z.Resize(n,0);
        return; 
    }

    const char uploChar = UpperOrLowerToChar( uplo );
    w.Resize( n, 1 );
    if( subset.indexSubset )
    {
        const Int numEigs = subset.upperIndex-subset.lowerIndex+1;
        Z.Resize( n, numEigs );
        lapack::HermitianEig
        ( uploChar, BlasInt(n), A.Buffer(), BlasInt(A.LDim()), w.Buffer(), 
          Z.Buffer(), BlasInt(Z.LDim()),
          BlasInt(subset.lowerIndex), BlasInt(subset.upperIndex) );
        w.Resize( numEigs, 1 );
    }
    else if( subset.rangeSubset )
    {
        Z.Resize( n, n );
        const Int numEigs = lapack::HermitianEig
          ( uploChar, BlasInt(n), A.Buffer(), BlasInt(A.LDim()), w.Buffer(), 
            Z.Buffer(), BlasInt(Z.LDim()),
            subset.lowerBound, subset.upperBound );
        w.Resize( numEigs, 1 );
        Z.Resize( n, numEigs );
    }
    else
    {
        Z.Resize( n, n );
        lapack::HermitianEig
        ( uploChar, BlasInt(n), A.Buffer(), BlasInt(A.LDim()), w.Buffer(), 
          Z.Buffer(), BlasInt(Z.LDim()) );
    }
    herm_eig::Sort( w, Z, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<Base<F>,STAR,STAR>& w,
  DistMatrix<F,STAR,STAR>& Z, 
  SortType sort,
  const HermitianEigSubset<Base<F>> subset,
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_ONLY(CSE cse("HermitianEig"))
    const Int n = A.Height();
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");
    if( ctrl.useSDC )
    {
        w.Resize(n,1);
        Z.Resize(n,n);
        herm_eig::SDC( uplo, A.Matrix(), w.Matrix(), Z.Matrix(), ctrl.sdcCtrl );
        herm_eig::Sort( w.Matrix(), Z.Matrix(), sort );
        return;
    }

    if( subset.indexSubset && subset.rangeSubset )
        LogicError("Cannot mix index and range subsets");
    if( (subset.rangeSubset && (subset.lowerBound >= subset.upperBound)) ||
        (subset.indexSubset && (subset.lowerIndex > subset.upperIndex)) )
    {
        w.Resize(0,1);
        Z.Resize(n,0);
        return; 
    }

    const char uploChar = UpperOrLowerToChar( uplo );
    w.Resize( n, 1 );
    if( subset.indexSubset )
    {
        const Int numEigs = subset.upperIndex-subset.lowerIndex+1;
        Z.Resize( n, numEigs );
        lapack::HermitianEig
        ( uploChar, BlasInt(n), A.Buffer(), A.LDim(), w.Buffer(), 
          Z.Buffer(), BlasInt(Z.LDim()),
          BlasInt(subset.lowerIndex), BlasInt(subset.upperIndex) );
        w.Resize( numEigs, 1 );
    }
    else if( subset.rangeSubset )
    {
        Z.Resize( n, n );
        const Int numEigs = lapack::HermitianEig
          ( uploChar, BlasInt(n), A.Buffer(), A.LDim(), w.Buffer(), 
            Z.Buffer(), BlasInt(Z.LDim()),
            subset.lowerBound, subset.upperBound );
        w.Resize( numEigs, 1 );
        Z.Resize( n, numEigs );
    }
    else
    {
        Z.Resize( n, n );
        lapack::HermitianEig
        ( uploChar, BlasInt(n), A.Buffer(), BlasInt(A.LDim()), w.Buffer(), 
          Z.Buffer(), BlasInt(Z.LDim()) );
    }
    herm_eig::Sort( w.Matrix(), Z.Matrix(), sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo,
  ElementalMatrix<F>& APre,
  ElementalMatrix<Base<F>>& w,
  ElementalMatrix<F>& ZPre,
  SortType sort,
  const HermitianEigSubset<Base<F>> subset,
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_ONLY(CSE cse("HermitianEig"))
    typedef Base<F> Real;
    const Int n = APre.Height();
    if( APre.Height() != APre.Width() )
        LogicError("Hermitian matrices must be square");

    if( ctrl.useSDC )
    {
        herm_eig::SDC( uplo, APre, w, ZPre, ctrl.sdcCtrl );
        herm_eig::Sort( w, ZPre, sort );
        return;
    }

    if( subset.indexSubset && subset.rangeSubset )
        LogicError("Cannot mix index and range subsets");
    if( (subset.rangeSubset && (subset.lowerBound >= subset.upperBound)) ||
        (subset.indexSubset && (subset.lowerIndex > subset.upperIndex)) )
    {
        w.Resize(0,1);
        ZPre.Resize(n,0);
        return; 
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
    DistMatrix<F,STAR,STAR> t(g);
    HermitianTridiag( uplo, A, t, ctrl.tridiagCtrl );

    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
        {
            cout << "  Condense time:      " << timer.Stop() << " secs" << endl;
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

    // We will use the same buffer for Z in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Z's dimensions slightly.
    const Int N = MaxLength(n,g.Height())*g.Height();
    const Int K = MaxLength(kEst,g.Size())*g.Size(); 

    ElementalProxyCtrl proxCtrl;
    proxCtrl.colConstrain = true;
    proxCtrl.rowConstrain = true;
    proxCtrl.colAlign = 0;
    proxCtrl.rowAlign = 0;
    
    DistMatrixWriteProxy<F,F,MC,MR> ZProx( ZPre, proxCtrl );
    auto& Z = ZProx.Get();

    Z.Resize( N, K );
    DistMatrix<Real,STAR,VR> Z_STAR_VR(g);
    {
        // Grab a slice of size Z_STAR_VR_BufferSize from the very end
        // of ZBuf so that we can later redistribute in place
        Real* ZBuf = (Real*)Z.Buffer();
        const Int ZBufSize =
            ( IsComplex<F>::value ? 2*Z.LDim()*Z.LocalWidth()
                                  :   Z.LDim()*Z.LocalWidth() );
        const Int Z_STAR_VR_LocalWidth = Length(kEst,g.VRRank(),g.Size());
        const Int Z_STAR_VR_BufSize = n*Z_STAR_VR_LocalWidth;
        Real* Z_STAR_VR_Buf = &ZBuf[ZBufSize-Z_STAR_VR_BufSize];
        Z_STAR_VR.Attach( n, kEst, g, 0, 0, Z_STAR_VR_Buf, n );
    }
    // NOTE: We should be guaranteeing that Z_STAR_VR does not need to
    //       reallocate a buffer
    if( subset.rangeSubset )
        HermitianTridiagEigPostEstimate
        ( d_STAR_STAR, e_STAR_STAR, w, Z_STAR_VR, UNSORTED,
          subset.lowerBound, subset.upperBound );
    else
        HermitianTridiagEig
        ( d_STAR_STAR, e_STAR_STAR, w, Z_STAR_VR, UNSORTED, subset );

    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
        {
            cout << "  TridiagEig time:    " << timer.Stop() << " secs" << endl;
            timer.Start();
        }
    }

    const Int k = w.Height();
    {
        // Redistribute Z piece-by-piece in place. This is to keep the 
        // send/recv buffer memory usage low.
        const Int p = g.Size();
        const Int numEqualPanels = K/p;
        const Int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const Int nbProp = numPanelsPerComm*p;

        // Manually maintain information about the implicit Z[* ,VR] stored 
        // at the end of the Z[MC,MR] buffers.
        Int alignment = 0;
        const Real* readBuffer = Z_STAR_VR.LockedBuffer();
        for( Int j=0; j<k; j+=nbProp )
        {
            const Int nb = Min(nbProp,k-j);
            auto Z1 = Z( IR(0,n), IR(j,j+nb) );

            // Redistribute Z1[MC,MR] <- Z1[* ,VR] in place.
            // NOTE: This assumes that Z_STAR_VR did not reallocate within
            //       HermitianTridiagEig[PostEstimate] above
            herm_eig::InPlaceRedist( Z1, alignment, readBuffer );

            // Update the Z1[* ,VR] information
            const Int localWidth = nb/p;
            readBuffer = &readBuffer[localWidth*n];
            alignment = (alignment+nb) % p;
        }
    }
    Z.Resize( n, k ); // We can simply shrink matrices

    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
        {
            cout << "  Redist time:        " << timer.Stop() << " secs" << endl;
            timer.Start();
        }
    }

    // Backtransform the tridiagonal eigenvectors, Z
    herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, t, Z );

    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
        {
            cout << "  Backtransform time: " << timer.Stop() << " secs" << endl;
            timer.Start();
        }
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        w *= 1/scale;

    herm_eig::Sort( w, Z, sort );

    if( ctrl.timeStages )
    {
        mpi::Barrier( A.DistComm() );
        if( A.Grid().Rank() == 0 )
            cout << "  Scale+sort time:    " << timer.Stop() << " secs" << endl;
    }
}

template<typename F>
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<F,MC,MR,BLOCK>& A,
  Matrix<Base<F>>& w,
  DistMatrix<F,MC,MR,BLOCK>& Z,
  const HermitianEigSubset<Base<F>> subset )
{
    DEBUG_ONLY(CSE cse("HermitianEig"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    const Int n = A.Height();
    if( subset.indexSubset && subset.rangeSubset )
        LogicError("Cannot mix index and range subsets");
    if( (subset.rangeSubset && (subset.lowerBound >= subset.upperBound)) ||
        (subset.indexSubset && (subset.lowerIndex > subset.upperIndex)) )
    {
        w.Resize(0,1);
        Z.Resize(n,0);
        return;
    }

    AssertScaLAPACKSupport();
#ifdef EL_HAVE_SCALAPACK
    w.Resize( n, 1 );
    Z.SetGrid( A.Grid() );
    Z.AlignWith( A );
    Z.Resize( n, n );

    const int bHandle = blacs::Handle( A );
    const int context = blacs::GridInit( bHandle, A );
    auto descA = FillDesc( A, context );
    auto descZ = FillDesc( Z, context );
    
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
          Z.Buffer(), descZ.data() );
    }
#endif
}

#define EIGVAL_PROTO(F) \
  template void HermitianEig\
  ( UpperOrLower uplo, \
    Matrix<F>& A, \
    Matrix<Base<F>>& w, \
    SortType sort, \
    const HermitianEigSubset<Base<F>> subset, \
    const HermitianEigCtrl<F>& ctrl ); \
  template void HermitianEig\
  ( UpperOrLower uplo, \
    DistMatrix<F,STAR,STAR>& A,\
    DistMatrix<Base<F>,STAR,STAR>& w, \
    SortType sort, \
    const HermitianEigSubset<Base<F>> subset, \
    const HermitianEigCtrl<F>& ctrl ); \
  template void HermitianEig\
  ( UpperOrLower uplo, \
    ElementalMatrix<F>& A, \
    ElementalMatrix<Base<F>>& w, \
    SortType sort, \
    const HermitianEigSubset<Base<F>> subset, \
    const HermitianEigCtrl<F>& ctrl ); \
  template void HermitianEig\
  ( UpperOrLower uplo, \
    DistMatrix<F,MC,MR,BLOCK>& A, \
    Matrix<Base<F>>& w, \
    const HermitianEigSubset<Base<F>> subset );

#define EIGPAIR_PROTO(F) \
  template void HermitianEig\
  ( UpperOrLower uplo, \
    Matrix<F>& A, \
    Matrix<Base<F>>& w, \
    Matrix<F>& Z,\
    SortType sort, \
    const HermitianEigSubset<Base<F>> subset, \
    const HermitianEigCtrl<F>& ctrl ); \
  template void HermitianEig\
  ( UpperOrLower uplo, \
    DistMatrix<F,STAR,STAR>& A,\
    DistMatrix<Base<F>,STAR,STAR>& w, \
    DistMatrix<F,STAR,STAR>& Z,\
    SortType sort, \
    const HermitianEigSubset<Base<F>> subset, \
    const HermitianEigCtrl<F>& ctrl ); \
  template void HermitianEig\
  ( UpperOrLower uplo, \
    ElementalMatrix<F>& A, \
    ElementalMatrix<Base<F>>& w, \
    ElementalMatrix<F>& Z, \
    SortType sort, \
    const HermitianEigSubset<Base<F>> subset, \
    const HermitianEigCtrl<F>& ctrl ); \
  template void HermitianEig\
  ( UpperOrLower uplo, \
    DistMatrix<F,MC,MR,BLOCK>& A, \
    Matrix<Base<F>>& w, \
    DistMatrix<F,MC,MR,BLOCK>& Z, \
    const HermitianEigSubset<Base<F>> subset );

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
#include "El/macros/Instantiate.h"

} // namespace El
