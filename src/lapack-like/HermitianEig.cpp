/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

#include ELEM_SCALE_INC
#include ELEM_SCALETRAPEZOID_INC

#include ELEM_APPLYPACKEDREFLECTORS_INC
#include ELEM_HERMITIANTRIDIAG_INC
#include ELEM_HERMITIANTRIDIAGEIG_INC
#include ELEM_HERMITIANEIG_INC
#include ELEM_MAXNORM_INC

// The targeted number of pieces to break the eigenvectors into during the
// redistribution from the [* ,VR] distribution after PMRRR to the [MC,MR]
// distribution needed for backtransformation.
#define TARGET_CHUNKS 20

namespace elem {
namespace herm_eig {

// We create specialized redistribution routines for redistributing the 
// real eigenvectors of the symmetric tridiagonal matrix at the core of our 
// eigensolver in order to minimize the temporary memory usage.
template<typename F>
void InPlaceRedist
( DistMatrix<F>& paddedZ, Int rowAlign, const Base<F>* readBuffer )
{
    typedef Base<F> Real;
    const Grid& g = paddedZ.Grid();
    const Int height = paddedZ.Height();
    const Int width = paddedZ.Width();

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = r * c;
    const Int row = g.Row();
    const Int col = g.Col();
    const Int rowShift = paddedZ.RowShift();
    const Int colAlign = paddedZ.ColAlign();
    const Int localWidth = Length(width,g.VRRank(),rowAlign,p);

    const Int maxHeight = MaxLength(height,r);
    const Int maxWidth = MaxLength(width,p);
    const Int portionSize = mpi::Pad( maxHeight*maxWidth );
    
    // Allocate our send/recv buffers
    std::vector<Real> buffer(2*r*portionSize);
    Real* sendBuffer = &buffer[0];
    Real* recvBuffer = &buffer[r*portionSize];

    // Pack
    ELEM_OUTER_PARALLEL_FOR
    for( Int k=0; k<r; ++k )
    {
        Real* data = &sendBuffer[k*portionSize];

        const Int thisColShift = Shift(k,colAlign,r);
        const Int thisLocalHeight = Length(height,thisColShift,r);

        ELEM_INNER_PARALLEL_FOR_COLLAPSE2
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
    ELEM_OUTER_PARALLEL_FOR
    for( Int k=0; k<r; ++k )
    {
        const Real* data = &recvBuffer[k*portionSize];

        const Int thisRank = col+k*c;
        const Int thisRowShift = Shift(thisRank,rowAlign,p);
        const Int thisRowOffset = (thisRowShift-rowShift) / c;
        const Int thisLocalWidth = Length(width,thisRowShift,p);

        ELEM_INNER_PARALLEL_FOR
        for( Int j=0; j<thisLocalWidth; ++j )
        {
            const Real* dataCol = &(data[j*localHeight]);
            Real* thisCol = (Real*)paddedZ.Buffer(0,thisRowOffset+j*r);
            if( IsComplex<F>::val )
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

// Compute the full set of eigenvalues
// ===================================

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, SortType sort,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEig"))
    if( ctrl.useSdc )
    {
        herm_eig::SDC( uplo, A, w, ctrl.sdcCtrl );
    } 
    else
    {
        const Int n = A.Height();
        const char uploChar = UpperOrLowerToChar( uplo );
        w.Resize( n, 1 );
        lapack::HermitianEig( uploChar, n, A.Buffer(), A.LDim(), w.Buffer() );
    }
    Sort( w, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<F,STAR,STAR>& A, DistMatrix<Base<F>,STAR,STAR>& w, SortType sort,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEig"))
    if( ctrl.useSdc )
    {
        herm_eig::SDC( uplo, A.Matrix(), w.Matrix(), ctrl.sdcCtrl );
    }
    else
    {
        const Int n = A.Height();
        const char uploChar = UpperOrLowerToChar( uplo );
        w.Resize( n, 1 );
        lapack::HermitianEig( uploChar, n, A.Buffer(), A.LDim(), w.Buffer() );
    }
    Sort( w, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A,
  DistMatrix<Base<F>,VR,STAR>& w, SortType sort,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEig"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    // Check if we need to rescale the matrix, and do so if necessary
    Base<F> scale;
    const bool needRescaling = herm_eig::CheckScale( uplo, A, scale );
    if( needRescaling )
        ScaleTrapezoid( F(scale), uplo, A );

    if( ctrl.useSdc )
    {
        herm_eig::SDC( uplo, A, w, ctrl.sdcCtrl );
    }
    else
    {
        // Tridiagonalize A
        HermitianTridiag( uplo, A, ctrl.tridiagCtrl );

        // Solve the symmetric tridiagonal EVP
        const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
        auto d = A.GetRealPartOfDiagonal();
        auto e = A.GetRealPartOfDiagonal( subdiagonal );
        HermitianTridiagEig( d, e, w, UNSORTED );
    }
    Sort( w, sort );

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );
}

template<>
void HermitianEig<float>
( UpperOrLower uplo, DistMatrix<float>& A,
  DistMatrix<float,VR,STAR>& w, SortType sort,
  const HermitianEigCtrl<float> ctrl )
{ LogicError("HermitianEig not yet implemented for float"); }

template<>
void HermitianEig<Complex<float>>
( UpperOrLower uplo, DistMatrix<Complex<float>>& A,
  DistMatrix<float,VR,STAR>& w, SortType sort,
  const HermitianEigCtrl<float> ctrl )
{ LogicError("HermitianEig not yet implemented for float"); }

// Compute the full set of eigenpairs
// ==================================

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, 
  Matrix<Base<F>>& w, Matrix<F>& Z, SortType sort,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEig"))
    if( ctrl.useSdc )
    {
        herm_eig::SDC( uplo, A, w, Z, ctrl.sdcCtrl );
    }
    else
    {
        const Int n = A.Height();
        const char uploChar = UpperOrLowerToChar( uplo );
        w.Resize( n, 1 );
        Z.Resize( n, n );
        lapack::HermitianEig
        ( uploChar, n, A.Buffer(), A.LDim(), w.Buffer(), Z.Buffer(), Z.LDim() );
    }
    herm_eig::Sort( w, Z, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<Base<F>,STAR,STAR>& w, DistMatrix<F,STAR,STAR>& Z, SortType sort,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEig"))
    if( ctrl.useSdc )
    {
        herm_eig::SDC( uplo, A.Matrix(), w.Matrix(), Z.Matrix(), ctrl.sdcCtrl );
    }
    else
    {
        const Int n = A.Height();
        const char uploChar = UpperOrLowerToChar( uplo );
        w.Resize( n, 1 );
        Z.Resize( n, n );
        lapack::HermitianEig
        ( uploChar, n, A.Buffer(), A.LDim(), w.Buffer(), Z.Buffer(), Z.LDim() );
    }
    herm_eig::Sort( w.Matrix(), Z.Matrix(), sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A,
  DistMatrix<Base<F>,VR,STAR>& w, DistMatrix<F>& paddedZ, SortType sort,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEig"))
    typedef Base<F> Real;
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    // Check if we need to rescale the matrix, and do so if necessary
    Real scale;
    const bool needRescaling = herm_eig::CheckScale( uplo, A, scale );
    if( needRescaling )
        ScaleTrapezoid( F(scale), uplo, A );

    if( ctrl.useSdc )
    {
        herm_eig::SDC( uplo, A, w, paddedZ, ctrl.sdcCtrl );
    }
    else
    {
        // Tridiagonalize A
        const Grid& g = A.Grid();
        DistMatrix<F,STAR,STAR> t(g);
        HermitianTridiag( uplo, A, t, ctrl.tridiagCtrl );

        // Solve the symmetric tridiagonal EVP. The main subtlety is that we are
        // using a padded [MC,MR] distributed matrix's local buffer to hold the
        // local buffer for a [* ,VR].
        const Int n = A.Height();
        const Int k = n; // full set of eigenpairs
        const Int N = MaxLength(n,g.Height())*g.Height();
        const Int K = MaxLength(k,g.Size())*g.Size(); 
        paddedZ.Align( 0, 0 );
        paddedZ.Resize( N, K );
        DistMatrix<Real,STAR,VR> Z_STAR_VR(g);
        {
            // Grab a slice of size Z_STAR_VR_BufferSize from the very end
            // of paddedZBuffer so that we can later redistribute in place
            Real* paddedZBuf = (Real*)paddedZ.Buffer();
            const Int paddedZBufSize =
                ( IsComplex<F>::val ? 2*paddedZ.LDim()*paddedZ.LocalWidth()
                                    :   paddedZ.LDim()*paddedZ.LocalWidth() );
            const Int Z_STAR_VR_LocalWidth = Length(k,g.VRRank(),g.Size());
            const Int Z_STAR_VR_BufSize = n*Z_STAR_VR_LocalWidth;
            Real* Z_STAR_VR_Buf = &paddedZBuf[paddedZBufSize-Z_STAR_VR_BufSize];
            Z_STAR_VR.Attach( n, k, g, 0, 0, Z_STAR_VR_Buf, n );
        }
        const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
        auto d = A.GetRealPartOfDiagonal();
        auto e = A.GetRealPartOfDiagonal( subdiagonal );
        HermitianTridiagEig( d, e, w, Z_STAR_VR, UNSORTED );
        {
            // Redistribute Z piece-by-piece in place. This is to keep the 
            // send/recv buffer memory usage low.
            const Int p = g.Size();
            const Int numEqualPanels = K/p;
            const Int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
            const Int nbProp = numPanelsPerComm*p;

            // Manually maintain information about the implicit Z[* ,VR] stored 
            // at the end of the paddedZ[MC,MR] buffers.
            Int alignment = 0;
            const Real* readBuffer = Z_STAR_VR.LockedBuffer();
            for( Int j=0; j<k; j+=nbProp )
            {
                const Int nb = Min(nbProp,k-j);
                auto paddedZ1 = View( paddedZ, 0, j, n, nb );

                // Redistribute Z1[MC,MR] <- Z1[* ,VR] in place.
                herm_eig::InPlaceRedist( paddedZ1, alignment, readBuffer );
            
                // Update the Z1[* ,VR] information
                const Int localWidth = nb/p;
                readBuffer = &readBuffer[localWidth*n];
                alignment = (alignment+nb) % p;
            }
        }
        paddedZ.Resize( n, k ); // We can simply shrink matrices

        // Backtransform the tridiagonal eigenvectors, Z
        herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, t, paddedZ );
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );

    herm_eig::Sort( w, paddedZ, sort );
}

template<>
void HermitianEig<float>
( UpperOrLower uplo, DistMatrix<float>& A,
  DistMatrix<float,VR,STAR>& w, DistMatrix<float>& paddedZ, SortType sort,
  const HermitianEigCtrl<float> ctrl )
{ LogicError("HermitianEig not yet implemented for float"); }

template<>
void HermitianEig<Complex<float>>
( UpperOrLower uplo, DistMatrix<Complex<float>>& A,
  DistMatrix<float,VR,STAR>& w, DistMatrix<Complex<float>>& paddedZ,
  SortType sort, const HermitianEigCtrl<float> ctrl )
{ LogicError("HermitianEig not yet implemented for float"); }

// Compute the eigenvalues within an inclusive index range
// =======================================================

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, 
  Matrix<Base<F>>& w, Int il, Int iu, SortType sort,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEig"))
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const Int numEigs = ( n==0 ? 0 : iu-il+1 );
    w.Resize( numEigs, 1 );
    lapack::HermitianEig
    ( uploChar, n, A.Buffer(), A.LDim(), w.Buffer(), il, iu );
    Sort( w, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<Base<F>,STAR,STAR>& w, Int il, Int iu, SortType sort,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEig"))
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const Int numEigs = ( n==0 ? 0 : iu-il+1 );
    w.Resize( numEigs, 1 );
    lapack::HermitianEig
    ( uploChar, n, A.Buffer(), A.LDim(), w.Buffer(), il, iu );
    Sort( w, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A,
  DistMatrix<Base<F>,VR,STAR>& w, Int il, Int iu, SortType sort,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEig"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    // Check if we need to rescale the matrix, and do so if necessary
    Base<F> scale;
    const bool needRescaling = herm_eig::CheckScale( uplo, A, scale );
    if( needRescaling )
        ScaleTrapezoid( F(scale), uplo, A );

    // Tridiagonalize A
    HermitianTridiag( uplo, A, ctrl.tridiagCtrl );

    // Solve the symmetric tridiagonal EVP
    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
    auto d = A.GetRealPartOfDiagonal();
    auto e = A.GetRealPartOfDiagonal( subdiagonal );
    HermitianTridiagEig( d, e, w, il, iu, sort );

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );
}

template<>
void HermitianEig<float>
( UpperOrLower uplo, DistMatrix<float >& A,
  DistMatrix<float,VR,STAR>& w, Int il, Int iu, SortType sort,
  const HermitianEigCtrl<float> ctrl )
{ LogicError("HermitianEig not yet implemented for float"); }

template<>
void HermitianEig<Complex<float>>
( UpperOrLower uplo, DistMatrix<Complex<float>>& A,
  DistMatrix<float,VR,STAR>& w, Int il, Int iu, SortType sort,
  const HermitianEigCtrl<float> ctrl )
{ LogicError("HermitianEig not yet implemented for float"); }

// Compute the eigenpairs within an inclusive index range
// ======================================================

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, 
  Matrix<Base<F>>& w, Matrix<F>& Z, Int il, Int iu, SortType sort,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEig"))
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const Int numEigs = ( n==0 ? 0 : iu-il+1 );
    w.Resize( numEigs, 1 );
    Z.Resize( n, numEigs );
    lapack::HermitianEig
    ( uploChar, n, A.Buffer(), A.LDim(), w.Buffer(), Z.Buffer(), Z.LDim(),
      il, iu );
    herm_eig::Sort( w, Z, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<Base<F>,STAR,STAR>& w, DistMatrix<F,STAR,STAR>& Z, 
  Int il, Int iu, SortType sort,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEig"))
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const Int numEigs = ( n==0 ? 0 : iu-il+1 );
    w.Resize( numEigs, 1 );
    Z.Resize( n, numEigs );
    lapack::HermitianEig
    ( uploChar, n, A.Buffer(), A.LDim(), w.Buffer(), Z.Buffer(), Z.LDim(),
      il, iu ); 
    herm_eig::Sort( w.Matrix(), Z.Matrix(), sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A,
  DistMatrix<Base<F>,VR,STAR>& w, DistMatrix<F>& paddedZ, 
  Int il, Int iu, SortType sort,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEig"))
    typedef Base<F> Real;
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    // Check if we need to rescale the matrix, and do so if necessary
    Real scale;
    const bool needRescaling = herm_eig::CheckScale( uplo, A, scale );
    if( needRescaling )
        ScaleTrapezoid( F(scale), uplo, A );

    // Tridiagonalize A
    const Grid& g = A.Grid();
    DistMatrix<F,STAR,STAR> t(g);
    HermitianTridiag( uplo, A, t, ctrl.tridiagCtrl );
  
    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR] in place, panel by panel.
    // We will use the same buffer for Z in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Z's dimensions slightly.
    const Int n = A.Height();
    const Int k = (iu-il) + 1;
    const Int N = MaxLength(n,g.Height())*g.Height();
    const Int K = MaxLength(k,g.Size())*g.Size(); 
    paddedZ.Align( 0, 0 );
    paddedZ.Resize( N, K );
    DistMatrix<Real,STAR,VR> Z_STAR_VR(g);
    {
        // Grab a slice of size Z_STAR_VR_BufferSize from the very end
        // of paddedZBuffer so that we can later redistribute in place
        Real* paddedZBuf = (Real*)paddedZ.Buffer();
        const Int paddedZBufSize =
            ( IsComplex<F>::val ? 2*paddedZ.LDim()*paddedZ.LocalWidth()
                                :   paddedZ.LDim()*paddedZ.LocalWidth() );
        const Int Z_STAR_VR_LocalWidth = Length(k,g.VRRank(),g.Size());
        const Int Z_STAR_VR_BufSize = n*Z_STAR_VR_LocalWidth;
        Real* Z_STAR_VR_Buf = &paddedZBuf[paddedZBufSize-Z_STAR_VR_BufSize];
        Z_STAR_VR.Attach( n, k, g, 0, 0, Z_STAR_VR_Buf, n );
    }
    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
    auto d = A.GetRealPartOfDiagonal();
    auto e = A.GetRealPartOfDiagonal( subdiagonal );
    HermitianTridiagEig( d, e, w, Z_STAR_VR, il, iu, UNSORTED );
    {
        // Redistribute Z piece-by-piece in place. This is to keep the 
        // send/recv buffer memory usage low.
        const Int p = g.Size();
        const Int numEqualPanels = K/p;
        const Int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const Int nbProp = numPanelsPerComm*p;

        // Manually maintain information about the implicit Z[* ,VR] stored 
        // at the end of the paddedZ[MC,MR] buffers.
        Int alignment = 0;
        const Real* readBuffer = Z_STAR_VR.LockedBuffer();
        for( Int j=0; j<k; j+=nbProp )
        {
            const Int nb = Min(nbProp,k-j);
            auto paddedZ1 = View( paddedZ, 0, j, n, nb );

            // Redistribute Z1[MC,MR] <- Z1[* ,VR] in place.
            herm_eig::InPlaceRedist( paddedZ1, alignment, readBuffer );

            // Update the Z1[* ,VR] information
            const Int localWidth = nb/p;
            readBuffer = &readBuffer[localWidth*n];
            alignment = (alignment+nb) % p;
        }
    }
    paddedZ.Resize( n, k ); // We can simply shrink matrices

    // Backtransform the tridiagonal eigenvectors, Z
    herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, t, paddedZ );

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );

    herm_eig::Sort( w, paddedZ, sort );
}

template<>
void HermitianEig<float>
( UpperOrLower uplo, DistMatrix<float>& A,
  DistMatrix<float,VR,STAR>& w, DistMatrix<float>& paddedZ,
  Int il, Int iu, SortType sort,
  const HermitianEigCtrl<float> ctrl )
{ LogicError("HermitianEig not yet implemented for float"); }

template<>
void HermitianEig<Complex<float>>
( UpperOrLower uplo, DistMatrix<Complex<float>>& A,
  DistMatrix<float,VR,STAR>& w, DistMatrix<Complex<float>>& paddedZ,
  Int il, Int iu, SortType sort,
  const HermitianEigCtrl<float> ctrl )
{ LogicError("HermitianEig not yet implemented for float"); }

// Compute the eigenvalues within a half-open, (a,b], interval
// ===========================================================

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, 
  Matrix<Base<F>>& w, Base<F> vl, Base<F> vu, SortType sort,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEig"))
    if( vl >= vu )
    {
        w.Resize(0,1);
        return; 
    }

    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    w.Resize( n, 1 );
    const Int numEigs = lapack::HermitianEig
    ( uploChar, n, A.Buffer(), A.LDim(), w.Buffer(), vl, vu );
    w.Resize( numEigs, 1 );
    Sort( w, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<Base<F>,STAR,STAR>& w, Base<F> vl, Base<F> vu, SortType sort,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEig"))
    if( vl >= vu )
    {
        w.Resize(0,1);
        return; 
    }

    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    w.Resize( n, 1 );
    const Int numEigs = lapack::HermitianEig
    ( uploChar, n, A.Buffer(), A.LDim(), w.Buffer(), vl, vu );
    w.Resize( numEigs, 1 );
    Sort( w, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A,
  DistMatrix<Base<F>,VR,STAR>& w, Base<F> vl, Base<F> vu,
  SortType sort, const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEig"))
    if( vl >= vu )
    {
        w.Resize(0,1);
        return; 
    }
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    // Check if we need to rescale the matrix, and do so if necessary
    Base<F> scale;
    const bool needRescaling = herm_eig::CheckScale( uplo, A, scale );
    if( needRescaling )
        ScaleTrapezoid( F(scale), uplo, A );

    // Tridiagonalize A
    HermitianTridiag( uplo, A, ctrl.tridiagCtrl );

    // Solve the symmetric tridiagonal EVP
    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
    auto d = A.GetRealPartOfDiagonal();
    auto e = A.GetRealPartOfDiagonal( subdiagonal );
    HermitianTridiagEig( d, e, w, vl, vu, sort );

    // Rescale the eigenvalues if necessary
    if( needRescaling ) 
        Scale( 1/scale, w );
}

template<>
void HermitianEig<float>
( UpperOrLower uplo, DistMatrix<float>& A,
  DistMatrix<float,VR,STAR>& w, float vl, float vu,
  SortType sort, const HermitianEigCtrl<float> ctrl )
{ LogicError("HermitianEig not yet implemented for float"); }

template<>
void HermitianEig<Complex<float>>
( UpperOrLower uplo, DistMatrix<Complex<float>>& A,
  DistMatrix<float,VR,STAR>& w, float vl, float vu,
  SortType sort, const HermitianEigCtrl<float> ctrl )
{ LogicError("HermitianEig not yet implemented for float"); }

// Compute the eigenpairs with eigenvalues in a half-open, (a,b], interval
// =======================================================================

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, Matrix<F>& Z, 
  Base<F> vl, Base<F> vu, SortType sort,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEig"))
    const Int n = A.Height();
    if( vl >= vu )
    {
        w.Resize(0,1);
        Z.Resize(n,0);
        return; 
    }

    const char uploChar = UpperOrLowerToChar( uplo );
    w.Resize( n, 1 );
    Z.Resize( n, n );
    const Int numEigs = lapack::HermitianEig
    ( uploChar, n, A.Buffer(), A.LDim(), w.Buffer(), Z.Buffer(), Z.LDim(),
      vl, vu );
    w.Resize( numEigs, 1 );
    Z.Resize( n, numEigs );
    herm_eig::Sort( w, Z, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<Base<F>,STAR,STAR>& w, 
  DistMatrix<F,STAR,STAR>& Z, 
  Base<F> vl, Base<F> vu, SortType sort,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEig"))
    const Int n = A.Height();
    if( vl >= vu )
    {
        w.Resize(0,1);
        Z.Resize(n,0);
        return; 
    }

    const char uploChar = UpperOrLowerToChar( uplo );
    w.Resize( n, 1 );
    Z.Resize( n, n );
    const Int numEigs = lapack::HermitianEig
    ( uploChar, n, A.Buffer(), A.LDim(), w.Buffer(), Z.Buffer(), Z.LDim(),
      vl, vu );
    w.Resize( numEigs, 1 );
    Z.Resize( n, numEigs );
    herm_eig::Sort( w.Matrix(), Z.Matrix(), sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A,
  DistMatrix<Base<F>,VR,STAR>& w, DistMatrix<F>& paddedZ,
  Base<F> vl, Base<F> vu, SortType sort, const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEig"))
    typedef Base<F> Real;
    const Int n = A.Height();
    if( vl >= vu )
    {
        w.Resize(0,1);
        paddedZ.Resize(n,0);
        return; 
    }
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    // Check if we need to rescale the matrix, and do so if necessary
    Real scale;
    const bool needRescaling = herm_eig::CheckScale( uplo, A, scale );
    if( needRescaling )
        ScaleTrapezoid( F(scale), uplo, A );

    // Tridiagonalize A
    const Grid& g = A.Grid();
    DistMatrix<F,STAR,STAR> t(g);
    HermitianTridiag( uplo, A, t, ctrl.tridiagCtrl );

    // Get an upper-bound on the number of local eigenvalues in the range
    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
    auto d = A.GetRealPartOfDiagonal();
    auto e = A.GetRealPartOfDiagonal( subdiagonal );
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<Real,STAR,STAR> e_STAR_STAR( g );
    e_STAR_STAR.Resize( n-1, 1, n );
    e_STAR_STAR = e;
    const Int kEst = HermitianTridiagEigEstimate( d, e, g.VRComm(), vl, vu );

    // We will use the same buffer for Z in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Z's dimensions slightly.
    const Int N = MaxLength(n,g.Height())*g.Height();
    const Int K = MaxLength(kEst,g.Size())*g.Size(); 
    paddedZ.Align( 0, 0 );
    paddedZ.Resize( N, K );
    DistMatrix<Real,STAR,VR> Z_STAR_VR(g);
    {
        // Grab a slice of size Z_STAR_VR_BufferSize from the very end
        // of paddedZBuffer so that we can later redistribute in place
        Real* paddedZBuf = (Real*)paddedZ.Buffer();
        const Int paddedZBufSize =
            ( IsComplex<F>::val ? 2*paddedZ.LDim()*paddedZ.LocalWidth()
                                :   paddedZ.LDim()*paddedZ.LocalWidth() );
        const Int Z_STAR_VR_LocalWidth = Length(kEst,g.VRRank(),g.Size());
        const Int Z_STAR_VR_BufSize = n*Z_STAR_VR_LocalWidth;
        Real* Z_STAR_VR_Buf = &paddedZBuf[paddedZBufSize-Z_STAR_VR_BufSize];
        Z_STAR_VR.Attach( n, kEst, g, 0, 0, Z_STAR_VR_Buf, n );
    }
    HermitianTridiagEigPostEstimate
    ( d_STAR_STAR, e_STAR_STAR, w, Z_STAR_VR, vl, vu, UNSORTED );
    const Int k = w.Height();
    {
        // Redistribute Z piece-by-piece in place. This is to keep the 
        // send/recv buffer memory usage low.
        const Int p = g.Size();
        const Int numEqualPanels = K/p;
        const Int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const Int nbProp = numPanelsPerComm*p;

        // Manually maintain information about the implicit Z[* ,VR] stored 
        // at the end of the paddedZ[MC,MR] buffers.
        Int alignment = 0;
        const Real* readBuffer = Z_STAR_VR.LockedBuffer();
        for( Int j=0; j<k; j+=nbProp )
        {
            const Int nb = Min(nbProp,k-j);
            auto paddedZ1 = View( paddedZ, 0, j, n, nb );

            // Redistribute Z1[MC,MR] <- Z1[* ,VR] in place.
            herm_eig::InPlaceRedist( paddedZ1, alignment, readBuffer );

            // Update the Z1[* ,VR] information
            const Int localWidth = nb/p;
            readBuffer = &readBuffer[localWidth*n];
            alignment = (alignment+nb) % p;
        }
    }
    paddedZ.Resize( n, k ); // We can simply shrink matrices

    // Backtransform the tridiagonal eigenvectors, Z
    herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, t, paddedZ );

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );

    herm_eig::Sort( w, paddedZ, sort );
}

template<>
void HermitianEig<float>
( UpperOrLower uplo, DistMatrix<float>& A,
  DistMatrix<float,VR,STAR>& w, DistMatrix<float>& paddedZ,
  float vl, float vu, SortType sort,
  const HermitianEigCtrl<float> ctrl )
{ LogicError("HermitianEig not yet implemented for float"); }

template<>
void HermitianEig<Complex<float>>
( UpperOrLower uplo, DistMatrix<Complex<float>>& A,
  DistMatrix<float,VR,STAR>& w, DistMatrix<Complex<float>>& paddedZ,
  float vl, float vu, SortType sort,
  const HermitianEigCtrl<float> ctrl )
{ LogicError("HermitianEig not yet implemented for float"); }

// Full set of eigenvalues
#define FULL_EIGVAL(F) \
  template void HermitianEig\
  ( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, SortType sort, \
    const HermitianEigCtrl<Base<F>> sdcCtrl ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,\
    DistMatrix<Base<F>,STAR,STAR>& w, SortType sort, \
    const HermitianEigCtrl<Base<F>> sdcCtrl ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w,\
    SortType sort, const HermitianEigCtrl<Base<F>> sdcCtrl );
// Full set of eigenpairs
#define FULL_EIGPAIR(F) \
  template void HermitianEig\
  ( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, Matrix<F>& Z,\
    SortType sort, const HermitianEigCtrl<Base<F>> sdcCtrl ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,\
    DistMatrix<Base<F>,STAR,STAR>& w, DistMatrix<F,STAR,STAR>& Z,\
    SortType sort, const HermitianEigCtrl<Base<F>> sdcCtrl ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w,\
    DistMatrix<F>& Z, SortType sort, \
    const HermitianEigCtrl<Base<F>> sdcCtrl );
// Integer range of eigenvalues
#define INT_EIGVAL(F) \
  template void HermitianEig\
  ( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, Int il, Int iu,\
    SortType sort, const HermitianEigCtrl<Base<F>> sdcCtrl ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,\
    DistMatrix<Base<F>,STAR,STAR>& w, Int il, Int iu, SortType sort, \
    const HermitianEigCtrl<Base<F>> sdcCtrl ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w,\
    Int il, Int iu, SortType sort, \
    const HermitianEigCtrl<Base<F>> sdcCtrl );
// Integer range of eigenpairs
#define INT_EIGPAIR(F) \
  template void HermitianEig\
  ( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, Matrix<F>& Z,\
    Int il, Int iu, SortType sort, const HermitianEigCtrl<Base<F>> sdcCtrl ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,\
    DistMatrix<Base<F>,STAR,STAR>& w, DistMatrix<F,STAR,STAR>& Z,\
    Int il, Int iu, SortType sort, const HermitianEigCtrl<Base<F>> sdcCtrl ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w,\
    DistMatrix<F>& Z, Int il, Int iu, SortType sort, \
    const HermitianEigCtrl<Base<F>> sdcCtrl );
// Floating-point range of eigenvalues
#define FLOAT_EIGVAL(F) \
  template void HermitianEig\
  ( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w,\
    Base<F> vl, Base<F> vu, SortType sort, \
    const HermitianEigCtrl<Base<F>> sdcCtrl ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,\
    DistMatrix<Base<F>,STAR,STAR>& w,\
    Base<F> vl, Base<F> vu, SortType sort, \
    const HermitianEigCtrl<Base<F>> sdcCtrl ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w,\
    Base<F> vl, Base<F> vu, SortType sort, \
    const HermitianEigCtrl<Base<F>> sdcCtrl );
// Floating-point range of eigenpairs
#define FLOAT_EIGPAIR(F) \
  template void HermitianEig\
  ( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, Matrix<F>& Z,\
    Base<F> vl, Base<F> vu, SortType sort, \
    const HermitianEigCtrl<Base<F>> sdcCtrl ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,\
    DistMatrix<Base<F>,STAR,STAR>& w, DistMatrix<F,STAR,STAR>& Z,\
    Base<F> vl, Base<F> vu, SortType sort, \
    const HermitianEigCtrl<Base<F>> sdcCtrl ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w,\
    DistMatrix<F>& Z, Base<F> vl, Base<F> vu, SortType sort, \
    const HermitianEigCtrl<Base<F>> sdcCtrl ); \
// All options
#define ALL_OPTS(F) \
  FULL_EIGVAL(F);\
  FULL_EIGPAIR(F);\
  INT_EIGVAL(F);\
  INT_EIGPAIR(F);\
  FLOAT_EIGVAL(F);\
  FLOAT_EIGPAIR(F);

#ifndef ELEM_DISABLE_FLOAT
ALL_OPTS(float);
#ifndef ELEM_DISABLE_COMPLEX
ALL_OPTS(Complex<float>);
#endif // ifndef ELEM_DISABLE_COMPLEX
#endif // ifndef ELEM_DISABLE_FLOAT
ALL_OPTS(double);
#ifndef ELEM_DISABLE_COMPLEX
ALL_OPTS(Complex<double>);
#endif // ifndef ELEM_DISABLE_COMPLEX

} // namespace elem
