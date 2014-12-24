/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T,Dist U,Dist V>
GeneralBlockDistMatrix<T,U,V>::GeneralBlockDistMatrix
( const El::Grid& g, Int root )
: AbstractBlockDistMatrix<T>(g,root)
{ }

template<typename T,Dist U,Dist V>
GeneralBlockDistMatrix<T,U,V>::GeneralBlockDistMatrix
( const El::Grid& g, Int blockHeight, Int blockWidth, Int root )
: AbstractBlockDistMatrix<T>(g,blockHeight,blockWidth,root)
{ }

template<typename T,Dist U,Dist V>
GeneralBlockDistMatrix<T,U,V>::GeneralBlockDistMatrix
( GeneralBlockDistMatrix<T,U,V>&& A ) EL_NOEXCEPT
: AbstractBlockDistMatrix<T>(std::move(A))
{ }

// Assignment and reconfiguration
// ==============================

template<typename T,Dist U,Dist V>
GeneralBlockDistMatrix<T,U,V>& 
GeneralBlockDistMatrix<T,U,V>::operator=( GeneralBlockDistMatrix<T,U,V>&& A )
{
    AbstractBlockDistMatrix<T>::operator=( std::move(A) );
    return *this;
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::RowSumScatterFrom
( const BlockDistMatrix<T,U,VGath>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::RowSumScatterFrom");
        AssertSameGrids( *this, A );
    )
    this->AlignColsAndResize
    ( A.BlockHeight(), A.ColAlign(), A.ColCut(), A.Height(), A.Width(), 
      false, false );
    // NOTE: This will be *slightly* slower than necessary due to the result
    //       of the MPI operations being added rather than just copied
    Zeros( this->Matrix(), this->LocalHeight(), this->LocalWidth() );
    this->RowSumScatterUpdate( T(1), A );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::ColSumScatterFrom
( const BlockDistMatrix<T,UGath,V>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::ColSumScatterFrom");
        AssertSameGrids( *this, A );
    )
    this->AlignRowsAndResize
    ( A.BlockWidth(), A.RowAlign(), A.RowCut(), A.Height(), A.Width(),
      false, false );
    // NOTE: This will be *slightly* slower than necessary due to the result
    //       of the MPI operations being added rather than just copied
    Zeros( this->Matrix(), this->LocalHeight(), this->LocalWidth() );
    this->ColSumScatterUpdate( T(1), A );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::SumScatterFrom
( const BlockDistMatrix<T,UGath,VGath>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::SumScatterFrom");
        AssertSameGrids( *this, A );
    )
    this->Resize( A.Height(), A.Width() );
    // NOTE: This will be *slightly* slower than necessary due to the result
    //       of the MPI operations being added rather than just copied
    Zeros( this->Matrix(), this->LocalHeight(), this->LocalWidth() );
    this->SumScatterUpdate( T(1), A );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::PartialRowSumScatterFrom
( const BlockDistMatrix<T,U,VPart>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::PartialRowSumScatterFrom");
        AssertSameGrids( *this, A );
    )
    this->AlignAndResize
    ( A.BlockHeight(), A.BlockWidth(), 
      A.ColAlign(), A.RowAlign(), A.ColCut(), A.RowCut(), 
      A.Height(), A.Width(), false, false );
    // NOTE: This will be *slightly* slower than necessary due to the result
    //       of the MPI operations being added rather than just copied
    Zeros( this->Matrix(), this->LocalHeight(), this->LocalWidth() );
    this->PartialRowSumScatterUpdate( T(1), A );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::PartialColSumScatterFrom
( const BlockDistMatrix<T,UPart,V>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::PartialColSumScatterFrom");
        AssertSameGrids( *this, A );
    )
    this->AlignAndResize
    ( A.BlockHeight(), A.BlockWidth(), 
      A.ColAlign(), A.RowAlign(), A.ColCut(), A.RowCut(), 
      A.Height(), A.Width(), false, false );
    // NOTE: This will be *slightly* slower than necessary due to the result
    //       of the MPI operations being added rather than just copied
    Zeros( this->Matrix(), this->LocalHeight(), this->LocalWidth() );
    this->PartialColSumScatterUpdate( T(1), A );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::RowSumScatterUpdate
( T alpha, const BlockDistMatrix<T,U,VGath>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::RowSumScatterUpdate");
        AssertSameGrids( *this, A );
        this->AssertNotLocked();
        this->AssertSameSize( A.Height(), A.Width() );
    )
    if( !this->Participating() )
        return;

    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::ColSumScatterUpdate
( T alpha, const BlockDistMatrix<T,UGath,V>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::ColSumScatterUpdate");
        AssertSameGrids( *this, A );
        this->AssertNotLocked();
        this->AssertSameSize( A.Height(), A.Width() );
    )
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::SumScatterUpdate
( T alpha, const BlockDistMatrix<T,UGath,VGath>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::SumScatterUpdate");
        AssertSameGrids( *this, A );
        this->AssertNotLocked();
        this->AssertSameSize( A.Height(), A.Width() );
    )
    if( !this->Participating() )
        return;
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::PartialRowSumScatterUpdate
( T alpha, const BlockDistMatrix<T,U,VPart>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::PartialRowSumScatterUpdate");
        AssertSameGrids( *this, A );
        this->AssertNotLocked();
        this->AssertSameSize( A.Height(), A.Width() );
    )
    if( !this->Participating() )
        return;

    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::PartialColSumScatterUpdate
( T alpha, const BlockDistMatrix<T,UPart,V>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::PartialColSumScatterUpdate");
        AssertSameGrids( *this, A );
        this->AssertNotLocked();
        this->AssertSameSize( A.Height(), A.Width() );
    )
    if( !this->Participating() )
        return;

    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposeColAllGather
( BlockDistMatrix<T,V,UGath>& A, bool conjugate ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposeColAllGather"))
    BlockDistMatrix<T,V,U> ATrans( this->Grid() );
    ATrans.AlignWith( *this );
    ATrans.Resize( this->Width(), this->Height() );
    Transpose( this->LockedMatrix(), ATrans.Matrix(), conjugate );
    copy::RowAllGather( ATrans, A );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposePartialColAllGather
( BlockDistMatrix<T,V,UPart>& A, bool conjugate ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposePartialColAllGather"))
    BlockDistMatrix<T,V,U> ATrans( this->Grid() );
    ATrans.AlignWith( *this );
    ATrans.Resize( this->Width(), this->Height() );
    Transpose( this->LockedMatrix(), ATrans.Matrix(), conjugate );
    copy::PartialRowAllGather( ATrans, A );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointColAllGather
( BlockDistMatrix<T,V,UGath>& A ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointRowAllGather"))
    this->TransposeColAllGather( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointPartialColAllGather
( BlockDistMatrix<T,V,UPart>& A ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointPartialColAllGather"))
    this->TransposePartialColAllGather( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposeColFilterFrom
( const BlockDistMatrix<T,V,UGath>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposeColFilterFrom"))
    BlockDistMatrix<T,V,U> AFilt( A.Grid() );
    if( this->ColConstrained() )
        AFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        AFilt.AlignColsWith( *this, false );
    copy::RowFilter( A, AFilt );
    if( !this->ColConstrained() )
        this->AlignColsWith( AFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( AFilt, false );
    this->Resize( A.Width(), A.Height() );
    Transpose( AFilt.LockedMatrix(), this->Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposeRowFilterFrom
( const BlockDistMatrix<T,VGath,U>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposeRowFilterFrom"))
    BlockDistMatrix<T,V,U> AFilt( A.Grid() );
    if( this->ColConstrained() )
        AFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        AFilt.AlignColsWith( *this, false );
    copy::ColFilter( A, AFilt );
    if( !this->ColConstrained() )
        this->AlignColsWith( AFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( AFilt, false );
    this->Resize( A.Width(), A.Height() );
    Transpose( AFilt.LockedMatrix(), this->Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposePartialColFilterFrom
( const BlockDistMatrix<T,V,UPart>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposePartialColFilterFrom"))
    BlockDistMatrix<T,V,U> AFilt( A.Grid() );
    if( this->ColConstrained() )
        AFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        AFilt.AlignColsWith( *this, false );
    copy::PartialRowFilter( A, AFilt );
    if( !this->ColConstrained() )
        this->AlignColsWith( AFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( AFilt, false );
    this->Resize( A.Width(), A.Height() );
    Transpose( AFilt.LockedMatrix(), this->Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposePartialRowFilterFrom
( const BlockDistMatrix<T,VPart,U>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposePartialRowFilterFrom"))
    BlockDistMatrix<T,V,U> AFilt( A.Grid() );
    if( this->ColConstrained() )
        AFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        AFilt.AlignColsWith( *this, false );
    copy::PartialColFilter( A, AFilt );
    if( !this->ColConstrained() )
        this->AlignColsWith( AFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( AFilt, false );
    this->Resize( A.Width(), A.Height() );
    Transpose( AFilt.LockedMatrix(), this->Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointColFilterFrom
( const BlockDistMatrix<T,V,UGath>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointColFilterFrom"))
    this->TransposeColFilterFrom( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointRowFilterFrom
( const BlockDistMatrix<T,VGath,U>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointRowFilterFrom"))
    this->TransposeRowFilterFrom( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointPartialColFilterFrom
( const BlockDistMatrix<T,V,UPart>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointPartialColFilterFrom"))
    this->TransposePartialColFilterFrom( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointPartialRowFilterFrom
( const BlockDistMatrix<T,VPart,U>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointPartialRowFilterFrom"))
    this->TransposePartialRowFilterFrom( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposeColSumScatterFrom
( const BlockDistMatrix<T,V,UGath>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposeColSumScatterFrom"))
    BlockDistMatrix<T,V,U> ASumFilt( A.Grid() );
    if( this->ColConstrained() )
        ASumFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        ASumFilt.AlignColsWith( *this, false );
    ASumFilt.RowSumScatterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( ASumFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( ASumFilt, false );
    this->Resize( A.Width(), A.Height() );
    Transpose( ASumFilt.LockedMatrix(), this->Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposePartialColSumScatterFrom
( const BlockDistMatrix<T,V,UPart>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposePartialColSumScatterFrom"))
    BlockDistMatrix<T,V,U> ASumFilt( A.Grid() );
    if( this->ColConstrained() )
        ASumFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        ASumFilt.AlignColsWith( *this, false );
    ASumFilt.PartialRowSumScatterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( ASumFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( ASumFilt, false );
    this->Resize( A.Width(), A.Height() );
    Transpose( ASumFilt.LockedMatrix(), this->Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointColSumScatterFrom
( const BlockDistMatrix<T,V,UGath>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointColSumScatterFrom"))
    this->TransposeColSumScatterFrom( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointPartialColSumScatterFrom
( const BlockDistMatrix<T,V,UPart>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointPartialColSumScatterFrom"))
    this->TransposePartialColSumScatterFrom( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposeColSumScatterUpdate
( T alpha, const BlockDistMatrix<T,V,UGath>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposeColSumScatterUpdate"))
    BlockDistMatrix<T,V,U> ASumFilt( A.Grid() );
    if( this->ColConstrained() )
        ASumFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        ASumFilt.AlignColsWith( *this, false );
    ASumFilt.RowSumScatterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( ASumFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( ASumFilt, false );
    // ALoc += alpha ASumFiltLoc'
    El::Matrix<T>& ALoc = this->Matrix();
    const El::Matrix<T>& BLoc = ASumFilt.LockedMatrix();
    const Int localHeight = ALoc.Height();
    const Int localWidth = ALoc.Width();
    if( conjugate )
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                ALoc.Update( iLoc, jLoc, alpha*Conj(BLoc.Get(jLoc,iLoc)) );
    }
    else
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                ALoc.Update( iLoc, jLoc, alpha*BLoc.Get(jLoc,iLoc) );
    }
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposePartialColSumScatterUpdate
( T alpha, const BlockDistMatrix<T,V,UPart>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposePartialColSumScatterUpdate"))
    BlockDistMatrix<T,V,U> ASumFilt( A.Grid() );
    if( this->ColConstrained() )
        ASumFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        ASumFilt.AlignColsWith( *this, false );
    ASumFilt.PartialRowSumScatterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( ASumFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( ASumFilt, false );
    // ALoc += alpha ASumFiltLoc'
    El::Matrix<T>& ALoc = this->Matrix();
    const El::Matrix<T>& BLoc = ASumFilt.LockedMatrix();
    const Int localHeight = ALoc.Height();
    const Int localWidth = ALoc.Width();
    if( conjugate )
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                ALoc.Update( iLoc, jLoc, alpha*Conj(BLoc.Get(jLoc,iLoc)) );
    }
    else
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                ALoc.Update( iLoc, jLoc, alpha*BLoc.Get(jLoc,iLoc) );
    }
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointColSumScatterUpdate
( T alpha, const BlockDistMatrix<T,V,UGath>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointColSumScatterUpdate"))
    this->TransposeColSumScatterUpdate( alpha, A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointPartialColSumScatterUpdate
( T alpha, const BlockDistMatrix<T,V,UPart>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointPartialColSumScatterUpdate"))
    this->TransposePartialColSumScatterUpdate( alpha, A, true );
}

// Basic queries
// =============

// Diagonal manipulation
// =====================
template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::GetDiagonal
( AbstractBlockDistMatrix<T>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::GetDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::GetRealPartOfDiagonal
( AbstractBlockDistMatrix<Base<T>>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::GetRealPartOfDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = RealPart(beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::GetImagPartOfDiagonal
( AbstractBlockDistMatrix<Base<T>>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::GetImagPartOfDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = ImagPart(beta); } );
}

template<typename T,Dist U,Dist V>
auto
GeneralBlockDistMatrix<T,U,V>::GetDiagonal( Int offset ) const
-> BlockDistMatrix<T,UDiag,VDiag>
{
    BlockDistMatrix<T,UDiag,VDiag> d( this->Grid() );
    GetDiagonal( d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
auto
GeneralBlockDistMatrix<T,U,V>::GetRealPartOfDiagonal( Int offset ) const
-> BlockDistMatrix<Base<T>,UDiag,VDiag>
{
    BlockDistMatrix<Base<T>,UDiag,VDiag> d( this->Grid() );
    GetRealPartOfDiagonal( d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
auto
GeneralBlockDistMatrix<T,U,V>::GetImagPartOfDiagonal( Int offset ) const
-> BlockDistMatrix<Base<T>,UDiag,VDiag>
{
    BlockDistMatrix<Base<T>,UDiag,VDiag> d( this->Grid() );
    GetImagPartOfDiagonal( d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::SetDiagonal
( const AbstractBlockDistMatrix<T>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::SetDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::SetRealPartOfDiagonal
( const AbstractBlockDistMatrix<Base<T>>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::SetRealPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { El::SetRealPart(alpha,beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::SetImagPartOfDiagonal
( const AbstractBlockDistMatrix<Base<T>>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::SetImagPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { El::SetImagPart(alpha,beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::UpdateDiagonal
( T gamma, const AbstractBlockDistMatrix<T>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::UpdateDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, [gamma]( T& alpha, T beta ) { alpha += gamma*beta; } );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::UpdateRealPartOfDiagonal
( Base<T> gamma, const AbstractBlockDistMatrix<Base<T>>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::UpdateRealPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      [gamma]( T& alpha, Base<T> beta ) 
      { El::UpdateRealPart(alpha,gamma*beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::UpdateImagPartOfDiagonal
( Base<T> gamma, const AbstractBlockDistMatrix<Base<T>>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::UpdateImagPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      [gamma]( T& alpha, Base<T> beta ) 
      { El::UpdateImagPart(alpha,gamma*beta); } );
}

// Private section
// ###############

// Diagonal helper functions
// =========================
template<typename T,Dist U,Dist V>
template<typename S,class Function>
void
GeneralBlockDistMatrix<T,U,V>::GetDiagonalHelper
( AbstractBlockDistMatrix<S>& d, Int offset, Function func ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::GetDiagonalHelper"))
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
template<typename S,class Function>
void
GeneralBlockDistMatrix<T,U,V>::SetDiagonalHelper
( const AbstractBlockDistMatrix<S>& d, Int offset, Function func ) 
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::SetDiagonalHelper");
        if( !this->DiagonalAlignedWith( d, offset ) )
            LogicError("Invalid diagonal alignment");
    )
    LogicError("This routine is not yet written");
}

// Instantiations for {Int,Real,Complex<Real>} for each Real in {float,double}
// ###########################################################################

#define DISTPROTO(T,U,V) template class GeneralBlockDistMatrix<T,U,V>
  
#define PROTO(T)\
  DISTPROTO(T,CIRC,CIRC);\
  DISTPROTO(T,MC,  MR  );\
  DISTPROTO(T,MC,  STAR);\
  DISTPROTO(T,MD,  STAR);\
  DISTPROTO(T,MR,  MC  );\
  DISTPROTO(T,MR,  STAR);\
  DISTPROTO(T,STAR,MC  );\
  DISTPROTO(T,STAR,MD  );\
  DISTPROTO(T,STAR,MR  );\
  DISTPROTO(T,STAR,STAR);\
  DISTPROTO(T,STAR,VC  );\
  DISTPROTO(T,STAR,VR  );\
  DISTPROTO(T,VC,  STAR);\
  DISTPROTO(T,VR,  STAR);

#include "El/macros/Instantiate.h"

} // namespace El
