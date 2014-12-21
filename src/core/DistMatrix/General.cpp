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
GeneralDistMatrix<T,U,V>::GeneralDistMatrix( const El::Grid& grid, Int root )
: AbstractDistMatrix<T>(grid,root)
{ }

template<typename T,Dist U,Dist V>
GeneralDistMatrix<T,U,V>::GeneralDistMatrix( GeneralDistMatrix<T,U,V>&& A ) 
EL_NOEXCEPT
: AbstractDistMatrix<T>(std::move(A))
{ }

// Assignment and reconfiguration
// ==============================

template<typename T,Dist U,Dist V>
GeneralDistMatrix<T,U,V>& 
GeneralDistMatrix<T,U,V>::operator=( GeneralDistMatrix<T,U,V>&& A )
{
    AbstractDistMatrix<T>::operator=( std::move(A) );
    return *this;
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::AlignColsWith
( const El::DistData& data, bool constrain, bool allowMismatch )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::AlignColsWith")) 
    this->SetGrid( *data.grid );
    this->SetRoot( data.root );
    if( data.colDist == U || data.colDist == UPart )
        this->AlignCols( data.colAlign, constrain );
    else if( data.rowDist == U || data.rowDist == UPart )
        this->AlignCols( data.rowAlign, constrain );
    else if( data.colDist == UScat )
        this->AlignCols( data.colAlign % this->ColStride(), constrain );
    else if( data.rowDist == UScat )
        this->AlignCols( data.rowAlign % this->ColStride(), constrain );
    else if( U != UGath && data.colDist != UGath && data.rowDist != UGath &&
            !allowMismatch ) 
        LogicError("Nonsensical alignment");
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::AlignRowsWith
( const El::DistData& data, bool constrain, bool allowMismatch )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::AlignRowsWith")) 
    this->SetGrid( *data.grid );
    this->SetRoot( data.root );
    if( data.colDist == V || data.colDist == VPart )
        this->AlignRows( data.colAlign, constrain );
    else if( data.rowDist == V || data.rowDist == VPart )
        this->AlignRows( data.rowAlign, constrain );
    else if( data.colDist == VScat )
        this->AlignRows( data.colAlign % this->RowStride(), constrain );
    else if( data.rowDist == VScat )
        this->AlignRows( data.rowAlign % this->RowStride(), constrain );
    else if( V != VGath && data.colDist != VGath && data.rowDist != VGath &&
             !allowMismatch ) 
        LogicError("Nonsensical alignment");
}

// Basic queries
// =============
// Distribution information
// ------------------------
template<typename T,Dist U,Dist V>
Dist GeneralDistMatrix<T,U,V>::ColDist() const { return U; }
template<typename T,Dist U,Dist V>
Dist GeneralDistMatrix<T,U,V>::RowDist() const { return V; }
template<typename T,Dist U,Dist V>
Dist GeneralDistMatrix<T,U,V>::PartialColDist() const 
{ return Partial<U>(); }
template<typename T,Dist U,Dist V>
Dist GeneralDistMatrix<T,U,V>::PartialRowDist() const 
{ return Partial<V>(); }
template<typename T,Dist U,Dist V>
Dist GeneralDistMatrix<T,U,V>::PartialUnionColDist() const 
{ return PartialUnionCol<U,V>(); }
template<typename T,Dist U,Dist V>
Dist GeneralDistMatrix<T,U,V>::PartialUnionRowDist() const 
{ return PartialUnionRow<U,V>(); }

// Diagonal manipulation
// =====================
template<typename T,Dist U,Dist V>
void GeneralDistMatrix<T,U,V>::GetDiagonal
( AbstractDistMatrix<T>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::GetDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::GetRealPartOfDiagonal
( AbstractDistMatrix<Base<T>>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::GetRealPartOfDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = RealPart(beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::GetImagPartOfDiagonal
( AbstractDistMatrix<Base<T>>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::GetImagPartOfDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = ImagPart(beta); } );
}

template<typename T,Dist U,Dist V>
auto
GeneralDistMatrix<T,U,V>::GetDiagonal( Int offset ) const
-> DistMatrix<T,UDiag,VDiag>
{
    DistMatrix<T,UDiag,VDiag> d( this->Grid() );
    GetDiagonal( d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
auto
GeneralDistMatrix<T,U,V>::GetRealPartOfDiagonal( Int offset ) const
-> DistMatrix<Base<T>,UDiag,VDiag>
{
    DistMatrix<Base<T>,UDiag,VDiag> d( this->Grid() );
    GetRealPartOfDiagonal( d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
auto
GeneralDistMatrix<T,U,V>::GetImagPartOfDiagonal( Int offset ) const
-> DistMatrix<Base<T>,UDiag,VDiag>
{
    DistMatrix<Base<T>,UDiag,VDiag> d( this->Grid() );
    GetImagPartOfDiagonal( d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::SetDiagonal
( const AbstractDistMatrix<T>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::SetDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::SetRealPartOfDiagonal
( const AbstractDistMatrix<Base<T>>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::SetRealPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { El::SetRealPart(alpha,beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::SetImagPartOfDiagonal
( const AbstractDistMatrix<Base<T>>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::SetImagPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { El::SetImagPart(alpha,beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::UpdateDiagonal
( T gamma, const AbstractDistMatrix<T>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::UpdateDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, [gamma]( T& alpha, T beta ) { alpha += gamma*beta; } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::UpdateRealPartOfDiagonal
( Base<T> gamma, const AbstractDistMatrix<Base<T>>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::UpdateRealPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      [gamma]( T& alpha, Base<T> beta ) 
      { El::UpdateRealPart(alpha,gamma*beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::UpdateImagPartOfDiagonal
( Base<T> gamma, const AbstractDistMatrix<Base<T>>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::UpdateImagPartOfDiagonal"))
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
GeneralDistMatrix<T,U,V>::GetDiagonalHelper
( AbstractDistMatrix<S>& dPre, Int offset, Function func ) const
{
    DEBUG_ONLY(
      CallStackEntry cse("GDM::GetDiagonalHelper");
      AssertSameGrids( *this, dPre );
    )
    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = this->DiagonalAlign(offset);
    ctrl.rootConstrain = true;
    ctrl.root = this->DiagonalRoot(offset);
    auto dPtr = WriteProxy<S,UDiag,VDiag>(&dPre,ctrl);
    auto& d = *dPtr;

    d.Resize( this->DiagonalLength(offset), 1 );
    if( !d.Participating() )
        return;

    const Int diagShift = d.ColShift();
    const Int iStart = diagShift + Max(-offset,0);
    const Int jStart = diagShift + Max( offset,0);

    const Int colStride = this->ColStride();
    const Int rowStride = this->RowStride();
    const Int iLocStart = (iStart-this->ColShift()) / colStride;
    const Int jLocStart = (jStart-this->RowShift()) / rowStride;
    const Int iLocStride = d.ColStride() / colStride;
    const Int jLocStride = d.ColStride() / rowStride;

    const Int localDiagLength = d.LocalHeight();
    S* dBuf = d.Buffer();
    const T* buffer = this->LockedBuffer();
    const Int ldim = this->LDim();
    EL_PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*iLocStride;
        const Int jLoc = jLocStart + k*jLocStride;
        func( dBuf[k], buffer[iLoc+jLoc*ldim] );
    }
}

template<typename T,Dist U,Dist V>
template<typename S,class Function>
void
GeneralDistMatrix<T,U,V>::SetDiagonalHelper
( const AbstractDistMatrix<S>& dPre, Int offset, Function func ) 
{
    DEBUG_ONLY(
      CallStackEntry cse("GDM::SetDiagonalHelper");
      AssertSameGrids( *this, dPre );
    )
    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = this->DiagonalAlign(offset);
    ctrl.rootConstrain = true;
    ctrl.root = this->DiagonalRoot(offset);
    auto dPtr = ReadProxy<S,UDiag,VDiag>(&dPre,ctrl); 
    const auto& d = *dPtr;

    if( !d.Participating() )
        return;

    const Int diagShift = d.ColShift();
    const Int iStart = diagShift + Max(-offset,0);
    const Int jStart = diagShift + Max( offset,0);

    const Int colStride = this->ColStride();
    const Int rowStride = this->RowStride();
    const Int iLocStart = (iStart-this->ColShift()) / colStride;
    const Int jLocStart = (jStart-this->RowShift()) / rowStride;
    const Int iLocStride = d.ColStride() / colStride;
    const Int jLocStride = d.ColStride() / rowStride;

    const Int localDiagLength = d.LocalHeight();
    const S* dBuf = d.LockedBuffer();
    T* buffer = this->Buffer();
    const Int ldim = this->LDim();
    EL_PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*iLocStride;
        const Int jLoc = jLocStart + k*jLocStride;
        func( buffer[iLoc+jLoc*ldim], dBuf[k] );
    }
}

// Instantiations for {Int,Real,Complex<Real>} for each Real in {float,double}
// ###########################################################################

#define DISTPROTO(T,U,V) template class GeneralDistMatrix<T,U,V>
  
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
