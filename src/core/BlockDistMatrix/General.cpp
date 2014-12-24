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
