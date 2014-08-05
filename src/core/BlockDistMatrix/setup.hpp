/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

// This file should be included into each of the BlockDistMatrix specializations
// as a workaround for the fact that C++11 constructor inheritance is not 
// yet widely supported.

namespace El {

#define DM DistMatrix<T,ColDist,RowDist>
#define BDM BlockDistMatrix<T,ColDist,RowDist>
#define GBDM GeneralBlockDistMatrix<T,ColDist,RowDist>

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T>
BDM::BlockDistMatrix( const El::Grid& g, Int root )
: GBDM(g,root)
{ 
    if( ColDist == CIRC && RowDist == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts(); 
}

template<typename T>
BDM::BlockDistMatrix
( const El::Grid& g, Int blockHeight, Int blockWidth, Int root )
: GBDM(g,blockHeight,blockWidth,root)
{ 
    if( ColDist == CIRC && RowDist == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts(); 
}

template<typename T>
BDM::BlockDistMatrix
( Int height, Int width, const El::Grid& g, Int root )
: GBDM(g,root)
{ 
    if( ColDist == CIRC && RowDist == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts(); this->Resize(height,width); 
}

template<typename T>
BDM::BlockDistMatrix
( Int height, Int width, const El::Grid& g,
  Int blockHeight, Int blockWidth, Int root )
: GBDM(g,blockHeight,blockWidth,root)
{ 
    if( ColDist == CIRC && RowDist == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts(); 
    this->Resize(height,width); 
}

template<typename T>
BDM::BlockDistMatrix( const BDM& A )
: GBDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("BlockDistMatrix::BlockDistMatrix"))
    if( ColDist == CIRC && RowDist == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct BlockDistMatrix with itself");
}

template<typename T>
template<Dist U,Dist V>
BDM::BlockDistMatrix( const BlockDistMatrix<T,U,V>& A )
: GBDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("BlockDistMatrix::BlockDistMatrix"))
    if( ColDist == CIRC && RowDist == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts();
    if( ColDist != U || RowDist != V ||
        reinterpret_cast<const BDM*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct BlockDistMatrix with itself");
}

template<typename T>
BDM::BlockDistMatrix( const AbstractBlockDistMatrix<T>& A )
: GBDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("BDM(ABDM)"))
    if( ColDist == CIRC && RowDist == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts();
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = dynamic_cast<const DistMatrix<T,CDIST,RDIST>&>(A); \
      if( ColDist != CDIST || RowDist != RDIST || \
          reinterpret_cast<const BDM*>(&A) != this ) \
          *this = ACast; \
      else \
          LogicError("Tried to construct DistMatrix with itself");
    #include "El/macros/GuardAndPayload.h"
}

template<typename T>
template<Dist U,Dist V>
BDM::BlockDistMatrix( const DistMatrix<T,U,V>& A )
: GBDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("BlockDistMatrix::BlockDistMatrix"))
    if( ColDist == CIRC && RowDist == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts();
    *this = A;
}

template<typename T>
BDM::BlockDistMatrix( BDM&& A ) EL_NOEXCEPT : GBDM(std::move(A)) { } 

template<typename T> BDM::~BlockDistMatrix() { }

template<typename T> 
BlockDistMatrix<T,ColDist,RowDist>* BDM::Construct
( const El::Grid& g, Int root ) const
{ return new BlockDistMatrix<T,ColDist,RowDist>(g,root); }

template<typename T> 
BlockDistMatrix<T,RowDist,ColDist>* BDM::ConstructTranspose
( const El::Grid& g, Int root ) const
{ return new BlockDistMatrix<T,RowDist,ColDist>(g,root); }

template<typename T> 
BlockDistMatrix<T,DiagColDist<ColDist,RowDist>(),
                  DiagRowDist<ColDist,RowDist>()>* BDM::ConstructDiagonal
( const El::Grid& g, Int root ) const
{ return new BlockDistMatrix<T,DiagColDist<ColDist,RowDist>(),
                               DiagRowDist<ColDist,RowDist>()>(g,root); }

template<typename T>
template<Dist U,Dist V>
BDM& BDM::operator=( const DistMatrix<T,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("BDM = DM[U,V]"))
    BlockDistMatrix<T,U,V> ABlock(A.Grid());
    LockedView( ABlock, A );
    *this = ABlock;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const AbstractBlockDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("BDM = ABDM"))
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = dynamic_cast<const BlockDistMatrix<T,CDIST,RDIST>&>(A); \
      *this = ACast;
    #include "El/macros/GuardAndPayload.h"
    return *this;
}

template<typename T>
BDM& BDM::operator=( BDM&& A )
{
    if( this->Viewing() && !A.Viewing() )
        this->operator=( (const BDM&)A );
    else
        GBDM::operator=( std::move(A) );
    return *this;
}

template<typename T>
El::BlockDistData BDM::DistData() const { return El::BlockDistData(*this); }

} // namespace El
