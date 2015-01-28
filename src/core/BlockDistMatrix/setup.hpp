/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

// This file should be included into each of the BlockDistMatrix specializations
// as a workaround for the fact that C++11 constructor inheritance is not 
// yet widely supported.

#include "El/blas_like/level1/copy_internal.hpp"

namespace El {

#define DM DistMatrix<T,COLDIST,ROWDIST>
#define BDM BlockDistMatrix<T,COLDIST,ROWDIST>
#define ABDM AbstractBlockDistMatrix<T>

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T>
BDM::BlockDistMatrix( const El::Grid& g, int root )
: ABDM(g,root)
{ 
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts(); 
}

template<typename T>
BDM::BlockDistMatrix
( const El::Grid& g, Int blockHeight, Int blockWidth, int root )
: ABDM(g,blockHeight,blockWidth,root)
{ 
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts(); 
}

template<typename T>
BDM::BlockDistMatrix
( Int height, Int width, const El::Grid& g, int root )
: ABDM(g,root)
{ 
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts(); this->Resize(height,width); 
}

template<typename T>
BDM::BlockDistMatrix
( Int height, Int width, const El::Grid& g,
  Int blockHeight, Int blockWidth, int root )
: ABDM(g,blockHeight,blockWidth,root)
{ 
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts(); 
    this->Resize(height,width); 
}

template<typename T>
BDM::BlockDistMatrix( const BDM& A )
: ABDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("BlockDistMatrix::BlockDistMatrix"))
    if( COLDIST == CIRC && ROWDIST == CIRC )
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
: ABDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("BlockDistMatrix::BlockDistMatrix"))
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts();
    if( COLDIST != U || ROWDIST != V ||
        reinterpret_cast<const BDM*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct BlockDistMatrix with itself");
}

template<typename T>
BDM::BlockDistMatrix( const AbstractBlockDistMatrix<T>& A )
: ABDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("BDM(ABDM)"))
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts();
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = dynamic_cast<const DistMatrix<T,CDIST,RDIST>&>(A); \
      if( COLDIST != CDIST || ROWDIST != RDIST || \
          reinterpret_cast<const BDM*>(&A) != this ) \
          *this = ACast; \
      else \
          LogicError("Tried to construct DistMatrix with itself");
    #include "El/macros/GuardAndPayload.h"
}

template<typename T>
template<Dist U,Dist V>
BDM::BlockDistMatrix( const DistMatrix<T,U,V>& A )
: ABDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("BlockDistMatrix::BlockDistMatrix"))
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts();
    *this = A;
}

template<typename T>
BDM::BlockDistMatrix( BDM&& A ) EL_NOEXCEPT : ABDM(std::move(A)) { } 

template<typename T> BDM::~BlockDistMatrix() { }

template<typename T> 
BlockDistMatrix<T,COLDIST,ROWDIST>* BDM::Construct
( const El::Grid& g, int root ) const
{ return new BlockDistMatrix<T,COLDIST,ROWDIST>(g,root); }

template<typename T> 
BlockDistMatrix<T,ROWDIST,COLDIST>* BDM::ConstructTranspose
( const El::Grid& g, int root ) const
{ return new BlockDistMatrix<T,ROWDIST,COLDIST>(g,root); }

template<typename T> 
BlockDistMatrix<T,DiagCol<COLDIST,ROWDIST>(),
                  DiagRow<COLDIST,ROWDIST>()>* 
BDM::ConstructDiagonal
( const El::Grid& g, int root ) const
{ return new BlockDistMatrix<T,DiagCol<COLDIST,ROWDIST>(),
                               DiagRow<COLDIST,ROWDIST>()>(g,root); }

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
BDM& BDM::operator=( BDM&& A )
{
    if( this->Viewing() || A.Viewing() )
        this->operator=( (const BDM&)A );
    else
        ABDM::operator=( std::move(A) );
    return *this;
}

// Distribution data
// =================

template<typename T>
El::BlockDistData BDM::DistData() const { return El::BlockDistData(*this); }

template<typename T>
Dist BDM::ColDist() const { return COLDIST; }
template<typename T>
Dist BDM::RowDist() const { return ROWDIST; }

template<typename T>
Dist BDM::PartialColDist() const { return Partial<COLDIST>(); }
template<typename T>
Dist BDM::PartialRowDist() const { return Partial<ROWDIST>(); }

template<typename T>
Dist BDM::PartialUnionColDist() const
{ return PartialUnionCol<COLDIST,ROWDIST>(); }
template<typename T>
Dist BDM::PartialUnionRowDist() const
{ return PartialUnionRow<COLDIST,ROWDIST>(); }

template<typename T>
Dist BDM::CollectedColDist() const { return Collect<COLDIST>(); }
template<typename T>
Dist BDM::CollectedRowDist() const { return Collect<ROWDIST>(); }

} // namespace El
