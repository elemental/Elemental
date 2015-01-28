/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

// This file should be included into each of the DistMatrix specializations
// as a workaround for the fact that C++11 constructor inheritance is not 
// yet widely supported.

#include "El/blas_like/level1/copy_internal.hpp"

namespace El {

#define DM DistMatrix<T,COLDIST,ROWDIST>
#define ADM AbstractDistMatrix<T>

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T>
DM::DistMatrix( const El::Grid& grid, int root )
: ADM(grid,root)
{ 
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts(); 
}

template<typename T>
DM::DistMatrix( Int height, Int width, const El::Grid& grid, int root )
: ADM(grid,root)
{ 
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts(); 
    this->Resize(height,width); 
}

template<typename T>
DM::DistMatrix( const DM& A )
: ADM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("DistMatrix::DistMatrix"))
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct DistMatrix with itself");
}

template<typename T>
template<Dist U,Dist V>
DM::DistMatrix( const DistMatrix<T,U,V>& A )
: ADM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("DistMatrix::DistMatrix"))
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts();
    if( COLDIST != U || ROWDIST != V ||
        reinterpret_cast<const DM*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct DistMatrix with itself");
}

template<typename T>
DM::DistMatrix( const AbstractDistMatrix<T>& A )
: ADM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("DM(ADM)"))
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts();
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = dynamic_cast<const DistMatrix<T,CDIST,RDIST>&>(A); \
      if( COLDIST != CDIST || ROWDIST != RDIST || \
          reinterpret_cast<const DM*>(&A) != this ) \
          *this = ACast; \
      else \
          LogicError("Tried to construct DistMatrix with itself");
    #include "El/macros/GuardAndPayload.h"
}

template<typename T>
template<Dist U,Dist V>
DM::DistMatrix( const BlockDistMatrix<T,U,V>& A )
: ADM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("DistMatrix::DistMatrix"))
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.viewType_ = OWNER;
    this->SetShifts();
    *this = A;
}

template<typename T>
DM::DistMatrix( DM&& A ) EL_NOEXCEPT : ADM(std::move(A)) { }

template<typename T> DM::~DistMatrix() { }

template<typename T> 
DistMatrix<T,COLDIST,ROWDIST>* DM::Construct
( const El::Grid& g, int root ) const
{ return new DistMatrix<T,COLDIST,ROWDIST>(g,root); }

template<typename T> 
DistMatrix<T,ROWDIST,COLDIST>* DM::ConstructTranspose
( const El::Grid& g, int root ) const
{ return new DistMatrix<T,ROWDIST,COLDIST>(g,root); }

template<typename T>
DistMatrix<T,DiagCol<COLDIST,ROWDIST>(),
             DiagRow<COLDIST,ROWDIST>()>* 
DM::ConstructDiagonal
( const El::Grid& g, int root ) const
{ return new DistMatrix<T,DiagCol<COLDIST,ROWDIST>(),
                          DiagRow<COLDIST,ROWDIST>()>(g,root); }

template<typename T>
template<Dist U,Dist V>
DM& DM::operator=( const BlockDistMatrix<T,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("DM = BDM[U,V]"))
    BlockDistMatrix<T,COLDIST,ROWDIST> AElem(A.Grid(),1,1);
    AElem = A;
    DistMatrix<T,COLDIST,ROWDIST> AElemView(A.Grid());
    LockedView( AElemView, AElem ); 
    *this = AElemView;
    return *this;
}

template<typename T>
DM& DM::operator=( DM&& A )
{
    if( this->Viewing() || A.Viewing() )
        this->operator=( (const DM&)A );
    else
        ADM::operator=( std::move(A) );
    return *this;
}

// Distribution data
// =================

template<typename T>
El::DistData DM::DistData() const { return El::DistData(*this); }

template<typename T>
Dist DM::ColDist() const { return COLDIST; }
template<typename T>
Dist DM::RowDist() const { return ROWDIST; }

template<typename T>
Dist DM::PartialColDist() const { return Partial<COLDIST>(); }
template<typename T>
Dist DM::PartialRowDist() const { return Partial<ROWDIST>(); }

template<typename T>
Dist DM::PartialUnionColDist() const 
{ return PartialUnionCol<COLDIST,ROWDIST>(); }
template<typename T>
Dist DM::PartialUnionRowDist() const
{ return PartialUnionRow<COLDIST,ROWDIST>(); }

template<typename T>
Dist DM::CollectedColDist() const { return Collect<COLDIST>(); }
template<typename T>
Dist DM::CollectedRowDist() const { return Collect<ROWDIST>(); }

} // namespace El
