/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "El/blas_like/level1/Copy/internal_impl.hpp"

namespace El {

#define DM DistMatrix<T,COLDIST,ROWDIST>
#define BDM DistMatrix<T,COLDIST,ROWDIST,BLOCK>
#define BCM BlockMatrix<T>
#define ADM AbstractDistMatrix<T>

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T>
BDM::DistMatrix( const El::Grid& g, int root )
: BCM(g,root)
{ 
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.SetViewType( OWNER );
    this->SetShifts(); 
}

template<typename T>
BDM::DistMatrix
( const El::Grid& g, Int blockHeight, Int blockWidth, int root )
: BCM(g,blockHeight,blockWidth,root)
{ 
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.SetViewType( OWNER );
    this->SetShifts(); 
}

template<typename T>
BDM::DistMatrix
( Int height, Int width, const El::Grid& g, int root )
: BCM(g,root)
{ 
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.SetViewType( OWNER );
    this->SetShifts(); this->Resize(height,width); 
}

template<typename T>
BDM::DistMatrix
( Int height, Int width, const El::Grid& g,
  Int blockHeight, Int blockWidth, int root )
: BCM(g,blockHeight,blockWidth,root)
{ 
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.SetViewType( OWNER );
    this->SetShifts(); 
    this->Resize(height,width); 
}

template<typename T>
BDM::DistMatrix( const BDM& A )
: BCM(A.Grid())
{
    DEBUG_ONLY(CSE cse("DistMatrix<T,U,V,BLOCK>::DistMatrix"))
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.SetViewType( OWNER );
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct block DistMatrix with itself");
}

template<typename T>
template<Dist U,Dist V>
BDM::DistMatrix( const DistMatrix<T,U,V,BLOCK>& A )
: BCM(A.Grid())
{
    DEBUG_ONLY(CSE cse("DistMatrix<T,U,V,BLOCK>::DistMatrix"))
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.SetViewType( OWNER );
    this->SetShifts();
    if( COLDIST != U || ROWDIST != V ||
        reinterpret_cast<const BDM*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct block DistMatrix with itself");
}

template<typename T>
BDM::DistMatrix( const AbstractDistMatrix<T>& A )
: BCM(A.Grid())
{
    DEBUG_ONLY(CSE cse("BDM(ADM)"))
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.SetViewType( OWNER );
    this->SetShifts();
    if( A.Wrap() == ELEMENT )
    {
        #define GUARD(CDIST,RDIST) A.ColDist() == CDIST && A.RowDist() == RDIST
        #define PAYLOAD(CDIST,RDIST) \
          auto& ACast = \
            static_cast<const DistMatrix<T,CDIST,RDIST>&>(A); \
          *this = ACast;
        #include "El/macros/GuardAndPayload.h"
    }
    else
    {
        #define GUARD(CDIST,RDIST) A.ColDist() == CDIST && A.RowDist() == RDIST
        #define PAYLOAD(CDIST,RDIST) \
          auto& ACast = \
            static_cast<const DistMatrix<T,CDIST,RDIST,BLOCK>&>(A); \
          *this = ACast;
        #include "El/macros/GuardAndPayload.h"
    }
}

template<typename T>
BDM::DistMatrix( const BlockMatrix<T>& A )
: BCM(A.Grid())
{
    DEBUG_ONLY(CSE cse("BDM(BCM)"))
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.SetViewType( OWNER );
    this->SetShifts();
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = \
        static_cast<const DistMatrix<T,CDIST,RDIST,BLOCK>&>(A); \
      if( COLDIST != CDIST || ROWDIST != RDIST || \
          reinterpret_cast<const BDM*>(&A) != this ) \
          *this = ACast; \
      else \
          LogicError("Tried to construct DistMatrix with itself");
    #include "El/macros/GuardAndPayload.h"
}

template<typename T>
template<Dist U,Dist V>
BDM::DistMatrix( const DistMatrix<T,U,V>& A )
: BCM(A.Grid())
{
    DEBUG_ONLY(CSE cse("DistMatrix<T,U,V,BLOCK>::DistMatrix"))
    if( COLDIST == CIRC && ROWDIST == CIRC )
        this->matrix_.SetViewType( OWNER );
    this->SetShifts();
    *this = A;
}

template<typename T>
BDM::DistMatrix( BDM&& A ) EL_NO_EXCEPT : BCM(std::move(A)) { } 

template<typename T> BDM::~DistMatrix() { }

template<typename T> 
BDM* BDM::Copy() const
{ return new DistMatrix<T,COLDIST,ROWDIST,BLOCK>(*this); }

template<typename T> 
BDM* BDM::Construct( const El::Grid& g, int root ) const
{ return new DistMatrix<T,COLDIST,ROWDIST,BLOCK>(g,root); }

template<typename T> 
DistMatrix<T,ROWDIST,COLDIST,BLOCK>* BDM::ConstructTranspose
( const El::Grid& g, int root ) const
{ return new DistMatrix<T,ROWDIST,COLDIST,BLOCK>(g,root); }

template<typename T> 
DistMatrix<T,DiagCol<COLDIST,ROWDIST>(),
             DiagRow<COLDIST,ROWDIST>(),BLOCK>* 
BDM::ConstructDiagonal
( const El::Grid& g, int root ) const
{ return new DistMatrix<T,DiagCol<COLDIST,ROWDIST>(),
                          DiagRow<COLDIST,ROWDIST>(),BLOCK>(g,root); }

// Operator overloading
// ====================

// Return a view
// -------------
template<typename T>
BDM BDM::operator()( Range<Int> I, Range<Int> J )
{
    DEBUG_ONLY(CSE cse("BDM( ind, ind )"))
    if( this->Locked() )
        return LockedView( *this, I, J );
    else
        return View( *this, I, J );
}

template<typename T>
const BDM BDM::operator()( Range<Int> I, Range<Int> J ) const
{
    DEBUG_ONLY(CSE cse("BDM( ind, ind ) const"))
    return LockedView( *this, I, J );
}

// Non-contiguous
// --------------
template<typename T>
BDM BDM::operator()( Range<Int> I, const vector<Int>& J ) const
{
    DEBUG_ONLY(CSE cse("DM( ind, vec ) const"))
    BDM ASub( this->Grid(), this->BlockHeight(), this->BlockWidth() );
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename T>
BDM BDM::operator()( const vector<Int>& I, Range<Int> J ) const
{
    DEBUG_ONLY(CSE cse("DM( vec, ind ) const"))
    BDM ASub( this->Grid(), this->BlockHeight(), this->BlockWidth() );
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename T>
BDM BDM::operator()( const vector<Int>& I, const vector<Int>& J ) const
{
    DEBUG_ONLY(CSE cse("DM( vec, vec ) const"))
    BDM ASub( this->Grid(), this->BlockHeight(), this->BlockWidth() );
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

// Copy
// ----

template<typename T>
BDM& BDM::operator=( const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("BDM = ADM"))
    if( A.Wrap() == ELEMENT )
    {
        #define GUARD(CDIST,RDIST) A.ColDist() == CDIST && A.RowDist() == RDIST
        #define PAYLOAD(CDIST,RDIST) \
          auto& ACast = \
            static_cast<const DistMatrix<T,CDIST,RDIST>&>(A); \
          *this = ACast;
        #include "El/macros/GuardAndPayload.h"
    }
    else
    {
        #define GUARD(CDIST,RDIST) A.ColDist() == CDIST && A.RowDist() == RDIST
        #define PAYLOAD(CDIST,RDIST) \
          auto& ACast = \
            static_cast<const DistMatrix<T,CDIST,RDIST,BLOCK>&>(A); \
          *this = ACast;
        #include "El/macros/GuardAndPayload.h"
    }
    return *this;
}

template<typename T>
template<Dist U,Dist V>
BDM& BDM::operator=( const DistMatrix<T,U,V>& A )
{
    DEBUG_ONLY(CSE cse("BDM = DM[U,V]"))
    // TODO: Use either AllGather or Gather if the distribution of this matrix
    //       is respectively either (STAR,STAR) or (CIRC,CIRC)
    // TODO: Specially handle cases where the block size is 1 x 1
    copy::GeneralPurpose( A, *this );
    return *this;
}

template<typename T>
BDM& BDM::operator=( BDM&& A )
{
    if( this->Viewing() || A.Viewing() )
        this->operator=( (const BDM&)A );
    else
        BCM::operator=( std::move(A) );
    return *this;
}

// Rescaling
// ---------
template<typename T>
const BDM& BDM::operator*=( T alpha )
{
    DEBUG_ONLY(CSE cse("BDM *= T"))
    Scale( alpha, *this );
    return *this;
}

// Addition/subtraction
// --------------------
template<typename T>
const BDM& BDM::operator+=( const BCM& A )
{
    DEBUG_ONLY(CSE cse("BDM += BCM&"))
    Axpy( T(1), A, *this );
    return *this;
}

template<typename T>
const BDM& BDM::operator+=( const ADM& A )
{
    DEBUG_ONLY(CSE cse("BDM += ADM&"))
    Axpy( T(1), A, *this );
    return *this;
}

template<typename T>
const BDM& BDM::operator-=( const BCM& A )
{
    DEBUG_ONLY(CSE cse("BDM += BCM&"))
    Axpy( T(-1), A, *this );
    return *this;
}

template<typename T>
const BDM& BDM::operator-=( const ADM& A )
{
    DEBUG_ONLY(CSE cse("BDM += ADM&"))
    Axpy( T(-1), A, *this );
    return *this;
}

// Distribution data
// =================

template<typename T>
El::DistData BDM::DistData() const { return El::DistData(*this); }

template<typename T>
Dist BDM::ColDist() const EL_NO_EXCEPT { return COLDIST; }
template<typename T>
Dist BDM::RowDist() const EL_NO_EXCEPT { return ROWDIST; }

template<typename T>
Dist BDM::PartialColDist() const EL_NO_EXCEPT { return Partial<COLDIST>(); }
template<typename T>
Dist BDM::PartialRowDist() const EL_NO_EXCEPT { return Partial<ROWDIST>(); }

template<typename T>
Dist BDM::PartialUnionColDist() const EL_NO_EXCEPT
{ return PartialUnionCol<COLDIST,ROWDIST>(); }
template<typename T>
Dist BDM::PartialUnionRowDist() const EL_NO_EXCEPT
{ return PartialUnionRow<COLDIST,ROWDIST>(); }

template<typename T>
Dist BDM::CollectedColDist() const EL_NO_EXCEPT { return Collect<COLDIST>(); }
template<typename T>
Dist BDM::CollectedRowDist() const EL_NO_EXCEPT { return Collect<ROWDIST>(); }

} // namespace El
