/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#define COLDIST MR
#define ROWDIST STAR

#include "./setup.hpp"

namespace El {

// Public section
// ##############

// Assignment and reconfiguration
// ==============================

// Return a view
// -------------
template<typename T>
DM DM::operator()( Range<Int> indVert, Range<Int> indHorz )
{
    DEBUG_ONLY(CallStackEntry cse("DM[MR,STAR]( ind, ind )"))
    if( this->Locked() )
        return LockedView( *this, indVert, indHorz );
    else
        return View( *this, indVert, indHorz );
}

template<typename T>
const DM DM::operator()( Range<Int> indVert, Range<Int> indHorz ) const
{
    DEBUG_ONLY(CallStackEntry cse("DM[MR,STAR]( ind, ind )"))
    return LockedView( *this, indVert, indHorz );
}

// Make a copy
// -----------
template<typename T>
DM& DM::operator=( const DM& A )
{
    DEBUG_ONLY(CallStackEntry cse("DM[MR,STAR] = DM[MR,STAR]"))
    copy::Translate( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MC,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MR,STAR] = [MC,MR]"))
    auto A_VC_STAR = MakeUnique<DistMatrix<T,VC,STAR>>( A );
    auto A_VR_STAR = MakeUnique<DistMatrix<T,VR,STAR>>( this->Grid() );
    A_VR_STAR->AlignColsWith(*this);
    *A_VR_STAR = *A_VC_STAR;
    A_VC_STAR.reset(); 

    *this = *A_VR_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
    DEBUG_ONLY(
      CallStackEntry cse("[MR,STAR] = [MC,STAR]");
      AssertSameGrids( *this, A );
    )
    const Grid& grid = A.Grid();
    if( grid.Height() == grid.Width() ) 
    {
        const int gridDim = grid.Height();
        const int transposeRank =
            A.RowOwner(this->ColShift()) + gridDim*this->RowOwner(A.ColShift());
        copy::Exchange( A, *this, transposeRank, transposeRank, grid.VCComm() );
    }
    else
    {
        auto A_VC_STAR = MakeUnique<DistMatrix<T,VC,STAR>>( A );
        auto A_VR_STAR = MakeUnique<DistMatrix<T,VR,STAR>>( grid );
        A_VR_STAR->AlignColsWith(*this);
        *A_VR_STAR = *A_VC_STAR;
        A_VC_STAR.reset(); 
        *this = *A_VR_STAR;
    }
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MR,STAR] = [STAR,MR]"))
    auto A_MC_MR = MakeUnique<DistMatrix<T,MC,MR>>( A );
    auto A_VC_STAR = MakeUnique<DistMatrix<T,VC,STAR>>( *A_MC_MR );
    A_MC_MR.reset(); 

    auto A_VR_STAR = MakeUnique<DistMatrix<T,VR,STAR>>( this->Grid() );
    A_VR_STAR->AlignColsWith(*this);
    *A_VR_STAR = *A_VC_STAR;
    A_VC_STAR.reset();

    *this = *A_VR_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MR,STAR] = [MD,STAR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MR,STAR] = [STAR,MD]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MR,STAR] = [MR,MC]"))
    copy::RowAllGather( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MR,STAR] = [STAR,MC]"))
    DistMatrix<T,MR,MC> A_MR_MC( A );
    *this = A_MR_MC;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MR,STAR] = [VC,STAR]"))
    DistMatrix<T,VR,STAR> A_VR_STAR(this->Grid());
    A_VR_STAR.AlignColsWith(*this);
    A_VR_STAR = A;
    *this = A_VR_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MR,STAR] = [STAR,VC]"))
    DistMatrix<T,MR,MC> A_MR_MC( A );
    *this = A_MR_MC;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MR,STAR] = [VR,STAR]"))
    copy::PartialColAllGather( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MR,STAR] = [STAR,VR]"))
    auto A_STAR_VC = MakeUnique<DistMatrix<T,STAR,VC>>( A );
    auto A_MR_MC = MakeUnique<DistMatrix<T,MR,MC>>( this->Grid() );
    A_MR_MC->AlignColsWith(*this);
    *A_MR_MC = *A_STAR_VC;
    A_STAR_VC.reset(); 

    *this = *A_MR_MC;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MR,STAR] = [STAR,STAR]"))
    copy::ColFilter( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MR,STAR] = [CIRC,CIRC]"))
    DistMatrix<T,MR,MC> A_MR_MC( this->Grid() );
    A_MR_MC.AlignWith( *this );
    A_MR_MC = A;
    *this = A_MR_MC;
    return *this;
}

template<typename T>
DM& DM::operator=( const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("DM = ADM"))
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = dynamic_cast<const DistMatrix<T,CDIST,RDIST>&>(A); \
      *this = ACast;
    #include "El/macros/GuardAndPayload.h"
    return *this;
}

// Basic queries
// =============

template<typename T>
mpi::Comm DM::DistComm() const { return this->grid_->MRComm(); }
template<typename T>
mpi::Comm DM::CrossComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::RedundantComm() const { return this->grid_->MCComm(); }
template<typename T>
mpi::Comm DM::ColComm() const { return this->grid_->MRComm(); }
template<typename T>
mpi::Comm DM::RowComm() const { return mpi::COMM_SELF; }

template<typename T>
int DM::ColStride() const { return this->grid_->MRSize(); }
template<typename T>
int DM::RowStride() const { return 1; }
template<typename T>
int DM::DistSize() const { return this->grid_->MRSize(); }
template<typename T>
int DM::CrossSize() const { return 1; }
template<typename T>
int DM::RedundantSize() const { return this->grid_->MCSize(); }

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define SELF(T,U,V) \
  template DistMatrix<T,COLDIST,ROWDIST>::DistMatrix \
  ( const DistMatrix<T,U,V>& A );
#define OTHER(T,U,V) \
  template DistMatrix<T,COLDIST,ROWDIST>::DistMatrix \
  ( const BlockDistMatrix<T,U,V>& A ); \
  template DistMatrix<T,COLDIST,ROWDIST>& \
           DistMatrix<T,COLDIST,ROWDIST>::operator= \
           ( const BlockDistMatrix<T,U,V>& A )
#define BOTH(T,U,V) \
  SELF(T,U,V); \
  OTHER(T,U,V)
#define PROTO(T) \
  template class DistMatrix<T,COLDIST,ROWDIST>; \
  BOTH( T,CIRC,CIRC); \
  BOTH( T,MC,  MR  ); \
  BOTH( T,MC,  STAR); \
  BOTH( T,MD,  STAR); \
  BOTH( T,MR,  MC  ); \
  OTHER(T,MR,  STAR); \
  BOTH( T,STAR,MC  ); \
  BOTH( T,STAR,MD  ); \
  BOTH( T,STAR,MR  ); \
  BOTH( T,STAR,STAR); \
  BOTH( T,STAR,VC  ); \
  BOTH( T,STAR,VR  ); \
  BOTH( T,VC,  STAR); \
  BOTH( T,VR,  STAR);

#include "El/macros/Instantiate.h"

} // namespace El
