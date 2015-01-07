/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#define COLDIST VC
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
    DEBUG_ONLY(CallStackEntry cse("DM[VC,STAR]( ind, ind )"))
    if( this->Locked() )
        return LockedView( *this, indVert, indHorz );
    else
        return View( *this, indVert, indHorz );
}

template<typename T>
const DM DM::operator()( Range<Int> indVert, Range<Int> indHorz ) const
{
    DEBUG_ONLY(CallStackEntry cse("DM[VC,STAR]( ind, ind )"))
    return LockedView( *this, indVert, indHorz );
}

// Make a copy
// -----------
template<typename T>
DM& DM::operator=( const DM& A )
{
    DEBUG_ONLY(CallStackEntry cse("DM[VC,STAR] = DM[VC,STAR]"))
    copy::Translate( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MC,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [MC,MR]"))
    copy::ColAllToAllDemote( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [MC,STAR]"))
    copy::PartialColFilter( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [STAR,MR]"))
    DistMatrix<T,MC,MR> A_MC_MR( A );
    *this = A_MC_MR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [MD,STAR]"))
    // TODO: Optimize this later if important
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [STAR,MD]"))
    // TODO: Optimize this later if important
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [MR,MC]"))
    DistMatrix<T,VR,STAR> A_VR_STAR( A );
    *this = A_VR_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [MR,STAR]"))
    DistMatrix<T,VR,STAR> A_VR_STAR( A );
    *this = A_VR_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [STAR,MC]"))
    auto A_MR_MC = MakeUnique<DistMatrix<T,MR,MC>>( A );
    auto A_VR_STAR = MakeUnique<DistMatrix<T,VR,STAR>>( *A_MR_MC );
    A_MR_MC.reset(); 
    *this = *A_VR_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [STAR,VC]"))
    auto A_MR_MC = MakeUnique<DistMatrix<T,MR,MC>>( A );
    auto A_VR_STAR = MakeUnique<DistMatrix<T,VR,STAR>>( *A_MR_MC );
    A_MR_MC.reset(); 
    *this = *A_VR_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [VR,STAR]"))
    copy::ColwiseVectorExchange<T,MR,MC>( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [STAR,VR]"))
    DistMatrix<T,MC,MR> A_MC_MR( A );
    *this = A_MC_MR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [STAR,STAR]"))
    copy::ColFilter( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [CIRC,CIRC]"))
    copy::Scatter( A, *this );
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
mpi::Comm DM::DistComm() const { return this->grid_->VCComm(); }
template<typename T>
mpi::Comm DM::CrossComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::RedundantComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::ColComm() const { return this->grid_->VCComm(); }
template<typename T>
mpi::Comm DM::RowComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::PartialColComm() const { return this->grid_->MCComm(); }
template<typename T>
mpi::Comm DM::PartialUnionColComm() const { return this->grid_->MRComm(); }

template<typename T>
int DM::ColStride() const { return this->grid_->VCSize(); }
template<typename T>
int DM::RowStride() const { return 1; }
template<typename T>
int DM::PartialColStride() const { return this->grid_->MCSize(); }
template<typename T>
int DM::PartialUnionColStride() const { return this->grid_->MRSize(); }
template<typename T>
int DM::DistSize() const { return this->grid_->VCSize(); }
template<typename T>
int DM::CrossSize() const { return 1; }
template<typename T>
int DM::RedundantSize() const { return 1; }

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
  BOTH( T,MR,  STAR); \
  BOTH( T,STAR,MC  ); \
  BOTH( T,STAR,MD  ); \
  BOTH( T,STAR,MR  ); \
  BOTH( T,STAR,STAR); \
  BOTH( T,STAR,VC  ); \
  BOTH( T,STAR,VR  ); \
  OTHER(T,VC,  STAR); \
  BOTH( T,VR,  STAR);

#include "El/macros/Instantiate.h"

} // namespace El
