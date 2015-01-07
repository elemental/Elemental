/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#define COLDIST STAR
#define ROWDIST VR

#include "./setup.hpp"

namespace El {

// Public section
// ##############

// Assignment and reconfiguration
// ==============================

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,MC,MR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,VR] = [MC,MR]"))
    copy::RowAllToAllDemote( A, *this );
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,MC,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,VR] = [MC,STAR]"))
    BlockDistMatrix<T,MC,MR> A_MC_MR( A );
    *this = A_MC_MR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,STAR,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,VR] = [STAR,MR]"))
    copy::PartialRowFilter( A, *this );
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,VR] = [MD,STAR]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,STAR,MD>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,VR] = [STAR,MD]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,MR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,VR] = [MR,MC]"))
    BlockDistMatrix<T,STAR,VC> A_STAR_VC( A );
    *this = A_STAR_VC;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,MR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,VR] = [MR,STAR]"))
    auto A_MR_MC = MakeUnique<BlockDistMatrix<T,MR,MC>>( A );
    auto A_STAR_VC = MakeUnique<BlockDistMatrix<T,STAR,VC>>( *A_MR_MC );
    A_MR_MC.reset();
    *this = *A_STAR_VC;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,STAR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,VR] = [STAR,MC]"))
    BlockDistMatrix<T,STAR,VC> A_STAR_VC( A );
    *this = A_STAR_VC;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,VC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,VR] = [VC,STAR]"))
    BlockDistMatrix<T,MC,MR> A_MC_MR( A );
    *this = A_MC_MR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,STAR,VC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,VR] = [STAR,VC]"))
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,VR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,VR] = [VR,STAR]"))
    auto A_MR_MC = MakeUnique<BlockDistMatrix<T,MR,MC>>( A );
    auto A_STAR_VC = MakeUnique<BlockDistMatrix<T,STAR,VC>>( *A_MR_MC );
    A_MR_MC.reset(); 
    *this = *A_STAR_VC;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BDM& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,VR] = [STAR,VR]"))
    copy::Translate( A, *this );
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,VR] = [STAR,STAR]"))
    copy::RowFilter( A, *this );
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,CIRC,CIRC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,VR] = [CIRC,CIRC]"))
    LogicError("This routine is not yet written");
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

// Basic queries
// =============

template<typename T>
mpi::Comm BDM::DistComm() const { return this->grid_->VRComm(); }
template<typename T>
mpi::Comm BDM::CrossComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::RedundantComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::ColComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::RowComm() const { return this->grid_->VRComm(); }
template<typename T>
mpi::Comm BDM::PartialRowComm() const { return this->grid_->MRComm(); }
template<typename T>
mpi::Comm BDM::PartialUnionRowComm() const { return this->grid_->MCComm(); }

template<typename T>
int BDM::ColStride() const { return 1; }
template<typename T>
int BDM::RowStride() const { return this->grid_->VRSize(); }
template<typename T>
int BDM::PartialRowStride() const { return this->grid_->MRSize(); }
template<typename T>
int BDM::PartialUnionRowStride() const { return this->grid_->MCSize(); }
template<typename T>
int BDM::DistSize() const { return this->grid_->VRSize(); }
template<typename T>
int BDM::CrossSize() const { return 1; }
template<typename T>
int BDM::RedundantSize() const { return 1; }

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define SELF(T,U,V) \
  template BlockDistMatrix<T,COLDIST,ROWDIST>::BlockDistMatrix \
  ( const BlockDistMatrix<T,U,V>& A );
#define OTHER(T,U,V) \
  template BlockDistMatrix<T,COLDIST,ROWDIST>::BlockDistMatrix \
  ( const DistMatrix<T,U,V>& A ); \
  template BlockDistMatrix<T,COLDIST,ROWDIST>& \
           BlockDistMatrix<T,COLDIST,ROWDIST>::operator= \
           ( const DistMatrix<T,U,V>& A )
#define BOTH(T,U,V) \
  SELF(T,U,V); \
  OTHER(T,U,V)
#define PROTO(T) \
  template class BlockDistMatrix<T,COLDIST,ROWDIST>; \
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
  OTHER(T,STAR,VR  ); \
  BOTH( T,VC,  STAR); \
  BOTH( T,VR,  STAR);

#include "El/macros/Instantiate.h"

} // namespace El
