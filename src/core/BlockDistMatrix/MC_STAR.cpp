/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#define COLDIST MC
#define ROWDIST STAR

#include "./setup.hpp"

namespace El {

// Public section
// ##############

// Assignment and reconfiguration
// ==============================

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,MC,MR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [MC,MR]"))
    copy::RowAllGather( A, *this );
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BDM& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [MC,STAR]"))
    copy::Translate( A, *this );
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,STAR,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [STAR,MR]"))
    BlockDistMatrix<T,MC,MR> A_MC_MR(this->Grid());
    A_MC_MR.AlignColsWith(*this);
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [MD,STAR]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,STAR,MD>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [STAR,MD]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,MR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [MR,MC]"))
    auto A_VR_STAR = MakeUnique<BlockDistMatrix<T,VR,STAR>>( A );
    auto A_VC_STAR = MakeUnique<BlockDistMatrix<T,VC,STAR>>( this->Grid() );
    A_VC_STAR->AlignColsWith(*this);
    *A_VC_STAR = *A_VR_STAR;
    A_VR_STAR.reset(); 
    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,MR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [MR,STAR]"))
    auto A_VR_STAR = MakeUnique<BlockDistMatrix<T,VR,STAR>>( A );
    auto A_VC_STAR = MakeUnique<BlockDistMatrix<T,VC,STAR>>( this->Grid() );
    A_VC_STAR->AlignColsWith(*this);
    *A_VC_STAR = *A_VR_STAR;
    A_VR_STAR.reset();
    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,STAR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [STAR,MC]"))
    auto A_MR_MC = MakeUnique<BlockDistMatrix<T,MR,MC>>( A );
    auto A_VR_STAR = MakeUnique<BlockDistMatrix<T,VR,STAR>>( *A_MR_MC );
    A_MR_MC.reset();

    auto A_VC_STAR = MakeUnique<BlockDistMatrix<T,VC,STAR>>( this->Grid() );
    A_VC_STAR->AlignColsWith(*this);
    *A_VC_STAR = *A_VR_STAR;
    A_VR_STAR.reset(); 

    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,VC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [VC,STAR]"))
    copy::PartialColAllGather( A, *this );
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,STAR,VC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [STAR,VC]"))
    auto A_STAR_VR = MakeUnique<BlockDistMatrix<T,STAR,VR>>( A );
    auto A_MC_MR = MakeUnique<BlockDistMatrix<T,MC,MR>>( this->Grid() );
    A_MC_MR->AlignColsWith(*this);
    *A_MC_MR = *A_STAR_VR;
    A_STAR_VR.reset();
    *this = *A_MC_MR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,VR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [VR,STAR]"))
    BlockDistMatrix<T,VC,STAR> A_VC_STAR(this->Grid());
    A_VC_STAR.AlignColsWith(*this);
    A_VC_STAR = A;
    *this = A_VC_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,STAR,VR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [STAR,VR]"))
    BlockDistMatrix<T,MC,MR> A_MC_MR(this->Grid());
    A_MC_MR.AlignColsWith( *this );
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [STAR,STAR]"))
    copy::ColFilter( A, *this );
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockDistMatrix<T,CIRC,CIRC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [CIRC,CIRC]"))
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
mpi::Comm BDM::DistComm() const { return this->grid_->MCComm(); }
template<typename T>
mpi::Comm BDM::RedundantComm() const { return this->grid_->MRComm(); }
template<typename T>
mpi::Comm BDM::CrossComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::ColComm() const { return this->grid_->MCComm(); }
template<typename T>
mpi::Comm BDM::RowComm() const { return mpi::COMM_SELF; }

template<typename T>
int BDM::ColStride() const { return this->grid_->MCSize(); }
template<typename T>
int BDM::RowStride() const { return 1; }
template<typename T>
int BDM::DistSize() const { return this->grid_->MCSize(); }
template<typename T>
int BDM::CrossSize() const { return 1; }
template<typename T>
int BDM::RedundantSize() const { return this->grid_->MRSize(); }

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
  OTHER(T,MC,  STAR); \
  BOTH( T,MD,  STAR); \
  BOTH( T,MR,  MC  ); \
  BOTH( T,MR,  STAR); \
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
