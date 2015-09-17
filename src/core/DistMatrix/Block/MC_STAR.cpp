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
BDM& BDM::operator=( const DistMatrix<T,MC,MR,BLOCK>& A )
{
    DEBUG_ONLY(CSE cse("[MC,STAR] = [MC,MR]"))
    copy::RowAllGather( A, *this );
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BDM& A )
{
    DEBUG_ONLY(CSE cse("[MC,STAR] = [MC,STAR]"))
    copy::Translate( A, *this );
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,STAR,MR,BLOCK>& A )
{ 
    DEBUG_ONLY(CSE cse("[MC,STAR] = [STAR,MR]"))
    DistMatrix<T,MC,MR,BLOCK> A_MC_MR(this->Grid());
    A_MC_MR.AlignColsWith(*this);
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,MD,STAR,BLOCK>& A )
{
    DEBUG_ONLY(CSE cse("[MC,STAR] = [MD,STAR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR,BLOCK> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,STAR,MD,BLOCK>& A )
{
    DEBUG_ONLY(CSE cse("[MC,STAR] = [STAR,MD]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR,BLOCK> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,MR,MC,BLOCK>& A )
{ 
    DEBUG_ONLY(CSE cse("[MC,STAR] = [MR,MC]"))
    auto A_VR_STAR = MakeUnique<DistMatrix<T,VR,STAR,BLOCK>>( A );
    auto A_VC_STAR =
      MakeUnique<DistMatrix<T,VC,STAR,BLOCK>>( this->Grid() );
    A_VC_STAR->AlignColsWith(*this);
    *A_VC_STAR = *A_VR_STAR;
    A_VR_STAR.reset(); 
    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,MR,STAR,BLOCK>& A )
{ 
    DEBUG_ONLY(CSE cse("[MC,STAR] = [MR,STAR]"))
    auto A_VR_STAR = MakeUnique<DistMatrix<T,VR,STAR,BLOCK>>( A );
    auto A_VC_STAR =
      MakeUnique<DistMatrix<T,VC,STAR,BLOCK>>( this->Grid() );
    A_VC_STAR->AlignColsWith(*this);
    *A_VC_STAR = *A_VR_STAR;
    A_VR_STAR.reset();
    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,STAR,MC,BLOCK>& A )
{ 
    DEBUG_ONLY(CSE cse("[MC,STAR] = [STAR,MC]"))
    auto A_MR_MC = MakeUnique<DistMatrix<T,MR,MC,BLOCK>>( A );
    auto A_VR_STAR = MakeUnique<DistMatrix<T,VR,STAR,BLOCK>>( *A_MR_MC );
    A_MR_MC.reset();

    auto A_VC_STAR =
      MakeUnique<DistMatrix<T,VC,STAR,BLOCK>>( this->Grid() );
    A_VC_STAR->AlignColsWith(*this);
    *A_VC_STAR = *A_VR_STAR;
    A_VR_STAR.reset(); 

    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,VC,STAR,BLOCK>& A )
{ 
    DEBUG_ONLY(CSE cse("[MC,STAR] = [VC,STAR]"))
    copy::PartialColAllGather( A, *this );
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,STAR,VC,BLOCK>& A )
{ 
    DEBUG_ONLY(CSE cse("[MC,STAR] = [STAR,VC]"))
    auto A_STAR_VR = MakeUnique<DistMatrix<T,STAR,VR,BLOCK>>( A );
    auto A_MC_MR = MakeUnique<DistMatrix<T,MC,MR,BLOCK>>( this->Grid() );
    A_MC_MR->AlignColsWith(*this);
    *A_MC_MR = *A_STAR_VR;
    A_STAR_VR.reset();
    *this = *A_MC_MR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,VR,STAR,BLOCK>& A )
{ 
    DEBUG_ONLY(CSE cse("[MC,STAR] = [VR,STAR]"))
    DistMatrix<T,VC,STAR,BLOCK> A_VC_STAR(this->Grid());
    A_VC_STAR.AlignColsWith(*this);
    A_VC_STAR = A;
    *this = A_VC_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,STAR,VR,BLOCK>& A )
{ 
    DEBUG_ONLY(CSE cse("[MC,STAR] = [STAR,VR]"))
    DistMatrix<T,MC,MR,BLOCK> A_MC_MR(this->Grid());
    A_MC_MR.AlignColsWith( *this );
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,STAR,STAR,BLOCK>& A )
{
    DEBUG_ONLY(CSE cse("[MC,STAR] = [STAR,STAR]"))
    copy::ColFilter( A, *this );
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,CIRC,CIRC,BLOCK>& A )
{
    DEBUG_ONLY(CSE cse("[MC,STAR] = [CIRC,CIRC]"))
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("BDM = ABDM"))
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = \
        dynamic_cast<const DistMatrix<T,CDIST,RDIST,BLOCK>&>(A); \
      *this = ACast;
    #include "El/macros/GuardAndPayload.h"
    return *this;
}

// Basic queries
// =============

template<typename T>
mpi::Comm BDM::DistComm() const EL_NO_EXCEPT
{ return this->grid_->MCComm(); }
template<typename T>
mpi::Comm BDM::RedundantComm() const EL_NO_EXCEPT
{ return this->grid_->MRComm(); }
template<typename T>
mpi::Comm BDM::CrossComm() const EL_NO_EXCEPT
{ return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::ColComm() const EL_NO_EXCEPT
{ return this->grid_->MCComm(); }
template<typename T>
mpi::Comm BDM::RowComm() const EL_NO_EXCEPT
{ return mpi::COMM_SELF; }

template<typename T>
int BDM::ColStride() const EL_NO_EXCEPT
{ return this->grid_->MCSize(); }
template<typename T>
int BDM::RowStride() const EL_NO_EXCEPT
{ return 1; }
template<typename T>
int BDM::DistSize() const EL_NO_EXCEPT
{ return this->grid_->MCSize(); }
template<typename T>
int BDM::CrossSize() const EL_NO_EXCEPT
{ return 1; }
template<typename T>
int BDM::RedundantSize() const EL_NO_EXCEPT
{ return this->grid_->MRSize(); }

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define SELF(T,U,V) \
  template DistMatrix<T,COLDIST,ROWDIST,BLOCK>::DistMatrix \
  ( const DistMatrix<T,U,V,BLOCK>& A );
#define OTHER(T,U,V) \
  template DistMatrix<T,COLDIST,ROWDIST,BLOCK>::DistMatrix \
  ( const DistMatrix<T,U,V>& A ); \
  template DistMatrix<T,COLDIST,ROWDIST,BLOCK>& \
           DistMatrix<T,COLDIST,ROWDIST,BLOCK>::operator= \
           ( const DistMatrix<T,U,V>& A )
#define BOTH(T,U,V) \
  SELF(T,U,V); \
  OTHER(T,U,V)
#define PROTO(T) \
  template class DistMatrix<T,COLDIST,ROWDIST,BLOCK>; \
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

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
