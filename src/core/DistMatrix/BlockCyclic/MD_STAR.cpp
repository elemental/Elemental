/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#define COLDIST MD
#define ROWDIST STAR

#include "./setup.hpp"

namespace El {

// Public section
// ##############

// Assignment and reconfiguration
// ==============================

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,MC,MR,BLOCK_CYCLIC>& A )
{
    DEBUG_ONLY(CSE cse("[MD,STAR] = [MC,MR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR,BLOCK_CYCLIC> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,MC,STAR,BLOCK_CYCLIC>& A )
{
    DEBUG_ONLY(CSE cse("[MD,STAR] = [MC,STAR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR,BLOCK_CYCLIC> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,STAR,MR,BLOCK_CYCLIC>& A )
{ 
    DEBUG_ONLY(CSE cse("[MD,STAR] = [STAR,MR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR,BLOCK_CYCLIC> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BDM& A )
{
    DEBUG_ONLY(CSE cse("[MD,STAR] = [MD,STAR]"))
    copy::Translate( A, *this );
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,STAR,MD,BLOCK_CYCLIC>& A )
{
    DEBUG_ONLY(CSE cse("[MD,STAR] = [STAR,MD]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR,BLOCK_CYCLIC> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,MR,MC,BLOCK_CYCLIC>& A )
{ 
    DEBUG_ONLY(CSE cse("[MD,STAR] = [MR,MC]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR,BLOCK_CYCLIC> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,MR,STAR,BLOCK_CYCLIC>& A )
{ 
    DEBUG_ONLY(CSE cse("[MD,STAR] = [MR,STAR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR,BLOCK_CYCLIC> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,STAR,MC,BLOCK_CYCLIC>& A )
{ 
    DEBUG_ONLY(CSE cse("[MD,STAR] = [STAR,MC]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR,BLOCK_CYCLIC> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,VC,STAR,BLOCK_CYCLIC>& A )
{ 
    DEBUG_ONLY(CSE cse("[MD,STAR] = [VC,STAR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR,BLOCK_CYCLIC> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,STAR,VC,BLOCK_CYCLIC>& A )
{ 
    DEBUG_ONLY(CSE cse("[MD,STAR] = [STAR,VC]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR,BLOCK_CYCLIC> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,VR,STAR,BLOCK_CYCLIC>& A )
{ 
    DEBUG_ONLY(CSE cse("[MD,STAR] = [VR,STAR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR,BLOCK_CYCLIC> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,STAR,VR,BLOCK_CYCLIC>& A )
{ 
    DEBUG_ONLY(CSE cse("[MD,STAR] = [STAR,VR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR,BLOCK_CYCLIC> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,STAR,STAR,BLOCK_CYCLIC>& A )
{
    DEBUG_ONLY(CSE cse("[MD,STAR] = [STAR,STAR]"))
    copy::ColFilter( A, *this );
    return *this;
}

template<typename T>
BDM& BDM::operator=( const DistMatrix<T,CIRC,CIRC,BLOCK_CYCLIC>& A )
{
    DEBUG_ONLY(CSE cse("[MD,STAR] = [CIRC,CIRC]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR,BLOCK_CYCLIC> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BlockCyclicMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("BDM = ABDM"))
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = \
        dynamic_cast<const DistMatrix<T,CDIST,RDIST,BLOCK_CYCLIC>&>(A); \
      *this = ACast;
    #include "El/macros/GuardAndPayload.h"
    return *this;
}

// Basic queries
// =============

template<typename T>
mpi::Comm BDM::DistComm() const EL_NO_EXCEPT
{ return this->grid_->MDComm(); }
template<typename T>
mpi::Comm BDM::CrossComm() const EL_NO_EXCEPT
{ return this->grid_->MDPerpComm(); }
template<typename T>
mpi::Comm BDM::RedundantComm() const EL_NO_EXCEPT
{ return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::ColComm() const EL_NO_EXCEPT
{ return this->grid_->MDComm(); }
template<typename T>
mpi::Comm BDM::RowComm() const EL_NO_EXCEPT
{ return mpi::COMM_SELF; }

template<typename T>
int BDM::ColStride() const EL_NO_EXCEPT
{ return this->grid_->LCM(); }
template<typename T>
int BDM::RowStride() const EL_NO_EXCEPT
{ return 1; }
template<typename T>
int BDM::DistSize() const EL_NO_EXCEPT
{ return this->grid_->LCM(); }
template<typename T>
int BDM::CrossSize() const EL_NO_EXCEPT
{ return this->grid_->GCD(); }
template<typename T>
int BDM::RedundantSize() const EL_NO_EXCEPT
{ return 1; }

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define SELF(T,U,V) \
  template DistMatrix<T,COLDIST,ROWDIST,BLOCK_CYCLIC>::DistMatrix \
  ( const DistMatrix<T,U,V,BLOCK_CYCLIC>& A );
#define OTHER(T,U,V) \
  template DistMatrix<T,COLDIST,ROWDIST,BLOCK_CYCLIC>::DistMatrix \
  ( const DistMatrix<T,U,V>& A ); \
  template DistMatrix<T,COLDIST,ROWDIST,BLOCK_CYCLIC>& \
           DistMatrix<T,COLDIST,ROWDIST,BLOCK_CYCLIC>::operator= \
           ( const DistMatrix<T,U,V>& A )
#define BOTH(T,U,V) \
  SELF(T,U,V); \
  OTHER(T,U,V)
#define PROTO(T) \
  template class DistMatrix<T,COLDIST,ROWDIST,BLOCK_CYCLIC>; \
  BOTH( T,CIRC,CIRC); \
  BOTH( T,MC,  MR  ); \
  BOTH( T,MC,  STAR); \
  OTHER(T,MD,  STAR); \
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
