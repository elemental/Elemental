/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#define COLDIST STAR
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
    DEBUG_ONLY(CallStackEntry cse("DM[STAR,STAR]( ind, ind )"))
    if( this->Locked() )
        return LockedView( *this, indVert, indHorz );
    else
        return View( *this, indVert, indHorz );
}

template<typename T>
const DM DM::operator()( Range<Int> indVert, Range<Int> indHorz ) const
{
    DEBUG_ONLY(CallStackEntry cse("DM[STAR,STAR]( ind, ind )"))
    return LockedView( *this, indVert, indHorz );
}

// Make a copy
// -----------
template<typename T>
DM& DM::operator=( const DM& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,STAR] = [STAR,STAR]"))
    copy::Translate( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MC,MR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,STAR] = [MC,MR]"))
    copy::AllGather( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MC,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,STAR] = [MC,STAR]"))
    copy::ColAllGather( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,MR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,STAR] = [STAR,MR]"))
    copy::RowAllGather( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,STAR] = [MD,STAR]"))
    copy::ColAllGather( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,STAR] = [STAR,MD]"))
    copy::RowAllGather( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MR,MC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,STAR] = [MR,MC]"))
    copy::AllGather( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,STAR] = [MR,STAR]"))
    copy::ColAllGather( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,MC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,STAR] = [STAR,MC]"))
    copy::RowAllGather( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,VC,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,STAR] = [VC,STAR]"))
    copy::ColAllGather( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,VC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,STAR] = [STAR,VC]"))
    copy::RowAllGather( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,VR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,STAR] = [VR,STAR]"))
    copy::ColAllGather( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,VR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,STAR] = [STAR,VR]"))
    copy::RowAllGather( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,STAR] = [CIRC,CIRC]"))
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
mpi::Comm DM::DistComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::CrossComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::RedundantComm() const { return this->grid_->VCComm(); }
template<typename T>
mpi::Comm DM::ColComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::RowComm() const { return mpi::COMM_SELF; }

template<typename T>
int DM::ColStride() const { return 1; }
template<typename T>
int DM::RowStride() const { return 1; }
template<typename T>
int DM::DistSize() const { return 1; }
template<typename T>
int DM::CrossSize() const { return 1; }
template<typename T>
int DM::RedundantSize() const { return this->grid_->VCSize(); }

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
  OTHER(T,STAR,STAR); \
  BOTH( T,STAR,VC  ); \
  BOTH( T,STAR,VR  ); \
  BOTH( T,VC,  STAR); \
  BOTH( T,VR,  STAR);

#include "El/macros/Instantiate.h"

} // namespace El
