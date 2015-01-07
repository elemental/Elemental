/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#define COLDIST CIRC
#define ROWDIST CIRC

#include "./setup.hpp"

namespace El {

// Public section
// ##############

// Assignment and reconfiguration
// ==============================

template<typename T>
BDM& BDM::operator=( const AbstractBlockDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = ABDM"))
    copy::Gather( A, *this );
    return *this;
}

template<typename T>
BDM& BDM::operator=( const BDM& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [CIRC,CIRC]"))
    copy::Translate( A, *this );
    return *this;
}

template<typename T>
void BDM::CopyFromRoot( const Matrix<T>& A, bool includingViewers )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC]::CopyFromRoot"))
    if( this->CrossRank() != this->Root() )
        LogicError("Called CopyFromRoot from non-root");

    this->Resize( A.Height(), A.Width() );
    this->MakeSizeConsistent( includingViewers );

    this->matrix_ = A;
}

template<typename T>
void BDM::CopyFromNonRoot( bool includingViewers )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC]::CopyFromNonRoot"))
    if( this->CrossRank() == this->Root() )
        LogicError("Called CopyFromNonRoot from root");

    this->MakeSizeConsistent( includingViewers );
}

// Basic queries
// =============

template<typename T>
mpi::Comm BDM::DistComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::RedundantComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::CrossComm() const { return this->grid_->VCComm(); }
template<typename T>
mpi::Comm BDM::ColComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::RowComm() const { return mpi::COMM_SELF; }

template<typename T>
int BDM::ColStride() const { return 1; }
template<typename T>
int BDM::RowStride() const { return 1; }
template<typename T>
int BDM::DistSize() const { return 1; }
template<typename T>
int BDM::CrossSize() const { return this->grid_->VCSize(); }
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
  OTHER(T,CIRC,CIRC); \
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
  BOTH( T,VC,  STAR); \
  BOTH( T,VR,  STAR);

#include "El/macros/Instantiate.h"

} // namespace El
