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

// Return a view
// -------------

template<typename T>
DM DM::operator()( Range<Int> indVert, Range<Int> indHorz )
{
    DEBUG_ONLY(CallStackEntry cse("DM[CIRC,CIRC]( ind, ind )"))
    if( this->Locked() )
        return LockedView( *this, indVert, indHorz );
    else
        return View( *this, indVert, indHorz );
}

template<typename T>
const DM DM::operator()( Range<Int> indVert, Range<Int> indHorz ) const
{
    DEBUG_ONLY(CallStackEntry cse("DM[CIRC,CIRC]( ind, ind )"))
    return LockedView( *this, indVert, indHorz );
}

// Make a copy
// -----------
template<typename T>
DM& DM::operator=( const DM& A )
{
    DEBUG_ONLY(CallStackEntry cse("DM[CIRC,CIRC] = DM[CIRC,CIRC]"))
    copy::Translate( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = ADM"));
    copy::Gather( A, *this );
    return *this;
}

// TODO: operator= for BlockDistMatrix

template<typename T>
void DM::CopyFromRoot( const Matrix<T>& A, bool includingViewers )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC]::CopyFromRoot"))
    if( this->CrossRank() != this->Root() )
        LogicError("Called CopyFromRoot from non-root");

    this->Resize( A.Height(), A.Width() );
    this->MakeSizeConsistent( includingViewers );

    this->matrix_ = A;
}

template<typename T>
void DM::CopyFromNonRoot( bool includingViewers )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC]::CopyFromNonRoot"))
    if( this->CrossRank() == this->Root() )
        LogicError("Called CopyFromNonRoot from root");

    this->MakeSizeConsistent( includingViewers );
}

// Basic queries
// =============

template<typename T>
mpi::Comm DM::DistComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::RedundantComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::CrossComm() const { return this->grid_->VCComm(); }
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
int DM::CrossSize() const { return this->grid_->VCSize(); }
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
