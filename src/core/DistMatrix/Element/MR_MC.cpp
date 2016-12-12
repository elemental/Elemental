/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like.hpp>

#define COLDIST MR
#define ROWDIST MC

#include "./setup.hpp"

namespace El {

// Public section
// ##############

// Assignment and reconfiguration
// ==============================

// Make a copy
// -----------
template<typename T>
DM& DM::operator=( const DistMatrix<T,MC,MR>& A )
{
    EL_DEBUG_CSE
    const Grid& grid = A.Grid();
    if( grid.Height() == grid.Width() )
    {
        const int gridDim = grid.Height();
        const int sendRank = this->ColOwner(A.RowShift()) +
                             this->RowOwner(A.ColShift())*gridDim;
        const int recvRank = A.RowOwner(this->ColShift()) +
                             A.ColOwner(this->RowShift())*gridDim;
        copy::Exchange( A, *this, sendRank, recvRank, grid.VCComm() );
    }
    else
    {
        copy::TransposeDist( A, *this );
    }
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MC,STAR>& A )
{
    EL_DEBUG_CSE
    DistMatrix<T,VC,STAR> A_VC_STAR( A );
    DistMatrix<T,VR,STAR> A_VR_STAR( this->Grid() );
    A_VR_STAR.AlignColsWith(*this);
    A_VR_STAR = A_VC_STAR;
    A_VC_STAR.Empty();
    *this = A_VR_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,MR>& A )
{
    EL_DEBUG_CSE
    DistMatrix<T,STAR,VR> A_STAR_VR( A );
    DistMatrix<T,STAR,VC> A_STAR_VC( this->Grid() );
    A_STAR_VC.AlignRowsWith(*this);
    A_STAR_VC = A_STAR_VR;
    A_STAR_VR.Empty();
    *this = A_STAR_VC;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MD,STAR>& A )
{
    EL_DEBUG_CSE
    // TODO: More efficient implementation
    copy::GeneralPurpose( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,MD>& A )
{
    EL_DEBUG_CSE
    // TODO: More efficient implementation
    copy::GeneralPurpose( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MR,STAR>& A )
{
    EL_DEBUG_CSE
    copy::RowFilter( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,MC>& A )
{
    EL_DEBUG_CSE
    copy::ColFilter( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,VC,STAR>& A )
{
    EL_DEBUG_CSE
    DistMatrix<T,VR,STAR> A_VR_STAR( A );
    *this = A_VR_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,VC>& A )
{
    EL_DEBUG_CSE
    copy::RowAllToAllPromote( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,VR,STAR>& A )
{
    EL_DEBUG_CSE
    copy::ColAllToAllPromote( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,VR>& A )
{
    EL_DEBUG_CSE
    DistMatrix<T,STAR,VC> A_STAR_VC( A );
    *this = A_STAR_VC;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,STAR>& A )
{
    EL_DEBUG_CSE
    copy::Filter( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
    EL_DEBUG_CSE
    copy::Scatter( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const ElementalMatrix<T>& A )
{
    EL_DEBUG_CSE
    #define GUARD(CDIST,RDIST,WRAP) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST && \
      ELEMENT == WRAP
    #define PAYLOAD(CDIST,RDIST,WRAP) \
      auto& ACast = static_cast<const DistMatrix<T,CDIST,RDIST>&>(A); \
      *this = ACast;
    #include "El/macros/GuardAndPayload.h"
    return *this;
}

// Basic queries
// =============
template<typename T>
mpi::Comm DM::DistComm() const EL_NO_EXCEPT { return this->grid_->VRComm(); }

template<typename T>
mpi::Comm DM::CrossComm() const EL_NO_EXCEPT
{ return ( this->Grid().InGrid() ? mpi::COMM_SELF : mpi::COMM_NULL ); }

template<typename T>
mpi::Comm DM::RedundantComm() const EL_NO_EXCEPT
{ return ( this->Grid().InGrid() ? mpi::COMM_SELF : mpi::COMM_NULL ); }

template<typename T>
mpi::Comm DM::ColComm() const EL_NO_EXCEPT { return this->grid_->MRComm(); }
template<typename T>
mpi::Comm DM::RowComm() const EL_NO_EXCEPT { return this->grid_->MCComm(); }

template<typename T>
mpi::Comm DM::PartialColComm() const EL_NO_EXCEPT { return this->ColComm(); }
template<typename T>
mpi::Comm DM::PartialRowComm() const EL_NO_EXCEPT { return this->RowComm(); }
template<typename T>
mpi::Comm DM::PartialUnionColComm() const EL_NO_EXCEPT
{ return ( this->Grid().InGrid() ? mpi::COMM_SELF : mpi::COMM_NULL ); }
template<typename T>
mpi::Comm DM::PartialUnionRowComm() const EL_NO_EXCEPT
{ return ( this->Grid().InGrid() ? mpi::COMM_SELF : mpi::COMM_NULL ); }

template<typename T>
int DM::ColStride() const EL_NO_EXCEPT { return this->grid_->MRSize(); }
template<typename T>
int DM::RowStride() const EL_NO_EXCEPT { return this->grid_->MCSize(); }
template<typename T>
int DM::DistSize() const EL_NO_EXCEPT { return this->grid_->VRSize(); }
template<typename T>
int DM::CrossSize() const EL_NO_EXCEPT { return 1; }
template<typename T>
int DM::RedundantSize() const EL_NO_EXCEPT { return 1; }
template<typename T>
int DM::PartialColStride() const EL_NO_EXCEPT { return this->ColStride(); }
template<typename T>
int DM::PartialRowStride() const EL_NO_EXCEPT { return this->RowStride(); }
template<typename T>
int DM::PartialUnionColStride() const EL_NO_EXCEPT { return 1; }
template<typename T>
int DM::PartialUnionRowStride() const EL_NO_EXCEPT { return 1; }

template<typename T>
int DM::ColRank() const EL_NO_EXCEPT { return this->grid_->MRRank(); }
template<typename T>
int DM::RowRank() const EL_NO_EXCEPT { return this->grid_->MCRank(); }
template<typename T>
int DM::DistRank() const EL_NO_EXCEPT { return this->grid_->VRRank(); }
template<typename T>
int DM::CrossRank() const EL_NO_EXCEPT
{ return ( this->Grid().InGrid() ? 0 : mpi::UNDEFINED ); }
template<typename T>
int DM::RedundantRank() const EL_NO_EXCEPT
{ return ( this->Grid().InGrid() ? 0 : mpi::UNDEFINED ); }
template<typename T>
int DM::PartialColRank() const EL_NO_EXCEPT { return this->ColRank(); }
template<typename T>
int DM::PartialRowRank() const EL_NO_EXCEPT { return this->RowRank(); }
template<typename T>
int DM::PartialUnionColRank() const EL_NO_EXCEPT
{ return ( this->Grid().InGrid() ? 0 : mpi::UNDEFINED ); }
template<typename T>
int DM::PartialUnionRowRank() const EL_NO_EXCEPT
{ return ( this->Grid().InGrid() ? 0 : mpi::UNDEFINED ); }

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define SELF(T,U,V) \
  template DistMatrix<T,COLDIST,ROWDIST>::DistMatrix \
  ( const DistMatrix<T,U,V>& A );
#define OTHER(T,U,V) \
  template DistMatrix<T,COLDIST,ROWDIST>::DistMatrix \
  ( const DistMatrix<T,U,V,BLOCK>& A ); \
  template DistMatrix<T,COLDIST,ROWDIST>& \
           DistMatrix<T,COLDIST,ROWDIST>::operator= \
           ( const DistMatrix<T,U,V,BLOCK>& A )
#define BOTH(T,U,V) \
  SELF(T,U,V); \
  OTHER(T,U,V)
#define PROTO(T) \
  template class DistMatrix<T,COLDIST,ROWDIST>; \
  BOTH( T,CIRC,CIRC); \
  BOTH( T,MC,  MR  ); \
  BOTH( T,MC,  STAR); \
  BOTH( T,MD,  STAR); \
  OTHER(T,MR,  MC  ); \
  BOTH( T,MR,  STAR); \
  BOTH( T,STAR,MC  ); \
  BOTH( T,STAR,MD  ); \
  BOTH( T,STAR,MR  ); \
  BOTH( T,STAR,STAR); \
  BOTH( T,STAR,VC  ); \
  BOTH( T,STAR,VR  ); \
  BOTH( T,VC,  STAR); \
  BOTH( T,VR,  STAR);

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
