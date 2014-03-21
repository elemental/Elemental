/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

#define ColDist STAR
#define RowDist MD

#include "./setup.hpp"

namespace elem {

// Public section
// ##############

// Assignment and reconfiguration
// ==============================

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MC,MR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,MD] = [MC,MR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MC,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,MD] = [MC,STAR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,MR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,MD] = [STAR,MR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,MD] = [MD,STAR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DM& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,MD] = [STAR,MD]"))
    A.Translate( *this );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MR,MC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,MD] = [MR,MC]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,MD] = [MR,STAR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,MC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,MD] = [STAR,MC]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,VC,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,MD] = [VC,STAR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,VC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,MD] = [STAR,VC]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,VR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,MD] = [VR,STAR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,VR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,MD] = [STAR,VR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,MD] = [STAR,STAR]"))
    this->RowFilterFrom( A );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,MD] = [CIRC,CIRC]"))
    DistMatrix<T,MC,MR> A_MC_MR( A.Grid() );
    A_MC_MR.AlignWith( *this );
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

// Basic queries
// =============

template<typename T>
mpi::Comm DM::DistComm() const { return this->grid_->MDComm(); }
template<typename T>
mpi::Comm DM::CrossComm() const { return this->grid_->MDPerpComm(); }
template<typename T>
mpi::Comm DM::RedundantComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::ColComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::RowComm() const { return this->grid_->MDComm(); }

template<typename T>
Int DM::ColStride() const { return 1; }
template<typename T>
Int DM::RowStride() const { return this->grid_->LCM(); }

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define PROTO(T) template class DistMatrix<T,ColDist,RowDist>
#define COPY(T,U,V) \
  template DistMatrix<T,ColDist,RowDist>::DistMatrix \
  ( const DistMatrix<T,U,V>& A )
#define FULL(T) \
  PROTO(T); \
  COPY(T,CIRC,CIRC); \
  COPY(T,MC,  MR  ); \
  COPY(T,MC,  STAR); \
  COPY(T,MD,  STAR); \
  COPY(T,MR,  MC  ); \
  COPY(T,MR,  STAR); \
  COPY(T,STAR,MC  ); \
  COPY(T,STAR,MR  ); \
  COPY(T,STAR,STAR); \
  COPY(T,STAR,VC  ); \
  COPY(T,STAR,VR  ); \
  COPY(T,VC,  STAR); \
  COPY(T,VR,  STAR);

FULL(Int);
#ifndef DISABLE_FLOAT
FULL(float);
#endif
FULL(double);

#ifndef DISABLE_COMPLEX
#ifndef DISABLE_FLOAT
FULL(Complex<float>);
#endif
FULL(Complex<double>);
#endif

} // namespace elem
