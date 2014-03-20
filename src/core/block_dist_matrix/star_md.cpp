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
BDM&
BDM::operator=( const BlockDistMatrix<T,MC,MR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MD,STAR] = [MC,MR]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MC,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MD,STAR] = [MC,STAR]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MD,STAR] = [STAR,MR]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MD,STAR] = [MD,STAR]"))
    LogicError("This routine is not yet implemented");
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BDM& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MD,STAR] = [STAR,MD]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MD,STAR] = [MR,MC]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MD,STAR] = [MR,STAR]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MD,STAR] = [STAR,MC]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,VC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MD,STAR] = [VC,STAR]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,VC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MD,STAR] = [STAR,VC]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,VR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MD,STAR] = [VR,STAR]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,VR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MD,STAR] = [STAR,VR]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MD,STAR] = [STAR,STAR]"))
    this->ColFilterFrom( A );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,CIRC,CIRC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MD,STAR] = [CIRC,CIRC]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

// Realignment
// -----------

template<typename T>
void
BDM::AlignWith( const elem::BlockDistData& data )
{
    DEBUG_ONLY(CallStackEntry cse("[MD,STAR]::AlignWith"))
    this->AlignRowsWith( data );
}

template<typename T>
void
BDM::AlignRowsWith( const elem::BlockDistData& data )
{
    DEBUG_ONLY(CallStackEntry cse("[MD,STAR]::AlignRowsWith"))
    this->SetGrid( *data.grid );
    if( data.colDist == MD && data.rowDist == STAR )
    {
        this->SetRoot( data.root );
        this->AlignRows( data.blockHeight, data.colAlign, data.colCut );
    }
    else if( data.colDist == STAR && data.rowDist == MD )
    {
        this->SetRoot( data.root );
        this->AlignRows( data.blockWidth, data.rowAlign, data.rowCut );
    }
    DEBUG_ONLY(else LogicError("Nonsensical alignment"))
}

// Basic queries
// =============

template<typename T>
mpi::Comm BDM::DistComm() const { return this->grid_->MDComm(); }
template<typename T>
mpi::Comm BDM::CrossComm() const { return this->grid_->MDPerpComm(); }
template<typename T>
mpi::Comm BDM::RedundantComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::ColComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::RowComm() const { return this->grid_->MDComm(); }

template<typename T>
Int BDM::ColStride() const { return 1; }
template<typename T>
Int BDM::RowStride() const { return this->grid_->LCM(); }

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define PROTO(T) template class BlockDistMatrix<T,ColDist,RowDist>
#define COPY(T,U,V) \
  template BlockDistMatrix<T,ColDist,RowDist>::BlockDistMatrix\
  ( const BlockDistMatrix<T,U,V>& A );\
  template BlockDistMatrix<T,ColDist,RowDist>::BlockDistMatrix\
  ( const DistMatrix<T,U,V>& A );
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
