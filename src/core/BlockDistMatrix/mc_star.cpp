/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

#define ColDist MC
#define RowDist STAR

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
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [MC,MR]"))
    A.RowAllGather( *this );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BDM& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [MC,STAR]"))
    A.Translate( *this );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [STAR,MR]"))
    BlockDistMatrix<T,MC,MR> A_MC_MR(this->Grid());
    A_MC_MR.AlignColsWith(*this);
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [MD,STAR]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,MD>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [STAR,MD]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [MR,MC]"))
    std::unique_ptr<BlockDistMatrix<T,VR,STAR>> A_VR_STAR
    ( new BlockDistMatrix<T,VR,STAR>(A) );
    std::unique_ptr<BlockDistMatrix<T,VC,STAR>> A_VC_STAR
    ( new BlockDistMatrix<T,VC,STAR>(this->Grid()) );
    A_VC_STAR->AlignColsWith(*this);
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater
    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [MR,STAR]"))
    std::unique_ptr<BlockDistMatrix<T,VR,STAR>> A_VR_STAR
    ( new BlockDistMatrix<T,VR,STAR>(A) );
    std::unique_ptr<BlockDistMatrix<T,VC,STAR>> A_VC_STAR
    ( new BlockDistMatrix<T,VC,STAR>(this->Grid()) );
    A_VC_STAR->AlignColsWith(*this);
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater
    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [STAR,MC]"))
    std::unique_ptr<BlockDistMatrix<T,MR,MC>>
        A_MR_MC( new BlockDistMatrix<T,MR,MC>(A) );
    std::unique_ptr<BlockDistMatrix<T,VR,STAR>>
        A_VR_STAR( new BlockDistMatrix<T,VR,STAR>(*A_MR_MC) );
    delete A_MR_MC.release(); // lowers memory highwater

    std::unique_ptr<BlockDistMatrix<T,VC,STAR>>
        A_VC_STAR( new BlockDistMatrix<T,VC,STAR>(this->Grid()) );
    A_VC_STAR->AlignColsWith(*this);
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater

    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,VC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [VC,STAR]"))
    A.PartialColAllGather( *this );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,VC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [STAR,VC]"))
    std::unique_ptr<BlockDistMatrix<T,STAR,VR>>
        A_STAR_VR( new BlockDistMatrix<T,STAR,VR>(A) );
    std::unique_ptr<BlockDistMatrix<T,MC,MR>>
        A_MC_MR
        ( new BlockDistMatrix<T,MC,MR>(this->Grid()) );
    A_MC_MR->AlignColsWith(*this);
    *A_MC_MR = *A_STAR_VR;
    delete A_STAR_VR.release(); // lowers memory highwater
    *this = *A_MC_MR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,VR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [VR,STAR]"))
    BlockDistMatrix<T,VC,STAR> A_VC_STAR(this->Grid());
    A_VC_STAR.AlignColsWith(*this);
    A_VC_STAR = A;
    *this = A_VC_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,VR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [STAR,VR]"))
    BlockDistMatrix<T,MC,MR> A_MC_MR(this->Grid());
    A_MC_MR.AlignColsWith( *this );
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [STAR,STAR]"))
    this->ColFilterFrom( A );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,CIRC,CIRC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [CIRC,CIRC]"))
    LogicError("This routine is not yet written");
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
Int BDM::ColStride() const { return this->grid_->MCSize(); }
template<typename T>
Int BDM::RowStride() const { return 1; }
template<typename T>
Int BDM::DistSize() const { return this->grid_->MCSize(); }
template<typename T>
Int BDM::CrossSize() const { return 1; }
template<typename T>
Int BDM::RedundantSize() const { return this->grid_->MRSize(); }

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define PROTO(T) template class BlockDistMatrix<T,ColDist,RowDist>
#define SELF(T,U,V) \
  template BlockDistMatrix<T,ColDist,RowDist>::BlockDistMatrix \
  ( const BlockDistMatrix<T,U,V>& A );
#define OTHER(T,U,V) \
  template BlockDistMatrix<T,ColDist,RowDist>::BlockDistMatrix \
  ( const DistMatrix<T,U,V>& A ); \
  template BlockDistMatrix<T,ColDist,RowDist>& \
           BlockDistMatrix<T,ColDist,RowDist>::operator= \
           ( const DistMatrix<T,U,V>& A )
#define BOTH(T,U,V) \
  SELF(T,U,V); \
  OTHER(T,U,V)
#define FULL(T) \
  PROTO(T); \
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

FULL(Int);
#ifndef ELEM_DISABLE_FLOAT
FULL(float);
#endif
FULL(double);

#ifndef ELEM_DISABLE_COMPLEX
#ifndef ELEM_DISABLE_FLOAT
FULL(Complex<float>);
#endif
FULL(Complex<double>);
#endif 

} // namespace elem
