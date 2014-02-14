/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

namespace elem {

#define DM DistMatrix<T,STAR,STAR>
#define GDM GeneralDistMatrix<T,STAR,STAR>

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T>
DM::DistMatrix( const elem::Grid& g )
: GDM(g)
{ }

template<typename T>
DM::DistMatrix( Int height, Int width, const elem::Grid& g )
: GDM(g)
{ this->Resize(height,width); }

template<typename T>
DM::DistMatrix( Int height, Int width, Int ldim, const elem::Grid& g )
: GDM(g)
{ this->Resize(height,width,ldim); }

template<typename T>
DM::DistMatrix
( Int height, Int width, const T* buffer, Int ldim, const elem::Grid& g )
: GDM(g)
{ this->LockedAttach(height,width,0,0,buffer,ldim,g); }

template<typename T>
DM::DistMatrix
( Int height, Int width, T* buffer, Int ldim, const elem::Grid& g )
: GDM(g)
{ this->Attach(height,width,0,0,buffer,ldim,g); }

template<typename T>
DM::DistMatrix( const DM& A )
: GDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[* ,* ]::DistMatrix"))
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct [* ,* ] with itself");
}

template<typename T>
template<Dist U,Dist V>
DM::DistMatrix( const DistMatrix<T,U,V>& A )
: GDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[* ,* ]::DistMatrix"))
    if( STAR != U || STAR != V || 
        reinterpret_cast<const DM*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct [* ,* ] with itself");
}

template<typename T>
DM::DistMatrix( DM&& A ) : GDM(std::move(A)) { }

template<typename T> DM::~DistMatrix() { }

// Assignment and reconfiguration
// ==============================

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,MC,MR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,* ] = [MC,MR]"))
    A.AllGather( *this );
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,MC,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,* ] = [MC,* ]"))
    A.ColAllGather( *this );
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,STAR,MR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,* ] = [* ,MR]"))
    A.RowAllGather( *this );
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,* ] = [MD,* ]"))
    A.ColAllGather( *this );
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,* ] = [* ,MD]"))
    A.RowAllGather( *this );
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,MR,MC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,* ] = [MR,MC]"))
    A.AllGather( *this );
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,MR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,* ] = [MR,* ]"))
    A.ColAllGather( *this );
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,STAR,MC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,* ] = [* ,MC]"))
    A.RowAllGather( *this );
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,VC,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,* ] = [VC,* ]"))
    A.ColAllGather( *this );
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,STAR,VC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,* ] = [* ,VC]"))
    A.RowAllGather( *this );
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,VR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,* ] = [VR,* ]"))
    A.ColAllGather( *this );
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,STAR,VR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,* ] = [* ,VR]"))
    A.RowAllGather( *this );
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DM& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[* ,* ] = [* ,* ]");
        this->AssertNotLocked();
    )
    this->Resize( A.Height(), A.Width() );

    if( this->Grid() == A.Grid() )
    {
        this->matrix_ = A.LockedMatrix();
    }
    else
    {
        // TODO: Remember why I wrote this...
        if( !mpi::CongruentComms( A.Grid().ViewingComm(),
                                  this->Grid().ViewingComm() ) )
            LogicError
            ("Redistributing between nonmatching grids currently requires"
             " the viewing communicators to match.");

        // Compute and allocate the amount of required memory
        Int requiredMemory = 0;
        if( A.Grid().VCRank() == 0 )
            requiredMemory += A.Height()*A.Width();
        if( this->Participating() )
            requiredMemory += A.Height()*A.Width();
        T* buffer = this->auxMemory_.Require( requiredMemory );
        Int offset = 0;
        T* sendBuf = &buffer[offset];
        if( A.Grid().VCRank() == 0 )
            offset += A.Height()*A.Width();
        T* bcastBuffer = &buffer[offset];

        // Send from the root of A to the root of this matrix's grid
        mpi::Request sendRequest;
        if( A.Grid().VCRank() == 0 )
        {
            for( Int j=0; j<A.Width(); ++j ) 
                for( Int i=0; i<A.Height(); ++i )
                    sendBuf[i+j*A.Height()] = A.GetLocal(i,j);
            const Int recvViewingRank = this->Grid().VCToViewingMap(0);
            mpi::ISend
            ( sendBuf, A.Height()*A.Width(), recvViewingRank,
              this->Grid().ViewingComm(), sendRequest );
        }

        // Receive on the root of this matrix's grid and then broadcast
        // over this matrix's owning communicator
        if( this->Participating() )
        {
            if( this->Grid().VCRank() == 0 )
            {
                const Int sendViewingRank = A.Grid().VCToViewingMap(0);
                mpi::Recv
                ( bcastBuffer, A.Height()*A.Width(), sendViewingRank,
                  this->Grid().ViewingComm() );
            }

            mpi::Broadcast
            ( bcastBuffer, A.Height()*A.Width(), 0, this->Grid().VCComm() );

            for( Int j=0; j<A.Width(); ++j )
                for( Int i=0; i<A.Height(); ++i )
                    this->SetLocal(i,j,bcastBuffer[i+j*A.Height()]);
        }

        if( A.Grid().VCRank() == 0 )
            mpi::Wait( sendRequest );
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,* ] = [o ,o ]"))
    const Grid& g = A.Grid();
    const Int m = A.Height(); 
    const Int n = A.Width();
    this->Resize( A.Height(), A.Width() );

    if( this->Participating() )
    {
        const Int pkgSize = mpi::Pad( m*n );
        T* commBuffer = this->auxMemory_.Require( pkgSize );

        if( A.Participating() )
        {
            // Pack            
            const Int ALDim = A.LDim();
            const T* ABuf = A.LockedBuffer();
            for( Int j=0; j<n; ++j )
                for( Int i=0; i<m; ++i )
                    commBuffer[i+j*m] = ABuf[i+j*ALDim];
        }

        // Broadcast from the process that packed
        mpi::Broadcast( commBuffer, pkgSize, A.Root(), g.VCComm() );

        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<m; ++i )
                buffer[i+j*ldim] = commBuffer[i+j*m];        
    }

    return *this;
}

template<typename T>
DM&
DM::operator=( DM&& A )
{
    GDM::operator=( std::move(A) );
    return *this;
}

// Basic queries
// =============

template<typename T>
elem::DistData DM::DistData() const { return elem::DistData(*this); }

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
Int DM::ColStride() const { return 1; }
template<typename T>
Int DM::RowStride() const { return 1; }

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define PROTO(T) template class DistMatrix<T,STAR,STAR>
#define COPY(T,U,V) \
  template DistMatrix<T,STAR,STAR>::DistMatrix( const DistMatrix<T,U,V>& A )
#define FULL(T) \
  PROTO(T); \
  COPY(T,CIRC,CIRC); \
  COPY(T,MC,  MR  ); \
  COPY(T,MC,  STAR); \
  COPY(T,MD,  STAR); \
  COPY(T,MR,  MC  ); \
  COPY(T,MR,  STAR); \
  COPY(T,STAR,MC  ); \
  COPY(T,STAR,MD  ); \
  COPY(T,STAR,MR  ); \
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
