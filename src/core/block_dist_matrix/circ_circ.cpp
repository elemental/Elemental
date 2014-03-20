/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

namespace elem {

#define DM DistMatrix<T,CIRC,CIRC>
#define BDM BlockDistMatrix<T,CIRC,CIRC>
#define GBDM GeneralBlockDistMatrix<T,CIRC,CIRC>

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T>
BDM::BlockDistMatrix
( const elem::Grid& g, Int blockHeight, Int blockWidth, Int root )
: GBDM(g,blockHeight,blockWidth,root)
{ this->SetShifts(); }

template<typename T>
BDM::BlockDistMatrix
( Int height, Int width, const elem::Grid& g,
  Int blockHeight, Int blockWidth, Int root )
: GBDM(g,blockHeight,blockWidth,root)
{ this->SetShifts(); this->Resize(height,width); }

template<typename T>
BDM::BlockDistMatrix
( Int height, Int width, const elem::Grid& g, 
  Int blockHeight, Int blockWidth, 
  Int colAlign, Int rowAlign, Int colCut, Int rowCut, Int root )
: GBDM(g,blockHeight,blockWidth,root)
{ 
    this->SetShifts(); 
    this->Align(blockHeight,blockWidth,colAlign,rowAlign,colCut,rowCut); 
    this->Resize(height,width); 
}

template<typename T>
BDM::BlockDistMatrix
( Int height, Int width, const elem::Grid& g,
  Int blockHeight, Int blockWidth, 
  Int colAlign, Int rowAlign, Int colCut, Int rowCut, Int ldim, Int root )
: GBDM(g,blockHeight,blockWidth,root)
{ 
    this->SetShifts();
    this->Align(blockHeight,blockWidth,colAlign,rowAlign,colCut,rowCut); 
    this->Resize(height,width,ldim); 
}

template<typename T>
BDM::BlockDistMatrix
( Int height, Int width, const elem::Grid& g,
  Int blockHeight, Int blockWidth, 
  Int colAlign, Int rowAlign, Int colCut, Int rowCut,
  const T* buffer, Int ldim, Int root )
: GBDM(g,blockHeight,blockWidth,root)
{ 
    this->LockedAttach
    (height,width,g,blockHeight,blockWidth,colAlign,rowAlign,colCut,rowCut,
     buffer,ldim,root); 
}

template<typename T>
BDM::BlockDistMatrix
( Int height, Int width, const elem::Grid& g,
  Int blockHeight, Int blockWidth,
  Int colAlign, Int rowAlign, Int colCut, Int rowCut,
  T* buffer, Int ldim, Int root )
: GBDM(g,blockHeight,blockWidth,root)
{ 
    this->Attach
    (height,width,g,blockHeight,blockWidth,colAlign,rowAlign,colCut,rowCut,
     buffer,ldim,root); 
}

template<typename T>
BDM::BlockDistMatrix( const BDM& A )
: GBDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC]::BlockDistMatrix"))
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct [CIRC,CIRC] with itself");
}

template<typename T>
template<Dist U,Dist V>
BDM::BlockDistMatrix( const DistMatrix<T,U,V>& A )
: GBDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC]::BlockDistMatrix"))
    this->SetShifts();
    *this = A;
}

template<typename T>
template<Dist U,Dist V>
BDM::BlockDistMatrix( const BlockDistMatrix<T,U,V>& A )
: GBDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC]::BlockDistMatrix"))
    this->SetShifts();
    if( CIRC != U || CIRC != V ||
        reinterpret_cast<const BDM*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct [CIRC,CIRC] with itself");
}

template<typename T>
BDM::BlockDistMatrix( BDM&& A ) noexcept : GBDM(std::move(A)) { }

template<typename T> BDM::~BlockDistMatrix() { }

// Assignment and reconfiguration
// ==============================

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MC,MR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [MC,MR]"))
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MC,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [MC,STAR]"))
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [* ,MR]"))
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [MD,* ]"))
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,MD>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [* ,MD]"))
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [MR,MC]"))
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [MR,* ]"))
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [* ,MC]"))
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,VC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [VC,* ]"))
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,VC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [* ,VC]"))
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,VR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [VR,* ]"))
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,VR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [* ,VR]"))
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [STAR,STAR]"))
    this->Resize( A.Height(), A.Width() );
    if( A.Grid().VCRank() == this->Root() )
        this->matrix_ = A.LockedMatrix();
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BDM& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[CIRC,CIRC] = [CIRC,CIRC]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
void
BDM::CopyFromRoot( const Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC]::CopyFromRoot"))
    const Grid& grid = this->Grid();
    if( grid.VCRank() != this->Root() )
        LogicError("Called CopyFromRoot from non-root");

    Int dims[2];
    dims[0] = A.Height();
    dims[1] = A.Width();
    mpi::Broadcast( dims, 2, this->Root(), grid.VCComm() );

    this->Resize( dims[0], dims[1] );
    this->matrix_ = A;
}

template<typename T>
void
BDM::CopyFromNonRoot()
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC]::CopyFromNonRoot"))
    const Grid& grid = this->Grid();
    if( grid.VCRank() == this->Root() )
        LogicError("Called CopyFromNonRoot from root");

    Int dims[2];
    mpi::Broadcast( dims, 2, this->Root(), grid.VCComm() );

    this->Resize( dims[0], dims[1] );
}

template<typename T>
BDM&
BDM::operator=( BDM&& A )
{
    if( this->Viewing() && !A.Viewing() )
    {
        const BDM& AConst = A;
        this->operator=( AConst );
    }
    else
    {
        GBDM::operator=( std::move(A) );
    }
    return *this;
}

// Basic queries
// =============

template<typename T>
elem::BlockDistData BDM::DistData() const { return elem::BlockDistData(*this); }

template<typename T>
mpi::Comm BDM::DistComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::CrossComm() const { return this->grid_->VCComm(); }
template<typename T>
mpi::Comm BDM::RedundantComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::ColComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::RowComm() const { return mpi::COMM_SELF; }

template<typename T>
Int BDM::ColStride() const { return 1; }
template<typename T>
Int BDM::RowStride() const { return 1; }

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define PROTO(T) template class BlockDistMatrix<T,CIRC,CIRC>
#define COPY(T,U,V) \
  template BlockDistMatrix<T,CIRC,CIRC>::BlockDistMatrix\
  ( const BlockDistMatrix<T,U,V>& A );\
  template BlockDistMatrix<T,CIRC,CIRC>::BlockDistMatrix\
  ( const DistMatrix<T,U,V>& A );
#define FULL(T) \
  PROTO(T); \
  COPY(T,MC,  MR  ); \
  COPY(T,MC,  STAR); \
  COPY(T,MD,  STAR); \
  COPY(T,MR,  MC  ); \
  COPY(T,MR,  STAR); \
  COPY(T,STAR,MC  ); \
  COPY(T,STAR,MD  ); \
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
