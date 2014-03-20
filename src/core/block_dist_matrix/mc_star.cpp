/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

namespace elem {

#define DM DistMatrix<T,MC,STAR>
#define BDM BlockDistMatrix<T,MC,STAR>
#define GBDM GeneralBlockDistMatrix<T,MC,STAR>

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
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR]::BlockDistMatrix"))
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct [MC,STAR] with itself");
}

template<typename T>
template<Dist U,Dist V>
BDM::BlockDistMatrix( const DistMatrix<T,U,V>& A )
: GBDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR]::BlockDistMatrix"))
    this->SetShifts();
    *this = A;
}

template<typename T>
template<Dist U,Dist V>
BDM::BlockDistMatrix( const BlockDistMatrix<T,U,V>& A )
: GBDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR]::BlockDistMatrix"))
    this->SetShifts();
    if( MC != U || STAR != V ||
        reinterpret_cast<const BDM*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct [MC,STAR] with itself");
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
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [MC,MR]"))
    A.RowAllGather( *this );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BDM& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [MC,STAR]"))
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [* ,MR]"))
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
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [MD,* ]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,MD>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [* ,MD]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MR,MC>& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[MC,STAR] = [MR,MC]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
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
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [MR,* ]"))
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
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [* ,MC]"))
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
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [VC,* ]"))
    A.PartialColAllGather( *this );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,VC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [* ,VC]"))
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
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [VR,* ]"))
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
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [* ,VR]"))
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
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR] = [* ,* ]"))
    this->ColFilterFrom( A );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,CIRC,CIRC>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[MC,STAR] = [o ,o ]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    LogicError("This routine is not yet written");
    return *this;
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

// Realignment
// -----------

template<typename T>
void
BDM::AlignWith( const elem::BlockDistData& data )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR]::AlignWith"))
    this->AlignColsWith( data );
}

template<typename T>
void
BDM::AlignColsWith( const elem::BlockDistData& data )
{
    DEBUG_ONLY(
        CallStackEntry cse("[MC,STAR]::AlignColsWith");
        // Consider the case where the row alignment is larger than that
        // permitted by the new grid
        if( *this->grid_ != *data.grid )
            LogicError("Grids do not match");
    )
    if( data.colDist == MC )
        this->AlignCols( data.blockHeight, data.colAlign, data.colCut );
    else if( data.rowDist == MC )
        this->AlignCols( data.blockWidth, data.rowAlign, data.rowCut );
    else if( data.colDist == VC )
        this->AlignCols
        ( data.blockHeight, data.colAlign % this->ColStride(), data.colCut );
    else if( data.rowDist == VC )
        this->AlignCols
        ( data.blockWidth, data.rowAlign % this->ColStride(), data.rowCut );
    DEBUG_ONLY(else LogicError("Nonsensical alignment"))
}

// Basic queries
// =============

template<typename T>
elem::BlockDistData BDM::DistData() const { return elem::BlockDistData(*this); }

template<typename T>
mpi::Comm BDM::DistComm() const { return this->grid_->MCComm(); }
template<typename T>
mpi::Comm BDM::CrossComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::RedundantComm() const { return this->grid_->MRComm(); }
template<typename T>
mpi::Comm BDM::ColComm() const { return this->grid_->MCComm(); }
template<typename T>
mpi::Comm BDM::RowComm() const { return mpi::COMM_SELF; }

template<typename T>
Int BDM::ColStride() const { return this->grid_->MCSize(); }
template<typename T>
Int BDM::RowStride() const { return 1; }

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define PROTO(T) template class BlockDistMatrix<T,MC,STAR>
#define COPY(T,U,V) \
  template BlockDistMatrix<T,MC,STAR>::BlockDistMatrix\
  ( const BlockDistMatrix<T,U,V>& A );\
  template BlockDistMatrix<T,MC,STAR>::BlockDistMatrix\
  ( const DistMatrix<T,U,V>& A );
#define FULL(T) \
  PROTO(T); \
  COPY(T,CIRC,CIRC); \
  COPY(T,MC,  MR  ); \
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
