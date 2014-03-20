/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

namespace elem {

#define DM DistMatrix<T,MC,MR>
#define BDM BlockDistMatrix<T,MC,MR>
#define GBDM GeneralBlockDistMatrix<T,MC,MR>

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
    DEBUG_ONLY(CallStackEntry cse("[MC,MR]::BlockDistMatrix"))
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct [MC,MR] with itself");
}

template<typename T>
template<Dist U,Dist V>
BDM::BlockDistMatrix( const DistMatrix<T,U,V>& A )
: GBDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[MC,MR]::BlockDistMatrix"))
    this->SetShifts();
    *this = A;
}

template<typename T>
template<Dist U,Dist V>
BDM::BlockDistMatrix( const BlockDistMatrix<T,U,V>& A )
: GBDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[MC,MR]::BlockDistMatrix"))
    this->SetShifts();
    if( MC != U || MR != V ||
        reinterpret_cast<const BDM*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct [MC,MR] with itself");
}

template<typename T>
BDM::BlockDistMatrix( BDM&& A ) noexcept : GBDM(std::move(A)) { }

template<typename T> BDM::~BlockDistMatrix() { }

// Assignment and reconfiguration
// ==============================

template<typename T>
BDM&
BDM::operator=( const BDM& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[MC,MR] = [MC,MR]");
        this->AssertNotLocked();
    )
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MC,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [MC,* ]"))
    this->RowFilterFrom( A );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [* ,MR]"))
    this->ColFilterFrom( A );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [MD,* ]"))
    // TODO: More efficient implementation?
    BlockDistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,MD>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [* ,MD]"))
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
        CallStackEntry cse("[MC,MR] = [MR,MC]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [MR,* ]"))
    const Grid& g = A.Grid();
    std::unique_ptr<BlockDistMatrix<T,VR,STAR>> A_VR_STAR
    ( new BlockDistMatrix<T,VR,STAR>(A) );
    std::unique_ptr<BlockDistMatrix<T,VC,STAR>> A_VC_STAR
    ( new BlockDistMatrix<T,VC,STAR>(g) );
    A_VC_STAR->AlignWith( *this );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater
    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [* ,MC]"))
    const Grid& g = A.Grid();
    std::unique_ptr<BlockDistMatrix<T,STAR,VC>> A_STAR_VC
    ( new BlockDistMatrix<T,STAR,VC>(A) );
    std::unique_ptr<BlockDistMatrix<T,STAR,VR>> A_STAR_VR
    ( new BlockDistMatrix<T,STAR,VR>(g) );
    A_STAR_VR->AlignWith( *this );
    *A_STAR_VR = *A_STAR_VC;
    delete A_STAR_VC.release(); // lowers memory highwater
    *this = *A_STAR_VR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,VC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [VC,* ]"))
    A.PartialColAllToAll( *this );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,VC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [* ,VC]"))
    const elem::Grid& g = this->Grid();
    BlockDistMatrix<T,STAR,VR> A_STAR_VR(g);
    A_STAR_VR.AlignWith( *this );
    A_STAR_VR = A;
    *this = A_STAR_VR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,VR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [VR,* ]"))
    const elem::Grid& g = this->Grid();
    BlockDistMatrix<T,VC,STAR> A_VC_STAR(g);
    A_VC_STAR.AlignWith( *this );
    A_VC_STAR = A;
    *this = A_VC_STAR;
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,VR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [* ,VR]"))
    A.PartialRowAllToAll( *this );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [* ,* ]"))
    this->FilterFrom( A );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,CIRC,CIRC>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[MC,MR] = [o ,o ]");
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
    DEBUG_ONLY(CallStackEntry cse("[MC,MR]::AlignWith"))
    this->SetGrid( *data.grid );
    if( data.colDist == MC && data.rowDist == MR )
        this->Align
        ( data.blockHeight, data.blockWidth, 
          data.colAlign, data.rowAlign, data.colCut, data.rowCut );
    else if( data.colDist == MC && data.rowDist == STAR )
        this->AlignCols( data.blockHeight, data.colAlign, data.colCut );
    else if( data.colDist == MR && data.rowDist == MC )
        this->Align
        ( data.blockWidth, data.blockHeight, 
          data.rowAlign, data.colAlign, data.rowCut, data.colCut );
    else if( data.colDist == MR && data.rowDist == STAR )
        this->AlignRows( data.blockHeight, data.colAlign, data.colCut );
    else if( data.colDist == STAR && data.rowDist == MC )
        this->AlignCols( data.blockWidth, data.rowAlign, data.rowCut );
    else if( data.colDist == STAR && data.rowDist == MR )
        this->AlignRows( data.blockWidth, data.rowAlign, data.rowCut );
    else if( data.colDist == STAR && data.rowDist == VC )
        this->AlignCols
        ( data.blockWidth, data.rowAlign % this->ColStride(), data.rowCut );
    else if( data.colDist == STAR && data.rowDist == VR )
        this->AlignRows
        ( data.blockWidth, data.rowAlign % this->RowStride(), data.rowCut );
    else if( data.colDist == VC && data.rowDist == STAR )
        this->AlignCols
        ( data.blockHeight, data.colAlign % this->ColStride(), data.colCut );
    else if( data.colDist == VR && data.rowDist == STAR )
        this->AlignRows
        ( data.blockHeight, data.colAlign % this->RowStride(), data.colCut );
    DEBUG_ONLY(else LogicError("Nonsensical alignment"))
}

template<typename T>
void
BDM::AlignColsWith( const elem::BlockDistData& data )
{
    DEBUG_ONLY(
        CallStackEntry cse("[MC,MR]::AlignColsWith");
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

template<typename T>
void
BDM::AlignRowsWith( const elem::BlockDistData& data )
{
    DEBUG_ONLY(
        CallStackEntry cse("[MC,MR]::AlignRowsWith");
        if( *this->grid_ != *data.grid )
            LogicError("Grids do not match");
    )
    if( data.colDist == MR )
        this->AlignRows( data.blockHeight, data.colAlign, data.colCut );
    else if( data.rowDist == MR )
        this->AlignRows( data.blockWidth, data.rowAlign, data.rowCut );
    else if( data.colDist == VR )
        this->AlignRows
        ( data.blockHeight, data.colAlign % this->RowStride(), data.colCut );
    else if( data.rowDist == VR )
        this->AlignRows
        ( data.blockWidth, data.rowAlign % this->RowStride(), data.rowCut );
    DEBUG_ONLY(else LogicError("Nonsensical alignment"))
}

// Basic queries
// =============

template<typename T>
elem::BlockDistData BDM::DistData() const { return elem::BlockDistData(*this); }

template<typename T>
mpi::Comm BDM::DistComm() const { return this->grid_->VCComm(); }
template<typename T>
mpi::Comm BDM::CrossComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::RedundantComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::ColComm() const { return this->grid_->MCComm(); }
template<typename T>
mpi::Comm BDM::RowComm() const { return this->grid_->MRComm(); }

template<typename T>
Int BDM::ColStride() const { return this->grid_->MCSize(); }
template<typename T>
Int BDM::RowStride() const { return this->grid_->MRSize(); }

// Private section
// ###############

// Redistribute from a different process grid
// ==========================================

template<typename T>
void BDM::CopyFromDifferentGrid( const BDM& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,MR]::CopyFromDifferentGrid"))
    LogicError("This routine is not yet written");
}

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define PROTO(T) template class BlockDistMatrix<T,MC,MR>
#define COPY(T,U,V) \
  template BlockDistMatrix<T,MC,MR>::BlockDistMatrix\
  ( const BlockDistMatrix<T,U,V>& A );\
  template BlockDistMatrix<T,MC,MR>::BlockDistMatrix\
  ( const DistMatrix<T,U,V>& A );
#define FULL(T) \
  PROTO(T); \
  COPY(T,CIRC,CIRC); \
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
