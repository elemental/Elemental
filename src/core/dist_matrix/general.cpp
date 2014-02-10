/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {

template<typename T>
using ADM = AbstractDistMatrix<T>;

template<typename T,Dist U,Dist V>
using GDM = GeneralDistMatrix<T,U,V>;

// NOTE: It seems that member functions cannot be defined using a 
//       fully-specified template alias, e.g., ADM<T,U,V>::AbstractDistMatrix(),
//       but DM<T> is okay if it is only partially specified, e.g., 
//       DM<T> = DistMatrix<T,MC,MR> and DM<T>::DistMatrix()

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T,Dist U,Dist V>
GeneralDistMatrix<T,U,V>::GeneralDistMatrix( GDM<T,U,V>&& A )
: ADM<T>(std::move(A))
{ }

// Assignment and reconfiguration
// ==============================

template<typename T,Dist U,Dist V>
GDM<T,U,V>& 
GeneralDistMatrix<T,U,V>::operator=( GDM<T,U,V>&& A )
{
    AbstractDistMatrix<T>::operator=( std::move(A) );
    return *this;
}

// Diagonal manipulation
// =====================
template<typename T,Dist U,Dist V>
template<typename S>
bool
GeneralDistMatrix<T,U,V>::DiagonalAligned
( const DistMatrix<S,UDiag,VDiag>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::DiagonalAligned"))
    if( this->Grid() != d.Grid() )
        return false;

    if( U == MC && V == MR )
    {
        // The result is an [MD,* ]
        const Int diagRow = d.ColAlign()            % this->ColStride();
        const Int diagCol = (d.ColAlign()+d.Root()) % this->RowStride();
        if( offset >= 0 )
        {
            const Int procRow = this->ColAlign();
            const Int procCol = (this->RowAlign()+offset) % this->RowStride();
            return procRow==diagRow && procCol==diagCol;
        }
        else
        {
            const Int procRow = (this->ColAlign()-offset) % this->ColStride(); 
            const Int procCol = this->RowAlign();
            return procRow==diagRow && procCol==diagCol;
        }
    }
    else if( U == MR && V == MC )
    {
        // The result is an [MD,* ]
        const Int diagRow = d.ColAlign()            % this->ColStride();
        const Int diagCol = (d.ColAlign()+d.Root()) % this->RowStride();
        if( offset >= 0 )
        {
            const Int procCol = this->ColAlign();
            const Int procRow = (this->RowAlign()+offset) % this->RowStride();
            return procRow==diagRow && procCol==diagCol;
        }
        else
        {
            const Int procCol = (this->ColAlign()-offset) % this->ColStride();
            const Int procRow = this->RowAlign();
            return procRow==diagRow && procCol==diagCol;
        }
    }
    else if( U == STAR )
    {
        // The result is a [V,* ]
        if( offset >= 0 )
            return this->Root()==d.Root() && 
                   ((this->RowAlign()+offset)%this->RowStride())==d.ColAlign();
        else
            return this->Root()==d.Root() && this->RowAlign()==d.ColAlign();
    }
    else
    {
        // The result is a [U,V], where V is either STAR or CIRC
        if( offset >= 0 )
            return this->Root()==d.Root() && this->ColAlign() == d.ColAlign();
        else
            return this->Root()==d.Root() && 
                   ((this->ColAlign()-offset)%this->ColStride())==d.ColAlign();
    }
}

template<typename T,Dist U,Dist V>
template<typename S>
void
GeneralDistMatrix<T,U,V>::ForceDiagonalAlign
( DistMatrix<S,UDiag,VDiag>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::ForceDiagonalAlign"))
    const elem::Grid& grid = this->Grid();
    d.SetGrid( grid );

    if( U == MC && V == MR )
    {
        // Result is an [MD,* ]
        Int owner;
        if( offset >= 0 )
        {
            const Int procRow = this->ColAlign();
            const Int procCol = (this->RowAlign()+offset) % this->RowStride();
            owner = procRow + this->ColStride()*procCol;
        }
        else
        {
            const Int procRow = (this->ColAlign()-offset) % this->ColStride();
            const Int procCol = this->RowAlign();
            owner = procRow + this->ColStride()*procCol;
        }
        d.SetRoot( grid.DiagPath(owner) );
        d.AlignCols( grid.DiagPathRank(owner) );
    }
    else if( U == MR && V == MC )
    {
        // Result is an [MD,* ]
        Int owner;
        if( offset >= 0 )
        {
            const Int procCol = this->ColAlign();
            const Int procRow = (this->RowAlign()+offset) % this->RowStride();
            owner = procRow + this->ColStride()*procCol;
        }
        else
        {
            const Int procCol = (this->ColAlign()-offset) % this->ColStride();
            const Int procRow = this->RowAlign();
            owner = procRow + this->ColStride()*procCol;
        }
        d.SetRoot( grid.DiagPath(owner) );
        d.AlignCols( grid.DiagPathRank(owner) );
    }
    else if( U == STAR )
    {
        // Result is a [V,* ] 
        Int colAlign;
        if( offset >= 0 )
            colAlign = (this->RowAlign()+offset) % this->RowStride();
        else
            colAlign = this->RowAlign();
        d.SetRoot( this->Root() );
        d.AlignCols( colAlign );
    }
    else
    {
        // Result is a [U,V], where V is either STAR or CIRC
        Int colAlign;
        if( offset >= 0 )
            colAlign = this->ColAlign();
        else
            colAlign = (this->ColAlign()-offset) % this->ColStride();
        d.SetRoot( this->Root() );
        d.AlignCols( colAlign );
    }
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::GetDiagonal
( DistMatrix<T,UDiag,VDiag>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::GetDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::GetRealPartOfDiagonal
( DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::GetRealPartOfDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = RealPart(beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::GetImagPartOfDiagonal
( DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::GetImagPartOfDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = ImagPart(beta); } );
}

template<typename T,Dist U,Dist V>
auto
GeneralDistMatrix<T,U,V>::GetDiagonal( Int offset ) const
-> DistMatrix<T,UDiag,VDiag>
{
    DistMatrix<T,UDiag,VDiag> d( this->Grid() );
    GetDiagonal( d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
auto
GeneralDistMatrix<T,U,V>::GetRealPartOfDiagonal( Int offset ) const
-> DistMatrix<Base<T>,UDiag,VDiag>
{
    DistMatrix<Base<T>,UDiag,VDiag> d( this->Grid() );
    GetRealPartOfDiagonal( d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
auto
GeneralDistMatrix<T,U,V>::GetImagPartOfDiagonal( Int offset ) const
-> DistMatrix<Base<T>,UDiag,VDiag>
{
    DistMatrix<Base<T>,UDiag,VDiag> d( this->Grid() );
    GetImagPartOfDiagonal( d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::SetDiagonal
( const DistMatrix<T,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::SetDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::SetRealPartOfDiagonal
( const DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::SetRealPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { elem::SetRealPart(alpha,beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::SetImagPartOfDiagonal
( const DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::SetImagPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { elem::SetImagPart(alpha,beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::UpdateDiagonal
( T gamma, const DistMatrix<T,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::UpdateDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, [gamma]( T& alpha, T beta ) { alpha += gamma*beta; } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::UpdateRealPartOfDiagonal
( Base<T> gamma, const DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::UpdateRealPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      [gamma]( T& alpha, Base<T> beta ) 
      { elem::UpdateRealPart(alpha,gamma*beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::UpdateImagPartOfDiagonal
( Base<T> gamma, const DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::UpdateImagPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      [gamma]( T& alpha, Base<T> beta ) 
      { elem::UpdateImagPart(alpha,gamma*beta); } );
}

// Private section
// ###############

// Construct using a particular process grid
// =========================================

template<typename T,Dist U,Dist V>
GeneralDistMatrix<T,U,V>::GeneralDistMatrix( const elem::Grid& grid )
: ADM<T>(grid)
{ }

// Helper functions
// ================
template<typename T,Dist U,Dist V>
template<typename S,class Function>
void
GeneralDistMatrix<T,U,V>::GetDiagonalHelper
( DistMatrix<S,UDiag,VDiag>& d, Int offset, Function func ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::GetDiagonalHelper"))
    d.SetGrid( this->Grid() );
    this->ForceDiagonalAlign( d, offset );
    d.Resize( this->DiagonalLength(offset), 1 );
    if( !d.Participating() )
        return;

    const Int diagShift = d.ColShift();
    const Int diagStride = d.ColStride();
    const Int iStart = ( offset>=0 ? diagShift        : diagShift-offset );
    const Int jStart = ( offset>=0 ? diagShift+offset : diagShift        );

    const Int colStride = this->ColStride();
    const Int rowStride = this->RowStride();
    const Int iLocStart = (iStart-this->ColShift()) / colStride;
    const Int jLocStart = (jStart-this->RowShift()) / rowStride;

    const Int localDiagLength = d.LocalHeight();
    S* dBuf = d.Buffer();
    const T* buffer = this->LockedBuffer();
    const Int ldim = this->LDim();

    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(diagStride/colStride);
        const Int jLoc = jLocStart + k*(diagStride/rowStride);
        func( dBuf[k], buffer[iLoc+jLoc*ldim] );
    }
}

template<typename T,Dist U,Dist V>
template<typename S,class Function>
void
GeneralDistMatrix<T,U,V>::SetDiagonalHelper
( const DistMatrix<S,UDiag,VDiag>& d, Int offset, Function func ) 
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::SetDiagonalHelper");
        if( !this->DiagonalAligned( d, offset ) )
            LogicError("Invalid diagonal alignment");
    )
    if( !d.Participating() )
        return;

    const Int diagShift = d.ColShift();
    const Int diagStride = d.ColStride();
    const Int iStart = ( offset>=0 ? diagShift        : diagShift-offset );
    const Int jStart = ( offset>=0 ? diagShift+offset : diagShift        );

    const Int colStride = this->ColStride();
    const Int rowStride = this->RowStride();
    const Int iLocStart = (iStart-this->ColShift()) / colStride;
    const Int jLocStart = (jStart-this->RowShift()) / rowStride;

    const Int localDiagLength = d.LocalHeight();
    const S* dBuf = d.LockedBuffer();
    T* buffer = this->Buffer();
    const Int ldim = this->LDim();

    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(diagStride/colStride);
        const Int jLoc = jLocStart + k*(diagStride/rowStride);
        func( buffer[iLoc+jLoc*ldim], dBuf[k] );
    }
}

// Instantiations for {Int,Real,Complex<Real>} for each Real in {float,double}
// ###########################################################################

#define DIAGALIGNED(S,T,U,V)\
  template bool GeneralDistMatrix<T,U,V>::DiagonalAligned\
  ( const DistMatrix<S,UDiag,VDiag>&, Int ) const;\
  template void GeneralDistMatrix<T,U,V>::ForceDiagonalAlign\
  ( DistMatrix<S,UDiag,VDiag>&, Int ) const;

#define DISTPROTO(T,U,V)\
  template class GeneralDistMatrix<T,U,V>;\
  DIAGALIGNED(Int,T,U,V);\
  DIAGALIGNED(float,T,U,V);\
  DIAGALIGNED(double,T,U,V);\
  DIAGALIGNED(Complex<float>,T,U,V);\
  DIAGALIGNED(Complex<double>,T,U,V);
  
#ifndef DISABLE_COMPLEX
#ifndef DISABLE_FLOAT

#define DISTPROTO(T,U,V)\
  template class GeneralDistMatrix<T,U,V>;\
  DIAGALIGNED(Int,T,U,V);\
  DIAGALIGNED(float,T,U,V);\
  DIAGALIGNED(double,T,U,V);\
  DIAGALIGNED(Complex<float>,T,U,V);\
  DIAGALIGNED(Complex<double>,T,U,V);
 
#else // ifndef DISABLE_FLOAT

#define DISTPROTO(T,U,V)\
  template class GeneralDistMatrix<T,U,V>;\
  DIAGALIGNED(Int,T,U,V);\
  DIAGALIGNED(double,T,U,V);\
  DIAGALIGNED(Complex<double>,T,U,V);

#endif // ifndef DISABLE_FLOAT
#else // ifndef DISABLE_COMPLEX
#ifndef DISABLE_FLOAT

#define DISTPROTO(T,U,V)\
  template class GeneralDistMatrix<T,U,V>;\
  DIAGALIGNED(Int,T,U,V);\
  DIAGALIGNED(float,T,U,V);\
  DIAGALIGNED(double,T,U,V);

#else // ifndef DISABLE_FLOAT

#define DISTPROTO(T,U,V)\
  template class GeneralDistMatrix<T,U,V>;\
  DIAGALIGNED(Int,T,U,V);\
  DIAGALIGNED(double,T,U,V);

#endif // ifndef DISABLE_FLOAT
#endif // ifndef DISABLE_COMPLEX

#define PROTO(T)\
  DISTPROTO(T,CIRC,CIRC);\
  DISTPROTO(T,MC,  MR  );\
  DISTPROTO(T,MC,  STAR);\
  DISTPROTO(T,MD,  STAR);\
  DISTPROTO(T,MR,  MC  );\
  DISTPROTO(T,MR,  STAR);\
  DISTPROTO(T,STAR,MC  );\
  DISTPROTO(T,STAR,MD  );\
  DISTPROTO(T,STAR,MR  );\
  DISTPROTO(T,STAR,STAR);\
  DISTPROTO(T,STAR,VC  );\
  DISTPROTO(T,STAR,VR  );\
  DISTPROTO(T,VC,  STAR);\
  DISTPROTO(T,VR,  STAR);

#ifndef DISABLE_COMPLEX
 #ifndef DISABLE_FLOAT
  PROTO(Int);
  PROTO(float);
  PROTO(double);
  PROTO(Complex<float>);
  PROTO(Complex<double>);
 #else // ifndef DISABLE_FLOAT
  PROTO(Int);
  PROTO(double);
  PROTO(Complex<double>);
 #endif // ifndef DISABLE_FLOAT
#else // ifndef DISABLE_COMPLEX
 #ifndef DISABLE_FLOAT
  PROTO(Int);
  PROTO(float);
  PROTO(double);
 #else // ifndef DISABLE_FLOAT
  PROTO(Int);
  PROTO(double);
 #endif // ifndef DISABLE_FLOAT
#endif // ifndef DISABLE_COMPLEX

} // namespace elem
