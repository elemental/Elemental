/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T,Dist U,Dist V>
std::shared_ptr<const DistMatrix<T,U,V>> ReadProxy
( const AbstractDistMatrix<T>* A, const ProxyCtrl& ctrl )
{
    typedef DistMatrix<T,U,V> DM;
    const DM* ACast = dynamic_cast<const DM*>(A);

    const bool haveDist = (ACast != nullptr);
    const bool haveColAlign = haveDist && 
        (!ctrl.colConstrain || A->ColAlign() == ctrl.colAlign);
    const bool haveRowAlign = haveDist &&
        (!ctrl.rowConstrain || A->RowAlign() == ctrl.rowAlign);
    const bool haveRoot = haveDist &&
        (!ctrl.rootConstrain || A->Root() == ctrl.root);

    if( haveColAlign && haveRowAlign && haveRoot )
    {
        return std::shared_ptr<const DM>( ACast, []( const DM* B ) { } );
    }
    else
    {
        auto AShared = std::make_shared<DM>( A->Grid() );
        if( ctrl.rootConstrain )
            AShared->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            AShared->AlignCols( ctrl.colAlign );
        if( ctrl.rowConstrain )
            AShared->AlignRows( ctrl.rowAlign );
        Copy( *A, *AShared );
        return AShared;
    }
}

template<typename Real>
std::shared_ptr<const Matrix<Complex<Real>>> ComplexReadProxy
( const Matrix<Real>* A )
{
    typedef Matrix<Complex<Real>> MCpx;
    auto AShared = std::make_shared<MCpx>();
    Copy( *A, *AShared );
    return AShared;
}

template<typename Real>
std::shared_ptr<const Matrix<Complex<Real>>> ComplexReadProxy
( const Matrix<Complex<Real>>* A )
{
    typedef Matrix<Complex<Real>> M;
    return std::shared_ptr<const M>( A, []( const M* B ) { } );
}

template<typename Real,Dist U,Dist V>
std::shared_ptr<const DistMatrix<Complex<Real>,U,V>> ComplexReadProxy
( const AbstractDistMatrix<Real>* A, const ProxyCtrl& ctrl )
{
    typedef DistMatrix<Complex<Real>,U,V> DMCpx;
    auto AShared = std::make_shared<DMCpx>( A->Grid() );
    if( ctrl.rootConstrain )
        AShared->SetRoot( ctrl.root );
    if( ctrl.colConstrain )
        AShared->AlignCols( ctrl.colAlign );
    if( ctrl.rowConstrain )
        AShared->AlignRows( ctrl.rowAlign );
    Copy( *A, *AShared );
    return AShared;
}

template<typename Real,Dist U,Dist V>
std::shared_ptr<const DistMatrix<Complex<Real>,U,V>> ComplexReadProxy
( const AbstractDistMatrix<Complex<Real>>* A, const ProxyCtrl& ctrl )
{ return ReadProxy<Complex<Real>,U,V>( A, ctrl ); }

template<typename T,Dist U,Dist V>
std::shared_ptr<DistMatrix<T,U,V>> ReadProxy
( AbstractDistMatrix<T>* A, const ProxyCtrl& ctrl )
{
    typedef DistMatrix<T,U,V> DM;
    DM* ACast = dynamic_cast<DM*>(A);

    const bool haveDist = (ACast != nullptr);
    const bool haveColAlign = haveDist && 
        (!ctrl.colConstrain || A->ColAlign() == ctrl.colAlign);
    const bool haveRowAlign = haveDist &&
        (!ctrl.rowConstrain || A->RowAlign() == ctrl.rowAlign);
    const bool haveRoot = haveDist &&
        (!ctrl.rootConstrain || A->Root() == ctrl.root);

    if( haveColAlign && haveRowAlign && haveRoot )
    {
        return std::shared_ptr<DM>( ACast, []( const DM* B ) { } );
    }
    else
    {
        auto AShared = std::make_shared<DM>( A->Grid() );
        if( ctrl.rootConstrain )
            AShared->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            AShared->AlignCols( ctrl.colAlign );
        if( ctrl.rowConstrain )
            AShared->AlignRows( ctrl.rowAlign );
        Copy( *A, *AShared );
        return AShared;
    }
}

template<typename Real>
std::shared_ptr<Matrix<Complex<Real>>> ComplexReadProxy( Matrix<Real>* A )
{
    typedef Matrix<Complex<Real>> MCpx;
    auto AShared = std::make_shared<MCpx>();
    Copy( *A, *AShared );
    return AShared;
}

template<typename Real>
std::shared_ptr<Matrix<Complex<Real>>> ComplexReadProxy
( Matrix<Complex<Real>>* A )
{
    typedef Matrix<Complex<Real>> M;
    return std::shared_ptr<M>( A, []( const M* B ) { } );
}

template<typename Real,Dist U,Dist V>
std::shared_ptr<DistMatrix<Complex<Real>,U,V>> ComplexReadProxy
( AbstractDistMatrix<Real>* A, const ProxyCtrl& ctrl )
{
    typedef DistMatrix<Complex<Real>,U,V> DMCpx;
    auto AShared = std::make_shared<DMCpx>( A->Grid() );
    if( ctrl.rootConstrain )
        AShared->SetRoot( ctrl.root );
    if( ctrl.colConstrain )
        AShared->AlignCols( ctrl.colAlign );
    if( ctrl.rowConstrain )
        AShared->AlignRows( ctrl.rowAlign );
    Copy( *A, *AShared );
    return AShared;
}

template<typename Real,Dist U,Dist V>
std::shared_ptr<DistMatrix<Complex<Real>,U,V>> ComplexReadProxy
( AbstractDistMatrix<Complex<Real>>* A, const ProxyCtrl& ctrl )
{ return ReadProxy<Complex<Real>,U,V>( A, ctrl ); }

template<typename T,Dist U,Dist V>
std::shared_ptr<DistMatrix<T,U,V>> ReadWriteProxy
( AbstractDistMatrix<T>* A, const ProxyCtrl& ctrl )
{
    typedef DistMatrix<T,U,V> DM;
    DM* ACast = dynamic_cast<DM*>(A);

    const bool haveDist = (ACast != nullptr);
    const bool haveColAlign = haveDist &&
        (!ctrl.colConstrain || A->ColAlign() == ctrl.colAlign);
    const bool haveRowAlign = haveDist &&
        (!ctrl.rowConstrain || A->RowAlign() == ctrl.rowAlign);
    const bool haveRoot = haveDist &&
        (!ctrl.rootConstrain || A->Root() == ctrl.root);

    if( haveColAlign && haveRowAlign && haveRoot )
    {
        return std::shared_ptr<DM>( ACast, []( const DM* B ) { } );
    }
    else
    {
        DM* ARaw = new DM( A->Grid() );

        try
        {
            if( ctrl.rootConstrain )
                ARaw->SetRoot( ctrl.root );
            if( ctrl.colConstrain )
                ARaw->AlignCols( ctrl.colAlign );
            if( ctrl.rowConstrain )
                ARaw->AlignRows( ctrl.rowAlign );
            Copy( *A, *ARaw );
        }
        catch( std::exception& e )
        {
            delete ARaw;
            throw e;
        }

        return std::shared_ptr<DM>
               ( ARaw, [=]( const DM* B ) { Copy( *B, *A ); delete B; } );
    }
}

template<typename T,Dist U,Dist V>
std::shared_ptr<DistMatrix<T,U,V>> WriteProxy
( AbstractDistMatrix<T>* A, const ProxyCtrl& ctrl )
{
    typedef DistMatrix<T,U,V> DM;
    DM* ACast = dynamic_cast<DM*>(A);

    const bool haveDist = (ACast != nullptr);
    const bool haveColAlign = haveDist &&
        (!ctrl.colConstrain || A->ColAlign() == ctrl.colAlign);
    const bool haveRowAlign = haveDist &&
        (!ctrl.rowConstrain || A->RowAlign() == ctrl.rowAlign);
    const bool haveRoot = haveDist &&
        (!ctrl.rootConstrain || A->Root() == ctrl.root);

    if( haveColAlign && haveRowAlign && haveRoot )
    {
        return std::shared_ptr<DM>( ACast, []( const DM* B ) { } );
    }
    else
    {
        DM* ARaw = new DM( A->Grid() );

        try
        {
            if( ctrl.rootConstrain )
                ARaw->SetRoot( ctrl.root );
            if( ctrl.colConstrain )
                ARaw->AlignCols( ctrl.colAlign );
            if( ctrl.rowConstrain )
                ARaw->AlignRows( ctrl.rowAlign );
            ARaw->Resize( A->Height(), A->Width() );
        }
        catch( std::exception& e )
        {
            delete ARaw;
            throw e;
        }

        return std::shared_ptr<DM>
               ( ARaw, [=]( const DM* B ) { Copy( *B, *A ); delete B; } );
    }
}

#define DIST_PROTO_BASE(T,U,V) \
  template std::shared_ptr<const DistMatrix<T,U,V>> \
  ReadProxy( const AbstractDistMatrix<T>* A, const ProxyCtrl& ctrl ); \
  template std::shared_ptr<DistMatrix<T,U,V>> \
  ReadProxy( AbstractDistMatrix<T>* A, const ProxyCtrl& ctrl ); \
  template std::shared_ptr<DistMatrix<T,U,V>> \
  ReadWriteProxy( AbstractDistMatrix<T>* A, const ProxyCtrl& ctrl ); \
  template std::shared_ptr<DistMatrix<T,U,V>> \
  WriteProxy( AbstractDistMatrix<T>* A, const ProxyCtrl& ctrl );

#define DIST_PROTO_FIELD(F,U,V) \
  DIST_PROTO_BASE(F,U,V) \
  template std::shared_ptr<const DistMatrix<Complex<Base<F>>,U,V>> \
  ComplexReadProxy( const AbstractDistMatrix<F>* A, const ProxyCtrl& ctrl ); \
  template std::shared_ptr<DistMatrix<Complex<Base<F>>,U,V>> \
  ComplexReadProxy( AbstractDistMatrix<F>* A, const ProxyCtrl& ctrl );

#define PROTO_BASE(T) \
  DIST_PROTO_BASE(T,CIRC,CIRC); \
  DIST_PROTO_BASE(T,MC,  MR  ); \
  DIST_PROTO_BASE(T,MC,  STAR); \
  DIST_PROTO_BASE(T,MD,  STAR); \
  DIST_PROTO_BASE(T,MR,  MC  ); \
  DIST_PROTO_BASE(T,MR,  STAR); \
  DIST_PROTO_BASE(T,STAR,MC  ); \
  DIST_PROTO_BASE(T,STAR,MD  ); \
  DIST_PROTO_BASE(T,STAR,MR  ); \
  DIST_PROTO_BASE(T,STAR,STAR); \
  DIST_PROTO_BASE(T,STAR,VC  ); \
  DIST_PROTO_BASE(T,STAR,VR  ); \
  DIST_PROTO_BASE(T,VC  ,STAR); \
  DIST_PROTO_BASE(T,VR  ,STAR);

#define PROTO_FIELD(F) \
  template std::shared_ptr<const Matrix<Complex<Base<F>>>> \
  ComplexReadProxy( const Matrix<F>* A ); \
  template std::shared_ptr<Matrix<Complex<Base<F>>>> \
  ComplexReadProxy( Matrix<F>* A ); \
  DIST_PROTO_FIELD(F,CIRC,CIRC); \
  DIST_PROTO_FIELD(F,MC,  MR  ); \
  DIST_PROTO_FIELD(F,MC,  STAR); \
  DIST_PROTO_FIELD(F,MD,  STAR); \
  DIST_PROTO_FIELD(F,MR,  MC  ); \
  DIST_PROTO_FIELD(F,MR,  STAR); \
  DIST_PROTO_FIELD(F,STAR,MC  ); \
  DIST_PROTO_FIELD(F,STAR,MD  ); \
  DIST_PROTO_FIELD(F,STAR,MR  ); \
  DIST_PROTO_FIELD(F,STAR,STAR); \
  DIST_PROTO_FIELD(F,STAR,VC  ); \
  DIST_PROTO_FIELD(F,STAR,VR  ); \
  DIST_PROTO_FIELD(F,VC  ,STAR); \
  DIST_PROTO_FIELD(F,VR  ,STAR);

#define PROTO_INT(T) PROTO_BASE(T)

#define PROTO(T) PROTO_FIELD(T)

#include "El/macros/Instantiate.h"

} // namespace El
