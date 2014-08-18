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
        auto APtr = std::make_shared<DM>( A->Grid());
        if( ctrl.rootConstrain )
            APtr->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            APtr->AlignCols( ctrl.colAlign );
        if( ctrl.rowConstrain )
            APtr->AlignRows( ctrl.rowAlign );
        Copy( *A, *APtr );
        return APtr;
    }
}

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
        auto APtr = std::make_shared<DM>( A->Grid());
        if( ctrl.rootConstrain )
            APtr->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            APtr->AlignCols( ctrl.colAlign );
        if( ctrl.rowConstrain )
            APtr->AlignRows( ctrl.rowAlign );
        Copy( *A, *APtr );
        return APtr;
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
        auto APtr = std::make_shared<DM>( A->Grid());
        if( ctrl.rootConstrain )
            APtr->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            APtr->AlignCols( ctrl.colAlign );
        if( ctrl.rowConstrain )
            APtr->AlignRows( ctrl.rowAlign );
        APtr->Resize( A->Height(), A->Width() );
        return APtr;
    }
}

template<typename T,Dist U,Dist V>
void RestoreReadWriteProxy
( std::shared_ptr<DistMatrix<T,U,V>>& AProx, AbstractDistMatrix<T>& A )
{
    typedef DistMatrix<T,U,V> DM;
    if( AProx.get() != dynamic_cast<DM*>(&A) )
        Copy( *AProx, A );
}

template<typename T,Dist U,Dist V>
void RestoreWriteProxy
( std::shared_ptr<DistMatrix<T,U,V>>& AProx, AbstractDistMatrix<T>& A )
{
    typedef DistMatrix<T,U,V> DM;
    if( AProx.get() != dynamic_cast<DM*>(&A) )
        Copy( *AProx, A );
}

#define DIST_PROTO(T,U,V) \
  template std::shared_ptr<const DistMatrix<T,U,V>> \
  ReadProxy( const AbstractDistMatrix<T>* A, const ProxyCtrl& ctrl ); \
  template std::shared_ptr<DistMatrix<T,U,V>> \
  ReadWriteProxy( AbstractDistMatrix<T>* A, const ProxyCtrl& ctrl ); \
  template std::shared_ptr<DistMatrix<T,U,V>> \
  WriteProxy( AbstractDistMatrix<T>* A, const ProxyCtrl& ctrl ); \
  template void RestoreReadWriteProxy \
  ( std::shared_ptr<DistMatrix<T,U,V>>& AProx, AbstractDistMatrix<T>& A ); \
  template void RestoreWriteProxy \
  ( std::shared_ptr<DistMatrix<T,U,V>>& AProx, AbstractDistMatrix<T>& A );

#define PROTO(T) \
  DIST_PROTO(T,CIRC,CIRC); \
  DIST_PROTO(T,MC,  MR  ); \
  DIST_PROTO(T,MC,  STAR); \
  DIST_PROTO(T,MD,  STAR); \
  DIST_PROTO(T,MR,  MC  ); \
  DIST_PROTO(T,MR,  STAR); \
  DIST_PROTO(T,STAR,MC  ); \
  DIST_PROTO(T,STAR,MD  ); \
  DIST_PROTO(T,STAR,MR  ); \
  DIST_PROTO(T,STAR,STAR); \
  DIST_PROTO(T,STAR,VC  ); \
  DIST_PROTO(T,STAR,VR  ); \
  DIST_PROTO(T,VC  ,STAR); \
  DIST_PROTO(T,VR  ,STAR);

#include "El/macros/Instantiate.h"

} // namespace El
