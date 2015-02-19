/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Read proxy
// ==========

// Sequential
// ----------

template<typename T,typename S>
shared_ptr<const Matrix<T>> ReadProxy( const Matrix<S>* A )
{
    if( std::is_same<S,T>::value )
    {
        auto ACast = reinterpret_cast<const Matrix<T>*>(A);
        return shared_ptr<const Matrix<T>>
               ( ACast, []( const Matrix<T>* B ) { } );
    }
    else
    {
        auto AShared = make_shared<Matrix<T>>();
        Copy( *A, *AShared );
        return AShared;
    }
}

template<typename T,typename S>
shared_ptr<Matrix<T>> ReadProxy( Matrix<S>* A )
{
    if( std::is_same<S,T>::value )
    {
        auto ACast = reinterpret_cast<Matrix<T>*>(A);
        return shared_ptr<Matrix<T>>
               ( ACast, []( const Matrix<T>* B ) { } );
    }
    else 
    {
        auto AShared = make_shared<Matrix<T>>();
        Copy( *A, *AShared );
        return AShared;
    }
}

// Distributed
// -----------

template<typename T,Dist U,Dist V,typename S>
shared_ptr<const DistMatrix<T,U,V>> 
ReadProxy( const AbstractDistMatrix<S>* A, const ProxyCtrl& ctrl )
{
    typedef DistMatrix<T,U,V> DM;
    if( std::is_same<S,T>::value )
    {
        const DM* ACast = dynamic_cast<const DM*>(A);

        const bool haveDist = (ACast != nullptr);
        const bool haveColAlign = haveDist && 
            (!ctrl.colConstrain || A->ColAlign() == ctrl.colAlign);
        const bool haveRowAlign = haveDist &&
            (!ctrl.rowConstrain || A->RowAlign() == ctrl.rowAlign);
        const bool haveRoot = haveDist &&
            (!ctrl.rootConstrain || A->Root() == ctrl.root);

        if( haveColAlign && haveRowAlign && haveRoot )
            return shared_ptr<const DM>( ACast, []( const DM* B ) { } );
    }

    auto AShared = make_shared<DM>( A->Grid() );
    if( ctrl.rootConstrain )
        AShared->SetRoot( ctrl.root );
    if( ctrl.colConstrain )
        AShared->AlignCols( ctrl.colAlign );
    if( ctrl.rowConstrain )
        AShared->AlignRows( ctrl.rowAlign );
    Copy( *A, *AShared );
    return AShared;
}

template<typename T,Dist U,Dist V,typename S>
shared_ptr<DistMatrix<T,U,V>> 
ReadProxy( AbstractDistMatrix<S>* A, const ProxyCtrl& ctrl )
{
    typedef DistMatrix<T,U,V> DM;
    if( std::is_same<S,T>::value )
    {
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
            // Constrain the proxy to have the forced alignemnts.
            // This is somewhat tricky since a subsequent write could otherwise
            // change the alignment.
            if( ctrl.colConstrain )
                A->AlignCols( ctrl.colAlign );
            if( ctrl.rowConstrain )
                A->AlignRows( ctrl.rowAlign );
            if( ctrl.rootConstrain )
                A->SetRoot( ctrl.root );
            return shared_ptr<DM>( ACast, []( const DM* B ) { } );
        }
    }

    auto AShared = make_shared<DM>( A->Grid() );
    if( ctrl.rootConstrain )
        AShared->SetRoot( ctrl.root );
    if( ctrl.colConstrain )
        AShared->AlignCols( ctrl.colAlign );
    if( ctrl.rowConstrain )
        AShared->AlignRows( ctrl.rowAlign );
    Copy( *A, *AShared );
    return AShared;
}

// Read-write proxy
// ================

// Sequential
// ----------

template<typename T,typename S>
shared_ptr<Matrix<T>> ReadWriteProxy( Matrix<S>* A )
{
    typedef Matrix<T> M;
    if( std::is_same<S,T>::value )
    {
        auto ACast = reinterpret_cast<Matrix<T>*>(A);
        return shared_ptr<Matrix<T>>
               ( ACast, []( const Matrix<T>* B ) { } );
    }
    else
    {
        auto B = new M;
        Copy( *A, *B );
        return shared_ptr<M>
               ( B, [=]( const M* C ) { Copy( *C, *A ); delete C; } );
    }
}

// Distributed
// -----------

template<typename T,Dist U,Dist V,typename S>
shared_ptr<DistMatrix<T,U,V>> 
ReadWriteProxy( AbstractDistMatrix<S>* A, const ProxyCtrl& ctrl )
{
    typedef DistMatrix<T,U,V> DM;
    if( std::is_same<S,T>::value )
    {
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
            // Constrain the proxy to have the forced alignemnts.
            // This is somewhat tricky since a subsequent write could otherwise
            // change the alignment.
            if( ctrl.colConstrain )
                A->AlignCols( ctrl.colAlign );
            if( ctrl.rowConstrain )
                A->AlignRows( ctrl.rowAlign );
            if( ctrl.rootConstrain )
                A->SetRoot( ctrl.root );
            return shared_ptr<DM>( ACast, []( const DM* B ) { } );
        }
    }

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

    return shared_ptr<DM>
           ( ARaw, [=]( const DM* B ) { Copy( *B, *A ); delete B; } );
}

// Write proxy
// ===========

// Sequential
// ----------

template<typename T,typename S>
shared_ptr<Matrix<T>> WriteProxy( Matrix<S>* A )
{
    typedef Matrix<T> M;
    if( std::is_same<S,T>::value )
    {
        auto ACast = reinterpret_cast<Matrix<T>*>(A);
        return shared_ptr<Matrix<T>>
               ( ACast, []( const Matrix<T>* B ) { } );
    }
    else
        return shared_ptr<M>
               ( new M(A->Height(),A->Width()), 
                 [=]( const M* B ) { Copy( *B, *A ); delete B; } );
}

// Distributed
// -----------

template<typename T,Dist U,Dist V,typename S>
shared_ptr<DistMatrix<T,U,V>> 
WriteProxy( AbstractDistMatrix<S>* A, const ProxyCtrl& ctrl )
{
    typedef DistMatrix<T,U,V> DM;
    if( std::is_same<S,T>::value )
    {
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
            // Constrain the proxy to have the forced alignemnts.
            // This is somewhat tricky since a subsequent write could otherwise
            // change the alignment.
            if( ctrl.colConstrain )
                A->AlignCols( ctrl.colAlign );
            if( ctrl.rowConstrain )
                A->AlignRows( ctrl.rowAlign );
            if( ctrl.rootConstrain )
                A->SetRoot( ctrl.root );
            return shared_ptr<DM>( ACast, []( const DM* B ) { } );
        }
    }

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

    return shared_ptr<DM>
           ( ARaw, [=]( const DM* B ) { Copy( *B, *A ); delete B; } );
}

// TODO: include guards so that certain datatypes can be properly disabled 

#define READ_CONVERT_DIST(S,T,U,V) \
  template shared_ptr<const DistMatrix<T,U,V>> \
  ReadProxy( const AbstractDistMatrix<S>* A, const ProxyCtrl& ctrl ); \
  template shared_ptr<DistMatrix<T,U,V>> \
  ReadProxy( AbstractDistMatrix<S>* A, const ProxyCtrl& ctrl );

#define READWRITE_CONVERT_DIST(S,T,U,V) \
  template shared_ptr<DistMatrix<T,U,V>> \
  ReadWriteProxy( AbstractDistMatrix<S>* A, const ProxyCtrl& ctrl );

#define WRITE_CONVERT_DIST(S,T,U,V) \
  template shared_ptr<DistMatrix<T,U,V>> \
  WriteProxy( AbstractDistMatrix<S>* A, const ProxyCtrl& ctrl );

#define CONVERT_DIST(S,T,U,V) \
  READ_CONVERT_DIST(S,T,U,V) \
  READWRITE_CONVERT_DIST(S,T,U,V) \
  WRITE_CONVERT_DIST(S,T,U,V)

#define READ_CONVERT(S,T) \
  template shared_ptr<const Matrix<T>> ReadProxy( const Matrix<S>* A ); \
  template shared_ptr<Matrix<T>> ReadProxy( Matrix<S>* A ); \
  READ_CONVERT_DIST(S,T,CIRC,CIRC) \
  READ_CONVERT_DIST(S,T,MC,  MR  ) \
  READ_CONVERT_DIST(S,T,MC,  STAR) \
  READ_CONVERT_DIST(S,T,MD,  STAR) \
  READ_CONVERT_DIST(S,T,MR,  MC  ) \
  READ_CONVERT_DIST(S,T,MR,  STAR) \
  READ_CONVERT_DIST(S,T,STAR,MC  ) \
  READ_CONVERT_DIST(S,T,STAR,MD  ) \
  READ_CONVERT_DIST(S,T,STAR,MR  ) \
  READ_CONVERT_DIST(S,T,STAR,STAR) \
  READ_CONVERT_DIST(S,T,STAR,VC  ) \
  READ_CONVERT_DIST(S,T,STAR,VR  ) \
  READ_CONVERT_DIST(S,T,VC,  STAR) \
  READ_CONVERT_DIST(S,T,VR,  STAR)

#define READWRITE_CONVERT(S,T) \
  template shared_ptr<Matrix<T>> ReadWriteProxy( Matrix<S>* A ); \
  READWRITE_CONVERT_DIST(S,T,CIRC,CIRC) \
  READWRITE_CONVERT_DIST(S,T,MC,  MR  ) \
  READWRITE_CONVERT_DIST(S,T,MC,  STAR) \
  READWRITE_CONVERT_DIST(S,T,MD,  STAR) \
  READWRITE_CONVERT_DIST(S,T,MR,  MC  ) \
  READWRITE_CONVERT_DIST(S,T,MR,  STAR) \
  READWRITE_CONVERT_DIST(S,T,STAR,MC  ) \
  READWRITE_CONVERT_DIST(S,T,STAR,MD  ) \
  READWRITE_CONVERT_DIST(S,T,STAR,MR  ) \
  READWRITE_CONVERT_DIST(S,T,STAR,STAR) \
  READWRITE_CONVERT_DIST(S,T,STAR,VC  ) \
  READWRITE_CONVERT_DIST(S,T,STAR,VR  ) \
  READWRITE_CONVERT_DIST(S,T,VC,  STAR) \
  READWRITE_CONVERT_DIST(S,T,VR,  STAR)

#define WRITE_CONVERT(S,T) \
  template shared_ptr<Matrix<T>> WriteProxy( Matrix<S>* A ); \
  WRITE_CONVERT_DIST(S,T,CIRC,CIRC) \
  WRITE_CONVERT_DIST(S,T,MC,  MR  ) \
  WRITE_CONVERT_DIST(S,T,MC,  STAR) \
  WRITE_CONVERT_DIST(S,T,MD,  STAR) \
  WRITE_CONVERT_DIST(S,T,MR,  MC  ) \
  WRITE_CONVERT_DIST(S,T,MR,  STAR) \
  WRITE_CONVERT_DIST(S,T,STAR,MC  ) \
  WRITE_CONVERT_DIST(S,T,STAR,MD  ) \
  WRITE_CONVERT_DIST(S,T,STAR,MR  ) \
  WRITE_CONVERT_DIST(S,T,STAR,STAR) \
  WRITE_CONVERT_DIST(S,T,STAR,VC  ) \
  WRITE_CONVERT_DIST(S,T,STAR,VR  ) \
  WRITE_CONVERT_DIST(S,T,VC,  STAR) \
  WRITE_CONVERT_DIST(S,T,VR,  STAR)

#define CONVERT(S,T) \
  READ_CONVERT(S,T) \
  READWRITE_CONVERT(S,T) \
  WRITE_CONVERT(S,T)

#define PROTO_INT(T) CONVERT(Int,Int)

#define PROTO_REAL(Real) \
  CONVERT(Real,Real) \
  /* Promotions up to Real */ \
  READ_CONVERT(Int,Real) \
  /* Promotions up from Real */ \
  READ_CONVERT(Real,Complex<Real>)

#define PROTO_COMPLEX(C) \
  CONVERT(C,C) \
  /* Promotions up to C */ \
  READ_CONVERT(Int,C)

#define PROTO_FLOAT \
  PROTO_REAL(float) \
  /* Promotions up from float */ \
  CONVERT(float,double) \
  READ_CONVERT(float,Complex<double>)

#define PROTO_DOUBLE \
  PROTO_REAL(double) \
  /* Promotions down to float */ \
  CONVERT(double,float) \
  /* Mixed conversion */ \
  READ_CONVERT(double,Complex<float>)

#define PROTO_COMPLEX_FLOAT \
  PROTO_COMPLEX(Complex<float>) \
  /* Promotions up from Complex<float> */ \
  CONVERT(Complex<float>,Complex<double>)

#define PROTO_COMPLEX_DOUBLE \
  PROTO_COMPLEX(Complex<double>) \
  /* Promotions down from Complex<double> */ \
  CONVERT(Complex<double>,Complex<float>)

#include "El/macros/Instantiate.h"

} // namespace El
