/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Sequential
// ==========

template<typename T,typename S>
shared_ptr<Matrix<T>> ReadWriteProxy( Matrix<S>* A )
{
    typedef Matrix<T> M;
    if( IsSame<S,T>::value )
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
// ===========

template<typename T,Dist U,Dist V,typename S>
shared_ptr<DistMatrix<T,U,V>> 
ReadWriteProxy( ElementalMatrix<S>* A, const ProxyCtrl& ctrl )
{
    typedef DistMatrix<T,U,V> DM;
    if( IsSame<S,T>::value )
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

#define CONVERT_DIST(S,T,U,V) \
  template shared_ptr<DistMatrix<T,U,V>> \
  ReadWriteProxy( ElementalMatrix<S>* A, const ProxyCtrl& ctrl );

#define CONVERT(S,T) \
  template shared_ptr<Matrix<T>> ReadWriteProxy( Matrix<S>* A ); \
  CONVERT_DIST(S,T,CIRC,CIRC) \
  CONVERT_DIST(S,T,MC,  MR  ) \
  CONVERT_DIST(S,T,MC,  STAR) \
  CONVERT_DIST(S,T,MD,  STAR) \
  CONVERT_DIST(S,T,MR,  MC  ) \
  CONVERT_DIST(S,T,MR,  STAR) \
  CONVERT_DIST(S,T,STAR,MC  ) \
  CONVERT_DIST(S,T,STAR,MD  ) \
  CONVERT_DIST(S,T,STAR,MR  ) \
  CONVERT_DIST(S,T,STAR,STAR) \
  CONVERT_DIST(S,T,STAR,VC  ) \
  CONVERT_DIST(S,T,STAR,VR  ) \
  CONVERT_DIST(S,T,VC,  STAR) \
  CONVERT_DIST(S,T,VR,  STAR)

#define PROTO_INT(T) CONVERT(Int,Int)

#define PROTO_REAL(Real) CONVERT(Real,Real)

#define PROTO_COMPLEX(C) CONVERT(C,C)

#ifdef EL_HAVE_QUAD

#define PROTO_FLOAT \
  PROTO_REAL(float) \
  CONVERT(float,double) \
  CONVERT(float,Quad)

#define PROTO_DOUBLE \
  PROTO_REAL(double) \
  CONVERT(double,float) \
  CONVERT(double,Quad)

#define PROTO_QUAD \
  PROTO_REAL(Quad) \
  CONVERT(Quad,float) \
  CONVERT(Quad,double)

#define PROTO_COMPLEX_FLOAT \
  PROTO_COMPLEX(Complex<float>) \
  CONVERT(Complex<float>,Complex<double>) \
  CONVERT(Complex<float>,Complex<Quad>)

#define PROTO_COMPLEX_DOUBLE \
  PROTO_COMPLEX(Complex<double>) \
  CONVERT(Complex<double>,Complex<float>) \
  CONVERT(Complex<double>,Complex<Quad>) 

#define PROTO_COMPLEX_QUAD \
  PROTO_COMPLEX(Complex<Quad>) \
  CONVERT(Complex<Quad>,Complex<float>) \
  CONVERT(Complex<Quad>,Complex<double>) 

#else

#define PROTO_FLOAT \
  PROTO_REAL(float) \
  CONVERT(float,double)

#define PROTO_DOUBLE \
  PROTO_REAL(double) \
  CONVERT(double,float)

#define PROTO_COMPLEX_FLOAT \
  PROTO_COMPLEX(Complex<float>) \
  CONVERT(Complex<float>,Complex<double>)

#define PROTO_COMPLEX_DOUBLE \
  PROTO_COMPLEX(Complex<double>) \
  CONVERT(Complex<double>,Complex<float>)

#endif

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
