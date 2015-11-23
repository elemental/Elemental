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
shared_ptr<const Matrix<T>> ReadProxy( const Matrix<S>* A )
{
    if( IsSame<S,T>::value )
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

// Distributed
// ===========

template<typename T,Dist U,Dist V,typename S>
shared_ptr<const DistMatrix<T,U,V>> 
ReadProxy( const ElementalMatrix<S>* A, const ProxyCtrl& ctrl )
{
    typedef DistMatrix<T,U,V> DM;
    if( IsSame<S,T>::value )
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

#define CONVERT_DIST(S,T,U,V) \
  template shared_ptr<const DistMatrix<T,U,V>> \
  ReadProxy( const ElementalMatrix<S>* A, const ProxyCtrl& ctrl );

#define CONVERT(S,T) \
  template shared_ptr<const Matrix<T>> ReadProxy( const Matrix<S>* A ); \
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

#define PROTO_REAL(Real) \
  CONVERT(Real,Real) \
  CONVERT(Int,Real) \
  CONVERT(Real,Complex<Real>)

#define PROTO_COMPLEX(C) \
  CONVERT(C,C) \
  CONVERT(Int,C)

#ifdef EL_HAVE_QUAD

#define PROTO_FLOAT \
  PROTO_REAL(float) \
  CONVERT(float,double) \
  CONVERT(float,Quad) \
  CONVERT(float,Complex<double>) \
  CONVERT(float,Complex<Quad>)

#define PROTO_DOUBLE \
  PROTO_REAL(double) \
  CONVERT(double,float) \
  CONVERT(double,Quad) \
  CONVERT(double,Complex<float>) \
  CONVERT(double,Complex<Quad>)

#define PROTO_QUAD \
  PROTO_REAL(Quad) \
  CONVERT(Quad,float) \
  CONVERT(Quad,double) \
  CONVERT(Quad,Complex<float>) \
  CONVERT(Quad,Complex<double>)

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
  CONVERT(float,double) \
  CONVERT(float,Complex<double>)

#define PROTO_DOUBLE \
  PROTO_REAL(double) \
  CONVERT(double,float) \
  CONVERT(double,Complex<float>)

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
