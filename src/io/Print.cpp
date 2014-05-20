/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include "El/io.hpp"

namespace El {

template<typename T>
void Print( const std::vector<T>& x, std::string title, std::ostream& os )
{
    DEBUG_ONLY(CallStackEntry cse("Print"))
    if( title != "" )
        os << title << std::endl;
    
    const Int length = x.size();
    for( Int i=0; i<length; ++i )
        os << x[i] << " ";
    os << std::endl;
}

template<typename T>
void Print( const Matrix<T>& A, std::string title, std::ostream& os )
{
    DEBUG_ONLY(CallStackEntry cse("Print"))
    if( title != "" )
        os << title << std::endl;
    
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int i=0; i<height; ++i )
    {
        for( Int j=0; j<width; ++j )
            os << A.Get(i,j) << " ";
        os << std::endl;
    }
    os << std::endl;
}

template<typename T,Dist U,Dist V>
void Print( const DistMatrix<T,U,V>& A, std::string title, std::ostream& os )
{
    DEBUG_ONLY(CallStackEntry cse("Print"))
    if( U == A.UGath && V == A.VGath )
    {
        if( A.CrossRank() == A.Root() && A.RedundantRank() == 0 )
            Print( A.LockedMatrix(), title, os );
    }
    else
    {
        DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
            Print( A_CIRC_CIRC.LockedMatrix(), title, os );
    }
}

template<typename T,Dist U,Dist V>
void Print
( const BlockDistMatrix<T,U,V>& A, std::string title, std::ostream& os )
{
    DEBUG_ONLY(CallStackEntry cse("Print"))
    if( U == A.UGath && V == A.VGath )
    {
        if( A.CrossRank() == A.Root() && A.RedundantRank() == 0 )
            Print( A.LockedMatrix(), title, os );
    }
    else
    {
        BlockDistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
            Print( A_CIRC_CIRC.LockedMatrix(), title, os );
    }
}

template<typename T>
void Print
( const AbstractDistMatrix<T>& A, std::string title, std::ostream& os )
{
    DEBUG_ONLY(CallStackEntry cse("Print"))
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = dynamic_cast<const DistMatrix<T,CDIST,RDIST>&>(A); \
      Print( ACast, title, os );
    #include "El/core/GuardAndPayload.h"
}

template<typename T>
void Print
( const AbstractBlockDistMatrix<T>& A, std::string title, std::ostream& os )
{
    DEBUG_ONLY(CallStackEntry cse("Print"))
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = dynamic_cast<const BlockDistMatrix<T,CDIST,RDIST>&>(A); \
      Print( ACast, title, os );
    #include "El/core/GuardAndPayload.h"
}

#define DISTPROTO(T,U,V) \
  template void Print \
  ( const DistMatrix<T,U,V>& A, std::string title, std::ostream& os ); \
  template void Print \
  ( const BlockDistMatrix<T,U,V>& A, std::string title, std::ostream& os );

#define PROTO(T) \
  template void Print \
  ( const std::vector<T>& x, std::string title, std::ostream& os ); \
  template void Print \
  ( const Matrix<T>& A, std::string title, std::ostream& os ); \
  DISTPROTO(T,CIRC,CIRC); \
  DISTPROTO(T,MC,  MR  ); \
  DISTPROTO(T,MC,  STAR); \
  DISTPROTO(T,MD,  STAR); \
  DISTPROTO(T,MR,  MC  ); \
  DISTPROTO(T,MR,  STAR); \
  DISTPROTO(T,STAR,MC  ); \
  DISTPROTO(T,STAR,MD  ); \
  DISTPROTO(T,STAR,MR  ); \
  DISTPROTO(T,STAR,STAR); \
  DISTPROTO(T,STAR,VC  ); \
  DISTPROTO(T,STAR,VR  ); \
  DISTPROTO(T,VC,  STAR); \
  DISTPROTO(T,VR,  STAR); \
  template void Print \
  ( const AbstractDistMatrix<T>& A, std::string title, std::ostream& os ); \
  template void Print \
  ( const AbstractBlockDistMatrix<T>& A, std::string title, std::ostream& os ); 

PROTO(Int);
#ifndef EL_DISABLE_FLOAT
PROTO(float);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<float>);
#endif // ifndef EL_DISABLE_COMPLEX
#endif // ifndef EL_DISABLE_FLOAT

PROTO(double);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<double>);
#endif // ifndef EL_DISABLE_COMPLEX

} // namespace El
