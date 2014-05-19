/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include "El/io.hpp"

#include "El/io/SpyWindow.hpp"

#ifdef EL_HAVE_QT5
# include <QApplication>
#endif

namespace El {

template<typename T>
void Spy( const Matrix<T>& A, std::string title, Base<T> tol )
{
    DEBUG_ONLY(CallStackEntry cse("Spy"))
#ifdef EL_HAVE_QT5
    if( GuiDisabled() )
        LogicError("GUI was disabled");

    // Convert A to double-precision since Qt's MOC does not support templates
    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<Int>* ASpy = new Matrix<Int>( m, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            ASpy->Set( i, j, ( Abs(A.Get(i,j))>tol ? 1 : 0 ) );

    QString qTitle = QString::fromStdString( title );
    SpyWindow* spyWindow = new SpyWindow;
    spyWindow->Spy( ASpy, qTitle );
    spyWindow->show();

    // Spend at most 200 milliseconds rendering
    ProcessEvents( 200 );
#else
    LogicError("Qt5 not available");
#endif // ifdef EL_HAVE_QT5
}

template<typename T,Distribution U,Distribution V>
void Spy( const DistMatrix<T,U,V>& A, std::string title, Base<T> tol )
{
    DEBUG_ONLY(CallStackEntry cse("Spy"))
#ifdef EL_HAVE_QT5
    if( GuiDisabled() )
        LogicError("GUI was disabled");
    if( U == A.UGath && V == A.VGath )
    {
        if( A.CrossRank() == A.Root() && A.RedundantRank() == 0 )
            Spy( A.LockedMatrix(), title, tol );
    }
    else
    {
        DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
            Spy( A_CIRC_CIRC.Matrix(), title, tol );
    }
#else
    LogicError("Qt5 not available");
#endif // ifdef EL_HAVE_QT5
}

template<typename T,Dist U,Dist V>
void Spy( const BlockDistMatrix<T,U,V>& A, std::string title, Base<T> tol )
{
    DEBUG_ONLY(CallStackEntry cse("Spy"))
#ifdef EL_HAVE_QT5
    if( GuiDisabled() )
        LogicError("GUI was disabled");
    if( U == A.UGath && V == A.VGath )
    {
        if( A.CrossRank() == A.Root() && A.RedundantRank() == 0 )
            Spy( A.LockedMatrix(), title, tol );
    }
    else
    {
        BlockDistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
            Spy( A_CIRC_CIRC.Matrix(), title, tol );
    }
#else
    LogicError("Qt5 not available");
#endif // ifdef EL_HAVE_QT5
}

template<typename T>
void Spy
( const AbstractDistMatrix<T>& A, std::string title, Base<T> tol )
{
    DEBUG_ONLY(CallStackEntry cse("Spy"))
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = dynamic_cast<const DistMatrix<T,CDIST,RDIST>&>(A); \
      Spy( ACast, title, tol );
    #include "El/core/GuardAndPayload.h"
}

template<typename T>
void Spy
( const AbstractBlockDistMatrix<T>& A, std::string title, Base<T> tol )
{
    DEBUG_ONLY(CallStackEntry cse("Spy"))
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = dynamic_cast<const BlockDistMatrix<T,CDIST,RDIST>&>(A); \
      Spy( ACast, title, tol );
    #include "El/core/GuardAndPayload.h"
}

#define PROTO(T) \
  template void Spy ( const Matrix<T>& A, std::string title, Base<T> tol ); \
  template void Spy \
  ( const AbstractDistMatrix<T>& A, std::string title, Base<T> tol ); \
  template void Spy \
  ( const AbstractBlockDistMatrix<T>& A, std::string title, Base<T> tol ); 

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
