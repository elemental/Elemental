/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include "El/io.hpp"

#ifdef EL_HAVE_QT5
# include "El/io/DisplayWindow-premoc.hpp"
# include "El/io/ComplexDisplayWindow-premoc.hpp"
# include <QApplication>
#endif

namespace El {

void ProcessEvents( int numMsecs )
{
#ifdef EL_HAVE_QT5
    QCoreApplication::instance()->processEvents
    ( QEventLoop::AllEvents, numMsecs );
#endif
}

template<typename T>
void Display( const Matrix<T>& A, std::string title )
{
    DEBUG_ONLY(CallStackEntry cse("Display"))
#ifdef EL_HAVE_QT5
    if( GuiDisabled() )
    {
        Print( A, title );
        return;
    }

    // Convert A to double-precision since Qt's MOC does not support templates
    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<double>* ADouble = new Matrix<double>( m, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            ADouble->Set( i, j, double(A.Get(i,j)) );

    QString qTitle = QString::fromStdString( title );
    DisplayWindow* displayWindow = new DisplayWindow;
    displayWindow->Display( ADouble, qTitle );
    displayWindow->show();

    // Spend at most 200 milliseconds rendering
    ProcessEvents( 200 );
#else
    Print( A, title );
#endif
}

template<typename T>
void Display( const Matrix<Complex<T>>& A, std::string title )
{
    DEBUG_ONLY(CallStackEntry cse("Display"))
#ifdef EL_HAVE_QT5
    if( GuiDisabled() )
    {
        Print( A, title );
        return;
    }

    // Convert A to double-precision since Qt's MOC does not support templates
    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<Complex<double>>* ADouble = new Matrix<Complex<double>>( m, n );
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            const Complex<T> alpha = A.Get(i,j);
            const Complex<double> alphaDouble = 
                Complex<double>(alpha.real(),alpha.imag()); 
            ADouble->Set( i, j, alphaDouble );
        }
    }

    QString qTitle = QString::fromStdString( title );
    ComplexDisplayWindow* displayWindow = new ComplexDisplayWindow;
    displayWindow->Display( ADouble, qTitle );
    displayWindow->show();

    // Spend at most 200 milliseconds rendering
    ProcessEvents( 200 );
#else
    Print( A, title );
#endif
}

template<typename T,Dist U,Dist V>
void Display( const DistMatrix<T,U,V>& A, std::string title )
{
    DEBUG_ONLY(CallStackEntry cse("Display"))
#ifdef EL_HAVE_QT5
    if( GuiDisabled() )
    {
        Print( A, title );
        return;
    }

    if( U == A.UGath && V == A.VGath )
    {
        if( A.CrossRank() == A.Root() && A.RedundantRank() == 0 )
            Display( A.LockedMatrix(), title );
    }
    else
    {
        DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
            Display( A_CIRC_CIRC.Matrix(), title );
    }
#else
    Print( A, title );
#endif
}

template<typename T,Dist U,Dist V>
void Display( const BlockDistMatrix<T,U,V>& A, std::string title )
{
    DEBUG_ONLY(CallStackEntry cse("Display"))
#ifdef EL_HAVE_QT5
    if( GuiDisabled() )
    {
        Print( A, title );
        return;
    }

    if( U == A.UGath && V == A.VGath )
    {
        if( A.CrossRank() == A.Root() && A.RedundantRank() == 0 )
            Display( A.LockedMatrix(), title );
    }
    else
    {
        BlockDistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
            Display( A_CIRC_CIRC.Matrix(), title );
    }
#else
    Print( A, title );
#endif
}

template<typename T>
void Display( const AbstractDistMatrix<T>& A, std::string title )
{
    DEBUG_ONLY(CallStackEntry cse("Display"))
    #define GUARD(CDIST,RDIST) \
        A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = dynamic_cast<const DistMatrix<T,CDIST,RDIST>&>(A); \
      Display( ACast, title );
    #include "El/core/GuardAndPayload.h"
}

template<typename T>
void Display( const AbstractBlockDistMatrix<T>& A, std::string title )
{
    DEBUG_ONLY(CallStackEntry cse("Display"))
    #define GUARD(CDIST,RDIST) \
        A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = dynamic_cast<const BlockDistMatrix<T,CDIST,RDIST>&>(A); \
      Display( ACast, title );
    #include "El/core/GuardAndPayload.h"
}

#define DISTPROTO(T,U,V) \
  template void Display( const DistMatrix<T,U,V>& A, std::string title ); \
  template void Display( const BlockDistMatrix<T,U,V>& A, std::string title );

#define PROTO(T) \
  template void Display( const Matrix<T>& A, std::string title ); \
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
  template void Display \
  ( const AbstractDistMatrix<T>& A, std::string title ); \
  template void Display \
  ( const AbstractBlockDistMatrix<T>& A, std::string title ); 

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
