/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DISPLAY_HPP
#define EL_DISPLAY_HPP

#ifdef EL_HAVE_QT5
# include "./DisplayWindow-premoc.hpp"
# include "./ComplexDisplayWindow-premoc.hpp"
# include "./DisplayWidget/impl.hpp"
# include <QApplication>
#endif

namespace El {

inline void
ProcessEvents( int numMsecs )
{
#ifdef EL_HAVE_QT5
    QCoreApplication::instance()->processEvents
    ( QEventLoop::AllEvents, numMsecs );
#endif
}

template<typename T>
inline void
Display( const Matrix<T>& A, std::string title="Default" )
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
inline void
Display( const Matrix<Complex<T>>& A, std::string title="Default" )
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
inline void
Display( const DistMatrix<T,U,V>& A, std::string title="Default" )
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
inline void
Display( const BlockDistMatrix<T,U,V>& A, std::string title="Default" )
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
inline void
Display( const DynamicDistMatrix<T>& ADyn, std::string title="" )
{
    DEBUG_ONLY(CallStackEntry cse("Display"))
    #define IF_CONVERT_AND_DISPLAY(ADYN,TITLE,CDIST,RDIST) \
      if( ADYN.U == CDIST && ADYN.V == RDIST ) \
      { \
          auto A = dynamic_cast<const DistMatrix<T,CDIST,RDIST>*>(ADYN.ADM); \
          if( A == nullptr ) \
              RuntimeError("Dynamic cast failed"); \
          Display( *A, TITLE ); \
      }
    #define ELSEIF_CONVERT_AND_DISPLAY(ADYN,TITLE,CDIST,RDIST) \
      else IF_CONVERT_AND_DISPLAY(ADYN,TITLE,CDIST,RDIST)

    IF_CONVERT_AND_DISPLAY(    ADyn,title,CIRC,CIRC)
    ELSEIF_CONVERT_AND_DISPLAY(ADyn,title,MC,  MR  )
    ELSEIF_CONVERT_AND_DISPLAY(ADyn,title,MC,  STAR)
    ELSEIF_CONVERT_AND_DISPLAY(ADyn,title,MD,  STAR)
    ELSEIF_CONVERT_AND_DISPLAY(ADyn,title,MR,  MC  )
    ELSEIF_CONVERT_AND_DISPLAY(ADyn,title,MR,  STAR)
    ELSEIF_CONVERT_AND_DISPLAY(ADyn,title,STAR,MC  )
    ELSEIF_CONVERT_AND_DISPLAY(ADyn,title,STAR,MD  )
    ELSEIF_CONVERT_AND_DISPLAY(ADyn,title,STAR,MR  )
    ELSEIF_CONVERT_AND_DISPLAY(ADyn,title,STAR,STAR)
    ELSEIF_CONVERT_AND_DISPLAY(ADyn,title,STAR,VC  )
    ELSEIF_CONVERT_AND_DISPLAY(ADyn,title,STAR,VR  )
    ELSEIF_CONVERT_AND_DISPLAY(ADyn,title,VC,  STAR)
    ELSEIF_CONVERT_AND_DISPLAY(ADyn,title,VR,  STAR)

    #undef ELSEIF_CONVERT_AND_DISPLAY
    #undef IF_CONVERT_AND_DISPLAY
}

} // namespace El

#endif // ifndef EL_DISPLAY_HPP
