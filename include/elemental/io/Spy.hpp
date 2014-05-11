/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SPY_HPP
#define ELEM_SPY_HPP

#include "./SpyWindow/decl.hpp"
#include "./SpyWidget/impl.hpp"
#include "./Display.hpp" // for ProcessEvents

#ifdef ELEM_HAVE_QT5
# include <QApplication>
#endif

namespace elem {

template<typename T>
inline void
Spy( const Matrix<T>& A, std::string title="Default", Base<T> tol=0 )
{
    DEBUG_ONLY(CallStackEntry cse("Spy"))
#ifdef ELEM_HAVE_QT5
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
#endif // ifdef ELEM_HAVE_QT5
}

template<typename T,Distribution U,Distribution V>
inline void
Spy( const DistMatrix<T,U,V>& A, std::string title="Default", Base<T> tol=0 )
{
    DEBUG_ONLY(CallStackEntry cse("Spy"))
#ifdef ELEM_HAVE_QT5
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
#endif // ifdef ELEM_HAVE_QT5
}

template<typename T,Dist U,Dist V>
inline void
Spy
( const BlockDistMatrix<T,U,V>& A, std::string title="Default", Base<T> tol=0 )
{
    DEBUG_ONLY(CallStackEntry cse("Spy"))
#ifdef ELEM_HAVE_QT5
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
#endif // ifdef ELEM_HAVE_QT5
}

} // namespace elem

#endif // ifndef ELEM_SPY_HPP
