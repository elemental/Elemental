/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_IO_SPY_HPP
#define ELEM_IO_SPY_HPP

#include "elemental/io/spy_window/decl.hpp"
#include "elemental/io/spy_widget/impl.hpp"
#include "elemental/io/Display.hpp" // for ProcessEvents

#ifdef HAVE_QT5
# include <QApplication>
#endif

namespace elem {

template<typename T>
inline void
Spy( const Matrix<T>& A, std::string title="Default", BASE(T) tol=0 )
{
    DEBUG_ONLY(CallStackEntry cse("Spy"))
#ifdef HAVE_QT5
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
#endif // ifdef HAVE_QT5
}

template<typename T,Distribution U,Distribution V>
inline void
Spy( const DistMatrix<T,U,V>& A, std::string title="Default", BASE(T) tol=0 )
{
    DEBUG_ONLY(CallStackEntry cse("Spy"))
#ifdef HAVE_QT5
    if( GuiDisabled() )
        LogicError("GUI was disabled");
    DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A );
    if( A.Grid().Rank() == A_CIRC_CIRC.Root() )
        Spy( A_CIRC_CIRC.Matrix(), title, tol );
#else
    LogicError("Qt5 not available");
#endif // ifdef HAVE_QT5
}

// If already in [* ,* ] or [o ,o ] distributions, no copy is needed
template<typename T>
inline void
Spy
( const DistMatrix<T,STAR,STAR>& A, std::string title="Default", BASE(T) tol=0 )
{
    DEBUG_ONLY(CallStackEntry cse("Spy"))
#ifdef HAVE_QT5
    if( GuiDisabled() )
        LogicError("GUI was disabled");
    if( A.Grid().Rank() == 0 )
        Spy( A.LockedMatrix(), title, tol );
#else
    LogicError("Qt5 not available");
#endif // ifdef HAVE_QT5
}
template<typename T>
inline void
Spy
( const DistMatrix<T,CIRC,CIRC>& A, std::string title="Default", BASE(T) tol=0 )
{
    DEBUG_ONLY(CallStackEntry cse("Spy"))
#ifdef HAVE_QT5
    if( GuiDisabled() )
        LogicError("GUI was disabled");
    if( A.Grid().Rank() == A.Root() )
        Spy( A.LockedMatrix(), title, tol );
#else
    LogicError("Qt5 not available");
#endif // ifdef HAVE_QT5
}

} // namespace elem

#endif // ifndef ELEM_IO_SPY_HPP
