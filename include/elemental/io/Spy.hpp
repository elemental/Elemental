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

#ifdef HAVE_QT5

#include "elemental/io/spy_window_decl.hpp"
#include "elemental/io/spy_widget_impl.hpp"
#include "elemental/io/Display.hpp" // for ProcessEvents

#include <QApplication>

namespace elem {

template<typename T>
inline void
Spy( const Matrix<T>& A, std::string title="Default", Base<T> tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("Spy");
#endif
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
}

template<typename T,Distribution U,Distribution V>
inline void
Spy( const DistMatrix<T,U,V>& A, std::string title="Default", Base<T> tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("Spy");
#endif
    DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A );
    if( A.Grid().Rank() == A_CIRC_CIRC.Root() )
        Spy( A_CIRC_CIRC.Matrix(), title, tol );
}

// If already in [* ,* ] or [o ,o ] distributions, no copy is needed
template<typename T>
inline void
Spy
( const DistMatrix<T,STAR,STAR>& A, std::string title="Default", Base<T> tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("Spy");
#endif
    if( A.Grid().Rank() == 0 )
        Spy( A.LockedMatrix(), title, tol );
}
template<typename T>
inline void
Spy
( const DistMatrix<T,CIRC,CIRC>& A, std::string title="Default", Base<T> tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("Spy");
#endif
    if( A.Grid().Rank() == A.Root() )
        Spy( A.LockedMatrix(), title, tol );
}

} // namespace elem

#endif // ifdef HAVE_QT5

#endif // ifndef ELEM_IO_SPY_HPP
