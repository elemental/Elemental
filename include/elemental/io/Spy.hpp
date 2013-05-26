/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef IO_SPY_HPP
#define IO_SPY_HPP

#ifdef HAVE_QT5

#include "elemental/io/spy_window_decl.hpp"
#include "elemental/io/spy_widget_impl.hpp"

#include <QApplication>

namespace elem {

template<typename T>
inline void
Spy( const Matrix<T>& A, std::string title="", BASE(T) tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("Spy");
#endif
    // Convert A to double-precision since Qt's MOC does not support templates
    const int m = A.Height();
    const int n = A.Width();
    Matrix<int>* ASpy = new Matrix<int>( m, n );
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            ASpy->Set( i, j, ( Abs(A.Get(i,j))>tol ? 1 : 0 ) );

    QString qTitle = QString::fromStdString( title );
    SpyWindow* spyWindow = new SpyWindow;
    spyWindow->Spy( ASpy, qTitle );
    spyWindow->show();

    // Spend at most 200 milliseconds rendering
    QCoreApplication::instance()->processEvents( QEventLoop::AllEvents, 200 );
}

template<typename T,Distribution U,Distribution V>
inline void
Spy( const DistMatrix<T,U,V>& A, std::string title="", BASE(T) tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("Spy");
#endif
    // TODO: Avoid giving every process a full copy and think about avoiding
    //       the extra copy needed when the underlying datatype is not 'double'
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    if( A.Grid().Rank() == 0 )
        Spy( A_STAR_STAR.Matrix(), title, tol );
}

} // namespace elem

#endif // ifdef HAVE_QT5

#endif // ifndef IO_SPY_HPP
