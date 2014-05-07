/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SPYWIDGET_IMPL_HPP
#define ELEM_SPYWIDGET_IMPL_HPP
#ifdef ELEM_HAVE_QT5

#include <QPainter>
#include <QPixmap>
#include <QStylePainter>

#include "elemental/io/DisplayWidget/decl.hpp"

namespace elem {

inline
SpyWidget::SpyWidget( QWidget* parent )
: QWidget(parent)
{ }

inline
SpyWidget::~SpyWidget()
{ }

inline void 
SpyWidget::paintEvent( QPaintEvent* event )
{
    DEBUG_ONLY(CallStackEntry cse("SpyWidget::paintEvent"))
    QStylePainter painter( this );
    painter.drawPixmap( 0, 0, pixmap_ );
}

inline void 
SpyWidget::Spy( const Matrix<Int>* A )
{
    DEBUG_ONLY(CallStackEntry cse("SpyWidget::Spy"))
    const Int m = A->Height();
    const Int n = A->Width();

    // TODO: Parameterize these instead
    const Int mPix = m;
    const Int nPix = n;
    const double mRatio = double(m) / double(mPix);
    const double nRatio = double(n) / double(nPix);
    pixmap_ = QPixmap( nPix, mPix );
    resize( nPix, mPix );

    // Paint the matrix
    QPainter painter( &pixmap_ );
    painter.initFrom( this );
    for( Int jPix=0; jPix<nPix; ++jPix ) 
    {
        const Int j = nRatio*jPix;
        for( Int iPix=0; iPix<mPix; ++iPix )
        {
            const Int i = mRatio*iPix;
            const Int z = A->Get(i,j);
            if( z == 0 )
                painter.setPen( qRgba(0,0,0,255) );
            else 
                painter.setPen( qRgba(255,255,255,255) );
            painter.drawPoint( jPix, iPix );
        }
    }

    // Refresh the widget
    update();
}

} // namespace elem

#endif // ifdef ELEM_HAVE_QT5
#endif // ifndef ELEM_SPYWIDGET_IMPL_HPP
