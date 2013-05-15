/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef GRAPHICS_SPYWIDGET_DECL_HPP
#define GRAPHICS_SPYWIDGET_DECL_HPP

#ifdef HAVE_QT5

#include <QWidget>

namespace elem {

class SpyWidget : public QWidget
{
public:
    SpyWidget( QWidget* parent=0 );
    ~SpyWidget();
    void Spy( const Matrix<int>* A );
    // TODO: Change style
protected:
    void paintEvent( QPaintEvent* event );

private:
    QPixmap pixmap_;
};

} // namespace elem

#endif // ifdef HAVE_QT5

#endif // ifndef GRAPHICS_SPYWIDGET_DECL_HPP
