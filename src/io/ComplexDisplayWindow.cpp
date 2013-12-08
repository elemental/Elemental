/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#include "elemental/io.hpp"

#ifdef HAVE_QT5

#include <QBoxLayout>
#include <QCheckBox>

namespace elem {

ComplexDisplayWindow::ComplexDisplayWindow( QWidget* parent )
: QWidget(parent)
{
#ifndef RELEASE
    CallStackEntry cse("ComplexDisplayWindow::ComplexDisplayWindow");
#endif
    matrix_ = 0;
    QVBoxLayout* mainLayout = new QVBoxLayout();

    QHBoxLayout* matrixLayout = new QHBoxLayout();
    // Real data
    realDisplay_ = new DisplayWidget<Complex<double>>();
    realScroll_ = new QScrollArea();
    realScroll_->setWidget( realDisplay_ );
    matrixLayout->addWidget( realScroll_ );
    // Imaginary data
    imagDisplay_ = new DisplayWidget<Complex<double>>();
    imagScroll_ = new QScrollArea();
    imagScroll_->setWidget( imagDisplay_ );
    matrixLayout->addWidget( imagScroll_ );
    // Push both
    mainLayout->addLayout( matrixLayout );

    // Two buttons for saving real and imaginary images
    QHBoxLayout* saveLayout = new QHBoxLayout();
    QPushButton* realSaveButton = new QPushButton("Save real");
    QPushButton* imagSaveButton = new QPushButton("Save imag");
    saveLayout->addWidget( realSaveButton );
    saveLayout->addWidget( imagSaveButton );
    mainLayout->addLayout( saveLayout );

    // Checkbox for switching to the global scale
    QCheckBox* scaleBox = new QCheckBox("Global scale");
    mainLayout->addWidget( scaleBox );

    setLayout( mainLayout );
    connect( realSaveButton, SIGNAL(clicked()), this, SLOT(SaveReal()) );
    connect( imagSaveButton, SIGNAL(clicked()), this, SLOT(SaveImag()) );
    connect( scaleBox, SIGNAL(clicked(bool)), this, SLOT(SetScale(bool)) );
    setAttribute( Qt::WA_DeleteOnClose );

    // Elemental needs to know if a window was opened for cleanup purposes
    OpenedWindow();
}

ComplexDisplayWindow::~ComplexDisplayWindow()
{ delete matrix_; } 

void 
ComplexDisplayWindow::Display
( const Matrix<Complex<double>>* matrix, QString title )
{
#ifndef RELEASE
    CallStackEntry cse("ComplexDisplayWindow::Display");
#endif
    if( matrix_ != 0 )
        delete matrix_; 
    matrix_ = matrix;

    setWindowTitle( title );
    realDisplay_->DisplayReal( matrix );
    imagDisplay_->DisplayImag( matrix );
}

void 
ComplexDisplayWindow::Display
( const Matrix<Complex<double>>* matrix, 
  double minRealVal, double maxRealVal, 
  double minImagVal, double maxImagVal,
  QString title )
{
#ifndef RELEASE
    CallStackEntry cse("ComplexDisplayWindow::Display");
#endif
    if( matrix_ != 0 )
        delete matrix_;
    matrix_ = matrix;

    setWindowTitle( title );
    realDisplay_->DisplayReal( matrix, minRealVal, maxRealVal );
    imagDisplay_->DisplayImag( matrix, minImagVal, maxImagVal );
}

void
ComplexDisplayWindow::SaveReal()
{
#ifndef RELEASE
    CallStackEntry cse("ComplexDisplayWindow::SaveReal");
#endif
    std::ostringstream os;
    os << windowTitle().toStdString() << " (real)";
    realDisplay_->SavePng( os.str() );
}

void
ComplexDisplayWindow::SaveImag()
{
#ifndef RELEASE
    CallStackEntry cse("ComplexDisplayWindow::SaveImag");
#endif
    std::ostringstream os;
    os << windowTitle().toStdString() << " (imag)";
    imagDisplay_->SavePng( os.str() );
}

void 
ComplexDisplayWindow::SetScale( bool global )
{
#ifndef RELEASE
    CallStackEntry cse("ComplexDisplayWindow::SetScale");
#endif
    if( global )
    {
        const double minRealVal = MinRealWindowVal();
        const double maxRealVal = MaxRealWindowVal();
        const double minImagVal = MinImagWindowVal();
        const double maxImagVal = MaxImagWindowVal();
        realDisplay_->DisplayReal( matrix_, minRealVal, maxRealVal );
        imagDisplay_->DisplayImag( matrix_, minImagVal, maxImagVal );
    }
    else
    {
        realDisplay_->DisplayReal( matrix_ );
        imagDisplay_->DisplayImag( matrix_ );
    }
}

} // namespace elem

#endif // ifdef HAVE_QT5
