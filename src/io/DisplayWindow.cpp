/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#ifdef EL_HAVE_QT5

#include <QBoxLayout>
#include <QCheckBox>

namespace El {

DisplayWindow::DisplayWindow( QWidget* parent )
: QWidget(parent)
{
    DEBUG_CSE
    matrix_ = 0;
    QVBoxLayout* mainLayout = new QVBoxLayout();

    QHBoxLayout* matrixLayout = new QHBoxLayout();
    // Real data
    display_ = new DisplayWidget<double>();
    scroll_ = new QScrollArea();
    scroll_->setWidget( display_ );
    matrixLayout->addWidget( scroll_ );
    // Push
    mainLayout->addLayout( matrixLayout );

    // Add a save button and a check box for the scale
    QHBoxLayout* optionsLayout = new QHBoxLayout();
    QPushButton* saveButton = new QPushButton("Save");
    QCheckBox* scaleBox = new QCheckBox("Global scale");
    optionsLayout->addWidget( saveButton );
    optionsLayout->addWidget( scaleBox );
    // Push
    mainLayout->addLayout( optionsLayout ); 

    setLayout( mainLayout );
    connect( saveButton, SIGNAL(clicked()), this, SLOT(Save()) );
    connect( scaleBox, SIGNAL(clicked(bool)), this, SLOT(SetScale(bool)) );
    setAttribute( Qt::WA_DeleteOnClose );

    // Elemental needs to know if a window was opened for cleanup purposes
    OpenedWindow();
}

DisplayWindow::~DisplayWindow()
{ delete matrix_; }

void
DisplayWindow::Display( const Matrix<double>* matrix, QString title )
{
    DEBUG_CSE
    if( matrix_ != 0 )
        delete matrix_;
    matrix_ = matrix;

    setWindowTitle( title );
    display_->DisplayReal( matrix );
}

void
DisplayWindow::Display
( const Matrix<double>* matrix, double minVal, double maxVal, QString title )
{
    DEBUG_CSE
    if( matrix_ != 0 )
        delete matrix_;
    matrix_ = matrix;

    setWindowTitle( title );
    display_->DisplayReal( matrix, minVal, maxVal );
}

void
DisplayWindow::Save()
{
    DEBUG_CSE
    display_->SavePng( windowTitle().toStdString() );
}

void
DisplayWindow::SetScale( bool global )
{
    DEBUG_CSE
    if( global )
    {
        const double minVal = MinRealWindowVal();
        const double maxVal = MaxRealWindowVal();
        display_->DisplayReal( matrix_, minVal, maxVal );
    }
    else
    {
        display_->DisplayReal( matrix_ );
    }
}

} // namespace El

#endif // ifdef EL_HAVE_QT5
