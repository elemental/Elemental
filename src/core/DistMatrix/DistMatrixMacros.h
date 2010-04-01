/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once

#ifndef RELEASE

#define REPORT_UNIMPLEMENTED_FEATURE \
{ \
    throw "Sorry, feature not yet implemented."; \
}
#else
#define REPORT_UNIMPLEMENTED_FEATURE \
{ \
    throw "Sorry, feature not yet implemented.\r\n" \
          "Run debug version to see which routine."; \
}
#endif

#ifndef RELEASE

#define CHECK_IF_ALIGNING_DIFF_GRID( A ) \
{ \
    if( GetGrid() != A.GetGrid() ) \
        throw "Tried to align with a matrix from diff. grid."; \
}

#define CHECK_IF_CONFORMING_1x2( AL, AR ) \
{ \
    if( AL.Height() != AR.Height() ) \
    { \
        std::ostringstream msg; \
        msg << "1x2 not conformant. Left is " << AL.Height() \
            << " x " << AL.Width() << ", right is " \
            << AR.Height() << " x " << AR.Width() << std::endl; \
        throw msg.str(); \
    } \
}
#define CHECK_IF_CONFORMING_2x1( AT, AB ) \
{ \
    if( AT.Width() != AB.Width() ) \
    { \
        std::ostringstream msg; \
        msg << "2x1 is not conformant. Top is " << AT.Height() \
            << " x " << AT.Width() << ", bottom is " \
            << AB.Height() << " x " << AB.Width() << std::endl; \
        throw msg.str(); \
    } \
}
#define CHECK_IF_CONFORMING_2x2( ATL, ATR, ABL, ABR ) \
{ \
    if( ATL.Width() != ABL.Width()   ||   \
        ATR.Width() != ABR.Width()   ||   \
        ATL.Height() != ATR.Height() ||   \
        ABL.Height() != ABR.Height()    ) \
    { \
        throw "2x2 is not conformant."; \
    } \
}

#define CHECK_IF_CONFORMING_DIFF_GRID( A ) \
{ \
    if( GetGrid() != A.GetGrid() ) \
        throw "Tried to conform with a matrix from diff. grid."; \
}

#define CHECK_IF_LOCKED_VIEW \
{ \
    if( _viewing && _lockedView ) \
        throw "Cannot alter data with locked view."; \
}

#define CHECK_IF_UNFREED_COL_CONSTRAINT \
{ \
    if( ConstrainedColDist() ) \
        throw "Forgot to free column constraint before reconstraining."; \
}

#define CHECK_IF_UNFREED_ROW_CONSTRAINT \
{ \
    if( ConstrainedRowDist() ) \
        throw "Forgot to free row constraint before reconstraining."; \
}

#define CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width ) \
{ \
    if( i < 0 || j < 0 ) \
        throw "Indices must be non-negative."; \
    if( height < 0 || width < 0 ) \
        throw "Height and width must be non-negative."; \
    if( (i+height) > A.Height() || (j+width) > A.Width() ) \
    { \
        ostringstream msg; \
        msg << "Out of bounds of distributed matrix: up to (" \
            << i+height-1 << "," << j+width-1 << ") of " \
            << A.Height() << "x" << A.Width() << " matrix." << std::endl;\
        throw msg.str(); \
    } \
}

#define CHECK_IF_REDIST_DIFF_GRID( A ) \
{ \
    if( GetGrid() != A.GetGrid() ) \
        throw "Tried to redist. btw. matrices from different grid."; \
}

#define CHECK_IF_DIFF_SIZE( A ) \
{ \
    if( Height() != A.Height() || Width() != A.Width() ) \
        throw "Distributed matrices are different sizes."; \
}

#define CHECK_IF_VIEWING_DIFF_SIZE( A ) \
{ \
    if( _viewing && (Height() != A.Height() || Width() != A.Width()) ) \
        throw "Cannot assign to different sized view."; \
}

#define CHECK_IF_VIEWING_AND_STORING \
{ \
    if( LockedLocalMatrix().MemorySize() > 0 ) \
        throw "Attempting to view and store with distributed matrix."; \
}

#define CHECK_IF_VIEWING_DIFF_GRID( A ) \
{ \
    if( GetGrid() != A.GetGrid() ) \
        throw "Tried to view a matrix from a different grid."; \
}

#endif

