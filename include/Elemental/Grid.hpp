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
#ifndef ELEMENTAL_PROCESSGRID_HPP
#define ELEMENTAL_PROCESSGRID_HPP 1

#include <iostream>
#include "mpi.h"

namespace Elemental {

class Grid
{
    int        _p;
    int        _c;
    int        _r;
    int        _gcd;
    int        _col;
    int        _row;
    int        _rank; 
    int        _matrixColRank;
    int        _matrixRowRank;
    int        _vectorColRank;
    int        _vectorRowRank;
    int*       _diagPathsAndRanks; 
    MPI_Comm   _comm;
    MPI_Comm   _matrixColComm;
    MPI_Comm   _matrixRowComm;
    MPI_Comm   _vectorColComm;
    MPI_Comm   _vectorRowComm;

    void Init( int r, int c );

    public:

    Grid( MPI_Comm comm );
    Grid( MPI_Comm comm, int r, int c );

    ~Grid();
        
    int Size() const;
    int Height() const;
    int Width() const;
    int GCD() const;
    int LCM() const;
    int DiagPath() const;
    int DiagPath( const int vectorColRank ) const;
    int DiagPathRank() const;
    int DiagPathRank( const int vectorColRank ) const;
    int MCRank() const;
    int MRRank() const;
    int VCRank() const;
    int VRRank() const;
    MPI_Comm MCComm() const;
    MPI_Comm MRComm() const;
    MPI_Comm VCComm() const;
    MPI_Comm VRComm() const;
};

bool operator== ( const Grid& A, const Grid& B );
bool operator!= ( const Grid& A, const Grid& B );

} // Elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline
Elemental::Grid::~Grid()
{
    MPI_Comm_free( &_matrixColComm );
    MPI_Comm_free( &_matrixRowComm );
    MPI_Comm_free( &_vectorColComm );
    MPI_Comm_free( &_vectorRowComm );
    MPI_Comm_free( &_comm );

    delete _diagPathsAndRanks;
}

inline int
Elemental::Grid::Size() const
{ return _p; }

inline int
Elemental::Grid::Height() const
{ return _r; }

inline int
Elemental::Grid::Width() const
{ return _c; }

inline int
Elemental::Grid::GCD() const
{ return _gcd; }

inline int
Elemental::Grid::LCM() const
{ return _p/_gcd; }

inline int
Elemental::Grid::DiagPath() const
{ return _diagPathsAndRanks[2*_vectorColRank]; }

inline int
Elemental::Grid::DiagPath( const int vectorColRank ) const
{ return _diagPathsAndRanks[2*vectorColRank]; }

inline int
Elemental::Grid::DiagPathRank() const
{ return _diagPathsAndRanks[2*_vectorColRank+1]; }

inline int
Elemental::Grid::DiagPathRank( const int vectorColRank ) const
{ return _diagPathsAndRanks[2*vectorColRank+1]; }

inline int
Elemental::Grid::MCRank() const
{ return _matrixColRank; }

inline int
Elemental::Grid::MRRank() const
{ return _matrixRowRank; }

inline int
Elemental::Grid::VCRank() const
{ return _vectorColRank; }

inline int
Elemental::Grid::VRRank() const
{ return _vectorRowRank; }

inline MPI_Comm
Elemental::Grid::MCComm() const
{ return _matrixColComm; }

inline MPI_Comm
Elemental::Grid::MRComm() const
{ return _matrixRowComm; }

inline MPI_Comm
Elemental::Grid::VCComm() const
{ return _vectorColComm; }

inline MPI_Comm
Elemental::Grid::VRComm() const
{ return _vectorRowComm; }

inline bool 
Elemental::operator== ( const Grid& A, const Grid& B )
{ return &A == &B; }

inline bool 
Elemental::operator!= ( const Grid& A, const Grid& B )
{ return &A != &B; }

#endif /* ELEMENTAL_PROCESSGRID_HPP */

