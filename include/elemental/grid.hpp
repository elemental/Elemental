/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#ifndef ELEMENTAL_PROCESSGRID_HPP
#define ELEMENTAL_PROCESSGRID_HPP 1

#include <iostream>
#include "mpi.h"

namespace elemental {

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
    int DiagPath( int vectorColRank ) const;
    int DiagPathRank() const;
    int DiagPathRank( int vectorColRank ) const;
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
elemental::Grid::~Grid()
{
    MPI_Comm_free( &_matrixColComm );
    MPI_Comm_free( &_matrixRowComm );
    MPI_Comm_free( &_vectorColComm );
    MPI_Comm_free( &_vectorRowComm );
    MPI_Comm_free( &_comm );

    delete _diagPathsAndRanks;
}

inline int
elemental::Grid::Size() const
{ return _p; }

inline int
elemental::Grid::Height() const
{ return _r; }

inline int
elemental::Grid::Width() const
{ return _c; }

inline int
elemental::Grid::GCD() const
{ return _gcd; }

inline int
elemental::Grid::LCM() const
{ return _p/_gcd; }

inline int
elemental::Grid::DiagPath() const
{ return _diagPathsAndRanks[2*_vectorColRank]; }

inline int
elemental::Grid::DiagPath( int vectorColRank ) const
{ return _diagPathsAndRanks[2*vectorColRank]; }

inline int
elemental::Grid::DiagPathRank() const
{ return _diagPathsAndRanks[2*_vectorColRank+1]; }

inline int
elemental::Grid::DiagPathRank( int vectorColRank ) const
{ return _diagPathsAndRanks[2*vectorColRank+1]; }

inline int
elemental::Grid::MCRank() const
{ return _matrixColRank; }

inline int
elemental::Grid::MRRank() const
{ return _matrixRowRank; }

inline int
elemental::Grid::VCRank() const
{ return _vectorColRank; }

inline int
elemental::Grid::VRRank() const
{ return _vectorRowRank; }

inline MPI_Comm
elemental::Grid::MCComm() const
{ return _matrixColComm; }

inline MPI_Comm
elemental::Grid::MRComm() const
{ return _matrixRowComm; }

inline MPI_Comm
elemental::Grid::VCComm() const
{ return _vectorColComm; }

inline MPI_Comm
elemental::Grid::VRComm() const
{ return _vectorRowComm; }

inline bool 
elemental::operator== ( const Grid& A, const Grid& B )
{ return &A == &B; }

inline bool 
elemental::operator!= ( const Grid& A, const Grid& B )
{ return &A != &B; }

#endif /* ELEMENTAL_PROCESSGRID_HPP */

