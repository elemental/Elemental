/*
   Copyright (c) 2009-2010, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#pragma once
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

