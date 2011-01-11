/*
   Copyright (c) 2009-2011, Jack Poulson
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
#include "elemental/dist_matrix.hpp"
#include "elemental/exports/c.h"

namespace {
std::vector<elemental::Grid*> 
    gridList;
std::vector<elemental::AbstractDistMatrix<double>*> 
    distMatrixDoubleList;
std::vector<elemental::AbstractDistMatrix<elemental::dcomplex>*> 
    distMatrixDComplexList;
}

extern "C" {

void ElementalInit( int* argc, char** argv[] )
{
    elemental::Init( argc, argv );    
}

void ElementalFinalize()
{
    ElementalClearDistMatrices();
    ElementalClearGrids();
    elemental::Finalize();
}

int ElementalBlocksize()
{
    return elemental::Blocksize();
}

void ElementalSetBlocksize( int blocksize )
{
    elemental::SetBlocksize( blocksize );
}

void ElementalPushBlocksizeStack( int blocksize )
{
    elemental::PushBlocksizeStack( blocksize );
}

void ElementalPopBlocksizeStack()
{
    elemental::PopBlocksizeStack();
}

void ElementalClearGrids()
{
    for( int j=0; j< ::gridList.size(); ++j )
        delete ::gridList[j];
    ::gridList.clear();
}

int ElementalDefaultGrid( MPI_Comm comm )
{
    int gridHandle = ::gridList.size();
    ::gridList.push_back( new elemental::Grid( comm ) );
    return gridHandle;
}

int ElementalGrid( MPI_Comm comm, int r, int c )
{
    int gridHandle = ::gridList.size();
    ::gridList.push_back( new elemental::Grid( comm, r, c ) );
    return gridHandle;
}

// Don't yet export the viewingComm/owningGroup constructors.

int ElementalGridHeight( int gridHandle )
{
    return gridList[gridHandle]->Height();
}

int ElementalGridWidth( int gridHandle )
{
    return gridList[gridHandle]->Width();
}

int ElementalGridSize( int gridHandle )
{
    return gridList[gridHandle]->Size();
}

int ElementalInGrid( int gridHandle )
{
    return gridList[gridHandle]->InGrid();
}

int ElementalGridVCRank( int gridHandle )
{
    return gridList[gridHandle]->VCRank();
}

int ElementalGridVRRank( int gridHandle )
{
    return gridList[gridHandle]->VRRank();
}

int ElementalGridMCRank( int gridHandle )
{
    return gridList[gridHandle]->MCRank();
}

int ElementalGridMRRank( int gridHandle )
{
    return gridList[gridHandle]->MRRank();
}

MPI_Comm ElementalGridVCComm( int gridHandle )
{
    return gridList[gridHandle]->VCComm();
}

MPI_Comm ElementalGridVRComm( int gridHandle )
{
    return gridList[gridHandle]->VRComm();
}

MPI_Comm ElementalGridMCComm( int gridHandle )
{
    return gridList[gridHandle]->MCComm();
}

MPI_Comm ElementalGridMRComm( int gridHandle )
{
    return gridList[gridHandle]->MRComm();
}

void ElementalClearDistMatrices()
{
    for( int j=0; j< ::distMatrixDoubleList.size(); ++j )
        delete ::distMatrixDoubleList[j];
    for( int j=0; j< ::distMatrixDComplexList.size(); ++j )
        delete ::distMatrixDComplexList[j];
    ::distMatrixDoubleList.clear();
    ::distMatrixDComplexList.clear();
}

int 
ElementalDistMatrix_MC_MR_Double
( int height, int width, int colAlignment, int rowAlignment,
  double* buffer, int ldim, int gridHandle )
{
    int distMatrixDoubleHandle = ::distMatrixDoubleList.size();
    ::distMatrixDoubleList.push_back
        ( new elemental::DistMatrix<double,elemental::MC,elemental::MR>
            ( height, width, colAlignment, rowAlignment, buffer, ldim, 
              *::gridList[gridHandle] ) );
    return distMatrixDoubleHandle;
}

int
ElementalDistMatrix_MC_MR_DComplex
( int height, int width, int colAlignment, int rowAlignment,
  ElementalDComplex* buffer, int ldim, int gridHandle )
{
    int distMatrixDComplexHandle = ::distMatrixDComplexList.size();
    ::distMatrixDComplexList.push_back
        ( new elemental::DistMatrix<elemental::dcomplex,
                                    elemental::MC,elemental::MR>
            ( height, width, colAlignment, rowAlignment, 
              (elemental::dcomplex*)buffer, ldim, *::gridList[gridHandle] ) );
    return distMatrixDComplexHandle;
}

void ElementalDistMatrixDoublePrint( char* msg, int distMatrixDoubleHandle )
{
    std::string s( msg );
    ::distMatrixDoubleList[distMatrixDoubleHandle]->Print( s );
}

void ElementalDistMatrixDComplexPrint( char* msg, int distMatrixDComplexHandle )
{
    std::string s( msg );
    ::distMatrixDComplexList[distMatrixDComplexHandle]->Print( s );
}

int
ElementalLocalLength
( int globalLength, int myIndex, int alignment, int modulus )
{ 
    return elemental::utilities::LocalLength
        ( globalLength, myIndex, alignment, modulus );
}

} // extern "C"

