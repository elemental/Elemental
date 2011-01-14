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
#include "elemental.hpp"
#include "elemental/exports/c.h"

namespace {
using elemental::MC;
using elemental::MD;
using elemental::MR;
using elemental::Star;
using elemental::VC;
using elemental::VR;
using elemental::scomplex;
using elemental::dcomplex;
using elemental::DistMatrix;

std::vector<elemental::Grid*> gridList;
std::vector<DistMatrix<float,MC,MR>*> MC_MR_SingleList;
std::vector<DistMatrix<double,MC,MR>*> MC_MR_DoubleList;
std::vector<DistMatrix<float,Star,VR>*> Star_VR_SingleList;
std::vector<DistMatrix<double,Star,VR>*> Star_VR_DoubleList;
#ifndef WITHOUT_COMPLEX
std::vector<DistMatrix<scomplex,MC,MR>*> MC_MR_SComplexList;
std::vector<DistMatrix<dcomplex,MC,MR>*> MC_MR_DComplexList;
std::vector<DistMatrix<scomplex,Star,VR>*> Star_VR_SComplexList;
std::vector<DistMatrix<dcomplex,Star,VR>*> Star_VR_DComplexList;
#endif // WITHOUT_COMPLEX
}

extern "C" {

/* TODO: Add wrapper macro for handling exceptions. */

void ElementalInit( int* argc, char** argv[] )
{ elemental::Init( argc, argv ); }

void ElementalFinalize()
{
    for( int j=0; j< ::MC_MR_SingleList.size(); ++j )
        delete ::MC_MR_SingleList[j];
    for( int j=0; j< ::MC_MR_DoubleList.size(); ++j )
        delete ::MC_MR_DoubleList[j];
    for( int j=0; j< ::Star_VR_SingleList.size(); ++j )
        delete ::Star_VR_SingleList[j];
    for( int j=0; j< ::Star_VR_DoubleList.size(); ++j )
        delete ::Star_VR_DoubleList[j];
    /* TODO: Add other distributions */

#ifndef WITHOUT_COMPLEX
    for( int j=0; j< ::MC_MR_SComplexList.size(); ++j )
        delete ::MC_MR_SComplexList[j];
    for( int j=0; j< ::MC_MR_DComplexList.size(); ++j )
        delete ::MC_MR_DComplexList[j];
    for( int j=0; j< ::Star_VR_SComplexList.size(); ++j )
        delete ::Star_VR_SComplexList[j];
    for( int j=0; j< ::Star_VR_DComplexList.size(); ++j )
        delete ::Star_VR_DComplexList[j];
    /* TODO: Add other distributions */
#endif

    for( int j=0; j< ::gridList.size(); ++j )
        delete ::gridList[j];

    elemental::Finalize();
}

int ElementalBlocksize()
{ return elemental::Blocksize(); }

void ElementalSetBlocksize( int blocksize )
{ elemental::SetBlocksize( blocksize ); }

void ElementalPushBlocksizeStack( int blocksize )
{ elemental::PushBlocksizeStack( blocksize ); }

void ElementalPopBlocksizeStack()
{ elemental::PopBlocksizeStack(); }

Grid ElementalDefaultGrid( MPI_Comm comm )
{
    Grid g = ::gridList.size();
    ::gridList.push_back( new elemental::Grid( comm ) );
    return g;
}

Grid ElementalGrid( MPI_Comm comm, int r, int c )
{
    Grid g = ::gridList.size();
    ::gridList.push_back( new elemental::Grid( comm, r, c ) );
    return g;
}

// Don't yet export the viewingComm/owningGroup constructors.

int ElementalGridHeight( Grid g )
{ return gridList[g]->Height(); }

int ElementalGridWidth( Grid g )
{ return gridList[g]->Width(); }

int ElementalGridSize( Grid g )
{ return gridList[g]->Size(); }

int ElementalInGrid( Grid g )
{ return gridList[g]->InGrid(); }

int ElementalGridVCRank( Grid g )
{ return gridList[g]->VCRank(); }

int ElementalGridVRRank( Grid g )
{ return gridList[g]->VRRank(); }

int ElementalGridMCRank( Grid g )
{ return gridList[g]->MCRank(); }

int ElementalGridMRRank( Grid g )
{ return gridList[g]->MRRank(); }

MPI_Comm ElementalGridVCComm( Grid g )
{ return gridList[g]->VCComm(); }

MPI_Comm ElementalGridVRComm( Grid g )
{ return gridList[g]->VRComm(); }

MPI_Comm ElementalGridMCComm( Grid g )
{ return gridList[g]->MCComm(); }

MPI_Comm ElementalGridMRComm( Grid g )
{ return gridList[g]->MRComm(); }

MC_MR_Single
ElementalCreateEmpty_MC_MR_Single( Grid g )
{
    MC_MR_Single handle = ::MC_MR_SingleList.size();
    ::MC_MR_SingleList.push_back
        ( new elemental::DistMatrix<float,elemental::MC,elemental::MR>
          ( *::gridList[g] ) );
    return handle;
}

MC_MR_Double
ElementalCreateEmpty_MC_MR_Double( Grid g )
{
    MC_MR_Double handle = ::MC_MR_DoubleList.size();
    ::MC_MR_DoubleList.push_back
        ( new elemental::DistMatrix<double,elemental::MC,elemental::MR>
          ( *::gridList[g] ) );
    return handle;
}

Star_VR_Single
ElementalCreateEmpty_Star_VR_Single( Grid g )
{
    Star_VR_Single handle = ::Star_VR_SingleList.size();
    ::Star_VR_SingleList.push_back
        ( new elemental::DistMatrix<float,elemental::Star,elemental::VR>
          ( *::gridList[g] ) );
    return handle;
}

Star_VR_Double
ElementalCreateEmpty_Star_VR_Double( Grid g )
{
    Star_VR_Double handle = ::Star_VR_DoubleList.size();
    ::Star_VR_DoubleList.push_back
        ( new elemental::DistMatrix<double,elemental::Star,elemental::VR>
          ( *::gridList[g] ) );
    return handle;
}

MC_MR_Single
ElementalRegister_MC_MR_Single
( int height, int width, int colAlignment, int rowAlignment,
  float* buffer, int ldim, Grid g )
{
    MC_MR_Single handle = ::MC_MR_SingleList.size();
    ::MC_MR_SingleList.push_back
        ( new elemental::DistMatrix<float,elemental::MC,elemental::MR>
            ( height, width, colAlignment, rowAlignment, buffer, ldim, 
              *::gridList[g] ) );
    return handle;
}

MC_MR_Double
ElementalRegister_MC_MR_Double
( int height, int width, int colAlignment, int rowAlignment,
  double* buffer, int ldim, Grid g )
{
    MC_MR_Double handle = ::MC_MR_DoubleList.size();
    ::MC_MR_DoubleList.push_back
        ( new elemental::DistMatrix<double,elemental::MC,elemental::MR>
            ( height, width, colAlignment, rowAlignment, buffer, ldim, 
              *::gridList[g] ) );
    return handle;
}

void 
ElementalPrint_MC_MR_Single( char* msg, MC_MR_Single handle )
{
    std::string s( msg );
    ::MC_MR_SingleList[handle]->Print( s );
}

void 
ElementalPrint_MC_MR_Double( char* msg, MC_MR_Double handle )
{
    std::string s( msg );
    ::MC_MR_DoubleList[handle]->Print( s );
}

void 
ElementalPrint_Star_VR_Single( char* msg, Star_VR_Single handle )
{
    std::string s( msg );
    ::Star_VR_SingleList[handle]->Print( s );
}

void 
ElementalPrint_Star_VR_Double( char* msg, Star_VR_Double handle )
{
    std::string s( msg );
    ::Star_VR_DoubleList[handle]->Print( s );
}

#ifndef WITHOUT_COMPLEX
MC_MR_SComplex
ElementalCreateEmpty_MC_MR_SComplex( Grid g )
{
    MC_MR_SComplex handle = ::MC_MR_SComplexList.size();
    ::MC_MR_SComplexList.push_back
        ( new elemental::DistMatrix<elemental::scomplex,
                                    elemental::MC,elemental::MR>
          ( *::gridList[g] ) );
    return handle;
}

MC_MR_DComplex
ElementalCreateEmpty_MC_MR_DComplex( Grid g )
{
    MC_MR_DComplex handle = ::MC_MR_DComplexList.size();
    ::MC_MR_DComplexList.push_back
        ( new elemental::DistMatrix<elemental::dcomplex,
                                    elemental::MC,elemental::MR>
          ( *::gridList[g] ) );
    return handle;
}

Star_VR_SComplex
ElementalCreateEmpty_Star_VR_SComplex( Grid g )
{
    Star_VR_SComplex handle = ::Star_VR_SComplexList.size();
    ::Star_VR_SComplexList.push_back
        ( new elemental::DistMatrix<elemental::scomplex,
                                    elemental::Star,elemental::VR>
          ( *::gridList[g] ) );
    return handle;
}

Star_VR_DComplex
ElementalCreateEmpty_Star_VR_DComplex( Grid g )
{
    Star_VR_DComplex handle = ::Star_VR_DComplexList.size();
    ::Star_VR_DComplexList.push_back
        ( new elemental::DistMatrix<elemental::dcomplex,
                                    elemental::Star,elemental::VR>
          ( *::gridList[g] ) );
    return handle;
}

MC_MR_SComplex
ElementalRegister_MC_MR_SComplex
( int height, int width, int colAlignment, int rowAlignment,
  SComplex* buffer, int ldim, Grid g )
{
    MC_MR_SComplex handle = ::MC_MR_SComplexList.size();
    ::MC_MR_SComplexList.push_back
        ( new elemental::DistMatrix<elemental::scomplex,
                                    elemental::MC,elemental::MR>
            ( height, width, colAlignment, rowAlignment, 
              (elemental::scomplex*)buffer, ldim, *::gridList[g] ) );
    return handle;
}

MC_MR_DComplex
ElementalRegister_MC_MR_DComplex
( int height, int width, int colAlignment, int rowAlignment,
  DComplex* buffer, int ldim, Grid g )
{
    MC_MR_DComplex handle = ::MC_MR_DComplexList.size();
    ::MC_MR_DComplexList.push_back
        ( new elemental::DistMatrix<elemental::dcomplex,
                                    elemental::MC,elemental::MR>
            ( height, width, colAlignment, rowAlignment, 
              (elemental::dcomplex*)buffer, ldim, *::gridList[g] ) );
    return handle;
}

Star_VR_SComplex
ElementalRegister_Star_VR_SComplex
( int height, int width, int rowAlignment,
  SComplex* buffer, int ldim, Grid g )
{
    Star_VR_SComplex handle = ::Star_VR_SComplexList.size();
    ::Star_VR_SComplexList.push_back
        ( new elemental::DistMatrix<elemental::scomplex,
                                    elemental::Star,elemental::VR>
            ( height, width, rowAlignment, 
              (elemental::scomplex*)buffer, ldim, *::gridList[g] ) );
    return handle;
}

Star_VR_DComplex
ElementalRegister_Star_VR_DComplex
( int height, int width, int rowAlignment,
  DComplex* buffer, int ldim, Grid g )
{
    Star_VR_DComplex handle = ::Star_VR_DComplexList.size();
    ::Star_VR_DComplexList.push_back
        ( new elemental::DistMatrix<elemental::dcomplex,
                                    elemental::Star,elemental::VR>
            ( height, width, rowAlignment, 
              (elemental::dcomplex*)buffer, ldim, *::gridList[g] ) );
    return handle;
}

void 
ElementalPrint_MC_MR_SComplex( char* msg, MC_MR_SComplex handle )
{
    std::string s( msg );
    ::MC_MR_SComplexList[handle]->Print( s );
}

void 
ElementalPrint_MC_MR_DComplex( char* msg, MC_MR_DComplex handle )
{
    std::string s( msg );
    ::MC_MR_DComplexList[handle]->Print( s );
}

void 
ElementalPrint_Star_VR_SComplex( char* msg, Star_VR_SComplex handle )
{
    std::string s( msg );
    ::Star_VR_SComplexList[handle]->Print( s );
}

void 
ElementalPrint_Star_VR_DComplex( char* msg, Star_VR_DComplex handle )
{
    std::string s( msg );
    ::Star_VR_DComplexList[handle]->Print( s );
}
#endif // WITHOUT_COMPLEX

int
ElementalLocalLength
( int globalLength, int myIndex, int alignment, int modulus )
{ 
    return elemental::utilities::LocalLength
        ( globalLength, myIndex, alignment, modulus );
}

/* LAPACK-level interface */
#ifndef WITHOUT_PMRRR
void
ElementalHermitianEigDouble
( char uplo,
  MC_MR_Double AHandle, Star_VR_Double wHandle, MC_MR_Double ZHandle,
  int tryForHighAccuracy )
{
    elemental::Shape shape = elemental::CharToShape( uplo );
    elemental::lapack::HermitianEig
    ( shape, *::MC_MR_DoubleList[AHandle], *::Star_VR_DoubleList[wHandle],
      *::MC_MR_DoubleList[ZHandle], tryForHighAccuracy );
}

void
ElementalHermitianEigDComplex
( char uplo,
  MC_MR_DComplex AHandle, Star_VR_Double wHandle, MC_MR_DComplex ZHandle,
  int tryForHighAccuracy )
{
    elemental::Shape shape = elemental::CharToShape( uplo );
    elemental::lapack::HermitianEig
    ( shape, *::MC_MR_DComplexList[AHandle], *::Star_VR_DoubleList[wHandle],
      *::MC_MR_DComplexList[ZHandle], tryForHighAccuracy );
}
#endif /* WITHOUT_PMRRR */

} // extern "C"

