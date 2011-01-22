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

#ifndef RELEASE
#define CATCH(expr) try { expr; } catch( std::exception& e ) \
{ \
    elemental::DumpCallStack(); \
    std::cerr << "Caught error message:\n" << e.what() << std::endl; \
    ElementalFinalize(); \
}
#else
#define CATCH(expr) try { expr; } catch( std::exception& e ) \
{ \
    std::cerr << "Caught error message:\n" << e.what() << std::endl; \
    ElementalFinalize(); \
}
#endif

extern "C" {

/* TODO: Add wrapper macro for handling exceptions. */

void ElementalInit( int* argc, char** argv[] )
{ elemental::Init( argc, argv ); }

void ElementalFinalize()
{
    for( unsigned j=0; j< ::MC_MR_SingleList.size(); ++j )
        delete ::MC_MR_SingleList[j];
    for( unsigned j=0; j< ::MC_MR_DoubleList.size(); ++j )
        delete ::MC_MR_DoubleList[j];
    for( unsigned j=0; j< ::Star_VR_SingleList.size(); ++j )
        delete ::Star_VR_SingleList[j];
    for( unsigned j=0; j< ::Star_VR_DoubleList.size(); ++j )
        delete ::Star_VR_DoubleList[j];
    /* TODO: Add other distributions */

#ifndef WITHOUT_COMPLEX
    for( unsigned j=0; j< ::MC_MR_SComplexList.size(); ++j )
        delete ::MC_MR_SComplexList[j];
    for( unsigned j=0; j< ::MC_MR_DComplexList.size(); ++j )
        delete ::MC_MR_DComplexList[j];
    for( unsigned j=0; j< ::Star_VR_SComplexList.size(); ++j )
        delete ::Star_VR_SComplexList[j];
    for( unsigned j=0; j< ::Star_VR_DComplexList.size(); ++j )
        delete ::Star_VR_DComplexList[j];
    /* TODO: Add other distributions */
#endif

    for( unsigned j=0; j< ::gridList.size(); ++j )
        delete ::gridList[j];

    elemental::Finalize();
}

int Blocksize()
{ return elemental::Blocksize(); }

void SetBlocksize( int blocksize )
{ elemental::SetBlocksize( blocksize ); }

void PushBlocksizeStack( int blocksize )
{ elemental::PushBlocksizeStack( blocksize ); }

void PopBlocksizeStack()
{ CATCH(elemental::PopBlocksizeStack()); }

Grid CreateDefaultGrid( MPI_Comm comm )
{
    Grid g;
    CATCH(
        g = ::gridList.size();
        ::gridList.push_back( new elemental::Grid( comm ) );
    );
    return g;
}

Grid CreateGrid( MPI_Comm comm, int r, int c )
{
    Grid g;
    CATCH(
        g = ::gridList.size();
        ::gridList.push_back( new elemental::Grid( comm, r, c ) );
    );
    return g;
}

// Don't yet export the viewingComm/owningGroup constructors.

int GridHeight( Grid g )
{ int r; CATCH(r=gridList[g]->Height()); return r; }

int GridWidth( Grid g )
{ int c; CATCH(c=gridList[g]->Width()); return c; }

int GridSize( Grid g )
{ int p; CATCH(p=gridList[g]->Size()); return p; }

int InGrid( Grid g )
{ int inGrid; CATCH(inGrid=gridList[g]->InGrid()); return inGrid; }

int GridVCRank( Grid g )
{ int VCRank; CATCH(VCRank=gridList[g]->VCRank()); return VCRank; }

int GridVRRank( Grid g )
{ int VRRank; CATCH(VRRank=gridList[g]->VRRank()); return VRRank; }

int GridMCRank( Grid g )
{ int MCRank; CATCH(MCRank=gridList[g]->MCRank()); return MCRank; }

int GridMRRank( Grid g )
{ int MRRank; CATCH(MRRank=gridList[g]->MRRank()); return MRRank; }

MPI_Comm GridVCComm( Grid g )
{ MPI_Comm VCComm; CATCH(VCComm=gridList[g]->VCComm()); return VCComm; }

MPI_Comm GridVRComm( Grid g )
{ MPI_Comm VRComm; CATCH(VRComm=gridList[g]->VRComm()); return VRComm; }

MPI_Comm GridMCComm( Grid g )
{ MPI_Comm MCComm; CATCH(MCComm=gridList[g]->MCComm()); return MCComm; }

MPI_Comm GridMRComm( Grid g )
{ MPI_Comm MRComm; CATCH(MRComm=gridList[g]->MRComm()); return MRComm; }

MC_MR_Single
CreateEmpty_MC_MR_Single( Grid g )
{
    MC_MR_Single A;
    CATCH(
        A = ::MC_MR_SingleList.size();
        ::MC_MR_SingleList.push_back
            ( new elemental::DistMatrix<float,elemental::MC,elemental::MR>
              ( *::gridList[g] ) );
    );
    return A;
}

MC_MR_Double
CreateEmpty_MC_MR_Double( Grid g )
{
    MC_MR_Double A;
    CATCH(
        A = ::MC_MR_DoubleList.size();
        ::MC_MR_DoubleList.push_back
            ( new elemental::DistMatrix<double,elemental::MC,elemental::MR>
              ( *::gridList[g] ) );
    );
    return A;
}

Star_VR_Single
CreateEmpty_Star_VR_Single( Grid g )
{
    Star_VR_Single A;
    CATCH(
        A = ::Star_VR_SingleList.size();
        ::Star_VR_SingleList.push_back
            ( new elemental::DistMatrix<float,elemental::Star,elemental::VR>
              ( *::gridList[g] ) );
    );
    return A;
}

Star_VR_Double
CreateEmpty_Star_VR_Double( Grid g )
{
    Star_VR_Double A;
    CATCH(
        A = ::Star_VR_DoubleList.size();
        ::Star_VR_DoubleList.push_back
            ( new elemental::DistMatrix<double,elemental::Star,elemental::VR>
              ( *::gridList[g] ) );
    );
    return A;
}

MC_MR_Single
Register_MC_MR_Single
( int height, int width, int colAlignment, int rowAlignment,
  float* buffer, int ldim, Grid g )
{
    MC_MR_Single A;
    CATCH(
        A = ::MC_MR_SingleList.size();
        ::MC_MR_SingleList.push_back
            ( new elemental::DistMatrix<float,elemental::MC,elemental::MR>
                ( height, width, colAlignment, rowAlignment, buffer, ldim, 
                  *::gridList[g] ) );
    );
    return A;
}

MC_MR_Double
Register_MC_MR_Double
( int height, int width, int colAlignment, int rowAlignment,
  double* buffer, int ldim, Grid g )
{
    MC_MR_Double A;
    CATCH(
        A = ::MC_MR_DoubleList.size();
        ::MC_MR_DoubleList.push_back
            ( new elemental::DistMatrix<double,elemental::MC,elemental::MR>
                ( height, width, colAlignment, rowAlignment, buffer, ldim, 
                  *::gridList[g] ) );
    );
    return A;
}

void 
Print_MC_MR_Single( char* msg, MC_MR_Single A )
{
    CATCH(
        std::string s( msg );
        ::MC_MR_SingleList[A]->Print( s );
    );
}

void 
Print_MC_MR_Double( char* msg, MC_MR_Double A )
{
    CATCH(
        std::string s( msg );
        ::MC_MR_DoubleList[A]->Print( s );
    );
}

void 
Print_Star_VR_Single( char* msg, Star_VR_Single A )
{
    CATCH(
        std::string s( msg );
        ::Star_VR_SingleList[A]->Print( s );
    );
}

void 
Print_Star_VR_Double( char* msg, Star_VR_Double A )
{
    CATCH(
        std::string s( msg );
        ::Star_VR_DoubleList[A]->Print( s );
    );
}

#ifndef WITHOUT_COMPLEX
MC_MR_SComplex
CreateEmpty_MC_MR_SComplex( Grid g )
{
    MC_MR_SComplex A;
    CATCH(
        A = ::MC_MR_SComplexList.size();
        ::MC_MR_SComplexList.push_back
            ( new elemental::DistMatrix<elemental::scomplex,
                                        elemental::MC,elemental::MR>
              ( *::gridList[g] ) );
    );
    return A;
}

MC_MR_DComplex
CreateEmpty_MC_MR_DComplex( Grid g )
{
    MC_MR_DComplex A;
    CATCH(
        A = ::MC_MR_DComplexList.size();
        ::MC_MR_DComplexList.push_back
            ( new elemental::DistMatrix<elemental::dcomplex,
                                        elemental::MC,elemental::MR>
              ( *::gridList[g] ) );
    );
    return A;
}

Star_VR_SComplex
CreateEmpty_Star_VR_SComplex( Grid g )
{
    Star_VR_SComplex A;
    CATCH(
        A = ::Star_VR_SComplexList.size();
        ::Star_VR_SComplexList.push_back
            ( new elemental::DistMatrix<elemental::scomplex,
                                        elemental::Star,elemental::VR>
              ( *::gridList[g] ) );
    );
    return A;
}

Star_VR_DComplex
CreateEmpty_Star_VR_DComplex( Grid g )
{
    Star_VR_DComplex A;
    CATCH(
        A = ::Star_VR_DComplexList.size();
        ::Star_VR_DComplexList.push_back
            ( new elemental::DistMatrix<elemental::dcomplex,
                                        elemental::Star,elemental::VR>
              ( *::gridList[g] ) );
    );
    return A;
}

MC_MR_SComplex
Register_MC_MR_SComplex
( int height, int width, int colAlignment, int rowAlignment,
  SComplex* buffer, int ldim, Grid g )
{
    MC_MR_SComplex A;
    CATCH(
        A = ::MC_MR_SComplexList.size();
        ::MC_MR_SComplexList.push_back
            ( new elemental::DistMatrix<elemental::scomplex,
                                        elemental::MC,elemental::MR>
                ( height, width, colAlignment, rowAlignment, 
                  (elemental::scomplex*)buffer, ldim, *::gridList[g] ) );
    );
    return A;
}

MC_MR_DComplex
Register_MC_MR_DComplex
( int height, int width, int colAlignment, int rowAlignment,
  DComplex* buffer, int ldim, Grid g )
{
    MC_MR_DComplex A;
    CATCH(
        A = ::MC_MR_DComplexList.size();
        ::MC_MR_DComplexList.push_back
            ( new elemental::DistMatrix<elemental::dcomplex,
                                        elemental::MC,elemental::MR>
                ( height, width, colAlignment, rowAlignment, 
                  (elemental::dcomplex*)buffer, ldim, *::gridList[g] ) );
    );
    return A;
}

Star_VR_SComplex
Register_Star_VR_SComplex
( int height, int width, int rowAlignment,
  SComplex* buffer, int ldim, Grid g )
{
    Star_VR_SComplex A;
    CATCH(
        A = ::Star_VR_SComplexList.size();
        ::Star_VR_SComplexList.push_back
            ( new elemental::DistMatrix<elemental::scomplex,
                                        elemental::Star,elemental::VR>
                ( height, width, rowAlignment, 
                  (elemental::scomplex*)buffer, ldim, *::gridList[g] ) );
    );
    return A;
}

Star_VR_DComplex
Register_Star_VR_DComplex
( int height, int width, int rowAlignment,
  DComplex* buffer, int ldim, Grid g )
{
    Star_VR_DComplex A;
    CATCH(
        A = ::Star_VR_DComplexList.size();
        ::Star_VR_DComplexList.push_back
            ( new elemental::DistMatrix<elemental::dcomplex,
                                        elemental::Star,elemental::VR>
                ( height, width, rowAlignment, 
                  (elemental::dcomplex*)buffer, ldim, *::gridList[g] ) );
    );
    return A;
}

void 
Print_MC_MR_SComplex( char* msg, MC_MR_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::MC_MR_SComplexList[A]->Print( s );
    );
}

void 
Print_MC_MR_DComplex( char* msg, MC_MR_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::MC_MR_DComplexList[A]->Print( s );
    );
}

void 
Print_Star_VR_SComplex( char* msg, Star_VR_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::Star_VR_SComplexList[A]->Print( s );
    );
}

void 
Print_Star_VR_DComplex( char* msg, Star_VR_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::Star_VR_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX

int
LocalLength
( int globalLength, int myIndex, int alignment, int modulus )
{ 
    int localLength;
    CATCH(
        localLength = elemental::utilities::LocalLength
        ( globalLength, myIndex, alignment, modulus );
    );
    return localLength;
}

/* LAPACK-level interface */
#ifndef WITHOUT_PMRRR
void
HermitianEigDouble
( char uplo, MC_MR_Double A, Star_VR_Double w, MC_MR_Double Z,
  int tryForHighAccuracy )
{
    CATCH(
        elemental::Shape shape = elemental::CharToShape( uplo );
        elemental::lapack::HermitianEig
        ( shape, *::MC_MR_DoubleList[A], *::Star_VR_DoubleList[w],
          *::MC_MR_DoubleList[Z], tryForHighAccuracy );
    );
}

void
HermitianEigDouble_IntegerSubset
( char uplo, MC_MR_Double A, Star_VR_Double w, MC_MR_Double Z,
  int a, int b, int tryForHighAccuracy )
{
    CATCH(
        elemental::Shape shape = elemental::CharToShape( uplo );
        elemental::lapack::HermitianEig
        ( shape, *::MC_MR_DoubleList[A], *::Star_VR_DoubleList[w],
          *::MC_MR_DoubleList[Z], a, b, tryForHighAccuracy );
    );
}

void
HermitianEigDouble_IntervalSubset
( char uplo, MC_MR_Double A, Star_VR_Double w, MC_MR_Double Z,
  double a, double b, int tryForHighAccuracy )
{
    CATCH(
        elemental::Shape shape = elemental::CharToShape( uplo );
        elemental::lapack::HermitianEig
        ( shape, *::MC_MR_DoubleList[A], *::Star_VR_DoubleList[w],
          *::MC_MR_DoubleList[Z], a, b, tryForHighAccuracy );
    );
}

void
HermitianEigDouble_OnlyEigvals
( char uplo, MC_MR_Double A, Star_VR_Double w, 
  int tryForHighAccuracy )
{
    CATCH(
        elemental::Shape shape = elemental::CharToShape( uplo );
        elemental::lapack::HermitianEig
        ( shape, *::MC_MR_DoubleList[A], *::Star_VR_DoubleList[w],
          tryForHighAccuracy );
    );
}

void
HermitianEigDouble_OnlyEigvals_IntegerSubset
( char uplo, MC_MR_Double A, Star_VR_Double w, 
  int a, int b, int tryForHighAccuracy )
{
    CATCH(
        elemental::Shape shape = elemental::CharToShape( uplo );
        elemental::lapack::HermitianEig
        ( shape, *::MC_MR_DoubleList[A], *::Star_VR_DoubleList[w],
          a, b, tryForHighAccuracy );
    );
}

void
HermitianEigDouble_OnlyEigvals_IntervalSubset
( char uplo, MC_MR_Double A, Star_VR_Double w, 
  double a, double b, int tryForHighAccuracy )
{
    CATCH(
        elemental::Shape shape = elemental::CharToShape( uplo );
        elemental::lapack::HermitianEig
        ( shape, *::MC_MR_DoubleList[A], *::Star_VR_DoubleList[w],
          a, b, tryForHighAccuracy );
    );
}

#ifndef WITHOUT_COMPLEX
void
HermitianEigDComplex
( char uplo,
  MC_MR_DComplex A, Star_VR_Double w, MC_MR_DComplex Z,
  int tryForHighAccuracy )
{
    CATCH(
        elemental::Shape shape = elemental::CharToShape( uplo );
        elemental::lapack::HermitianEig
        ( shape, *::MC_MR_DComplexList[A], *::Star_VR_DoubleList[w],
          *::MC_MR_DComplexList[Z], tryForHighAccuracy );
    );
}

void
HermitianEigDComplex_IntegerSubset
( char uplo,
  MC_MR_DComplex A, Star_VR_Double w, MC_MR_DComplex Z,
  int a, int b, int tryForHighAccuracy )
{
    CATCH(
        elemental::Shape shape = elemental::CharToShape( uplo );
        elemental::lapack::HermitianEig
        ( shape, *::MC_MR_DComplexList[A], *::Star_VR_DoubleList[w],
          *::MC_MR_DComplexList[Z], a, b, tryForHighAccuracy );
    );
}

void
HermitianEigDComplex_IntervalSubset
( char uplo,
  MC_MR_DComplex A, Star_VR_Double w, MC_MR_DComplex Z,
  double a, double b, int tryForHighAccuracy )
{
    CATCH(
        elemental::Shape shape = elemental::CharToShape( uplo );
        elemental::lapack::HermitianEig
        ( shape, *::MC_MR_DComplexList[A], *::Star_VR_DoubleList[w],
          *::MC_MR_DComplexList[Z], a, b, tryForHighAccuracy );
    );
}

void
HermitianEigDComplex_OnlyEigvals
( char uplo,
  MC_MR_DComplex A, Star_VR_Double w, 
  int tryForHighAccuracy )
{
    CATCH(
        elemental::Shape shape = elemental::CharToShape( uplo );
        elemental::lapack::HermitianEig
        ( shape, *::MC_MR_DComplexList[A], *::Star_VR_DoubleList[w],
          tryForHighAccuracy );
    );
}

void
HermitianEigDComplex_OnlyEigvals_IntegerSubset
( char uplo,
  MC_MR_DComplex A, Star_VR_Double w, 
  int a, int b, int tryForHighAccuracy )
{
    CATCH(
        elemental::Shape shape = elemental::CharToShape( uplo );
        elemental::lapack::HermitianEig
        ( shape, *::MC_MR_DComplexList[A], *::Star_VR_DoubleList[w],
          tryForHighAccuracy );
    );
}

void
HermitianEigDComplex_OnlyEigvals_IntervalSubset
( char uplo,
  MC_MR_DComplex A, Star_VR_Double w, 
  double a, double b, int tryForHighAccuracy )
{
    CATCH(
        elemental::Shape shape = elemental::CharToShape( uplo );
        elemental::lapack::HermitianEig
        ( shape, *::MC_MR_DComplexList[A], *::Star_VR_DoubleList[w],
          tryForHighAccuracy );
    );
}
#endif /* WITHOUT_COMPLEX */
#endif /* WITHOUT_PMRRR */

} // extern "C"

