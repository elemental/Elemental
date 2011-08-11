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

using elemental::MC;
using elemental::MD;
using elemental::MR;
using elemental::STAR;
using elemental::VC;
using elemental::VR;
using elemental::scomplex;
using elemental::dcomplex;
using elemental::DistMatrix;

namespace {
std::vector<elemental::Grid*> gridList;
std::vector<DistMatrix<float, MC,  MR  >*> MC_MR_SingleList;
std::vector<DistMatrix<double,MC,  MR  >*> MC_MR_DoubleList;
std::vector<DistMatrix<float, MC,  STAR>*> MC_STAR_SingleList;
std::vector<DistMatrix<double,MC,  STAR>*> MC_STAR_DoubleList;
std::vector<DistMatrix<float, MD,  STAR>*> MD_STAR_SingleList;
std::vector<DistMatrix<double,MD,  STAR>*> MD_STAR_DoubleList;
std::vector<DistMatrix<float, MR,  MC  >*> MR_MC_SingleList;
std::vector<DistMatrix<double,MR,  MC  >*> MR_MC_DoubleList;
std::vector<DistMatrix<float, MR,  STAR>*> MR_STAR_SingleList;
std::vector<DistMatrix<double,MR,  STAR>*> MR_STAR_DoubleList;
std::vector<DistMatrix<float, STAR,MC  >*> STAR_MC_SingleList;
std::vector<DistMatrix<double,STAR,MC  >*> STAR_MC_DoubleList;
std::vector<DistMatrix<float, STAR,MD  >*> STAR_MD_SingleList;
std::vector<DistMatrix<double,STAR,MD  >*> STAR_MD_DoubleList;
std::vector<DistMatrix<float, STAR,MR  >*> STAR_MR_SingleList;
std::vector<DistMatrix<double,STAR,MR  >*> STAR_MR_DoubleList;
std::vector<DistMatrix<float, STAR,STAR>*> STAR_STAR_SingleList;
std::vector<DistMatrix<double,STAR,STAR>*> STAR_STAR_DoubleList;
std::vector<DistMatrix<float, STAR,VC  >*> STAR_VC_SingleList;
std::vector<DistMatrix<double,STAR,VC  >*> STAR_VC_DoubleList;
std::vector<DistMatrix<float, STAR,VR  >*> STAR_VR_SingleList;
std::vector<DistMatrix<double,STAR,VR  >*> STAR_VR_DoubleList;
std::vector<DistMatrix<float, VC,  STAR>*> VC_STAR_SingleList;
std::vector<DistMatrix<double,VC,  STAR>*> VC_STAR_DoubleList;
std::vector<DistMatrix<float, VR,  STAR>*> VR_STAR_SingleList;
std::vector<DistMatrix<double,VR,  STAR>*> VR_STAR_DoubleList;
#ifndef WITHOUT_COMPLEX
std::vector<DistMatrix<scomplex,MC,  MR  >*> MC_MR_SComplexList;
std::vector<DistMatrix<dcomplex,MC,  MR  >*> MC_MR_DComplexList;
std::vector<DistMatrix<scomplex,MC,  STAR>*> MC_STAR_SComplexList;
std::vector<DistMatrix<dcomplex,MC,  STAR>*> MC_STAR_DComplexList;
std::vector<DistMatrix<scomplex,MD,  STAR>*> MD_STAR_SComplexList;
std::vector<DistMatrix<dcomplex,MD,  STAR>*> MD_STAR_DComplexList;
std::vector<DistMatrix<scomplex,MR,  MC  >*> MR_MC_SComplexList;
std::vector<DistMatrix<dcomplex,MR,  MC  >*> MR_MC_DComplexList;
std::vector<DistMatrix<scomplex,MR,  STAR>*> MR_STAR_SComplexList;
std::vector<DistMatrix<dcomplex,MR,  STAR>*> MR_STAR_DComplexList;
std::vector<DistMatrix<scomplex,STAR,MC  >*> STAR_MC_SComplexList;
std::vector<DistMatrix<dcomplex,STAR,MC  >*> STAR_MC_DComplexList;
std::vector<DistMatrix<scomplex,STAR,MD  >*> STAR_MD_SComplexList;
std::vector<DistMatrix<dcomplex,STAR,MD  >*> STAR_MD_DComplexList;
std::vector<DistMatrix<scomplex,STAR,MR  >*> STAR_MR_SComplexList;
std::vector<DistMatrix<dcomplex,STAR,MR  >*> STAR_MR_DComplexList;
std::vector<DistMatrix<scomplex,STAR,STAR>*> STAR_STAR_SComplexList;
std::vector<DistMatrix<dcomplex,STAR,STAR>*> STAR_STAR_DComplexList;
std::vector<DistMatrix<scomplex,STAR,VC  >*> STAR_VC_SComplexList;
std::vector<DistMatrix<dcomplex,STAR,VC  >*> STAR_VC_DComplexList;
std::vector<DistMatrix<scomplex,STAR,VR  >*> STAR_VR_SComplexList;
std::vector<DistMatrix<dcomplex,STAR,VR  >*> STAR_VR_DComplexList;
std::vector<DistMatrix<scomplex,VC,  STAR>*> VC_STAR_SComplexList;
std::vector<DistMatrix<dcomplex,VC,  STAR>*> VC_STAR_DComplexList;
std::vector<DistMatrix<scomplex,VR,  STAR>*> VR_STAR_SComplexList;
std::vector<DistMatrix<dcomplex,VR,  STAR>*> VR_STAR_DComplexList;
#endif // WITHOUT_COMPLEX
} // anonymous namespace

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

void ElementalInit( int* argc, char** argv[] )
{ elemental::Init( *argc, *argv ); }

void ElementalFinalize()
{
    // Real [MC,MR]
    for( unsigned j=0; j< ::MC_MR_SingleList.size(); ++j )
        delete ::MC_MR_SingleList[j];
    for( unsigned j=0; j< ::MC_MR_DoubleList.size(); ++j )
        delete ::MC_MR_DoubleList[j];
    ::MC_MR_SingleList.resize(0);
    ::MC_MR_DoubleList.resize(0);
    // Real [MC,* ]
    for( unsigned j=0; j< ::MC_STAR_SingleList.size(); ++j )
        delete ::MC_STAR_SingleList[j];
    for( unsigned j=0; j< ::MC_STAR_DoubleList.size(); ++j )
        delete ::MC_STAR_DoubleList[j];
    ::MC_STAR_SingleList.resize(0);
    ::MC_STAR_DoubleList.resize(0);
    // Real [MD,* ]
    for( unsigned j=0; j< ::MD_STAR_SingleList.size(); ++j )
        delete ::MD_STAR_SingleList[j];
    for( unsigned j=0; j< ::MD_STAR_DoubleList.size(); ++j )
        delete ::MD_STAR_DoubleList[j];
    ::MD_STAR_SingleList.resize(0);
    ::MD_STAR_DoubleList.resize(0);
    // Real [MR,MC]
    for( unsigned j=0; j< ::MR_MC_SingleList.size(); ++j )
        delete ::MR_MC_SingleList[j];
    for( unsigned j=0; j< ::MR_MC_DoubleList.size(); ++j )
        delete ::MR_MC_DoubleList[j];
    ::MR_MC_SingleList.resize(0);
    ::MR_MC_DoubleList.resize(0);
    // Real [MR,* ]
    for( unsigned j=0; j< ::MR_STAR_SingleList.size(); ++j )
        delete ::MR_STAR_SingleList[j];
    for( unsigned j=0; j< ::MR_STAR_DoubleList.size(); ++j )
        delete ::MR_STAR_DoubleList[j];
    ::MR_STAR_SingleList.resize(0);
    ::MR_STAR_DoubleList.resize(0);
    // Real [* ,MC]
    for( unsigned j=0; j< ::STAR_MC_SingleList.size(); ++j )
        delete ::STAR_MC_SingleList[j];
    for( unsigned j=0; j< ::STAR_MC_DoubleList.size(); ++j )
        delete ::STAR_MC_DoubleList[j];
    ::STAR_MC_SingleList.resize(0);
    ::STAR_MC_DoubleList.resize(0);
    // Real [* ,MD]
    for( unsigned j=0; j< ::STAR_MD_SingleList.size(); ++j )
        delete ::STAR_MD_SingleList[j];
    for( unsigned j=0; j< ::STAR_MD_DoubleList.size(); ++j )
        delete ::STAR_MD_DoubleList[j];
    ::STAR_MD_SingleList.resize(0);
    ::STAR_MD_DoubleList.resize(0);
    // Real [* ,MR]
    for( unsigned j=0; j< ::STAR_MR_SingleList.size(); ++j )
        delete ::STAR_MR_SingleList[j];
    for( unsigned j=0; j< ::STAR_MR_DoubleList.size(); ++j )
        delete ::STAR_MR_DoubleList[j];
    ::STAR_MR_SingleList.resize(0);
    ::STAR_MR_DoubleList.resize(0);
    // Real [* ,* ]
    for( unsigned j=0; j< ::STAR_STAR_SingleList.size(); ++j )
        delete ::STAR_STAR_SingleList[j];
    for( unsigned j=0; j< ::STAR_STAR_DoubleList.size(); ++j )
        delete ::STAR_STAR_DoubleList[j];
    ::STAR_STAR_SingleList.resize(0);
    ::STAR_STAR_DoubleList.resize(0);
    // Real [* ,VC]
    for( unsigned j=0; j< ::STAR_VC_SingleList.size(); ++j )
        delete ::STAR_VC_SingleList[j];
    for( unsigned j=0; j< ::STAR_VC_DoubleList.size(); ++j )
        delete ::STAR_VC_DoubleList[j];
    ::STAR_VC_SingleList.resize(0);
    ::STAR_VC_DoubleList.resize(0);
    // Real [* ,VR]
    for( unsigned j=0; j< ::STAR_VR_SingleList.size(); ++j )
        delete ::STAR_VR_SingleList[j];
    for( unsigned j=0; j< ::STAR_VR_DoubleList.size(); ++j )
        delete ::STAR_VR_DoubleList[j];
    ::STAR_VR_SingleList.resize(0);
    ::STAR_VR_DoubleList.resize(0);
    // Real [VC,* ]
    for( unsigned j=0; j< ::VC_STAR_SingleList.size(); ++j )
        delete ::VC_STAR_SingleList[j];
    for( unsigned j=0; j< ::VC_STAR_DoubleList.size(); ++j )
        delete ::VC_STAR_DoubleList[j];
    ::VC_STAR_SingleList.resize(0);
    ::VC_STAR_DoubleList.resize(0);
    // Real [VR,* ]
    for( unsigned j=0; j< ::VR_STAR_SingleList.size(); ++j )
        delete ::VR_STAR_SingleList[j];
    for( unsigned j=0; j< ::VR_STAR_DoubleList.size(); ++j )
        delete ::VR_STAR_DoubleList[j];
    ::VR_STAR_SingleList.resize(0);
    ::VR_STAR_DoubleList.resize(0);
#ifndef WITHOUT_COMPLEX
    // Complex [MC,MR]
    for( unsigned j=0; j< ::MC_MR_SComplexList.size(); ++j )
        delete ::MC_MR_SComplexList[j];
    for( unsigned j=0; j< ::MC_MR_DComplexList.size(); ++j )
        delete ::MC_MR_DComplexList[j];
    ::MC_MR_SComplexList.resize(0);
    ::MC_MR_DComplexList.resize(0);
    // Complex [MC,* ]
    for( unsigned j=0; j< ::MC_STAR_SComplexList.size(); ++j )
        delete ::MC_STAR_SComplexList[j];
    for( unsigned j=0; j< ::MC_STAR_DComplexList.size(); ++j )
        delete ::MC_STAR_DComplexList[j];
    ::MC_STAR_SComplexList.resize(0);
    ::MC_STAR_DComplexList.resize(0);
    // Complex [MD,* ]
    for( unsigned j=0; j< ::MD_STAR_SComplexList.size(); ++j )
        delete ::MD_STAR_SComplexList[j];
    for( unsigned j=0; j< ::MD_STAR_DComplexList.size(); ++j )
        delete ::MD_STAR_DComplexList[j];
    ::MD_STAR_SComplexList.resize(0);
    ::MD_STAR_DComplexList.resize(0);
    // Complex [MR,MC]
    for( unsigned j=0; j< ::MR_MC_SComplexList.size(); ++j )
        delete ::MR_MC_SComplexList[j];
    for( unsigned j=0; j< ::MR_MC_DComplexList.size(); ++j )
        delete ::MR_MC_DComplexList[j];
    ::MR_MC_SComplexList.resize(0);
    ::MR_MC_DComplexList.resize(0);
    // Complex [MR,* ]
    for( unsigned j=0; j< ::MR_STAR_SComplexList.size(); ++j )
        delete ::MR_STAR_SComplexList[j];
    for( unsigned j=0; j< ::MR_STAR_DComplexList.size(); ++j )
        delete ::MR_STAR_DComplexList[j];
    ::MR_STAR_SComplexList.resize(0);
    ::MR_STAR_DComplexList.resize(0);
    // Complex [* ,MC]
    for( unsigned j=0; j< ::STAR_MC_SComplexList.size(); ++j )
        delete ::STAR_MC_SComplexList[j];
    for( unsigned j=0; j< ::STAR_MC_DComplexList.size(); ++j )
        delete ::STAR_MC_DComplexList[j];
    ::STAR_MC_SComplexList.resize(0);
    ::STAR_MC_DComplexList.resize(0);
    // Complex [* ,MD]
    for( unsigned j=0; j< ::STAR_MD_SComplexList.size(); ++j )
        delete ::STAR_MD_SComplexList[j];
    for( unsigned j=0; j< ::STAR_MD_DComplexList.size(); ++j )
        delete ::STAR_MD_DComplexList[j];
    ::STAR_MD_SComplexList.resize(0);
    ::STAR_MD_DComplexList.resize(0);
    // Complex [* ,MR]
    for( unsigned j=0; j< ::STAR_MR_SComplexList.size(); ++j )
        delete ::STAR_MR_SComplexList[j];
    for( unsigned j=0; j< ::STAR_MR_DComplexList.size(); ++j )
        delete ::STAR_MR_DComplexList[j];
    ::STAR_MR_SComplexList.resize(0);
    ::STAR_MR_DComplexList.resize(0);
    // Complex [* ,* ]
    for( unsigned j=0; j< ::STAR_STAR_SComplexList.size(); ++j )
        delete ::STAR_STAR_SComplexList[j];
    for( unsigned j=0; j< ::STAR_STAR_DComplexList.size(); ++j )
        delete ::STAR_STAR_DComplexList[j];
    ::STAR_STAR_SComplexList.resize(0);
    ::STAR_STAR_DComplexList.resize(0);
    // Complex [* ,VC]
    for( unsigned j=0; j< ::STAR_VC_SComplexList.size(); ++j )
        delete ::STAR_VC_SComplexList[j];
    for( unsigned j=0; j< ::STAR_VC_DComplexList.size(); ++j )
        delete ::STAR_VC_DComplexList[j];
    ::STAR_VC_SComplexList.resize(0);
    ::STAR_VC_DComplexList.resize(0);
    // Complex [* ,VR]
    for( unsigned j=0; j< ::STAR_VR_SComplexList.size(); ++j )
        delete ::STAR_VR_SComplexList[j];
    for( unsigned j=0; j< ::STAR_VR_DComplexList.size(); ++j )
        delete ::STAR_VR_DComplexList[j];
    ::STAR_VR_SComplexList.resize(0);
    ::STAR_VR_DComplexList.resize(0);
    // Complex [VC,* ]
    for( unsigned j=0; j< ::VC_STAR_SComplexList.size(); ++j )
        delete ::VC_STAR_SComplexList[j];
    for( unsigned j=0; j< ::VC_STAR_DComplexList.size(); ++j )
        delete ::VC_STAR_DComplexList[j];
    ::VC_STAR_SComplexList.resize(0);
    ::VC_STAR_DComplexList.resize(0);
    // Complex [VR,* ]
    for( unsigned j=0; j< ::VR_STAR_SComplexList.size(); ++j )
        delete ::VR_STAR_SComplexList[j];
    for( unsigned j=0; j< ::VR_STAR_DComplexList.size(); ++j )
        delete ::VR_STAR_DComplexList[j];
    ::VR_STAR_SComplexList.resize(0);
    ::VR_STAR_DComplexList.resize(0);
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

Grid CreateDefaultGrid( ElementalComm comm )
{
    Grid g;
    CATCH(
        g = ::gridList.size();
        ::gridList.push_back( new elemental::Grid( comm ) );
    );
    return g;
}

Grid CreateGrid( ElementalComm comm, int r, int c )
{
    Grid g;
    CATCH(
        g = ::gridList.size();
        ::gridList.push_back( new elemental::Grid( comm, r, c ) );
    );
    return g;
}

void DestroyGrid( Grid g )
{
    delete ::gridList[g];
    ::gridList[g] = 0;
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

ElementalComm GridVCComm( Grid g )
{ ElementalComm VCComm; CATCH(VCComm=gridList[g]->VCComm()); return VCComm; }

ElementalComm GridVRComm( Grid g )
{ ElementalComm VRComm; CATCH(VRComm=gridList[g]->VRComm()); return VRComm; }

ElementalComm GridMCComm( Grid g )
{ ElementalComm MCComm; CATCH(MCComm=gridList[g]->MCComm()); return MCComm; }

ElementalComm GridMRComm( Grid g )
{ ElementalComm MRComm; CATCH(MRComm=gridList[g]->MRComm()); return MRComm; }

//----------------------------------------------------------------------------//
// [MC,MR] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
MC_MR_Single
CreateEmpty_MC_MR_Single( Grid g )
{
    MC_MR_Single A;
    CATCH(
        A = ::MC_MR_SingleList.size();
        ::MC_MR_SingleList.push_back
        ( new DistMatrix<float,MC,MR>( *::gridList[g] ) );
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
        ( new DistMatrix<double,MC,MR>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
MC_MR_SComplex
CreateEmpty_MC_MR_SComplex( Grid g )
{
    MC_MR_SComplex A;
    CATCH(
        A = ::MC_MR_SComplexList.size();
        ::MC_MR_SComplexList.push_back
        ( new DistMatrix<scomplex,MC,MR>( *::gridList[g] ) );
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
        ( new DistMatrix<dcomplex,MC,MR>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
MC_MR_Single
Register_MC_MR_Single
( int height, int width, int colAlignment, int rowAlignment,
  float* buffer, int ldim, Grid g )
{
    MC_MR_Single A;
    CATCH(
        A = ::MC_MR_SingleList.size();
        ::MC_MR_SingleList.push_back
        ( new DistMatrix<float,MC,MR>
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
        ( new DistMatrix<double,MC,MR>
          ( height, width, colAlignment, rowAlignment, buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
MC_MR_SComplex
Register_MC_MR_SComplex
( int height, int width, int colAlignment, int rowAlignment,
  SComplex* buffer, int ldim, Grid g )
{
    MC_MR_SComplex A;
    CATCH(
        A = ::MC_MR_SComplexList.size();
        ::MC_MR_SComplexList.push_back
        ( new DistMatrix<scomplex,MC,MR>
          ( height, width, colAlignment, rowAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
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
        ( new DistMatrix<dcomplex,MC,MR>
          ( height, width, colAlignment, rowAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the row/col alignments
int
ColAlignment_MC_MR_Single( MC_MR_Single A )
{ return ::MC_MR_SingleList[A]->ColAlignment(); }
int
RowAlignment_MC_MR_Single( MC_MR_Single A )
{ return ::MC_MR_SingleList[A]->RowAlignment(); }
int
ColAlignment_MC_MR_Double( MC_MR_Double A )
{ return ::MC_MR_DoubleList[A]->ColAlignment(); }
int
RowAlignment_MC_MR_Double( MC_MR_Double A )
{ return ::MC_MR_DoubleList[A]->RowAlignment(); }
#ifndef WITHOUT_COMPLEX
int
ColAlignment_MC_MR_SComplex( MC_MR_SComplex A )
{ return ::MC_MR_SComplexList[A]->ColAlignment(); }
int
RowAlignment_MC_MR_SComplex( MC_MR_SComplex A )
{ return ::MC_MR_SComplexList[A]->RowAlignment(); }
int
ColAlignment_MC_MR_DComplex( MC_MR_DComplex A )
{ return ::MC_MR_DComplexList[A]->ColAlignment(); }
int
RowAlignment_MC_MR_DComplex( MC_MR_DComplex A )
{ return ::MC_MR_DComplexList[A]->RowAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
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
#ifndef WITHOUT_COMPLEX
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
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [MC,* ] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
MC_STAR_Single
CreateEmpty_MC_STAR_Single( Grid g )
{
    MC_STAR_Single A;
    CATCH(
        A = ::MC_STAR_SingleList.size();
        ::MC_STAR_SingleList.push_back
        ( new DistMatrix<float,MC,STAR>( *::gridList[g] ) );
    );
    return A;
}
MC_STAR_Double
CreateEmpty_MC_STAR_Double( Grid g )
{
    MC_STAR_Double A;
    CATCH(
        A = ::MC_STAR_DoubleList.size();
        ::MC_STAR_DoubleList.push_back
        ( new DistMatrix<double,MC,STAR>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
MC_STAR_SComplex
CreateEmpty_MC_STAR_SComplex( Grid g )
{
    MC_STAR_SComplex A;
    CATCH(
        A = ::MC_STAR_SComplexList.size();
        ::MC_STAR_SComplexList.push_back
        ( new DistMatrix<scomplex,MC,STAR>( *::gridList[g] ) );
    );
    return A;
}
MC_STAR_DComplex
CreateEmpty_MC_STAR_DComplex( Grid g )
{
    MC_STAR_DComplex A;
    CATCH(
        A = ::MC_STAR_DComplexList.size();
        ::MC_STAR_DComplexList.push_back
        ( new DistMatrix<dcomplex,MC,STAR>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
MC_STAR_Single
Register_MC_STAR_Single
( int height, int width, int colAlignment, float* buffer, int ldim, Grid g )
{
    MC_STAR_Single A;
    CATCH(
        A = ::MC_STAR_SingleList.size();
        ::MC_STAR_SingleList.push_back
        ( new DistMatrix<float,MC,STAR>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
MC_STAR_Double
Register_MC_STAR_Double
( int height, int width, int colAlignment, double* buffer, int ldim, Grid g )
{
    MC_STAR_Double A;
    CATCH(
        A = ::MC_STAR_DoubleList.size();
        ::MC_STAR_DoubleList.push_back
        ( new DistMatrix<double,MC,STAR>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
MC_STAR_SComplex
Register_MC_STAR_SComplex
( int height, int width, int colAlignment, SComplex* buffer, int ldim, Grid g )
{
    MC_STAR_SComplex A;
    CATCH(
        A = ::MC_STAR_SComplexList.size();
        ::MC_STAR_SComplexList.push_back
        ( new DistMatrix<scomplex,MC,STAR>
          ( height, width, colAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
MC_STAR_DComplex
Register_MC_STAR_DComplex
( int height, int width, int colAlignment, DComplex* buffer, int ldim, Grid g )
{
    MC_STAR_DComplex A;
    CATCH(
        A = ::MC_STAR_DComplexList.size();
        ::MC_STAR_DComplexList.push_back
        ( new DistMatrix<dcomplex,MC,STAR>
          ( height, width, colAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the col alignments
int
ColAlignment_MC_STAR_Single( MC_STAR_Single A )
{ return ::MC_STAR_SingleList[A]->ColAlignment(); }
int
ColAlignment_MC_STAR_Double( MC_STAR_Double A )
{ return ::MC_STAR_DoubleList[A]->ColAlignment(); }
#ifndef WITHOUT_COMPLEX
int
ColAlignment_MC_STAR_SComplex( MC_STAR_SComplex A )
{ return ::MC_STAR_SComplexList[A]->ColAlignment(); }
int
ColAlignment_MC_STAR_DComplex( MC_STAR_DComplex A )
{ return ::MC_STAR_DComplexList[A]->ColAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_MC_STAR_Single( char* msg, MC_STAR_Single A )
{
    CATCH(
        std::string s( msg );
        ::MC_STAR_SingleList[A]->Print( s );
    );
}
void 
Print_MC_STAR_Double( char* msg, MC_STAR_Double A )
{
    CATCH(
        std::string s( msg );
        ::MC_STAR_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_MC_STAR_SComplex( char* msg, MC_STAR_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::MC_STAR_SComplexList[A]->Print( s );
    );
}
void 
Print_MC_STAR_DComplex( char* msg, MC_STAR_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::MC_STAR_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [MD,* ] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
MD_STAR_Single
CreateEmpty_MD_STAR_Single( Grid g )
{
    MD_STAR_Single A;
    CATCH(
        A = ::MD_STAR_SingleList.size();
        ::MD_STAR_SingleList.push_back
        ( new DistMatrix<float,MD,STAR>( *::gridList[g] ) );
    );
    return A;
}
MD_STAR_Double
CreateEmpty_MD_STAR_Double( Grid g )
{
    MD_STAR_Double A;
    CATCH(
        A = ::MD_STAR_DoubleList.size();
        ::MD_STAR_DoubleList.push_back
        ( new DistMatrix<double,MD,STAR>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
MD_STAR_SComplex
CreateEmpty_MD_STAR_SComplex( Grid g )
{
    MD_STAR_SComplex A;
    CATCH(
        A = ::MD_STAR_SComplexList.size();
        ::MD_STAR_SComplexList.push_back
        ( new DistMatrix<scomplex,MD,STAR>( *::gridList[g] ) );
    );
    return A;
}
MD_STAR_DComplex
CreateEmpty_MD_STAR_DComplex( Grid g )
{
    MD_STAR_DComplex A;
    CATCH(
        A = ::MD_STAR_DComplexList.size();
        ::MD_STAR_DComplexList.push_back
        ( new DistMatrix<dcomplex,MD,STAR>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
MD_STAR_Single
Register_MD_STAR_Single
( int height, int width, int colAlignment, float* buffer, int ldim, Grid g )
{
    MD_STAR_Single A;
    CATCH(
        A = ::MD_STAR_SingleList.size();
        ::MD_STAR_SingleList.push_back
        ( new DistMatrix<float,MD,STAR>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
MD_STAR_Double
Register_MD_STAR_Double
( int height, int width, int colAlignment, double* buffer, int ldim, Grid g )
{
    MD_STAR_Double A;
    CATCH(
        A = ::MD_STAR_DoubleList.size();
        ::MD_STAR_DoubleList.push_back
        ( new DistMatrix<double,MD,STAR>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
MD_STAR_SComplex
Register_MD_STAR_SComplex
( int height, int width, int colAlignment, SComplex* buffer, int ldim, Grid g )
{
    MD_STAR_SComplex A;
    CATCH(
        A = ::MD_STAR_SComplexList.size();
        ::MD_STAR_SComplexList.push_back
        ( new DistMatrix<scomplex,MD,STAR>
          ( height, width, colAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
MD_STAR_DComplex
Register_MD_STAR_DComplex
( int height, int width, int colAlignment, DComplex* buffer, int ldim, Grid g )
{
    MD_STAR_DComplex A;
    CATCH(
        A = ::MD_STAR_DComplexList.size();
        ::MD_STAR_DComplexList.push_back
        ( new DistMatrix<dcomplex,MD,STAR>
          ( height, width, colAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the col alignments
int
ColAlignment_MD_STAR_Single( MD_STAR_Single A )
{ return ::MD_STAR_SingleList[A]->ColAlignment(); }
int
ColAlignment_MD_STAR_Double( MD_STAR_Double A )
{ return ::MD_STAR_DoubleList[A]->ColAlignment(); }
#ifndef WITHOUT_COMPLEX
int
ColAlignment_MD_STAR_SComplex( MD_STAR_SComplex A )
{ return ::MD_STAR_SComplexList[A]->ColAlignment(); }
int
ColAlignment_MD_STAR_DComplex( MD_STAR_DComplex A )
{ return ::MD_STAR_DComplexList[A]->ColAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_MD_STAR_Single( char* msg, MD_STAR_Single A )
{
    CATCH(
        std::string s( msg );
        ::MD_STAR_SingleList[A]->Print( s );
    );
}
void 
Print_MD_STAR_Double( char* msg, MD_STAR_Double A )
{
    CATCH(
        std::string s( msg );
        ::MD_STAR_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_MD_STAR_SComplex( char* msg, MD_STAR_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::MD_STAR_SComplexList[A]->Print( s );
    );
}
void 
Print_MD_STAR_DComplex( char* msg, MD_STAR_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::MD_STAR_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [MR,MC] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
MR_MC_Single
CreateEmpty_MR_MC_Single( Grid g )
{
    MR_MC_Single A;
    CATCH(
        A = ::MR_MC_SingleList.size();
        ::MR_MC_SingleList.push_back
        ( new DistMatrix<float,MR,MC>( *::gridList[g] ) );
    );
    return A;
}
MR_MC_Double
CreateEmpty_MR_MC_Double( Grid g )
{
    MR_MC_Double A;
    CATCH(
        A = ::MR_MC_DoubleList.size();
        ::MR_MC_DoubleList.push_back
        ( new DistMatrix<double,MR,MC>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
MR_MC_SComplex
CreateEmpty_MR_MC_SComplex( Grid g )
{
    MR_MC_SComplex A;
    CATCH(
        A = ::MR_MC_SComplexList.size();
        ::MR_MC_SComplexList.push_back
        ( new DistMatrix<scomplex,MR,MC>( *::gridList[g] ) );
    );
    return A;
}
MR_MC_DComplex
CreateEmpty_MR_MC_DComplex( Grid g )
{
    MR_MC_DComplex A;
    CATCH(
        A = ::MR_MC_DComplexList.size();
        ::MR_MC_DComplexList.push_back
        ( new DistMatrix<dcomplex,MR,MC>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
MR_MC_Single
Register_MR_MC_Single
( int height, int width, int colAlignment, int rowAlignment,
  float* buffer, int ldim, Grid g )
{
    MR_MC_Single A;
    CATCH(
        A = ::MR_MC_SingleList.size();
        ::MR_MC_SingleList.push_back
        ( new DistMatrix<float,MR,MC>
          ( height, width, colAlignment, rowAlignment, buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
MR_MC_Double
Register_MR_MC_Double
( int height, int width, int colAlignment, int rowAlignment,
  double* buffer, int ldim, Grid g )
{
    MR_MC_Double A;
    CATCH(
        A = ::MR_MC_DoubleList.size();
        ::MR_MC_DoubleList.push_back
        ( new DistMatrix<double,MR,MC>
          ( height, width, colAlignment, rowAlignment, buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
MR_MC_SComplex
Register_MR_MC_SComplex
( int height, int width, int colAlignment, int rowAlignment,
  SComplex* buffer, int ldim, Grid g )
{
    MR_MC_SComplex A;
    CATCH(
        A = ::MR_MC_SComplexList.size();
        ::MR_MC_SComplexList.push_back
        ( new DistMatrix<scomplex,MR,MC>
          ( height, width, colAlignment, rowAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
MR_MC_DComplex
Register_MR_MC_DComplex
( int height, int width, int colAlignment, int rowAlignment,
  DComplex* buffer, int ldim, Grid g )
{
    MR_MC_DComplex A;
    CATCH(
        A = ::MR_MC_DComplexList.size();
        ::MR_MC_DComplexList.push_back
        ( new DistMatrix<dcomplex,MR,MC>
          ( height, width, colAlignment, rowAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the row/col alignments
int
ColAlignment_MR_MC_Single( MR_MC_Single A )
{ return ::MR_MC_SingleList[A]->ColAlignment(); }
int
RowAlignment_MR_MC_Single( MR_MC_Single A )
{ return ::MR_MC_SingleList[A]->RowAlignment(); }
int
ColAlignment_MR_MC_Double( MR_MC_Double A )
{ return ::MR_MC_DoubleList[A]->ColAlignment(); }
int
RowAlignment_MR_MC_Double( MR_MC_Double A )
{ return ::MR_MC_DoubleList[A]->RowAlignment(); }
#ifndef WITHOUT_COMPLEX
int
ColAlignment_MR_MC_SComplex( MR_MC_SComplex A )
{ return ::MR_MC_SComplexList[A]->ColAlignment(); }
int
RowAlignment_MR_MC_SComplex( MR_MC_SComplex A )
{ return ::MR_MC_SComplexList[A]->RowAlignment(); }
int
ColAlignment_MR_MC_DComplex( MR_MC_DComplex A )
{ return ::MR_MC_DComplexList[A]->ColAlignment(); }
int
RowAlignment_MR_MC_DComplex( MR_MC_DComplex A )
{ return ::MR_MC_DComplexList[A]->RowAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_MR_MC_Single( char* msg, MR_MC_Single A )
{
    CATCH(
        std::string s( msg );
        ::MR_MC_SingleList[A]->Print( s );
    );
}
void 
Print_MR_MC_Double( char* msg, MR_MC_Double A )
{
    CATCH(
        std::string s( msg );
        ::MR_MC_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_MR_MC_SComplex( char* msg, MR_MC_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::MR_MC_SComplexList[A]->Print( s );
    );
}
void 
Print_MR_MC_DComplex( char* msg, MR_MC_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::MR_MC_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [MR,* ] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
MR_STAR_Single
CreateEmpty_MR_STAR_Single( Grid g )
{
    MR_STAR_Single A;
    CATCH(
        A = ::MR_STAR_SingleList.size();
        ::MR_STAR_SingleList.push_back
        ( new DistMatrix<float,MR,STAR>( *::gridList[g] ) );
    );
    return A;
}
MR_STAR_Double
CreateEmpty_MR_STAR_Double( Grid g )
{
    MR_STAR_Double A;
    CATCH(
        A = ::MR_STAR_DoubleList.size();
        ::MR_STAR_DoubleList.push_back
        ( new DistMatrix<double,MR,STAR>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
MR_STAR_SComplex
CreateEmpty_MR_STAR_SComplex( Grid g )
{
    MR_STAR_SComplex A;
    CATCH(
        A = ::MR_STAR_SComplexList.size();
        ::MR_STAR_SComplexList.push_back
        ( new DistMatrix<scomplex,MR,STAR>( *::gridList[g] ) );
    );
    return A;
}
MR_STAR_DComplex
CreateEmpty_MR_STAR_DComplex( Grid g )
{
    MR_STAR_DComplex A;
    CATCH(
        A = ::MR_STAR_DComplexList.size();
        ::MR_STAR_DComplexList.push_back
        ( new DistMatrix<dcomplex,MR,STAR>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
MR_STAR_Single
Register_MR_STAR_Single
( int height, int width, int colAlignment, float* buffer, int ldim, Grid g )
{
    MR_STAR_Single A;
    CATCH(
        A = ::MR_STAR_SingleList.size();
        ::MR_STAR_SingleList.push_back
        ( new DistMatrix<float,MR,STAR>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
MR_STAR_Double
Register_MR_STAR_Double
( int height, int width, int colAlignment, double* buffer, int ldim, Grid g )
{
    MR_STAR_Double A;
    CATCH(
        A = ::MR_STAR_DoubleList.size();
        ::MR_STAR_DoubleList.push_back
        ( new DistMatrix<double,MR,STAR>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
MR_STAR_SComplex
Register_MR_STAR_SComplex
( int height, int width, int colAlignment, SComplex* buffer, int ldim, Grid g )
{
    MR_STAR_SComplex A;
    CATCH(
        A = ::MR_STAR_SComplexList.size();
        ::MR_STAR_SComplexList.push_back
        ( new DistMatrix<scomplex,MR,STAR>
          ( height, width, colAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
MR_STAR_DComplex
Register_MR_STAR_DComplex
( int height, int width, int colAlignment, DComplex* buffer, int ldim, Grid g )
{
    MR_STAR_DComplex A;
    CATCH(
        A = ::MR_STAR_DComplexList.size();
        ::MR_STAR_DComplexList.push_back
        ( new DistMatrix<dcomplex,MR,STAR>
          ( height, width, colAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the col alignments
int
ColAlignment_MR_STAR_Single( MR_STAR_Single A )
{ return ::MR_STAR_SingleList[A]->ColAlignment(); }
int
ColAlignment_MR_STAR_Double( MR_STAR_Double A )
{ return ::MR_STAR_DoubleList[A]->ColAlignment(); }
#ifndef WITHOUT_COMPLEX
int
ColAlignment_MR_STAR_SComplex( MR_STAR_SComplex A )
{ return ::MR_STAR_SComplexList[A]->ColAlignment(); }
int
ColAlignment_MR_STAR_DComplex( MR_STAR_DComplex A )
{ return ::MR_STAR_DComplexList[A]->ColAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_MR_STAR_Single( char* msg, MR_STAR_Single A )
{
    CATCH(
        std::string s( msg );
        ::MR_STAR_SingleList[A]->Print( s );
    );
}
void 
Print_MR_STAR_Double( char* msg, MR_STAR_Double A )
{
    CATCH(
        std::string s( msg );
        ::MR_STAR_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_MR_STAR_SComplex( char* msg, MR_STAR_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::MR_STAR_SComplexList[A]->Print( s );
    );
}
void 
Print_MR_STAR_DComplex( char* msg, MR_STAR_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::MR_STAR_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [* ,MC] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
STAR_MC_Single
CreateEmpty_STAR_MC_Single( Grid g )
{
    STAR_MC_Single A;
    CATCH(
        A = ::STAR_MC_SingleList.size();
        ::STAR_MC_SingleList.push_back
        ( new DistMatrix<float,STAR,MC>( *::gridList[g] ) );
    );
    return A;
}
STAR_MC_Double
CreateEmpty_STAR_MC_Double( Grid g )
{
    STAR_MC_Double A;
    CATCH(
        A = ::STAR_MC_DoubleList.size();
        ::STAR_MC_DoubleList.push_back
        ( new DistMatrix<double,STAR,MC>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
STAR_MC_SComplex
CreateEmpty_STAR_MC_SComplex( Grid g )
{
    STAR_MC_SComplex A;
    CATCH(
        A = ::STAR_MC_SComplexList.size();
        ::STAR_MC_SComplexList.push_back
        ( new DistMatrix<scomplex,STAR,MC>( *::gridList[g] ) );
    );
    return A;
}
STAR_MC_DComplex
CreateEmpty_STAR_MC_DComplex( Grid g )
{
    STAR_MC_DComplex A;
    CATCH(
        A = ::STAR_MC_DComplexList.size();
        ::STAR_MC_DComplexList.push_back
        ( new DistMatrix<dcomplex,STAR,MC>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
STAR_MC_Single
Register_STAR_MC_Single
( int height, int width, int rowAlignment, float* buffer, int ldim, Grid g )
{
    STAR_MC_Single A;
    CATCH(
        A = ::STAR_MC_SingleList.size();
        ::STAR_MC_SingleList.push_back
        ( new DistMatrix<float,STAR,MC>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
STAR_MC_Double
Register_STAR_MC_Double
( int height, int width, int rowAlignment, double* buffer, int ldim, Grid g )
{
    STAR_MC_Double A;
    CATCH(
        A = ::STAR_MC_DoubleList.size();
        ::STAR_MC_DoubleList.push_back
        ( new DistMatrix<double,STAR,MC>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
STAR_MC_SComplex
Register_STAR_MC_SComplex
( int height, int width, int rowAlignment, SComplex* buffer, int ldim, Grid g )
{
    STAR_MC_SComplex A;
    CATCH(
        A = ::STAR_MC_SComplexList.size();
        ::STAR_MC_SComplexList.push_back
        ( new DistMatrix<scomplex,STAR,MC>
          ( height, width, rowAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
STAR_MC_DComplex
Register_STAR_MC_DComplex
( int height, int width, int rowAlignment, DComplex* buffer, int ldim, Grid g )
{
    STAR_MC_DComplex A;
    CATCH(
        A = ::STAR_MC_DComplexList.size();
        ::STAR_MC_DComplexList.push_back
        ( new DistMatrix<dcomplex,STAR,MC>
          ( height, width, rowAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the row alignments
int
RowAlignment_STAR_MC_Single( STAR_MC_Single A )
{ return ::STAR_MC_SingleList[A]->RowAlignment(); }
int
RowAlignment_STAR_MC_Double( STAR_MC_Double A )
{ return ::STAR_MC_DoubleList[A]->RowAlignment(); }
#ifndef WITHOUT_COMPLEX
int
RowAlignment_STAR_MC_SComplex( STAR_MC_SComplex A )
{ return ::STAR_MC_SComplexList[A]->RowAlignment(); }
int
RowAlignment_STAR_MC_DComplex( STAR_MC_DComplex A )
{ return ::STAR_MC_DComplexList[A]->RowAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_STAR_MC_Single( char* msg, STAR_MC_Single A )
{
    CATCH(
        std::string s( msg );
        ::STAR_MC_SingleList[A]->Print( s );
    );
}
void 
Print_STAR_MC_Double( char* msg, STAR_MC_Double A )
{
    CATCH(
        std::string s( msg );
        ::STAR_MC_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_STAR_MC_SComplex( char* msg, STAR_MC_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::STAR_MC_SComplexList[A]->Print( s );
    );
}
void 
Print_STAR_MC_DComplex( char* msg, STAR_MC_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::STAR_MC_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [* ,MD] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
STAR_MD_Single
CreateEmpty_STAR_MD_Single( Grid g )
{
    STAR_MD_Single A;
    CATCH(
        A = ::STAR_MD_SingleList.size();
        ::STAR_MD_SingleList.push_back
        ( new DistMatrix<float,STAR,MD>( *::gridList[g] ) );
    );
    return A;
}
STAR_MD_Double
CreateEmpty_STAR_MD_Double( Grid g )
{
    STAR_MD_Double A;
    CATCH(
        A = ::STAR_MD_DoubleList.size();
        ::STAR_MD_DoubleList.push_back
        ( new DistMatrix<double,STAR,MD>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
STAR_MD_SComplex
CreateEmpty_STAR_MD_SComplex( Grid g )
{
    STAR_MD_SComplex A;
    CATCH(
        A = ::STAR_MD_SComplexList.size();
        ::STAR_MD_SComplexList.push_back
        ( new DistMatrix<scomplex,STAR,MD>( *::gridList[g] ) );
    );
    return A;
}
STAR_MD_DComplex
CreateEmpty_STAR_MD_DComplex( Grid g )
{
    STAR_MD_DComplex A;
    CATCH(
        A = ::STAR_MD_DComplexList.size();
        ::STAR_MD_DComplexList.push_back
        ( new DistMatrix<dcomplex,STAR,MD>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
STAR_MD_Single
Register_STAR_MD_Single
( int height, int width, int rowAlignment, float* buffer, int ldim, Grid g )
{
    STAR_MD_Single A;
    CATCH(
        A = ::STAR_MD_SingleList.size();
        ::STAR_MD_SingleList.push_back
        ( new DistMatrix<float,STAR,MD>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
STAR_MD_Double
Register_STAR_MD_Double
( int height, int width, int rowAlignment, double* buffer, int ldim, Grid g )
{
    STAR_MD_Double A;
    CATCH(
        A = ::STAR_MD_DoubleList.size();
        ::STAR_MD_DoubleList.push_back
        ( new DistMatrix<double,STAR,MD>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
STAR_MD_SComplex
Register_STAR_MD_SComplex
( int height, int width, int rowAlignment, SComplex* buffer, int ldim, Grid g )
{
    STAR_MD_SComplex A;
    CATCH(
        A = ::STAR_MD_SComplexList.size();
        ::STAR_MD_SComplexList.push_back
        ( new DistMatrix<scomplex,STAR,MD>
          ( height, width, rowAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
STAR_MD_DComplex
Register_STAR_MD_DComplex
( int height, int width, int rowAlignment, DComplex* buffer, int ldim, Grid g )
{
    STAR_MD_DComplex A;
    CATCH(
        A = ::STAR_MD_DComplexList.size();
        ::STAR_MD_DComplexList.push_back
        ( new DistMatrix<dcomplex,STAR,MD>
          ( height, width, rowAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the row alignments
int
RowAlignment_MD_STAR_Single( MD_STAR_Single A )
{ return ::MD_STAR_SingleList[A]->RowAlignment(); }
int
RowAlignment_MD_STAR_Double( MD_STAR_Double A )
{ return ::MD_STAR_DoubleList[A]->RowAlignment(); }
#ifndef WITHOUT_COMPLEX
int
RowAlignment_MD_STAR_SComplex( MD_STAR_SComplex A )
{ return ::MD_STAR_SComplexList[A]->RowAlignment(); }
int
RowAlignment_MD_STAR_DComplex( MD_STAR_DComplex A )
{ return ::MD_STAR_DComplexList[A]->RowAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_STAR_MD_Single( char* msg, STAR_MD_Single A )
{
    CATCH(
        std::string s( msg );
        ::STAR_MD_SingleList[A]->Print( s );
    );
}
void 
Print_STAR_MD_Double( char* msg, STAR_MD_Double A )
{
    CATCH(
        std::string s( msg );
        ::STAR_MD_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_STAR_MD_SComplex( char* msg, STAR_MD_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::STAR_MD_SComplexList[A]->Print( s );
    );
}
void 
Print_STAR_MD_DComplex( char* msg, STAR_MD_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::STAR_MD_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [* ,MR] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
STAR_MR_Single
CreateEmpty_STAR_MR_Single( Grid g )
{
    STAR_MR_Single A;
    CATCH(
        A = ::STAR_MR_SingleList.size();
        ::STAR_MR_SingleList.push_back
        ( new DistMatrix<float,STAR,MR>( *::gridList[g] ) );
    );
    return A;
}
STAR_MR_Double
CreateEmpty_STAR_MR_Double( Grid g )
{
    STAR_MR_Double A;
    CATCH(
        A = ::STAR_MR_DoubleList.size();
        ::STAR_MR_DoubleList.push_back
        ( new DistMatrix<double,STAR,MR>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
STAR_MR_SComplex
CreateEmpty_STAR_MR_SComplex( Grid g )
{
    STAR_MR_SComplex A;
    CATCH(
        A = ::STAR_MR_SComplexList.size();
        ::STAR_MR_SComplexList.push_back
        ( new DistMatrix<scomplex,STAR,MR>( *::gridList[g] ) );
    );
    return A;
}
STAR_MR_DComplex
CreateEmpty_STAR_MR_DComplex( Grid g )
{
    STAR_MR_DComplex A;
    CATCH(
        A = ::STAR_MR_DComplexList.size();
        ::STAR_MR_DComplexList.push_back
        ( new DistMatrix<dcomplex,STAR,MR>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
STAR_MR_Single
Register_STAR_MR_Single
( int height, int width, int rowAlignment, float* buffer, int ldim, Grid g )
{
    STAR_MR_Single A;
    CATCH(
        A = ::STAR_MR_SingleList.size();
        ::STAR_MR_SingleList.push_back
        ( new DistMatrix<float,STAR,MR>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
STAR_MR_Double
Register_STAR_MR_Double
( int height, int width, int rowAlignment, double* buffer, int ldim, Grid g )
{
    STAR_MR_Double A;
    CATCH(
        A = ::STAR_MR_DoubleList.size();
        ::STAR_MR_DoubleList.push_back
        ( new DistMatrix<double,STAR,MR>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
STAR_MR_SComplex
Register_STAR_MR_SComplex
( int height, int width, int rowAlignment, SComplex* buffer, int ldim, Grid g )
{
    STAR_MR_SComplex A;
    CATCH(
        A = ::STAR_MR_SComplexList.size();
        ::STAR_MR_SComplexList.push_back
        ( new DistMatrix<scomplex,STAR,MR>
          ( height, width, rowAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
STAR_MR_DComplex
Register_STAR_MR_DComplex
( int height, int width, int rowAlignment, DComplex* buffer, int ldim, Grid g )
{
    STAR_MR_DComplex A;
    CATCH(
        A = ::STAR_MR_DComplexList.size();
        ::STAR_MR_DComplexList.push_back
        ( new DistMatrix<dcomplex,STAR,MR>
          ( height, width, rowAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the row alignments
int
RowAlignment_STAR_MR_Single( STAR_MR_Single A )
{ return ::STAR_MR_SingleList[A]->RowAlignment(); }
int
RowAlignment_STAR_MR_Double( STAR_MR_Double A )
{ return ::STAR_MR_DoubleList[A]->RowAlignment(); }
#ifndef WITHOUT_COMPLEX
int
RowAlignment_STAR_MR_SComplex( STAR_MR_SComplex A )
{ return ::STAR_MR_SComplexList[A]->RowAlignment(); }
int
RowAlignment_STAR_MR_DComplex( STAR_MR_DComplex A )
{ return ::STAR_MR_DComplexList[A]->RowAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_STAR_MR_Single( char* msg, STAR_MR_Single A )
{
    CATCH(
        std::string s( msg );
        ::STAR_MR_SingleList[A]->Print( s );
    );
}
void 
Print_STAR_MR_Double( char* msg, STAR_MR_Double A )
{
    CATCH(
        std::string s( msg );
        ::STAR_MR_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_STAR_MR_SComplex( char* msg, STAR_MR_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::STAR_MR_SComplexList[A]->Print( s );
    );
}
void 
Print_STAR_MR_DComplex( char* msg, STAR_MR_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::STAR_MR_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [* ,* ] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
STAR_STAR_Single
CreateEmpty_STAR_STAR_Single( Grid g )
{
    STAR_STAR_Single A;
    CATCH(
        A = ::STAR_STAR_SingleList.size();
        ::STAR_STAR_SingleList.push_back
        ( new DistMatrix<float,STAR,STAR>( *::gridList[g] ) );
    );
    return A;
}
STAR_STAR_Double
CreateEmpty_STAR_STAR_Double( Grid g )
{
    STAR_STAR_Double A;
    CATCH(
        A = ::STAR_STAR_DoubleList.size();
        ::STAR_STAR_DoubleList.push_back
        ( new DistMatrix<double,STAR,STAR>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
STAR_STAR_SComplex
CreateEmpty_STAR_STAR_SComplex( Grid g )
{
    STAR_STAR_SComplex A;
    CATCH(
        A = ::STAR_STAR_SComplexList.size();
        ::STAR_STAR_SComplexList.push_back
        ( new DistMatrix<scomplex,STAR,STAR>( *::gridList[g] ) );
    );
    return A;
}
STAR_STAR_DComplex
CreateEmpty_STAR_STAR_DComplex( Grid g )
{
    STAR_STAR_DComplex A;
    CATCH(
        A = ::STAR_STAR_DComplexList.size();
        ::STAR_STAR_DComplexList.push_back
        ( new DistMatrix<dcomplex,STAR,STAR>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
STAR_STAR_Single
Register_STAR_STAR_Single
( int height, int width, float* buffer, int ldim, Grid g )
{
    STAR_STAR_Single A;
    CATCH(
        A = ::STAR_STAR_SingleList.size();
        ::STAR_STAR_SingleList.push_back
        ( new DistMatrix<float,STAR,STAR>
          ( height, width, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
STAR_STAR_Double
Register_STAR_STAR_Double
( int height, int width, double* buffer, int ldim, Grid g )
{
    STAR_STAR_Double A;
    CATCH(
        A = ::STAR_STAR_DoubleList.size();
        ::STAR_STAR_DoubleList.push_back
        ( new DistMatrix<double,STAR,STAR>
          ( height, width, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
STAR_STAR_SComplex
Register_STAR_STAR_SComplex
( int height, int width, SComplex* buffer, int ldim, Grid g )
{
    STAR_STAR_SComplex A;
    CATCH(
        A = ::STAR_STAR_SComplexList.size();
        ::STAR_STAR_SComplexList.push_back
        ( new DistMatrix<scomplex,STAR,STAR>
          ( height, width, (scomplex*)buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
STAR_STAR_DComplex
Register_STAR_STAR_DComplex
( int height, int width, DComplex* buffer, int ldim, Grid g )
{
    STAR_STAR_DComplex A;
    CATCH(
        A = ::STAR_STAR_DComplexList.size();
        ::STAR_STAR_DComplexList.push_back
        ( new DistMatrix<dcomplex,STAR,STAR>
          ( height, width, (dcomplex*)buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Print
void 
Print_STAR_STAR_Single( char* msg, STAR_STAR_Single A )
{
    CATCH(
        std::string s( msg );
        ::STAR_STAR_SingleList[A]->Print( s );
    );
}
void 
Print_STAR_STAR_Double( char* msg, STAR_STAR_Double A )
{
    CATCH(
        std::string s( msg );
        ::STAR_STAR_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_STAR_STAR_SComplex( char* msg, STAR_STAR_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::STAR_STAR_SComplexList[A]->Print( s );
    );
}
void 
Print_STAR_STAR_DComplex( char* msg, STAR_STAR_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::STAR_STAR_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [* ,VC] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
STAR_VC_Single
CreateEmpty_STAR_VC_Single( Grid g )
{
    STAR_VC_Single A;
    CATCH(
        A = ::STAR_VC_SingleList.size();
        ::STAR_VC_SingleList.push_back
        ( new DistMatrix<float,STAR,VC>( *::gridList[g] ) );
    );
    return A;
}
STAR_VC_Double
CreateEmpty_STAR_VC_Double( Grid g )
{
    STAR_VC_Double A;
    CATCH(
        A = ::STAR_VC_DoubleList.size();
        ::STAR_VC_DoubleList.push_back
        ( new DistMatrix<double,STAR,VC>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
STAR_VC_SComplex
CreateEmpty_STAR_VC_SComplex( Grid g )
{
    STAR_VC_SComplex A;
    CATCH(
        A = ::STAR_VC_SComplexList.size();
        ::STAR_VC_SComplexList.push_back
        ( new DistMatrix<scomplex,STAR,VC>( *::gridList[g] ) );
    );
    return A;
}
STAR_VC_DComplex
CreateEmpty_STAR_VC_DComplex( Grid g )
{
    STAR_VC_DComplex A;
    CATCH(
        A = ::STAR_VC_DComplexList.size();
        ::STAR_VC_DComplexList.push_back
        ( new DistMatrix<dcomplex,STAR,VC>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
STAR_VC_Single
Register_STAR_VC_Single
( int height, int width, int rowAlignment, float* buffer, int ldim, Grid g )
{
    STAR_VC_Single A;
    CATCH(
        A = ::STAR_VC_SingleList.size();
        ::STAR_VC_SingleList.push_back
        ( new DistMatrix<float,STAR,VC>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
STAR_VC_Double
Register_STAR_VC_Double
( int height, int width, int rowAlignment, double* buffer, int ldim, Grid g )
{
    STAR_VC_Double A;
    CATCH(
        A = ::STAR_VC_DoubleList.size();
        ::STAR_VC_DoubleList.push_back
        ( new DistMatrix<double,STAR,VC>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
STAR_VC_SComplex
Register_STAR_VC_SComplex
( int height, int width, int rowAlignment, SComplex* buffer, int ldim, Grid g )
{
    STAR_VC_SComplex A;
    CATCH(
        A = ::STAR_VC_SComplexList.size();
        ::STAR_VC_SComplexList.push_back
        ( new DistMatrix<scomplex,STAR,VC>
          ( height, width, rowAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
STAR_VC_DComplex
Register_STAR_VC_DComplex
( int height, int width, int rowAlignment, DComplex* buffer, int ldim, Grid g )
{
    STAR_VC_DComplex A;
    CATCH(
        A = ::STAR_VC_DComplexList.size();
        ::STAR_VC_DComplexList.push_back
        ( new DistMatrix<dcomplex,STAR,VC>
          ( height, width, rowAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the row alignments
int
RowAlignment_STAR_VC_Single( STAR_VC_Single A )
{ return ::STAR_VC_SingleList[A]->RowAlignment(); }
int
RowAlignment_STAR_VC_Double( STAR_VC_Double A )
{ return ::STAR_VC_DoubleList[A]->RowAlignment(); }
#ifndef WITHOUT_COMPLEX
int
RowAlignment_STAR_VC_SComplex( STAR_VC_SComplex A )
{ return ::STAR_VC_SComplexList[A]->RowAlignment(); }
int
RowAlignment_STAR_VC_DComplex( STAR_VC_DComplex A )
{ return ::STAR_VC_DComplexList[A]->RowAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_STAR_VC_Single( char* msg, STAR_VC_Single A )
{
    CATCH(
        std::string s( msg );
        ::STAR_VC_SingleList[A]->Print( s );
    );
}
void 
Print_STAR_VC_Double( char* msg, STAR_VC_Double A )
{
    CATCH(
        std::string s( msg );
        ::STAR_VC_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_STAR_VC_SComplex( char* msg, STAR_VC_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::STAR_VC_SComplexList[A]->Print( s );
    );
}
void 
Print_STAR_VC_DComplex( char* msg, STAR_VC_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::STAR_VC_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [* ,VR] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
STAR_VR_Single
CreateEmpty_STAR_VR_Single( Grid g )
{
    STAR_VR_Single A;
    CATCH(
        A = ::STAR_VR_SingleList.size();
        ::STAR_VR_SingleList.push_back
        ( new DistMatrix<float,STAR,VR>( *::gridList[g] ) );
    );
    return A;
}
STAR_VR_Double
CreateEmpty_STAR_VR_Double( Grid g )
{
    STAR_VR_Double A;
    CATCH(
        A = ::STAR_VR_DoubleList.size();
        ::STAR_VR_DoubleList.push_back
        ( new DistMatrix<double,STAR,VR>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
STAR_VR_SComplex
CreateEmpty_STAR_VR_SComplex( Grid g )
{
    STAR_VR_SComplex A;
    CATCH(
        A = ::STAR_VR_SComplexList.size();
        ::STAR_VR_SComplexList.push_back
        ( new DistMatrix<scomplex,STAR,VR>( *::gridList[g] ) );
    );
    return A;
}
STAR_VR_DComplex
CreateEmpty_STAR_VR_DComplex( Grid g )
{
    STAR_VR_DComplex A;
    CATCH(
        A = ::STAR_VR_DComplexList.size();
        ::STAR_VR_DComplexList.push_back
        ( new DistMatrix<dcomplex,STAR,VR>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
STAR_VR_Single
Register_STAR_VR_Single
( int height, int width, int rowAlignment, float* buffer, int ldim, Grid g )
{
    STAR_VR_Single A;
    CATCH(
        A = ::STAR_VR_SingleList.size();
        ::STAR_VR_SingleList.push_back
        ( new DistMatrix<float,STAR,VR>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
STAR_VR_Double
Register_STAR_VR_Double
( int height, int width, int rowAlignment, double* buffer, int ldim, Grid g )
{
    STAR_VR_Double A;
    CATCH(
        A = ::STAR_VR_DoubleList.size();
        ::STAR_VR_DoubleList.push_back
        ( new DistMatrix<double,STAR,VR>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
STAR_VR_SComplex
Register_STAR_VR_SComplex
( int height, int width, int rowAlignment, SComplex* buffer, int ldim, Grid g )
{
    STAR_VR_SComplex A;
    CATCH(
        A = ::STAR_VR_SComplexList.size();
        ::STAR_VR_SComplexList.push_back
        ( new DistMatrix<scomplex,STAR,VR>
          ( height, width, rowAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
STAR_VR_DComplex
Register_STAR_VR_DComplex
( int height, int width, int rowAlignment, DComplex* buffer, int ldim, Grid g )
{
    STAR_VR_DComplex A;
    CATCH(
        A = ::STAR_VR_DComplexList.size();
        ::STAR_VR_DComplexList.push_back
        ( new DistMatrix<dcomplex,STAR,VR>
          ( height, width, rowAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the row alignments
int
RowAlignment_STAR_VR_Single( STAR_VR_Single A )
{ return ::STAR_VR_SingleList[A]->RowAlignment(); }
int
RowAlignment_STAR_VR_Double( STAR_VR_Double A )
{ return ::STAR_VR_DoubleList[A]->RowAlignment(); }
#ifndef WITHOUT_COMPLEX
int
RowAlignment_STAR_VR_SComplex( STAR_VR_SComplex A )
{ return ::STAR_VR_SComplexList[A]->RowAlignment(); }
int
RowAlignment_STAR_VR_DComplex( STAR_VR_DComplex A )
{ return ::STAR_VR_DComplexList[A]->RowAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_STAR_VR_Single( char* msg, STAR_VR_Single A )
{
    CATCH(
        std::string s( msg );
        ::STAR_VR_SingleList[A]->Print( s );
    );
}
void 
Print_STAR_VR_Double( char* msg, STAR_VR_Double A )
{
    CATCH(
        std::string s( msg );
        ::STAR_VR_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_STAR_VR_SComplex( char* msg, STAR_VR_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::STAR_VR_SComplexList[A]->Print( s );
    );
}
void 
Print_STAR_VR_DComplex( char* msg, STAR_VR_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::STAR_VR_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [VC,* ] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
VC_STAR_Single
CreateEmpty_VC_STAR_Single( Grid g )
{
    VC_STAR_Single A;
    CATCH(
        A = ::VC_STAR_SingleList.size();
        ::VC_STAR_SingleList.push_back
        ( new DistMatrix<float,VC,STAR>( *::gridList[g] ) );
    );
    return A;
}
VC_STAR_Double
CreateEmpty_VC_STAR_Double( Grid g )
{
    VC_STAR_Double A;
    CATCH(
        A = ::VC_STAR_DoubleList.size();
        ::VC_STAR_DoubleList.push_back
        ( new DistMatrix<double,VC,STAR>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
VC_STAR_SComplex
CreateEmpty_VC_STAR_SComplex( Grid g )
{
    VC_STAR_SComplex A;
    CATCH(
        A = ::VC_STAR_SComplexList.size();
        ::VC_STAR_SComplexList.push_back
        ( new DistMatrix<scomplex,VC,STAR>( *::gridList[g] ) );
    );
    return A;
}
VC_STAR_DComplex
CreateEmpty_VC_STAR_DComplex( Grid g )
{
    VC_STAR_DComplex A;
    CATCH(
        A = ::VC_STAR_DComplexList.size();
        ::VC_STAR_DComplexList.push_back
        ( new DistMatrix<dcomplex,VC,STAR>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
VC_STAR_Single
Register_VC_STAR_Single
( int height, int width, int colAlignment, float* buffer, int ldim, Grid g )
{
    VC_STAR_Single A;
    CATCH(
        A = ::VC_STAR_SingleList.size();
        ::VC_STAR_SingleList.push_back
        ( new DistMatrix<float,VC,STAR>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
VC_STAR_Double
Register_VC_STAR_Double
( int height, int width, int colAlignment, double* buffer, int ldim, Grid g )
{
    VC_STAR_Double A;
    CATCH(
        A = ::VC_STAR_DoubleList.size();
        ::VC_STAR_DoubleList.push_back
        ( new DistMatrix<double,VC,STAR>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
VC_STAR_SComplex
Register_VC_STAR_SComplex
( int height, int width, int colAlignment, SComplex* buffer, int ldim, Grid g )
{
    VC_STAR_SComplex A;
    CATCH(
        A = ::VC_STAR_SComplexList.size();
        ::VC_STAR_SComplexList.push_back
        ( new DistMatrix<scomplex,VC,STAR>
          ( height, width, colAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
VC_STAR_DComplex
Register_VC_STAR_DComplex
( int height, int width, int colAlignment, DComplex* buffer, int ldim, Grid g )
{
    VC_STAR_DComplex A;
    CATCH(
        A = ::VC_STAR_DComplexList.size();
        ::VC_STAR_DComplexList.push_back
        ( new DistMatrix<dcomplex,VC,STAR>
          ( height, width, colAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the col alignments
int
ColAlignment_VC_STAR_Single( VC_STAR_Single A )
{ return ::VC_STAR_SingleList[A]->ColAlignment(); }
int
ColAlignment_VC_STAR_Double( VC_STAR_Double A )
{ return ::VC_STAR_DoubleList[A]->ColAlignment(); }
#ifndef WITHOUT_COMPLEX
int
ColAlignment_VC_STAR_SComplex( VC_STAR_SComplex A )
{ return ::VC_STAR_SComplexList[A]->ColAlignment(); }
int
ColAlignment_VC_STAR_DComplex( VC_STAR_DComplex A )
{ return ::VC_STAR_DComplexList[A]->ColAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_VC_STAR_Single( char* msg, VC_STAR_Single A )
{
    CATCH(
        std::string s( msg );
        ::VC_STAR_SingleList[A]->Print( s );
    );
}
void 
Print_VC_STAR_Double( char* msg, VC_STAR_Double A )
{
    CATCH(
        std::string s( msg );
        ::VC_STAR_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_VC_STAR_SComplex( char* msg, VC_STAR_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::VC_STAR_SComplexList[A]->Print( s );
    );
}
void 
Print_VC_STAR_DComplex( char* msg, VC_STAR_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::VC_STAR_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [VR,* ] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
VR_STAR_Single
CreateEmpty_VR_STAR_Single( Grid g )
{
    VR_STAR_Single A;
    CATCH(
        A = ::VR_STAR_SingleList.size();
        ::VR_STAR_SingleList.push_back
        ( new DistMatrix<float,VR,STAR>( *::gridList[g] ) );
    );
    return A;
}
VR_STAR_Double
CreateEmpty_VR_STAR_Double( Grid g )
{
    VR_STAR_Double A;
    CATCH(
        A = ::VR_STAR_DoubleList.size();
        ::VR_STAR_DoubleList.push_back
        ( new DistMatrix<double,VR,STAR>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
VR_STAR_SComplex
CreateEmpty_VR_STAR_SComplex( Grid g )
{
    VR_STAR_SComplex A;
    CATCH(
        A = ::VR_STAR_SComplexList.size();
        ::VR_STAR_SComplexList.push_back
        ( new DistMatrix<scomplex,VR,STAR>( *::gridList[g] ) );
    );
    return A;
}
VR_STAR_DComplex
CreateEmpty_VR_STAR_DComplex( Grid g )
{
    VR_STAR_DComplex A;
    CATCH(
        A = ::VR_STAR_DComplexList.size();
        ::VR_STAR_DComplexList.push_back
        ( new DistMatrix<dcomplex,VR,STAR>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
VR_STAR_Single
Register_VR_STAR_Single
( int height, int width, int colAlignment, float* buffer, int ldim, Grid g )
{
    VR_STAR_Single A;
    CATCH(
        A = ::VR_STAR_SingleList.size();
        ::VR_STAR_SingleList.push_back
        ( new DistMatrix<float,VR,STAR>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
VR_STAR_Double
Register_VR_STAR_Double
( int height, int width, int colAlignment, double* buffer, int ldim, Grid g )
{
    VR_STAR_Double A;
    CATCH(
        A = ::VR_STAR_DoubleList.size();
        ::VR_STAR_DoubleList.push_back
        ( new DistMatrix<double,VR,STAR>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
VR_STAR_SComplex
Register_VR_STAR_SComplex
( int height, int width, int colAlignment, SComplex* buffer, int ldim, Grid g )
{
    VR_STAR_SComplex A;
    CATCH(
        A = ::VR_STAR_SComplexList.size();
        ::VR_STAR_SComplexList.push_back
        ( new DistMatrix<scomplex,VR,STAR>
          ( height, width, colAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
VR_STAR_DComplex
Register_VR_STAR_DComplex
( int height, int width, int colAlignment, DComplex* buffer, int ldim, Grid g )
{
    VR_STAR_DComplex A;
    CATCH(
        A = ::VR_STAR_DComplexList.size();
        ::VR_STAR_DComplexList.push_back
        ( new DistMatrix<dcomplex,VR,STAR>
          ( height, width, colAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the col alignments
int
ColAlignment_VR_STAR_Single( VR_STAR_Single A )
{ return ::VR_STAR_SingleList[A]->ColAlignment(); }
int
ColAlignment_VR_STAR_Double( VR_STAR_Double A )
{ return ::VR_STAR_DoubleList[A]->ColAlignment(); }
#ifndef WITHOUT_COMPLEX
int
ColAlignment_VR_STAR_SComplex( VR_STAR_SComplex A )
{ return ::VR_STAR_SComplexList[A]->ColAlignment(); }
int
ColAlignment_VR_STAR_DComplex( VR_STAR_DComplex A )
{ return ::VR_STAR_DComplexList[A]->ColAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_VR_STAR_Single( char* msg, VR_STAR_Single A )
{
    CATCH(
        std::string s( msg );
        ::VR_STAR_SingleList[A]->Print( s );
    );
}
void 
Print_VR_STAR_Double( char* msg, VR_STAR_Double A )
{
    CATCH(
        std::string s( msg );
        ::VR_STAR_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_VR_STAR_SComplex( char* msg, VR_STAR_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::VR_STAR_SComplexList[A]->Print( s );
    );
}
void 
Print_VR_STAR_DComplex( char* msg, VR_STAR_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::VR_STAR_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX

// Utilities
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

//----------------------------------------------------------------------------//
// LAPACK-level interface                                                     //
//----------------------------------------------------------------------------//
// Cholesky
void CholSingle( char uplo, MC_MR_Single A )
{
    CATCH(
        elemental::Shape shape = elemental::CharToShape( uplo );
        elemental::advanced::Chol( shape, *::MC_MR_SingleList[A] );
    );
}
void CholDouble( char uplo, MC_MR_Double A )
{
    CATCH(
        elemental::Shape shape = elemental::CharToShape( uplo );
        elemental::advanced::Chol( shape, *::MC_MR_DoubleList[A] );
    );
}
#ifndef WITHOUT_COMPLEX
void CholSComplex( char uplo, MC_MR_SComplex A )
{
    CATCH(
        elemental::Shape shape = elemental::CharToShape( uplo );
        elemental::advanced::Chol( shape, *::MC_MR_SComplexList[A] );
    );
}
void CholDComplex( char uplo, MC_MR_DComplex A )
{
    CATCH(
        elemental::Shape shape = elemental::CharToShape( uplo );
        elemental::advanced::Chol( shape, *::MC_MR_DComplexList[A] );
    );
}
#endif

// HermitianEig
#ifndef WITHOUT_PMRRR
void
HermitianEigDouble
( char uplo, int tryForHighAccuracy, 
  char job, int a, int b, double u, double v,
  MC_MR_Double A, STAR_VR_Double w, MC_MR_Double Z )
{
    elemental::Shape shape;
    CATCH(shape = elemental::CharToShape( uplo ););
    if( job == 'a' || job == 'A' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DoubleList[A], *::STAR_VR_DoubleList[w],
              *::MC_MR_DoubleList[Z], tryForHighAccuracy );
        );
    }
    else if( job == 'i' || job == 'I' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DoubleList[A], *::STAR_VR_DoubleList[w],
              *::MC_MR_DoubleList[Z], a, b, tryForHighAccuracy );
        );
    }
    else if( job == 'v' || job == 'V' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DoubleList[A], *::STAR_VR_DoubleList[w],
              *::MC_MR_DoubleList[Z], u, v, tryForHighAccuracy );
        );
    }
    else
    {
        std::cerr << "job must be in { a, A, i, I, v, V }" << std::endl;
    }
}

void
HermitianEigDouble_OnlyEigvals
( char uplo, int tryForHighAccuracy,
  char job, int a, int b, double u, double v,
  MC_MR_Double A, STAR_VR_Double w )
{
    elemental::Shape shape;
    CATCH(shape = elemental::CharToShape( uplo ););
    if( job == 'a' || job == 'A' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DoubleList[A], *::STAR_VR_DoubleList[w],
              tryForHighAccuracy );
        );
    }
    else if( job == 'i' || job == 'I' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DoubleList[A], *::STAR_VR_DoubleList[w],
              a, b, tryForHighAccuracy );
        );
    }
    else if( job == 'v' || job == 'V' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DoubleList[A], *::STAR_VR_DoubleList[w],
              u, v, tryForHighAccuracy );
        );
    }
    else
    {
        std::cerr << "job must be in { a, A, i, I, v, V }" << std::endl;
    }
}

#ifndef WITHOUT_COMPLEX
void
HermitianEigDComplex
( char uplo, int tryForHighAccuracy,
  char job, int a, int b, double u, double v,
  MC_MR_DComplex A, STAR_VR_Double w, MC_MR_DComplex Z )
{
    elemental::Shape shape;
    CATCH(shape = elemental::CharToShape( uplo ););
    if( job == 'a' || job == 'A' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DComplexList[A], *::STAR_VR_DoubleList[w],
              *::MC_MR_DComplexList[Z], tryForHighAccuracy );
        );
    }
    else if( job == 'i' || job == 'I' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DComplexList[A], *::STAR_VR_DoubleList[w],
              *::MC_MR_DComplexList[Z], a, b, tryForHighAccuracy );
        );
    }
    else if( job == 'v' || job == 'V' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DComplexList[A], *::STAR_VR_DoubleList[w],
              *::MC_MR_DComplexList[Z], u, v, tryForHighAccuracy );
        );
    }
    else
    {
        std::cerr << "job must be in { a, A, i, I, v, V }" << std::endl;
    }
}

void
HermitianEigDComplex_OnlyEigvals
( char uplo, int tryForHighAccuracy, 
  char job, int a, int b, double u, double v,
  MC_MR_DComplex A, STAR_VR_Double w )
{
    elemental::Shape shape;
    CATCH(shape = elemental::CharToShape( uplo ););
    if( job == 'a' || job == 'A' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DComplexList[A], *::STAR_VR_DoubleList[w],
              tryForHighAccuracy );
        );
    }
    else if( job == 'i' || job == 'I' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DComplexList[A], *::STAR_VR_DoubleList[w],
              a, b, tryForHighAccuracy );
        );
    }
    else if( job == 'v' || job == 'V' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DComplexList[A], *::STAR_VR_DoubleList[w],
              u, v, tryForHighAccuracy );
        );
    }
    else
    {
        std::cerr << "job must be in { a, A, i, I, v, V }" << std::endl;
    }
}
#endif // WITHOUT_COMPLEX 

void
HermitianGenDefiniteEigDouble
( int genEigInt, char uplo, int tryForHighAccuracy, 
  char job, int a, int b, double u, double v,
  MC_MR_Double A, MC_MR_Double B, STAR_VR_Double w, MC_MR_Double Z )
{
    elemental::HermitianGenDefiniteEigType genEigType;
    if( genEigInt == 1 )
        genEigType = elemental::AXBX;
    else if( genEigInt == 2 )
        genEigType = elemental::ABX;
    else if( genEigInt == 3 )
        genEigType = elemental::BAX;
    else
    {
        std::cerr << "Invalid genEigType, choose from {1,2,3}" << std::endl;
        return;
    }
    elemental::Shape shape;
    CATCH(shape = elemental::CharToShape( uplo ););
    if( job == 'a' || job == 'A' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DoubleList[A], *::MC_MR_DoubleList[B], 
              *::STAR_VR_DoubleList[w], *::MC_MR_DoubleList[Z], 
              tryForHighAccuracy );
        );
    }
    else if( job == 'i' || job == 'I' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DoubleList[A], *::MC_MR_DoubleList[B], 
              *::STAR_VR_DoubleList[w], *::MC_MR_DoubleList[Z], 
              a, b, tryForHighAccuracy );
        );
    }
    else if( job == 'v' || job == 'V' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DoubleList[A], *::MC_MR_DoubleList[B], 
              *::STAR_VR_DoubleList[w], *::MC_MR_DoubleList[Z], 
              u, v, tryForHighAccuracy );
        );
    }
    else
    {
        std::cerr << "job must be in { a, A, i, I, v, V }" << std::endl;
    }
}

void
HermitianGenDefiniteEigDouble_OnlyEigvals
( int genEigInt, char uplo, int tryForHighAccuracy, 
  char job, int a, int b, double u, double v,
  MC_MR_Double A, MC_MR_Double B, STAR_VR_Double w )
{
    elemental::advanced::GenEigType genEigType;
    if( genEigInt == 1 )
    {
        genEigType = elemental::advanced::AXBX;
    }
    else if( genEigInt == 2 )
    {
        genEigType = elemental::advanced::ABX;
    }
    else if( genEigInt == 3 )
    {
        genEigType = elemental::advanced::BAX;
    }
    else
    {
        std::cerr << "Invalid genEigType, choose from {1,2,3}" << std::endl;
        return;
    }
    elemental::Shape shape;
    CATCH(shape = elemental::CharToShape( uplo ););
    if( job == 'a' || job == 'A' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DoubleList[A], *::MC_MR_DoubleList[B], 
              *::STAR_VR_DoubleList[w], tryForHighAccuracy );
        );
    }
    else if( job == 'i' || job == 'I' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DoubleList[A], *::MC_MR_DoubleList[B], 
              *::STAR_VR_DoubleList[w], a, b, tryForHighAccuracy );
        );
    }
    else if( job == 'v' || job == 'V' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DoubleList[A], *::MC_MR_DoubleList[B], 
              *::STAR_VR_DoubleList[w], u, v, tryForHighAccuracy );
        );
    }
    else
    {
        std::cerr << "job must be in { a, A, i, I, v, V }" << std::endl;
    }
}

#ifndef WITHOUT_COMPLEX
void
HermitianGenDefiniteEigDComplex
( int genEigInt, char uplo, int tryForHighAccuracy, 
  char job, int a, int b, double u, double v,
  MC_MR_DComplex A, MC_MR_DComplex B, STAR_VR_Double w, MC_MR_DComplex Z )
{
    elemental::advanced::GenEigType genEigType;
    if( genEigInt == 1 )
    {
        genEigType = elemental::advanced::AXBX;
    }
    else if( genEigInt == 2 )
    {
        genEigType = elemental::advanced::ABX;
    }
    else if( genEigInt == 3 )
    {
        genEigType = elemental::advanced::BAX;
    }
    else
    {
        std::cerr << "Invalid genEigType, choose from {1,2,3}" << std::endl;
        return;
    }
    elemental::Shape shape;
    CATCH(shape = elemental::CharToShape( uplo ););
    if( job == 'a' || job == 'A' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DComplexList[A], *::MC_MR_DComplexList[B], 
              *::STAR_VR_DoubleList[w], *::MC_MR_DComplexList[Z], 
              tryForHighAccuracy );
        );
    }
    else if( job == 'i' || job == 'I' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DComplexList[A], *::MC_MR_DComplexList[B], 
              *::STAR_VR_DoubleList[w], *::MC_MR_DComplexList[Z], 
              a, b, tryForHighAccuracy );
        );
    }
    else if( job == 'v' || job == 'V' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DComplexList[A], *::MC_MR_DComplexList[B], 
              *::STAR_VR_DoubleList[w], *::MC_MR_DComplexList[Z], 
              u, v, tryForHighAccuracy );
        );
    }
    else
    {
        std::cerr << "job must be in { a, A, i, I, v, V }" << std::endl;
    }
}

void
HermitianGenDefiniteEigDComplex_OnlyEigvals
( int genEigInt, char uplo, int tryForHighAccuracy, 
  char job, int a, int b, double u, double v,
  MC_MR_DComplex A, MC_MR_DComplex B, STAR_VR_Double w )
{
    elemental::advanced::GenEigType genEigType;
    if( genEigInt == 1 )
    {
        genEigType = elemental::advanced::AXBX;
    }
    else if( genEigInt == 2 )
    {
        genEigType = elemental::advanced::ABX;
    }
    else if( genEigInt == 3 )
    {
        genEigType = elemental::advanced::BAX;
    }
    else
    {
        std::cerr << "Invalid genEigType, choose from {1,2,3}" << std::endl;
        return;
    }
    elemental::Shape shape;
    CATCH(shape = elemental::CharToShape( uplo ););
    if( job == 'a' || job == 'A' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DComplexList[A], *::MC_MR_DComplexList[B], 
              *::STAR_VR_DoubleList[w], tryForHighAccuracy );
        );
    }
    else if( job == 'i' || job == 'I' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DComplexList[A], *::MC_MR_DComplexList[B], 
              *::STAR_VR_DoubleList[w], a, b, tryForHighAccuracy );
        );
    }
    else if( job == 'v' || job == 'V' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DComplexList[A], *::MC_MR_DComplexList[B], 
              *::STAR_VR_DoubleList[w], u, v, tryForHighAccuracy );
        );
    }
    else
    {
        std::cerr << "job must be in { a, A, i, I, v, V }" << std::endl;
    }
}
#endif // WITHOUT_COMPLEX 
#endif // WITHOUT_PMRRR 

} // extern "C"

