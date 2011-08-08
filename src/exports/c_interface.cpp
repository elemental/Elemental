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
using elemental::Star;
using elemental::VC;
using elemental::VR;
using elemental::scomplex;
using elemental::dcomplex;
using elemental::DistMatrix;

namespace {
std::vector<elemental::Grid*> gridList;
std::vector<DistMatrix<float, MC,  MR  >*> MC_MR_SingleList;
std::vector<DistMatrix<double,MC,  MR  >*> MC_MR_DoubleList;
std::vector<DistMatrix<float, MC,  Star>*> MC_Star_SingleList;
std::vector<DistMatrix<double,MC,  Star>*> MC_Star_DoubleList;
std::vector<DistMatrix<float, MD,  Star>*> MD_Star_SingleList;
std::vector<DistMatrix<double,MD,  Star>*> MD_Star_DoubleList;
std::vector<DistMatrix<float, MR,  MC  >*> MR_MC_SingleList;
std::vector<DistMatrix<double,MR,  MC  >*> MR_MC_DoubleList;
std::vector<DistMatrix<float, MR,  Star>*> MR_Star_SingleList;
std::vector<DistMatrix<double,MR,  Star>*> MR_Star_DoubleList;
std::vector<DistMatrix<float, Star,MC  >*> Star_MC_SingleList;
std::vector<DistMatrix<double,Star,MC  >*> Star_MC_DoubleList;
std::vector<DistMatrix<float, Star,MD  >*> Star_MD_SingleList;
std::vector<DistMatrix<double,Star,MD  >*> Star_MD_DoubleList;
std::vector<DistMatrix<float, Star,MR  >*> Star_MR_SingleList;
std::vector<DistMatrix<double,Star,MR  >*> Star_MR_DoubleList;
std::vector<DistMatrix<float, Star,Star>*> Star_Star_SingleList;
std::vector<DistMatrix<double,Star,Star>*> Star_Star_DoubleList;
std::vector<DistMatrix<float, Star,VC  >*> Star_VC_SingleList;
std::vector<DistMatrix<double,Star,VC  >*> Star_VC_DoubleList;
std::vector<DistMatrix<float, Star,VR  >*> Star_VR_SingleList;
std::vector<DistMatrix<double,Star,VR  >*> Star_VR_DoubleList;
std::vector<DistMatrix<float, VC,  Star>*> VC_Star_SingleList;
std::vector<DistMatrix<double,VC,  Star>*> VC_Star_DoubleList;
std::vector<DistMatrix<float, VR,  Star>*> VR_Star_SingleList;
std::vector<DistMatrix<double,VR,  Star>*> VR_Star_DoubleList;
#ifndef WITHOUT_COMPLEX
std::vector<DistMatrix<scomplex,MC,  MR  >*> MC_MR_SComplexList;
std::vector<DistMatrix<dcomplex,MC,  MR  >*> MC_MR_DComplexList;
std::vector<DistMatrix<scomplex,MC,  Star>*> MC_Star_SComplexList;
std::vector<DistMatrix<dcomplex,MC,  Star>*> MC_Star_DComplexList;
std::vector<DistMatrix<scomplex,MD,  Star>*> MD_Star_SComplexList;
std::vector<DistMatrix<dcomplex,MD,  Star>*> MD_Star_DComplexList;
std::vector<DistMatrix<scomplex,MR,  MC  >*> MR_MC_SComplexList;
std::vector<DistMatrix<dcomplex,MR,  MC  >*> MR_MC_DComplexList;
std::vector<DistMatrix<scomplex,MR,  Star>*> MR_Star_SComplexList;
std::vector<DistMatrix<dcomplex,MR,  Star>*> MR_Star_DComplexList;
std::vector<DistMatrix<scomplex,Star,MC  >*> Star_MC_SComplexList;
std::vector<DistMatrix<dcomplex,Star,MC  >*> Star_MC_DComplexList;
std::vector<DistMatrix<scomplex,Star,MD  >*> Star_MD_SComplexList;
std::vector<DistMatrix<dcomplex,Star,MD  >*> Star_MD_DComplexList;
std::vector<DistMatrix<scomplex,Star,MR  >*> Star_MR_SComplexList;
std::vector<DistMatrix<dcomplex,Star,MR  >*> Star_MR_DComplexList;
std::vector<DistMatrix<scomplex,Star,Star>*> Star_Star_SComplexList;
std::vector<DistMatrix<dcomplex,Star,Star>*> Star_Star_DComplexList;
std::vector<DistMatrix<scomplex,Star,VC  >*> Star_VC_SComplexList;
std::vector<DistMatrix<dcomplex,Star,VC  >*> Star_VC_DComplexList;
std::vector<DistMatrix<scomplex,Star,VR  >*> Star_VR_SComplexList;
std::vector<DistMatrix<dcomplex,Star,VR  >*> Star_VR_DComplexList;
std::vector<DistMatrix<scomplex,VC,  Star>*> VC_Star_SComplexList;
std::vector<DistMatrix<dcomplex,VC,  Star>*> VC_Star_DComplexList;
std::vector<DistMatrix<scomplex,VR,  Star>*> VR_Star_SComplexList;
std::vector<DistMatrix<dcomplex,VR,  Star>*> VR_Star_DComplexList;
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
    for( unsigned j=0; j< ::MC_Star_SingleList.size(); ++j )
        delete ::MC_Star_SingleList[j];
    for( unsigned j=0; j< ::MC_Star_DoubleList.size(); ++j )
        delete ::MC_Star_DoubleList[j];
    ::MC_Star_SingleList.resize(0);
    ::MC_Star_DoubleList.resize(0);
    // Real [MD,* ]
    for( unsigned j=0; j< ::MD_Star_SingleList.size(); ++j )
        delete ::MD_Star_SingleList[j];
    for( unsigned j=0; j< ::MD_Star_DoubleList.size(); ++j )
        delete ::MD_Star_DoubleList[j];
    ::MD_Star_SingleList.resize(0);
    ::MD_Star_DoubleList.resize(0);
    // Real [MR,MC]
    for( unsigned j=0; j< ::MR_MC_SingleList.size(); ++j )
        delete ::MR_MC_SingleList[j];
    for( unsigned j=0; j< ::MR_MC_DoubleList.size(); ++j )
        delete ::MR_MC_DoubleList[j];
    ::MR_MC_SingleList.resize(0);
    ::MR_MC_DoubleList.resize(0);
    // Real [MR,* ]
    for( unsigned j=0; j< ::MR_Star_SingleList.size(); ++j )
        delete ::MR_Star_SingleList[j];
    for( unsigned j=0; j< ::MR_Star_DoubleList.size(); ++j )
        delete ::MR_Star_DoubleList[j];
    ::MR_Star_SingleList.resize(0);
    ::MR_Star_DoubleList.resize(0);
    // Real [* ,MC]
    for( unsigned j=0; j< ::Star_MC_SingleList.size(); ++j )
        delete ::Star_MC_SingleList[j];
    for( unsigned j=0; j< ::Star_MC_DoubleList.size(); ++j )
        delete ::Star_MC_DoubleList[j];
    ::Star_MC_SingleList.resize(0);
    ::Star_MC_DoubleList.resize(0);
    // Real [* ,MD]
    for( unsigned j=0; j< ::Star_MD_SingleList.size(); ++j )
        delete ::Star_MD_SingleList[j];
    for( unsigned j=0; j< ::Star_MD_DoubleList.size(); ++j )
        delete ::Star_MD_DoubleList[j];
    ::Star_MD_SingleList.resize(0);
    ::Star_MD_DoubleList.resize(0);
    // Real [* ,MR]
    for( unsigned j=0; j< ::Star_MR_SingleList.size(); ++j )
        delete ::Star_MR_SingleList[j];
    for( unsigned j=0; j< ::Star_MR_DoubleList.size(); ++j )
        delete ::Star_MR_DoubleList[j];
    ::Star_MR_SingleList.resize(0);
    ::Star_MR_DoubleList.resize(0);
    // Real [* ,* ]
    for( unsigned j=0; j< ::Star_Star_SingleList.size(); ++j )
        delete ::Star_Star_SingleList[j];
    for( unsigned j=0; j< ::Star_Star_DoubleList.size(); ++j )
        delete ::Star_Star_DoubleList[j];
    ::Star_Star_SingleList.resize(0);
    ::Star_Star_DoubleList.resize(0);
    // Real [* ,VC]
    for( unsigned j=0; j< ::Star_VC_SingleList.size(); ++j )
        delete ::Star_VC_SingleList[j];
    for( unsigned j=0; j< ::Star_VC_DoubleList.size(); ++j )
        delete ::Star_VC_DoubleList[j];
    ::Star_VC_SingleList.resize(0);
    ::Star_VC_DoubleList.resize(0);
    // Real [* ,VR]
    for( unsigned j=0; j< ::Star_VR_SingleList.size(); ++j )
        delete ::Star_VR_SingleList[j];
    for( unsigned j=0; j< ::Star_VR_DoubleList.size(); ++j )
        delete ::Star_VR_DoubleList[j];
    ::Star_VR_SingleList.resize(0);
    ::Star_VR_DoubleList.resize(0);
    // Real [VC,* ]
    for( unsigned j=0; j< ::VC_Star_SingleList.size(); ++j )
        delete ::VC_Star_SingleList[j];
    for( unsigned j=0; j< ::VC_Star_DoubleList.size(); ++j )
        delete ::VC_Star_DoubleList[j];
    ::VC_Star_SingleList.resize(0);
    ::VC_Star_DoubleList.resize(0);
    // Real [VR,* ]
    for( unsigned j=0; j< ::VR_Star_SingleList.size(); ++j )
        delete ::VR_Star_SingleList[j];
    for( unsigned j=0; j< ::VR_Star_DoubleList.size(); ++j )
        delete ::VR_Star_DoubleList[j];
    ::VR_Star_SingleList.resize(0);
    ::VR_Star_DoubleList.resize(0);
#ifndef WITHOUT_COMPLEX
    // Complex [MC,MR]
    for( unsigned j=0; j< ::MC_MR_SComplexList.size(); ++j )
        delete ::MC_MR_SComplexList[j];
    for( unsigned j=0; j< ::MC_MR_DComplexList.size(); ++j )
        delete ::MC_MR_DComplexList[j];
    ::MC_MR_SComplexList.resize(0);
    ::MC_MR_DComplexList.resize(0);
    // Complex [MC,* ]
    for( unsigned j=0; j< ::MC_Star_SComplexList.size(); ++j )
        delete ::MC_Star_SComplexList[j];
    for( unsigned j=0; j< ::MC_Star_DComplexList.size(); ++j )
        delete ::MC_Star_DComplexList[j];
    ::MC_Star_SComplexList.resize(0);
    ::MC_Star_DComplexList.resize(0);
    // Complex [MD,* ]
    for( unsigned j=0; j< ::MD_Star_SComplexList.size(); ++j )
        delete ::MD_Star_SComplexList[j];
    for( unsigned j=0; j< ::MD_Star_DComplexList.size(); ++j )
        delete ::MD_Star_DComplexList[j];
    ::MD_Star_SComplexList.resize(0);
    ::MD_Star_DComplexList.resize(0);
    // Complex [MR,MC]
    for( unsigned j=0; j< ::MR_MC_SComplexList.size(); ++j )
        delete ::MR_MC_SComplexList[j];
    for( unsigned j=0; j< ::MR_MC_DComplexList.size(); ++j )
        delete ::MR_MC_DComplexList[j];
    ::MR_MC_SComplexList.resize(0);
    ::MR_MC_DComplexList.resize(0);
    // Complex [MR,* ]
    for( unsigned j=0; j< ::MR_Star_SComplexList.size(); ++j )
        delete ::MR_Star_SComplexList[j];
    for( unsigned j=0; j< ::MR_Star_DComplexList.size(); ++j )
        delete ::MR_Star_DComplexList[j];
    ::MR_Star_SComplexList.resize(0);
    ::MR_Star_DComplexList.resize(0);
    // Complex [* ,MC]
    for( unsigned j=0; j< ::Star_MC_SComplexList.size(); ++j )
        delete ::Star_MC_SComplexList[j];
    for( unsigned j=0; j< ::Star_MC_DComplexList.size(); ++j )
        delete ::Star_MC_DComplexList[j];
    ::Star_MC_SComplexList.resize(0);
    ::Star_MC_DComplexList.resize(0);
    // Complex [* ,MD]
    for( unsigned j=0; j< ::Star_MD_SComplexList.size(); ++j )
        delete ::Star_MD_SComplexList[j];
    for( unsigned j=0; j< ::Star_MD_DComplexList.size(); ++j )
        delete ::Star_MD_DComplexList[j];
    ::Star_MD_SComplexList.resize(0);
    ::Star_MD_DComplexList.resize(0);
    // Complex [* ,MR]
    for( unsigned j=0; j< ::Star_MR_SComplexList.size(); ++j )
        delete ::Star_MR_SComplexList[j];
    for( unsigned j=0; j< ::Star_MR_DComplexList.size(); ++j )
        delete ::Star_MR_DComplexList[j];
    ::Star_MR_SComplexList.resize(0);
    ::Star_MR_DComplexList.resize(0);
    // Complex [* ,* ]
    for( unsigned j=0; j< ::Star_Star_SComplexList.size(); ++j )
        delete ::Star_Star_SComplexList[j];
    for( unsigned j=0; j< ::Star_Star_DComplexList.size(); ++j )
        delete ::Star_Star_DComplexList[j];
    ::Star_Star_SComplexList.resize(0);
    ::Star_Star_DComplexList.resize(0);
    // Complex [* ,VC]
    for( unsigned j=0; j< ::Star_VC_SComplexList.size(); ++j )
        delete ::Star_VC_SComplexList[j];
    for( unsigned j=0; j< ::Star_VC_DComplexList.size(); ++j )
        delete ::Star_VC_DComplexList[j];
    ::Star_VC_SComplexList.resize(0);
    ::Star_VC_DComplexList.resize(0);
    // Complex [* ,VR]
    for( unsigned j=0; j< ::Star_VR_SComplexList.size(); ++j )
        delete ::Star_VR_SComplexList[j];
    for( unsigned j=0; j< ::Star_VR_DComplexList.size(); ++j )
        delete ::Star_VR_DComplexList[j];
    ::Star_VR_SComplexList.resize(0);
    ::Star_VR_DComplexList.resize(0);
    // Complex [VC,* ]
    for( unsigned j=0; j< ::VC_Star_SComplexList.size(); ++j )
        delete ::VC_Star_SComplexList[j];
    for( unsigned j=0; j< ::VC_Star_DComplexList.size(); ++j )
        delete ::VC_Star_DComplexList[j];
    ::VC_Star_SComplexList.resize(0);
    ::VC_Star_DComplexList.resize(0);
    // Complex [VR,* ]
    for( unsigned j=0; j< ::VR_Star_SComplexList.size(); ++j )
        delete ::VR_Star_SComplexList[j];
    for( unsigned j=0; j< ::VR_Star_DComplexList.size(); ++j )
        delete ::VR_Star_DComplexList[j];
    ::VR_Star_SComplexList.resize(0);
    ::VR_Star_DComplexList.resize(0);
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
MC_Star_Single
CreateEmpty_MC_Star_Single( Grid g )
{
    MC_Star_Single A;
    CATCH(
        A = ::MC_Star_SingleList.size();
        ::MC_Star_SingleList.push_back
        ( new DistMatrix<float,MC,Star>( *::gridList[g] ) );
    );
    return A;
}
MC_Star_Double
CreateEmpty_MC_Star_Double( Grid g )
{
    MC_Star_Double A;
    CATCH(
        A = ::MC_Star_DoubleList.size();
        ::MC_Star_DoubleList.push_back
        ( new DistMatrix<double,MC,Star>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
MC_Star_SComplex
CreateEmpty_MC_Star_SComplex( Grid g )
{
    MC_Star_SComplex A;
    CATCH(
        A = ::MC_Star_SComplexList.size();
        ::MC_Star_SComplexList.push_back
        ( new DistMatrix<scomplex,MC,Star>( *::gridList[g] ) );
    );
    return A;
}
MC_Star_DComplex
CreateEmpty_MC_Star_DComplex( Grid g )
{
    MC_Star_DComplex A;
    CATCH(
        A = ::MC_Star_DComplexList.size();
        ::MC_Star_DComplexList.push_back
        ( new DistMatrix<dcomplex,MC,Star>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
MC_Star_Single
Register_MC_Star_Single
( int height, int width, int colAlignment, float* buffer, int ldim, Grid g )
{
    MC_Star_Single A;
    CATCH(
        A = ::MC_Star_SingleList.size();
        ::MC_Star_SingleList.push_back
        ( new DistMatrix<float,MC,Star>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
MC_Star_Double
Register_MC_Star_Double
( int height, int width, int colAlignment, double* buffer, int ldim, Grid g )
{
    MC_Star_Double A;
    CATCH(
        A = ::MC_Star_DoubleList.size();
        ::MC_Star_DoubleList.push_back
        ( new DistMatrix<double,MC,Star>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
MC_Star_SComplex
Register_MC_Star_SComplex
( int height, int width, int colAlignment, SComplex* buffer, int ldim, Grid g )
{
    MC_Star_SComplex A;
    CATCH(
        A = ::MC_Star_SComplexList.size();
        ::MC_Star_SComplexList.push_back
        ( new DistMatrix<scomplex,MC,Star>
          ( height, width, colAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
MC_Star_DComplex
Register_MC_Star_DComplex
( int height, int width, int colAlignment, DComplex* buffer, int ldim, Grid g )
{
    MC_Star_DComplex A;
    CATCH(
        A = ::MC_Star_DComplexList.size();
        ::MC_Star_DComplexList.push_back
        ( new DistMatrix<dcomplex,MC,Star>
          ( height, width, colAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the col alignments
int
ColAlignment_MC_Star_Single( MC_Star_Single A )
{ return ::MC_Star_SingleList[A]->ColAlignment(); }
int
ColAlignment_MC_Star_Double( MC_Star_Double A )
{ return ::MC_Star_DoubleList[A]->ColAlignment(); }
#ifndef WITHOUT_COMPLEX
int
ColAlignment_MC_Star_SComplex( MC_Star_SComplex A )
{ return ::MC_Star_SComplexList[A]->ColAlignment(); }
int
ColAlignment_MC_Star_DComplex( MC_Star_DComplex A )
{ return ::MC_Star_DComplexList[A]->ColAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_MC_Star_Single( char* msg, MC_Star_Single A )
{
    CATCH(
        std::string s( msg );
        ::MC_Star_SingleList[A]->Print( s );
    );
}
void 
Print_MC_Star_Double( char* msg, MC_Star_Double A )
{
    CATCH(
        std::string s( msg );
        ::MC_Star_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_MC_Star_SComplex( char* msg, MC_Star_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::MC_Star_SComplexList[A]->Print( s );
    );
}
void 
Print_MC_Star_DComplex( char* msg, MC_Star_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::MC_Star_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [MD,* ] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
MD_Star_Single
CreateEmpty_MD_Star_Single( Grid g )
{
    MD_Star_Single A;
    CATCH(
        A = ::MD_Star_SingleList.size();
        ::MD_Star_SingleList.push_back
        ( new DistMatrix<float,MD,Star>( *::gridList[g] ) );
    );
    return A;
}
MD_Star_Double
CreateEmpty_MD_Star_Double( Grid g )
{
    MD_Star_Double A;
    CATCH(
        A = ::MD_Star_DoubleList.size();
        ::MD_Star_DoubleList.push_back
        ( new DistMatrix<double,MD,Star>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
MD_Star_SComplex
CreateEmpty_MD_Star_SComplex( Grid g )
{
    MD_Star_SComplex A;
    CATCH(
        A = ::MD_Star_SComplexList.size();
        ::MD_Star_SComplexList.push_back
        ( new DistMatrix<scomplex,MD,Star>( *::gridList[g] ) );
    );
    return A;
}
MD_Star_DComplex
CreateEmpty_MD_Star_DComplex( Grid g )
{
    MD_Star_DComplex A;
    CATCH(
        A = ::MD_Star_DComplexList.size();
        ::MD_Star_DComplexList.push_back
        ( new DistMatrix<dcomplex,MD,Star>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
MD_Star_Single
Register_MD_Star_Single
( int height, int width, int colAlignment, float* buffer, int ldim, Grid g )
{
    MD_Star_Single A;
    CATCH(
        A = ::MD_Star_SingleList.size();
        ::MD_Star_SingleList.push_back
        ( new DistMatrix<float,MD,Star>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
MD_Star_Double
Register_MD_Star_Double
( int height, int width, int colAlignment, double* buffer, int ldim, Grid g )
{
    MD_Star_Double A;
    CATCH(
        A = ::MD_Star_DoubleList.size();
        ::MD_Star_DoubleList.push_back
        ( new DistMatrix<double,MD,Star>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
MD_Star_SComplex
Register_MD_Star_SComplex
( int height, int width, int colAlignment, SComplex* buffer, int ldim, Grid g )
{
    MD_Star_SComplex A;
    CATCH(
        A = ::MD_Star_SComplexList.size();
        ::MD_Star_SComplexList.push_back
        ( new DistMatrix<scomplex,MD,Star>
          ( height, width, colAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
MD_Star_DComplex
Register_MD_Star_DComplex
( int height, int width, int colAlignment, DComplex* buffer, int ldim, Grid g )
{
    MD_Star_DComplex A;
    CATCH(
        A = ::MD_Star_DComplexList.size();
        ::MD_Star_DComplexList.push_back
        ( new DistMatrix<dcomplex,MD,Star>
          ( height, width, colAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the col alignments
int
ColAlignment_MD_Star_Single( MD_Star_Single A )
{ return ::MD_Star_SingleList[A]->ColAlignment(); }
int
ColAlignment_MD_Star_Double( MD_Star_Double A )
{ return ::MD_Star_DoubleList[A]->ColAlignment(); }
#ifndef WITHOUT_COMPLEX
int
ColAlignment_MD_Star_SComplex( MD_Star_SComplex A )
{ return ::MD_Star_SComplexList[A]->ColAlignment(); }
int
ColAlignment_MD_Star_DComplex( MD_Star_DComplex A )
{ return ::MD_Star_DComplexList[A]->ColAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_MD_Star_Single( char* msg, MD_Star_Single A )
{
    CATCH(
        std::string s( msg );
        ::MD_Star_SingleList[A]->Print( s );
    );
}
void 
Print_MD_Star_Double( char* msg, MD_Star_Double A )
{
    CATCH(
        std::string s( msg );
        ::MD_Star_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_MD_Star_SComplex( char* msg, MD_Star_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::MD_Star_SComplexList[A]->Print( s );
    );
}
void 
Print_MD_Star_DComplex( char* msg, MD_Star_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::MD_Star_DComplexList[A]->Print( s );
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
MR_Star_Single
CreateEmpty_MR_Star_Single( Grid g )
{
    MR_Star_Single A;
    CATCH(
        A = ::MR_Star_SingleList.size();
        ::MR_Star_SingleList.push_back
        ( new DistMatrix<float,MR,Star>( *::gridList[g] ) );
    );
    return A;
}
MR_Star_Double
CreateEmpty_MR_Star_Double( Grid g )
{
    MR_Star_Double A;
    CATCH(
        A = ::MR_Star_DoubleList.size();
        ::MR_Star_DoubleList.push_back
        ( new DistMatrix<double,MR,Star>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
MR_Star_SComplex
CreateEmpty_MR_Star_SComplex( Grid g )
{
    MR_Star_SComplex A;
    CATCH(
        A = ::MR_Star_SComplexList.size();
        ::MR_Star_SComplexList.push_back
        ( new DistMatrix<scomplex,MR,Star>( *::gridList[g] ) );
    );
    return A;
}
MR_Star_DComplex
CreateEmpty_MR_Star_DComplex( Grid g )
{
    MR_Star_DComplex A;
    CATCH(
        A = ::MR_Star_DComplexList.size();
        ::MR_Star_DComplexList.push_back
        ( new DistMatrix<dcomplex,MR,Star>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
MR_Star_Single
Register_MR_Star_Single
( int height, int width, int colAlignment, float* buffer, int ldim, Grid g )
{
    MR_Star_Single A;
    CATCH(
        A = ::MR_Star_SingleList.size();
        ::MR_Star_SingleList.push_back
        ( new DistMatrix<float,MR,Star>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
MR_Star_Double
Register_MR_Star_Double
( int height, int width, int colAlignment, double* buffer, int ldim, Grid g )
{
    MR_Star_Double A;
    CATCH(
        A = ::MR_Star_DoubleList.size();
        ::MR_Star_DoubleList.push_back
        ( new DistMatrix<double,MR,Star>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
MR_Star_SComplex
Register_MR_Star_SComplex
( int height, int width, int colAlignment, SComplex* buffer, int ldim, Grid g )
{
    MR_Star_SComplex A;
    CATCH(
        A = ::MR_Star_SComplexList.size();
        ::MR_Star_SComplexList.push_back
        ( new DistMatrix<scomplex,MR,Star>
          ( height, width, colAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
MR_Star_DComplex
Register_MR_Star_DComplex
( int height, int width, int colAlignment, DComplex* buffer, int ldim, Grid g )
{
    MR_Star_DComplex A;
    CATCH(
        A = ::MR_Star_DComplexList.size();
        ::MR_Star_DComplexList.push_back
        ( new DistMatrix<dcomplex,MR,Star>
          ( height, width, colAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the col alignments
int
ColAlignment_MR_Star_Single( MR_Star_Single A )
{ return ::MR_Star_SingleList[A]->ColAlignment(); }
int
ColAlignment_MR_Star_Double( MR_Star_Double A )
{ return ::MR_Star_DoubleList[A]->ColAlignment(); }
#ifndef WITHOUT_COMPLEX
int
ColAlignment_MR_Star_SComplex( MR_Star_SComplex A )
{ return ::MR_Star_SComplexList[A]->ColAlignment(); }
int
ColAlignment_MR_Star_DComplex( MR_Star_DComplex A )
{ return ::MR_Star_DComplexList[A]->ColAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_MR_Star_Single( char* msg, MR_Star_Single A )
{
    CATCH(
        std::string s( msg );
        ::MR_Star_SingleList[A]->Print( s );
    );
}
void 
Print_MR_Star_Double( char* msg, MR_Star_Double A )
{
    CATCH(
        std::string s( msg );
        ::MR_Star_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_MR_Star_SComplex( char* msg, MR_Star_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::MR_Star_SComplexList[A]->Print( s );
    );
}
void 
Print_MR_Star_DComplex( char* msg, MR_Star_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::MR_Star_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [* ,MC] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
Star_MC_Single
CreateEmpty_Star_MC_Single( Grid g )
{
    Star_MC_Single A;
    CATCH(
        A = ::Star_MC_SingleList.size();
        ::Star_MC_SingleList.push_back
        ( new DistMatrix<float,Star,MC>( *::gridList[g] ) );
    );
    return A;
}
Star_MC_Double
CreateEmpty_Star_MC_Double( Grid g )
{
    Star_MC_Double A;
    CATCH(
        A = ::Star_MC_DoubleList.size();
        ::Star_MC_DoubleList.push_back
        ( new DistMatrix<double,Star,MC>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
Star_MC_SComplex
CreateEmpty_Star_MC_SComplex( Grid g )
{
    Star_MC_SComplex A;
    CATCH(
        A = ::Star_MC_SComplexList.size();
        ::Star_MC_SComplexList.push_back
        ( new DistMatrix<scomplex,Star,MC>( *::gridList[g] ) );
    );
    return A;
}
Star_MC_DComplex
CreateEmpty_Star_MC_DComplex( Grid g )
{
    Star_MC_DComplex A;
    CATCH(
        A = ::Star_MC_DComplexList.size();
        ::Star_MC_DComplexList.push_back
        ( new DistMatrix<dcomplex,Star,MC>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
Star_MC_Single
Register_Star_MC_Single
( int height, int width, int rowAlignment, float* buffer, int ldim, Grid g )
{
    Star_MC_Single A;
    CATCH(
        A = ::Star_MC_SingleList.size();
        ::Star_MC_SingleList.push_back
        ( new DistMatrix<float,Star,MC>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
Star_MC_Double
Register_Star_MC_Double
( int height, int width, int rowAlignment, double* buffer, int ldim, Grid g )
{
    Star_MC_Double A;
    CATCH(
        A = ::Star_MC_DoubleList.size();
        ::Star_MC_DoubleList.push_back
        ( new DistMatrix<double,Star,MC>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
Star_MC_SComplex
Register_Star_MC_SComplex
( int height, int width, int rowAlignment, SComplex* buffer, int ldim, Grid g )
{
    Star_MC_SComplex A;
    CATCH(
        A = ::Star_MC_SComplexList.size();
        ::Star_MC_SComplexList.push_back
        ( new DistMatrix<scomplex,Star,MC>
          ( height, width, rowAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
Star_MC_DComplex
Register_Star_MC_DComplex
( int height, int width, int rowAlignment, DComplex* buffer, int ldim, Grid g )
{
    Star_MC_DComplex A;
    CATCH(
        A = ::Star_MC_DComplexList.size();
        ::Star_MC_DComplexList.push_back
        ( new DistMatrix<dcomplex,Star,MC>
          ( height, width, rowAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the row alignments
int
RowAlignment_Star_MC_Single( Star_MC_Single A )
{ return ::Star_MC_SingleList[A]->RowAlignment(); }
int
RowAlignment_Star_MC_Double( Star_MC_Double A )
{ return ::Star_MC_DoubleList[A]->RowAlignment(); }
#ifndef WITHOUT_COMPLEX
int
RowAlignment_Star_MC_SComplex( Star_MC_SComplex A )
{ return ::Star_MC_SComplexList[A]->RowAlignment(); }
int
RowAlignment_Star_MC_DComplex( Star_MC_DComplex A )
{ return ::Star_MC_DComplexList[A]->RowAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_Star_MC_Single( char* msg, Star_MC_Single A )
{
    CATCH(
        std::string s( msg );
        ::Star_MC_SingleList[A]->Print( s );
    );
}
void 
Print_Star_MC_Double( char* msg, Star_MC_Double A )
{
    CATCH(
        std::string s( msg );
        ::Star_MC_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_Star_MC_SComplex( char* msg, Star_MC_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::Star_MC_SComplexList[A]->Print( s );
    );
}
void 
Print_Star_MC_DComplex( char* msg, Star_MC_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::Star_MC_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [* ,MD] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
Star_MD_Single
CreateEmpty_Star_MD_Single( Grid g )
{
    Star_MD_Single A;
    CATCH(
        A = ::Star_MD_SingleList.size();
        ::Star_MD_SingleList.push_back
        ( new DistMatrix<float,Star,MD>( *::gridList[g] ) );
    );
    return A;
}
Star_MD_Double
CreateEmpty_Star_MD_Double( Grid g )
{
    Star_MD_Double A;
    CATCH(
        A = ::Star_MD_DoubleList.size();
        ::Star_MD_DoubleList.push_back
        ( new DistMatrix<double,Star,MD>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
Star_MD_SComplex
CreateEmpty_Star_MD_SComplex( Grid g )
{
    Star_MD_SComplex A;
    CATCH(
        A = ::Star_MD_SComplexList.size();
        ::Star_MD_SComplexList.push_back
        ( new DistMatrix<scomplex,Star,MD>( *::gridList[g] ) );
    );
    return A;
}
Star_MD_DComplex
CreateEmpty_Star_MD_DComplex( Grid g )
{
    Star_MD_DComplex A;
    CATCH(
        A = ::Star_MD_DComplexList.size();
        ::Star_MD_DComplexList.push_back
        ( new DistMatrix<dcomplex,Star,MD>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
Star_MD_Single
Register_Star_MD_Single
( int height, int width, int rowAlignment, float* buffer, int ldim, Grid g )
{
    Star_MD_Single A;
    CATCH(
        A = ::Star_MD_SingleList.size();
        ::Star_MD_SingleList.push_back
        ( new DistMatrix<float,Star,MD>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
Star_MD_Double
Register_Star_MD_Double
( int height, int width, int rowAlignment, double* buffer, int ldim, Grid g )
{
    Star_MD_Double A;
    CATCH(
        A = ::Star_MD_DoubleList.size();
        ::Star_MD_DoubleList.push_back
        ( new DistMatrix<double,Star,MD>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
Star_MD_SComplex
Register_Star_MD_SComplex
( int height, int width, int rowAlignment, SComplex* buffer, int ldim, Grid g )
{
    Star_MD_SComplex A;
    CATCH(
        A = ::Star_MD_SComplexList.size();
        ::Star_MD_SComplexList.push_back
        ( new DistMatrix<scomplex,Star,MD>
          ( height, width, rowAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
Star_MD_DComplex
Register_Star_MD_DComplex
( int height, int width, int rowAlignment, DComplex* buffer, int ldim, Grid g )
{
    Star_MD_DComplex A;
    CATCH(
        A = ::Star_MD_DComplexList.size();
        ::Star_MD_DComplexList.push_back
        ( new DistMatrix<dcomplex,Star,MD>
          ( height, width, rowAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the row alignments
int
RowAlignment_MD_Star_Single( MD_Star_Single A )
{ return ::MD_Star_SingleList[A]->RowAlignment(); }
int
RowAlignment_MD_Star_Double( MD_Star_Double A )
{ return ::MD_Star_DoubleList[A]->RowAlignment(); }
#ifndef WITHOUT_COMPLEX
int
RowAlignment_MD_Star_SComplex( MD_Star_SComplex A )
{ return ::MD_Star_SComplexList[A]->RowAlignment(); }
int
RowAlignment_MD_Star_DComplex( MD_Star_DComplex A )
{ return ::MD_Star_DComplexList[A]->RowAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_Star_MD_Single( char* msg, Star_MD_Single A )
{
    CATCH(
        std::string s( msg );
        ::Star_MD_SingleList[A]->Print( s );
    );
}
void 
Print_Star_MD_Double( char* msg, Star_MD_Double A )
{
    CATCH(
        std::string s( msg );
        ::Star_MD_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_Star_MD_SComplex( char* msg, Star_MD_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::Star_MD_SComplexList[A]->Print( s );
    );
}
void 
Print_Star_MD_DComplex( char* msg, Star_MD_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::Star_MD_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [* ,MR] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
Star_MR_Single
CreateEmpty_Star_MR_Single( Grid g )
{
    Star_MR_Single A;
    CATCH(
        A = ::Star_MR_SingleList.size();
        ::Star_MR_SingleList.push_back
        ( new DistMatrix<float,Star,MR>( *::gridList[g] ) );
    );
    return A;
}
Star_MR_Double
CreateEmpty_Star_MR_Double( Grid g )
{
    Star_MR_Double A;
    CATCH(
        A = ::Star_MR_DoubleList.size();
        ::Star_MR_DoubleList.push_back
        ( new DistMatrix<double,Star,MR>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
Star_MR_SComplex
CreateEmpty_Star_MR_SComplex( Grid g )
{
    Star_MR_SComplex A;
    CATCH(
        A = ::Star_MR_SComplexList.size();
        ::Star_MR_SComplexList.push_back
        ( new DistMatrix<scomplex,Star,MR>( *::gridList[g] ) );
    );
    return A;
}
Star_MR_DComplex
CreateEmpty_Star_MR_DComplex( Grid g )
{
    Star_MR_DComplex A;
    CATCH(
        A = ::Star_MR_DComplexList.size();
        ::Star_MR_DComplexList.push_back
        ( new DistMatrix<dcomplex,Star,MR>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
Star_MR_Single
Register_Star_MR_Single
( int height, int width, int rowAlignment, float* buffer, int ldim, Grid g )
{
    Star_MR_Single A;
    CATCH(
        A = ::Star_MR_SingleList.size();
        ::Star_MR_SingleList.push_back
        ( new DistMatrix<float,Star,MR>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
Star_MR_Double
Register_Star_MR_Double
( int height, int width, int rowAlignment, double* buffer, int ldim, Grid g )
{
    Star_MR_Double A;
    CATCH(
        A = ::Star_MR_DoubleList.size();
        ::Star_MR_DoubleList.push_back
        ( new DistMatrix<double,Star,MR>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
Star_MR_SComplex
Register_Star_MR_SComplex
( int height, int width, int rowAlignment, SComplex* buffer, int ldim, Grid g )
{
    Star_MR_SComplex A;
    CATCH(
        A = ::Star_MR_SComplexList.size();
        ::Star_MR_SComplexList.push_back
        ( new DistMatrix<scomplex,Star,MR>
          ( height, width, rowAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
Star_MR_DComplex
Register_Star_MR_DComplex
( int height, int width, int rowAlignment, DComplex* buffer, int ldim, Grid g )
{
    Star_MR_DComplex A;
    CATCH(
        A = ::Star_MR_DComplexList.size();
        ::Star_MR_DComplexList.push_back
        ( new DistMatrix<dcomplex,Star,MR>
          ( height, width, rowAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the row alignments
int
RowAlignment_Star_MR_Single( Star_MR_Single A )
{ return ::Star_MR_SingleList[A]->RowAlignment(); }
int
RowAlignment_Star_MR_Double( Star_MR_Double A )
{ return ::Star_MR_DoubleList[A]->RowAlignment(); }
#ifndef WITHOUT_COMPLEX
int
RowAlignment_Star_MR_SComplex( Star_MR_SComplex A )
{ return ::Star_MR_SComplexList[A]->RowAlignment(); }
int
RowAlignment_Star_MR_DComplex( Star_MR_DComplex A )
{ return ::Star_MR_DComplexList[A]->RowAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_Star_MR_Single( char* msg, Star_MR_Single A )
{
    CATCH(
        std::string s( msg );
        ::Star_MR_SingleList[A]->Print( s );
    );
}
void 
Print_Star_MR_Double( char* msg, Star_MR_Double A )
{
    CATCH(
        std::string s( msg );
        ::Star_MR_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_Star_MR_SComplex( char* msg, Star_MR_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::Star_MR_SComplexList[A]->Print( s );
    );
}
void 
Print_Star_MR_DComplex( char* msg, Star_MR_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::Star_MR_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [* ,* ] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
Star_Star_Single
CreateEmpty_Star_Star_Single( Grid g )
{
    Star_Star_Single A;
    CATCH(
        A = ::Star_Star_SingleList.size();
        ::Star_Star_SingleList.push_back
        ( new DistMatrix<float,Star,Star>( *::gridList[g] ) );
    );
    return A;
}
Star_Star_Double
CreateEmpty_Star_Star_Double( Grid g )
{
    Star_Star_Double A;
    CATCH(
        A = ::Star_Star_DoubleList.size();
        ::Star_Star_DoubleList.push_back
        ( new DistMatrix<double,Star,Star>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
Star_Star_SComplex
CreateEmpty_Star_Star_SComplex( Grid g )
{
    Star_Star_SComplex A;
    CATCH(
        A = ::Star_Star_SComplexList.size();
        ::Star_Star_SComplexList.push_back
        ( new DistMatrix<scomplex,Star,Star>( *::gridList[g] ) );
    );
    return A;
}
Star_Star_DComplex
CreateEmpty_Star_Star_DComplex( Grid g )
{
    Star_Star_DComplex A;
    CATCH(
        A = ::Star_Star_DComplexList.size();
        ::Star_Star_DComplexList.push_back
        ( new DistMatrix<dcomplex,Star,Star>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
Star_Star_Single
Register_Star_Star_Single
( int height, int width, float* buffer, int ldim, Grid g )
{
    Star_Star_Single A;
    CATCH(
        A = ::Star_Star_SingleList.size();
        ::Star_Star_SingleList.push_back
        ( new DistMatrix<float,Star,Star>
          ( height, width, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
Star_Star_Double
Register_Star_Star_Double
( int height, int width, double* buffer, int ldim, Grid g )
{
    Star_Star_Double A;
    CATCH(
        A = ::Star_Star_DoubleList.size();
        ::Star_Star_DoubleList.push_back
        ( new DistMatrix<double,Star,Star>
          ( height, width, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
Star_Star_SComplex
Register_Star_Star_SComplex
( int height, int width, SComplex* buffer, int ldim, Grid g )
{
    Star_Star_SComplex A;
    CATCH(
        A = ::Star_Star_SComplexList.size();
        ::Star_Star_SComplexList.push_back
        ( new DistMatrix<scomplex,Star,Star>
          ( height, width, (scomplex*)buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
Star_Star_DComplex
Register_Star_Star_DComplex
( int height, int width, DComplex* buffer, int ldim, Grid g )
{
    Star_Star_DComplex A;
    CATCH(
        A = ::Star_Star_DComplexList.size();
        ::Star_Star_DComplexList.push_back
        ( new DistMatrix<dcomplex,Star,Star>
          ( height, width, (dcomplex*)buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Print
void 
Print_Star_Star_Single( char* msg, Star_Star_Single A )
{
    CATCH(
        std::string s( msg );
        ::Star_Star_SingleList[A]->Print( s );
    );
}
void 
Print_Star_Star_Double( char* msg, Star_Star_Double A )
{
    CATCH(
        std::string s( msg );
        ::Star_Star_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_Star_Star_SComplex( char* msg, Star_Star_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::Star_Star_SComplexList[A]->Print( s );
    );
}
void 
Print_Star_Star_DComplex( char* msg, Star_Star_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::Star_Star_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [* ,VC] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
Star_VC_Single
CreateEmpty_Star_VC_Single( Grid g )
{
    Star_VC_Single A;
    CATCH(
        A = ::Star_VC_SingleList.size();
        ::Star_VC_SingleList.push_back
        ( new DistMatrix<float,Star,VC>( *::gridList[g] ) );
    );
    return A;
}
Star_VC_Double
CreateEmpty_Star_VC_Double( Grid g )
{
    Star_VC_Double A;
    CATCH(
        A = ::Star_VC_DoubleList.size();
        ::Star_VC_DoubleList.push_back
        ( new DistMatrix<double,Star,VC>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
Star_VC_SComplex
CreateEmpty_Star_VC_SComplex( Grid g )
{
    Star_VC_SComplex A;
    CATCH(
        A = ::Star_VC_SComplexList.size();
        ::Star_VC_SComplexList.push_back
        ( new DistMatrix<scomplex,Star,VC>( *::gridList[g] ) );
    );
    return A;
}
Star_VC_DComplex
CreateEmpty_Star_VC_DComplex( Grid g )
{
    Star_VC_DComplex A;
    CATCH(
        A = ::Star_VC_DComplexList.size();
        ::Star_VC_DComplexList.push_back
        ( new DistMatrix<dcomplex,Star,VC>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
Star_VC_Single
Register_Star_VC_Single
( int height, int width, int rowAlignment, float* buffer, int ldim, Grid g )
{
    Star_VC_Single A;
    CATCH(
        A = ::Star_VC_SingleList.size();
        ::Star_VC_SingleList.push_back
        ( new DistMatrix<float,Star,VC>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
Star_VC_Double
Register_Star_VC_Double
( int height, int width, int rowAlignment, double* buffer, int ldim, Grid g )
{
    Star_VC_Double A;
    CATCH(
        A = ::Star_VC_DoubleList.size();
        ::Star_VC_DoubleList.push_back
        ( new DistMatrix<double,Star,VC>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
Star_VC_SComplex
Register_Star_VC_SComplex
( int height, int width, int rowAlignment, SComplex* buffer, int ldim, Grid g )
{
    Star_VC_SComplex A;
    CATCH(
        A = ::Star_VC_SComplexList.size();
        ::Star_VC_SComplexList.push_back
        ( new DistMatrix<scomplex,Star,VC>
          ( height, width, rowAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
Star_VC_DComplex
Register_Star_VC_DComplex
( int height, int width, int rowAlignment, DComplex* buffer, int ldim, Grid g )
{
    Star_VC_DComplex A;
    CATCH(
        A = ::Star_VC_DComplexList.size();
        ::Star_VC_DComplexList.push_back
        ( new DistMatrix<dcomplex,Star,VC>
          ( height, width, rowAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the row alignments
int
RowAlignment_Star_VC_Single( Star_VC_Single A )
{ return ::Star_VC_SingleList[A]->RowAlignment(); }
int
RowAlignment_Star_VC_Double( Star_VC_Double A )
{ return ::Star_VC_DoubleList[A]->RowAlignment(); }
#ifndef WITHOUT_COMPLEX
int
RowAlignment_Star_VC_SComplex( Star_VC_SComplex A )
{ return ::Star_VC_SComplexList[A]->RowAlignment(); }
int
RowAlignment_Star_VC_DComplex( Star_VC_DComplex A )
{ return ::Star_VC_DComplexList[A]->RowAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_Star_VC_Single( char* msg, Star_VC_Single A )
{
    CATCH(
        std::string s( msg );
        ::Star_VC_SingleList[A]->Print( s );
    );
}
void 
Print_Star_VC_Double( char* msg, Star_VC_Double A )
{
    CATCH(
        std::string s( msg );
        ::Star_VC_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_Star_VC_SComplex( char* msg, Star_VC_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::Star_VC_SComplexList[A]->Print( s );
    );
}
void 
Print_Star_VC_DComplex( char* msg, Star_VC_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::Star_VC_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [* ,VR] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
Star_VR_Single
CreateEmpty_Star_VR_Single( Grid g )
{
    Star_VR_Single A;
    CATCH(
        A = ::Star_VR_SingleList.size();
        ::Star_VR_SingleList.push_back
        ( new DistMatrix<float,Star,VR>( *::gridList[g] ) );
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
        ( new DistMatrix<double,Star,VR>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
Star_VR_SComplex
CreateEmpty_Star_VR_SComplex( Grid g )
{
    Star_VR_SComplex A;
    CATCH(
        A = ::Star_VR_SComplexList.size();
        ::Star_VR_SComplexList.push_back
        ( new DistMatrix<scomplex,Star,VR>( *::gridList[g] ) );
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
        ( new DistMatrix<dcomplex,Star,VR>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
Star_VR_Single
Register_Star_VR_Single
( int height, int width, int rowAlignment, float* buffer, int ldim, Grid g )
{
    Star_VR_Single A;
    CATCH(
        A = ::Star_VR_SingleList.size();
        ::Star_VR_SingleList.push_back
        ( new DistMatrix<float,Star,VR>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
Star_VR_Double
Register_Star_VR_Double
( int height, int width, int rowAlignment, double* buffer, int ldim, Grid g )
{
    Star_VR_Double A;
    CATCH(
        A = ::Star_VR_DoubleList.size();
        ::Star_VR_DoubleList.push_back
        ( new DistMatrix<double,Star,VR>
          ( height, width, rowAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
Star_VR_SComplex
Register_Star_VR_SComplex
( int height, int width, int rowAlignment, SComplex* buffer, int ldim, Grid g )
{
    Star_VR_SComplex A;
    CATCH(
        A = ::Star_VR_SComplexList.size();
        ::Star_VR_SComplexList.push_back
        ( new DistMatrix<scomplex,Star,VR>
          ( height, width, rowAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
Star_VR_DComplex
Register_Star_VR_DComplex
( int height, int width, int rowAlignment, DComplex* buffer, int ldim, Grid g )
{
    Star_VR_DComplex A;
    CATCH(
        A = ::Star_VR_DComplexList.size();
        ::Star_VR_DComplexList.push_back
        ( new DistMatrix<dcomplex,Star,VR>
          ( height, width, rowAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the row alignments
int
RowAlignment_Star_VR_Single( Star_VR_Single A )
{ return ::Star_VR_SingleList[A]->RowAlignment(); }
int
RowAlignment_Star_VR_Double( Star_VR_Double A )
{ return ::Star_VR_DoubleList[A]->RowAlignment(); }
#ifndef WITHOUT_COMPLEX
int
RowAlignment_Star_VR_SComplex( Star_VR_SComplex A )
{ return ::Star_VR_SComplexList[A]->RowAlignment(); }
int
RowAlignment_Star_VR_DComplex( Star_VR_DComplex A )
{ return ::Star_VR_DComplexList[A]->RowAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
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
//----------------------------------------------------------------------------//
// [VC,* ] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
VC_Star_Single
CreateEmpty_VC_Star_Single( Grid g )
{
    VC_Star_Single A;
    CATCH(
        A = ::VC_Star_SingleList.size();
        ::VC_Star_SingleList.push_back
        ( new DistMatrix<float,VC,Star>( *::gridList[g] ) );
    );
    return A;
}
VC_Star_Double
CreateEmpty_VC_Star_Double( Grid g )
{
    VC_Star_Double A;
    CATCH(
        A = ::VC_Star_DoubleList.size();
        ::VC_Star_DoubleList.push_back
        ( new DistMatrix<double,VC,Star>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
VC_Star_SComplex
CreateEmpty_VC_Star_SComplex( Grid g )
{
    VC_Star_SComplex A;
    CATCH(
        A = ::VC_Star_SComplexList.size();
        ::VC_Star_SComplexList.push_back
        ( new DistMatrix<scomplex,VC,Star>( *::gridList[g] ) );
    );
    return A;
}
VC_Star_DComplex
CreateEmpty_VC_Star_DComplex( Grid g )
{
    VC_Star_DComplex A;
    CATCH(
        A = ::VC_Star_DComplexList.size();
        ::VC_Star_DComplexList.push_back
        ( new DistMatrix<dcomplex,VC,Star>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
VC_Star_Single
Register_VC_Star_Single
( int height, int width, int colAlignment, float* buffer, int ldim, Grid g )
{
    VC_Star_Single A;
    CATCH(
        A = ::VC_Star_SingleList.size();
        ::VC_Star_SingleList.push_back
        ( new DistMatrix<float,VC,Star>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
VC_Star_Double
Register_VC_Star_Double
( int height, int width, int colAlignment, double* buffer, int ldim, Grid g )
{
    VC_Star_Double A;
    CATCH(
        A = ::VC_Star_DoubleList.size();
        ::VC_Star_DoubleList.push_back
        ( new DistMatrix<double,VC,Star>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
VC_Star_SComplex
Register_VC_Star_SComplex
( int height, int width, int colAlignment, SComplex* buffer, int ldim, Grid g )
{
    VC_Star_SComplex A;
    CATCH(
        A = ::VC_Star_SComplexList.size();
        ::VC_Star_SComplexList.push_back
        ( new DistMatrix<scomplex,VC,Star>
          ( height, width, colAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
VC_Star_DComplex
Register_VC_Star_DComplex
( int height, int width, int colAlignment, DComplex* buffer, int ldim, Grid g )
{
    VC_Star_DComplex A;
    CATCH(
        A = ::VC_Star_DComplexList.size();
        ::VC_Star_DComplexList.push_back
        ( new DistMatrix<dcomplex,VC,Star>
          ( height, width, colAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the col alignments
int
ColAlignment_VC_Star_Single( VC_Star_Single A )
{ return ::VC_Star_SingleList[A]->ColAlignment(); }
int
ColAlignment_VC_Star_Double( VC_Star_Double A )
{ return ::VC_Star_DoubleList[A]->ColAlignment(); }
#ifndef WITHOUT_COMPLEX
int
ColAlignment_VC_Star_SComplex( VC_Star_SComplex A )
{ return ::VC_Star_SComplexList[A]->ColAlignment(); }
int
ColAlignment_VC_Star_DComplex( VC_Star_DComplex A )
{ return ::VC_Star_DComplexList[A]->ColAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_VC_Star_Single( char* msg, VC_Star_Single A )
{
    CATCH(
        std::string s( msg );
        ::VC_Star_SingleList[A]->Print( s );
    );
}
void 
Print_VC_Star_Double( char* msg, VC_Star_Double A )
{
    CATCH(
        std::string s( msg );
        ::VC_Star_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_VC_Star_SComplex( char* msg, VC_Star_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::VC_Star_SComplexList[A]->Print( s );
    );
}
void 
Print_VC_Star_DComplex( char* msg, VC_Star_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::VC_Star_DComplexList[A]->Print( s );
    );
}
#endif // WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// [VR,* ] manipulation routines                                              //
//----------------------------------------------------------------------------//
// Create empty
VR_Star_Single
CreateEmpty_VR_Star_Single( Grid g )
{
    VR_Star_Single A;
    CATCH(
        A = ::VR_Star_SingleList.size();
        ::VR_Star_SingleList.push_back
        ( new DistMatrix<float,VR,Star>( *::gridList[g] ) );
    );
    return A;
}
VR_Star_Double
CreateEmpty_VR_Star_Double( Grid g )
{
    VR_Star_Double A;
    CATCH(
        A = ::VR_Star_DoubleList.size();
        ::VR_Star_DoubleList.push_back
        ( new DistMatrix<double,VR,Star>( *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
VR_Star_SComplex
CreateEmpty_VR_Star_SComplex( Grid g )
{
    VR_Star_SComplex A;
    CATCH(
        A = ::VR_Star_SComplexList.size();
        ::VR_Star_SComplexList.push_back
        ( new DistMatrix<scomplex,VR,Star>( *::gridList[g] ) );
    );
    return A;
}
VR_Star_DComplex
CreateEmpty_VR_Star_DComplex( Grid g )
{
    VR_Star_DComplex A;
    CATCH(
        A = ::VR_Star_DComplexList.size();
        ::VR_Star_DComplexList.push_back
        ( new DistMatrix<dcomplex,VR,Star>( *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Create from existing buffer
VR_Star_Single
Register_VR_Star_Single
( int height, int width, int colAlignment, float* buffer, int ldim, Grid g )
{
    VR_Star_Single A;
    CATCH(
        A = ::VR_Star_SingleList.size();
        ::VR_Star_SingleList.push_back
        ( new DistMatrix<float,VR,Star>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
VR_Star_Double
Register_VR_Star_Double
( int height, int width, int colAlignment, double* buffer, int ldim, Grid g )
{
    VR_Star_Double A;
    CATCH(
        A = ::VR_Star_DoubleList.size();
        ::VR_Star_DoubleList.push_back
        ( new DistMatrix<double,VR,Star>
          ( height, width, colAlignment, buffer, ldim, *::gridList[g] ) );
    );
    return A;
}
#ifndef WITHOUT_COMPLEX
VR_Star_SComplex
Register_VR_Star_SComplex
( int height, int width, int colAlignment, SComplex* buffer, int ldim, Grid g )
{
    VR_Star_SComplex A;
    CATCH(
        A = ::VR_Star_SComplexList.size();
        ::VR_Star_SComplexList.push_back
        ( new DistMatrix<scomplex,VR,Star>
          ( height, width, colAlignment, (scomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
VR_Star_DComplex
Register_VR_Star_DComplex
( int height, int width, int colAlignment, DComplex* buffer, int ldim, Grid g )
{
    VR_Star_DComplex A;
    CATCH(
        A = ::VR_Star_DComplexList.size();
        ::VR_Star_DComplexList.push_back
        ( new DistMatrix<dcomplex,VR,Star>
          ( height, width, colAlignment, (dcomplex*)buffer, ldim, 
            *::gridList[g] ) );
    );
    return A;
}
#endif // WITHOUT_COMPLEX
// Get the col alignments
int
ColAlignment_VR_Star_Single( VR_Star_Single A )
{ return ::VR_Star_SingleList[A]->ColAlignment(); }
int
ColAlignment_VR_Star_Double( VR_Star_Double A )
{ return ::VR_Star_DoubleList[A]->ColAlignment(); }
#ifndef WITHOUT_COMPLEX
int
ColAlignment_VR_Star_SComplex( VR_Star_SComplex A )
{ return ::VR_Star_SComplexList[A]->ColAlignment(); }
int
ColAlignment_VR_Star_DComplex( VR_Star_DComplex A )
{ return ::VR_Star_DComplexList[A]->ColAlignment(); }
#endif // WITHOUT_COMPLEX
// Print
void 
Print_VR_Star_Single( char* msg, VR_Star_Single A )
{
    CATCH(
        std::string s( msg );
        ::VR_Star_SingleList[A]->Print( s );
    );
}
void 
Print_VR_Star_Double( char* msg, VR_Star_Double A )
{
    CATCH(
        std::string s( msg );
        ::VR_Star_DoubleList[A]->Print( s );
    );
}
#ifndef WITHOUT_COMPLEX
void 
Print_VR_Star_SComplex( char* msg, VR_Star_SComplex A )
{
    CATCH(
        std::string s( msg );
        ::VR_Star_SComplexList[A]->Print( s );
    );
}
void 
Print_VR_Star_DComplex( char* msg, VR_Star_DComplex A )
{
    CATCH(
        std::string s( msg );
        ::VR_Star_DComplexList[A]->Print( s );
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
  MC_MR_Double A, Star_VR_Double w, MC_MR_Double Z )
{
    elemental::Shape shape;
    CATCH(shape = elemental::CharToShape( uplo ););
    if( job == 'a' || job == 'A' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DoubleList[A], *::Star_VR_DoubleList[w],
              *::MC_MR_DoubleList[Z], tryForHighAccuracy );
        );
    }
    else if( job == 'i' || job == 'I' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DoubleList[A], *::Star_VR_DoubleList[w],
              *::MC_MR_DoubleList[Z], a, b, tryForHighAccuracy );
        );
    }
    else if( job == 'v' || job == 'V' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DoubleList[A], *::Star_VR_DoubleList[w],
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
  MC_MR_Double A, Star_VR_Double w )
{
    elemental::Shape shape;
    CATCH(shape = elemental::CharToShape( uplo ););
    if( job == 'a' || job == 'A' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DoubleList[A], *::Star_VR_DoubleList[w],
              tryForHighAccuracy );
        );
    }
    else if( job == 'i' || job == 'I' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DoubleList[A], *::Star_VR_DoubleList[w],
              a, b, tryForHighAccuracy );
        );
    }
    else if( job == 'v' || job == 'V' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DoubleList[A], *::Star_VR_DoubleList[w],
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
  MC_MR_DComplex A, Star_VR_Double w, MC_MR_DComplex Z )
{
    elemental::Shape shape;
    CATCH(shape = elemental::CharToShape( uplo ););
    if( job == 'a' || job == 'A' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DComplexList[A], *::Star_VR_DoubleList[w],
              *::MC_MR_DComplexList[Z], tryForHighAccuracy );
        );
    }
    else if( job == 'i' || job == 'I' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DComplexList[A], *::Star_VR_DoubleList[w],
              *::MC_MR_DComplexList[Z], a, b, tryForHighAccuracy );
        );
    }
    else if( job == 'v' || job == 'V' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DComplexList[A], *::Star_VR_DoubleList[w],
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
  MC_MR_DComplex A, Star_VR_Double w )
{
    elemental::Shape shape;
    CATCH(shape = elemental::CharToShape( uplo ););
    if( job == 'a' || job == 'A' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DComplexList[A], *::Star_VR_DoubleList[w],
              tryForHighAccuracy );
        );
    }
    else if( job == 'i' || job == 'I' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DComplexList[A], *::Star_VR_DoubleList[w],
              a, b, tryForHighAccuracy );
        );
    }
    else if( job == 'v' || job == 'V' )
    {
        CATCH(
            elemental::advanced::HermitianEig
            ( shape, *::MC_MR_DComplexList[A], *::Star_VR_DoubleList[w],
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
  MC_MR_Double A, MC_MR_Double B, Star_VR_Double w, MC_MR_Double Z )
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
              *::Star_VR_DoubleList[w], *::MC_MR_DoubleList[Z], 
              tryForHighAccuracy );
        );
    }
    else if( job == 'i' || job == 'I' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DoubleList[A], *::MC_MR_DoubleList[B], 
              *::Star_VR_DoubleList[w], *::MC_MR_DoubleList[Z], 
              a, b, tryForHighAccuracy );
        );
    }
    else if( job == 'v' || job == 'V' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DoubleList[A], *::MC_MR_DoubleList[B], 
              *::Star_VR_DoubleList[w], *::MC_MR_DoubleList[Z], 
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
  MC_MR_Double A, MC_MR_Double B, Star_VR_Double w )
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
              *::Star_VR_DoubleList[w], tryForHighAccuracy );
        );
    }
    else if( job == 'i' || job == 'I' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DoubleList[A], *::MC_MR_DoubleList[B], 
              *::Star_VR_DoubleList[w], a, b, tryForHighAccuracy );
        );
    }
    else if( job == 'v' || job == 'V' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DoubleList[A], *::MC_MR_DoubleList[B], 
              *::Star_VR_DoubleList[w], u, v, tryForHighAccuracy );
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
  MC_MR_DComplex A, MC_MR_DComplex B, Star_VR_Double w, MC_MR_DComplex Z )
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
              *::Star_VR_DoubleList[w], *::MC_MR_DComplexList[Z], 
              tryForHighAccuracy );
        );
    }
    else if( job == 'i' || job == 'I' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DComplexList[A], *::MC_MR_DComplexList[B], 
              *::Star_VR_DoubleList[w], *::MC_MR_DComplexList[Z], 
              a, b, tryForHighAccuracy );
        );
    }
    else if( job == 'v' || job == 'V' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DComplexList[A], *::MC_MR_DComplexList[B], 
              *::Star_VR_DoubleList[w], *::MC_MR_DComplexList[Z], 
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
  MC_MR_DComplex A, MC_MR_DComplex B, Star_VR_Double w )
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
              *::Star_VR_DoubleList[w], tryForHighAccuracy );
        );
    }
    else if( job == 'i' || job == 'I' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DComplexList[A], *::MC_MR_DComplexList[B], 
              *::Star_VR_DoubleList[w], a, b, tryForHighAccuracy );
        );
    }
    else if( job == 'v' || job == 'V' )
    {
        CATCH(
            elemental::advanced::HermitianGenDefiniteEig
            ( genEigType, shape, 
              *::MC_MR_DComplexList[A], *::MC_MR_DComplexList[B], 
              *::Star_VR_DoubleList[w], u, v, tryForHighAccuracy );
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

