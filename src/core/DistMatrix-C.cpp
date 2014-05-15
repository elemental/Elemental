/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include EL_COPY_INC 

#include "El-C.h"
using namespace El;

#define RC(TYPE,INPUT) reinterpret_cast<TYPE>(INPUT)

#define RCG(gridHandle) RC(Grid*,gridHandle)
#define RCG_const(gridHandle) RC(const Grid*,gridHandle)

#define RCDDM_s(AHandle) RC(DynamicDistMatrix<float          >*,AHandle)
#define RCDDM_d(AHandle) RC(DynamicDistMatrix<double         >*,AHandle)
#define RCDDM_c(AHandle) RC(DynamicDistMatrix<Complex<float >>*,AHandle)
#define RCDDM_z(AHandle) RC(DynamicDistMatrix<Complex<double>>*,AHandle)

#define RCDDM_s_const(AHandle) \
  RC(const DynamicDistMatrix<float          >*,AHandle)
#define RCDDM_d_const(AHandle) \
  RC(const DynamicDistMatrix<double         >*,AHandle)
#define RCDDM_c_const(AHandle) \
  RC(const DynamicDistMatrix<Complex<float >>*,AHandle)
#define RCDDM_z_const(AHandle) \
  RC(const DynamicDistMatrix<Complex<double>>*,AHandle)

#define RCB_c(buffer) RC(Complex<float>*,buffer)
#define RCB_z(buffer) RC(Complex<double>*,buffer)

#define RCB_c_const(buffer) RC(const Complex<float >*,buffer)
#define RCB_z_const(buffer) RC(const Complex<double>*,buffer)

#define CATCH catch( std::exception& e ) { ReportException(e); }

extern "C" {

// Simple contructor for [MC,MR] option for DynamicDistMatrix
// ----------------------------------------------------------
ElDistMatrix_s* ElDistMatrixCreate_s( const ElGrid* gridHandle )
{
    ElDistMatrix_s* AHandle = 0;
    try 
    {
        auto A = new DynamicDistMatrix<float>;
        A->ADM = new DistMatrix<float>(*RCG_const(gridHandle));
        AHandle = RC(ElDistMatrix_s*,A);
    }
    CATCH
    return AHandle;
}

ElDistMatrix_d* ElDistMatrixCreate_d( const ElGrid* gridHandle )
{
    ElDistMatrix_d* AHandle = 0;
    try 
    {
        auto A = new DynamicDistMatrix<double>;
        A->ADM = new DistMatrix<double>(*RCG_const(gridHandle));
        AHandle = RC(ElDistMatrix_d*,A);
    }
    CATCH
    return AHandle;
}

ElDistMatrix_c* ElDistMatrixCreate_c( const ElGrid* gridHandle )
{
    ElDistMatrix_c* AHandle = 0;
    try 
    {
        auto A = new DynamicDistMatrix<Complex<float>>;
        A->ADM = new DistMatrix<Complex<float>>(*RCG_const(gridHandle));
        AHandle = RC(ElDistMatrix_c*,A);
    }
    CATCH
    return AHandle;
}

ElDistMatrix_z* ElDistMatrixCreate_z( const ElGrid* gridHandle )
{
    ElDistMatrix_z* AHandle = 0;
    try 
    {
        auto A = new DynamicDistMatrix<Complex<double>>;
        A->ADM = new DistMatrix<Complex<double>>(*RCG_const(gridHandle));
        AHandle = RC(ElDistMatrix_z*,A);
    }
    CATCH
    return AHandle;
}

// Simple contructor for DynamicDistMatrix
// ---------------------------------------
ElDistMatrix_s* ElDistMatrixCreateSpecific_s
( ElDist UC, ElDist VC, const ElGrid* gridHandle )
{
    ElDistMatrix_s* AHandle = 0;
    try 
    {
        auto A = new DynamicDistMatrix<float>;
        Dist U = static_cast<Dist>(UC);
        Dist V = static_cast<Dist>(VC);
        A->SetDistribution( U, V, *RCG_const(gridHandle) );
        AHandle = RC(ElDistMatrix_s*,A);
    }
    CATCH
    return AHandle;
}

ElDistMatrix_d* ElDistMatrixCreateSpecific_d
( ElDist UC, ElDist VC, const ElGrid* gridHandle )
{
    ElDistMatrix_d* AHandle = 0;
    try 
    {
        auto A = new DynamicDistMatrix<double>;
        Dist U = static_cast<Dist>(UC);
        Dist V = static_cast<Dist>(VC);
        A->SetDistribution( U, V, *RCG_const(gridHandle) );
        AHandle = RC(ElDistMatrix_d*,A);
    }
    CATCH
    return AHandle;
}

ElDistMatrix_c* ElDistMatrixCreateSpecific_c
( ElDist UC, ElDist VC, const ElGrid* gridHandle )
{
    ElDistMatrix_c* AHandle = 0;
    try 
    {
        auto A = new DynamicDistMatrix<Complex<float>>;
        Dist U = static_cast<Dist>(UC);
        Dist V = static_cast<Dist>(VC);
        A->SetDistribution( U, V, *RCG_const(gridHandle) );
        AHandle = RC(ElDistMatrix_c*,A);
    }
    CATCH
    return AHandle;
}

ElDistMatrix_z* ElDistMatrixCreateSpecific_z
( ElDist UC, ElDist VC, const ElGrid* gridHandle )
{
    ElDistMatrix_z* AHandle = 0;
    try 
    {
        auto A = new DynamicDistMatrix<Complex<double>>;
        Dist U = static_cast<Dist>(UC);
        Dist V = static_cast<Dist>(VC);
        A->SetDistribution( U, V, *RCG_const(gridHandle) );
        AHandle = RC(ElDistMatrix_z*,A);
    }
    CATCH
    return AHandle;
}

// DynamicDistMatrix::~DynamicDistMatrix()
// ---------------------------------------
void ElDistMatrixDestroy_s( const ElDistMatrix_s* AHandle )
{ delete RCDDM_s_const(AHandle); }

void ElDistMatrixDestroy_d( const ElDistMatrix_d* AHandle )
{ delete RCDDM_d_const(AHandle); }

void ElDistMatrixDestroy_c( const ElDistMatrix_c* AHandle )
{ delete RCDDM_c_const(AHandle); }

void ElDistMatrixDestroy_z( const ElDistMatrix_z* AHandle )
{ delete RCDDM_z_const(AHandle); }

// AbstractDistMatrix<T>::Empty()
// ------------------------------
void ElDistMatrixEmpty_s( ElDistMatrix_s* AHandle )
{
    try { RCDDM_s(AHandle)->ADM->Empty(); }
    CATCH
}

void ElDistMatrixEmpty_d( ElDistMatrix_d* AHandle )
{
    try { RCDDM_d(AHandle)->ADM->Empty(); }
    CATCH
}

void ElDistMatrixEmpty_c( ElDistMatrix_c* AHandle )
{
    try { RCDDM_c(AHandle)->ADM->Empty(); }
    CATCH
}

void ElDistMatrixEmpty_z( ElDistMatrix_z* AHandle )
{
    try { RCDDM_z(AHandle)->ADM->Empty(); }
    CATCH
}

// AbstractDistMatrix<T>::EmptyData()
// ----------------------------------
void ElDistMatrixEmptyData_s( ElDistMatrix_s* AHandle )
{
    try { RCDDM_s(AHandle)->ADM->EmptyData(); }
    CATCH
}

void ElDistMatrixEmptyData_d( ElDistMatrix_d* AHandle )
{
    try { RCDDM_d(AHandle)->ADM->EmptyData(); }
    CATCH
}

void ElDistMatrixEmptyData_c( ElDistMatrix_c* AHandle )
{
    try { RCDDM_c(AHandle)->ADM->EmptyData(); }
    CATCH
}

void ElDistMatrixEmptyData_z( ElDistMatrix_z* AHandle )
{
    try { RCDDM_z(AHandle)->ADM->EmptyData(); }
    CATCH
}

// AbstractDistMatrix<T>::SetGrid()
// --------------------------------
void ElDistMatrixSetGrid_s( ElDistMatrix_s* AHandle, const ElGrid* gridHandle )
{
    try { RCDDM_s(AHandle)->ADM->SetGrid(*RCG_const(gridHandle)); }
    CATCH
}

void ElDistMatrixSetGrid_d( ElDistMatrix_d* AHandle, const ElGrid* gridHandle )
{
    try { RCDDM_d(AHandle)->ADM->SetGrid(*RCG_const(gridHandle)); }
    CATCH
}

void ElDistMatrixSetGrid_c( ElDistMatrix_c* AHandle, const ElGrid* gridHandle )
{
    try { RCDDM_c(AHandle)->ADM->SetGrid(*RCG_const(gridHandle)); }
    CATCH
}

void ElDistMatrixSetGrid_z( ElDistMatrix_z* AHandle, const ElGrid* gridHandle )
{
    try { RCDDM_z(AHandle)->ADM->SetGrid(*RCG_const(gridHandle)); }
    CATCH
}

// AbstractDistMatrix<T>::Resize( Int height, Int width )
// ------------------------------------------------------
void ElDistMatrixResize_s( ElDistMatrix_s* AHandle, ElInt height, ElInt width )
{
    try { RCDDM_s(AHandle)->ADM->Resize(height,width); }
    CATCH
}

void ElDistMatrixResize_d( ElDistMatrix_d* AHandle, ElInt height, ElInt width )
{
    try { RCDDM_d(AHandle)->ADM->Resize(height,width); }
    CATCH
}

void ElDistMatrixResize_c( ElDistMatrix_c* AHandle, ElInt height, ElInt width )
{
    try { RCDDM_c(AHandle)->ADM->Resize(height,width); }
    CATCH
}

void ElDistMatrixResize_z( ElDistMatrix_z* AHandle, ElInt height, ElInt width )
{
    try { RCDDM_z(AHandle)->ADM->Resize(height,width); }
    CATCH
}

// AbstractDistMatrix<T>::Resize( Int height, Int width, Int ldim )
// ----------------------------------------------------------------
void ElDistMatrixResizeWithLDim_s
( ElDistMatrix_s* AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { RCDDM_s(AHandle)->ADM->Resize(height,width,ldim); }
    CATCH
}

void ElDistMatrixResizeWithLDim_d
( ElDistMatrix_d* AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { RCDDM_d(AHandle)->ADM->Resize(height,width,ldim); }
    CATCH
}

void ElDistMatrixResizeWithLDim_c
( ElDistMatrix_c* AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { RCDDM_c(AHandle)->ADM->Resize(height,width,ldim); }
    CATCH
}

void ElDistMatrixResizeWithLDim_z
( ElDistMatrix_z* AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { RCDDM_z(AHandle)->ADM->Resize(height,width,ldim); }
    CATCH
}

// AbstractDistMatrix<T>::MakeConsistent()
// ---------------------------------------
void ElDistMatrixMakeConsistent_s
( ElDistMatrix_s* AHandle, bool includeViewers )
{
    try { RCDDM_s(AHandle)->ADM->MakeConsistent(includeViewers); }
    CATCH
}

void ElDistMatrixMakeConsistent_d
( ElDistMatrix_d* AHandle, bool includeViewers )
{
    try { RCDDM_d(AHandle)->ADM->MakeConsistent(includeViewers); }
    CATCH
}

void ElDistMatrixMakeConsistent_c
( ElDistMatrix_c* AHandle, bool includeViewers )
{
    try { RCDDM_c(AHandle)->ADM->MakeConsistent(includeViewers); }
    CATCH
}

void ElDistMatrixMakeConsistent_z
( ElDistMatrix_z* AHandle, bool includeViewers )
{
    try { RCDDM_z(AHandle)->ADM->MakeConsistent(includeViewers); }
    CATCH
}

// AbstractDistMatrix<T>::MakeSizeConsistent()
// -------------------------------------------
void ElDistMatrixMakeSizeConsistent_s
( ElDistMatrix_s* AHandle, bool includeViewers )
{
    try { RCDDM_s(AHandle)->ADM->MakeSizeConsistent(includeViewers); }
    CATCH
}

void ElDistMatrixMakeSizeConsistent_d
( ElDistMatrix_d* AHandle, bool includeViewers )
{
    try { RCDDM_d(AHandle)->ADM->MakeSizeConsistent(includeViewers); }
    CATCH
}

void ElDistMatrixMakeSizeConsistent_c
( ElDistMatrix_c* AHandle, bool includeViewers )
{
    try { RCDDM_c(AHandle)->ADM->MakeSizeConsistent(includeViewers); }
    CATCH
}

void ElDistMatrixMakeSizeConsistent_z
( ElDistMatrix_z* AHandle, bool includeViewers )
{
    try { RCDDM_z(AHandle)->ADM->MakeSizeConsistent(includeViewers); }
    CATCH
}

// AbstractDistMatrix<T>::Get( Int i, Int j ) const
// ------------------------------------------------
float ElDistMatrixGet_s( const ElDistMatrix_s* AHandle, ElInt i, ElInt j )
{
    float alpha = -1;
    try { alpha = RCDDM_s_const(AHandle)->ADM->Get(i,j); }
    CATCH
    return alpha;
}

double ElDistMatrixGet_d( const ElDistMatrix_d* AHandle, ElInt i, ElInt j )
{
    double alpha = -1;
    try { alpha = RCDDM_d_const(AHandle)->ADM->Get(i,j); }
    CATCH
    return alpha;
}

void ElDistMatrixGet_c
( const ElDistMatrix_c* AHandle, ElInt i, ElInt j, void* alpha )
{
    try { *RCB_c(alpha) = RCDDM_c_const(AHandle)->ADM->Get(i,j); }
    CATCH
}

void ElDistMatrixGet_z
( const ElDistMatrix_z* AHandle, ElInt i, ElInt j, void* alpha )
{
    try { *RCB_z(alpha) = RCDDM_z_const(AHandle)->ADM->Get(i,j); }
    CATCH
}

// AbstractDistMatrix<T>::GetRealPart( Int i, Int j ) const
// --------------------------------------------------------
float ElDistMatrixGetRealPart_c
( const ElDistMatrix_c* AHandle, ElInt i, ElInt j )
{
    float alpha = -1;
    try { alpha = RCDDM_c_const(AHandle)->ADM->GetRealPart(i,j); }
    CATCH
    return alpha;
}

double ElDistMatrixGetRealPart_z
( const ElDistMatrix_z* AHandle, ElInt i, ElInt j )
{
    double alpha = -1;
    try { alpha = RCDDM_c_const(AHandle)->ADM->GetRealPart(i,j); }
    CATCH
    return alpha;
}

// AbstractDistMatrix<T>::GetImagPart( Int i, Int j ) const
// --------------------------------------------------------
float ElDistMatrixGetImagPart_c
( const ElDistMatrix_c* AHandle, ElInt i, ElInt j )
{
    float alpha = -1;
    try { alpha = RCDDM_c_const(AHandle)->ADM->GetImagPart(i,j); }
    CATCH
    return alpha;
}

double ElDistMatrixGetImagPart_z
( const ElDistMatrix_z* AHandle, ElInt i, ElInt j )
{
    double alpha = -1;
    try { alpha = RCDDM_c_const(AHandle)->ADM->GetImagPart(i,j); }
    CATCH
    return alpha;
}

// AbstractDistMatrix<T>::Set( Int i, Int j, T alpha )
// ---------------------------------------------------
void ElDistMatrixSet_s( ElDistMatrix_s* AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCDDM_s(AHandle)->ADM->Set(i,j,alpha); }
    CATCH
}

void ElDistMatrixSet_d
( ElDistMatrix_d* AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCDDM_d(AHandle)->ADM->Set(i,j,alpha); }
    CATCH
}

void ElDistMatrixSet_c( ElDistMatrix_c* AHandle, ElInt i, ElInt j, void* alpha )
{
    try { RCDDM_c(AHandle)->ADM->Set(i,j,*RCB_c(alpha)); }
    CATCH
}

void ElDistMatrixSet_z( ElDistMatrix_z* AHandle, ElInt i, ElInt j, void* alpha )
{
    try { RCDDM_z(AHandle)->ADM->Set(i,j,*RCB_z(alpha)); }
    CATCH
}

// B = A
// -----
void ElDistMatrixCopy_s
( const ElDistMatrix_s* AHandle, ElDistMatrix_s* BHandle )
{ 
    try { Copy( *RCDDM_s_const(AHandle), *RCDDM_s(BHandle) ); }
    CATCH
}

void ElDistMatrixCopy_d
( const ElDistMatrix_d* AHandle, ElDistMatrix_d* BHandle )
{
    try { Copy( *RCDDM_d_const(AHandle), *RCDDM_d(BHandle) ); }
    CATCH
}

void ElDistMatrixCopy_c
( const ElDistMatrix_c* AHandle, ElDistMatrix_c* BHandle )
{
    try { Copy( *RCDDM_c_const(AHandle), *RCDDM_c(BHandle) ); }
    CATCH
}

void ElDistMatrixCopy_z
( const ElDistMatrix_z* AHandle, ElDistMatrix_z* BHandle )
{
    try { Copy( *RCDDM_z_const(AHandle), *RCDDM_z(BHandle) ); }
    CATCH
}

/// TODO: Many more member functions

} // extern "C"
