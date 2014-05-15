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

// DistMatrix<T,MC,MR>::DistMatrix( const Grid& g )
// ------------------------------------------------
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

// DistMatrix<T,U,V>::DistMatrix( const Grid& g )
// ----------------------------------------------
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

// DistMatrix<T,U,V>::~DistMatrix()
// --------------------------------
void ElDistMatrixDestroy_s( const ElDistMatrix_s* AHandle )
{ delete RCDDM_s_const(AHandle); }

void ElDistMatrixDestroy_d( const ElDistMatrix_d* AHandle )
{ delete RCDDM_d_const(AHandle); }

void ElDistMatrixDestroy_c( const ElDistMatrix_c* AHandle )
{ delete RCDDM_c_const(AHandle); }

void ElDistMatrixDestroy_z( const ElDistMatrix_z* AHandle )
{ delete RCDDM_z_const(AHandle); }

// void DistMatrix<T,U,V>::Empty()
// -------------------------------
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

// void DistMatrix<T,U,V>::EmptyData()
// -----------------------------------
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

// void DistMatrix<T,U,V>::SetGrid( const Grid& g )
// ------------------------------------------------
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

// void DistMatrix<T,U,V>::Resize( Int height, Int width )
// -------------------------------------------------------
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

// void DistMatrix<T,U,V>::Resize( Int height, Int width, Int ldim )
// -----------------------------------------------------------------
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

// void DistMatrix<T,U,V>::MakeConsistent()
// ----------------------------------------
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

// void DistMatrix<T,U,V>::MakeSizeConsistent()
// --------------------------------------------
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

// TODO: Fill in a large number of missing routines here
// =====================================================

// const Grid& DistMatrix<T,U,V>::Grid() const
// -------------------------------------------
const ElGrid* ElDistMatrixGrid_s( const ElDistMatrix_s* AHandle )
{ return reinterpret_cast<const ElGrid*>
         (&RCDDM_s_const(AHandle)->ADM->Grid()); }

const ElGrid* ElDistMatrixGrid_d( const ElDistMatrix_d* AHandle )
{ return reinterpret_cast<const ElGrid*>
         (&RCDDM_d_const(AHandle)->ADM->Grid()); }

const ElGrid* ElDistMatrixGrid_c( const ElDistMatrix_c* AHandle )
{ return reinterpret_cast<const ElGrid*>
         (&RCDDM_c_const(AHandle)->ADM->Grid()); }

const ElGrid* ElDistMatrixGrid_z( const ElDistMatrix_z* AHandle )
{ return reinterpret_cast<const ElGrid*>
         (&RCDDM_z_const(AHandle)->ADM->Grid()); }

// void DistMatrix<T,U,V>::Get( Int i, Int j ) const
// -------------------------------------------------
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

// void DistMatrix<T,U,V>::GetRealPart( Int i, Int j ) const
// ---------------------------------------------------------
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

// void DistMatrix<T,U,V>::GetImagPart( Int i, Int j ) const
// ---------------------------------------------------------
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

// void DistMatrix<T,U,V>::Set( Int i, Int j, T alpha )
// ----------------------------------------------------
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

// void DistMatrix<T,U,V>::SetRealPart( Int i, Int j, Base<T> alpha )
// ------------------------------------------------------------------
void ElDistMatrixSetRealPart_c
( ElDistMatrix_c* AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCDDM_c(AHandle)->ADM->SetRealPart(i,j,alpha); }
    CATCH
}

void ElDistMatrixSetRealPart_z
( ElDistMatrix_z* AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCDDM_z(AHandle)->ADM->SetRealPart(i,j,alpha); }
    CATCH
}

// void DistMatrix<T,U,V>::SetImagPart( Int i, Int j, Base<T> alpha )
// ------------------------------------------------------------------
void ElDistMatrixSetImagPart_c
( ElDistMatrix_c* AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCDDM_c(AHandle)->ADM->SetImagPart(i,j,alpha); }
    CATCH
}

void ElDistMatrixSetImagPart_z
( ElDistMatrix_z* AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCDDM_z(AHandle)->ADM->SetImagPart(i,j,alpha); }
    CATCH
}

// void DistMatrix<T,U,V>::Update( Int i, Int j, T alpha )
// -------------------------------------------------------
void ElDistMatrixUpdate_s
( ElDistMatrix_s* AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCDDM_s(AHandle)->ADM->Update(i,j,alpha); }
    CATCH
}

void ElDistMatrixUpdate_d
( ElDistMatrix_d* AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCDDM_d(AHandle)->ADM->Update(i,j,alpha); }
    CATCH
}

void ElDistMatrixUpdate_c
( ElDistMatrix_c* AHandle, ElInt i, ElInt j, void* alpha )
{
    try { RCDDM_c(AHandle)->ADM->Update(i,j,*RCB_c(alpha)); }
    CATCH
}

void ElDistMatrixUpdate_z
( ElDistMatrix_z* AHandle, ElInt i, ElInt j, void* alpha )
{
    try { RCDDM_z(AHandle)->ADM->Update(i,j,*RCB_z(alpha)); }
    CATCH
}

// void DistMatrix<T,U,V>::UpdateRealPart( Int i, Int j, Base<T> alpha )
// ---------------------------------------------------------------------
void ElDistMatrixUpdateRealPart_c
( ElDistMatrix_c* AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCDDM_c(AHandle)->ADM->UpdateRealPart(i,j,alpha); }
    CATCH
}

void ElDistMatrixUpdateRealPart_z
( ElDistMatrix_z* AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCDDM_z(AHandle)->ADM->UpdateRealPart(i,j,alpha); }
    CATCH
}

// void DistMatrix<T,U,V>::UpdateImagPart( Int i, Int j, Base<T> alpha )
// ---------------------------------------------------------------------
void ElDistMatrixUpdateImagPart_c
( ElDistMatrix_c* AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCDDM_c(AHandle)->ADM->UpdateImagPart(i,j,alpha); }
    CATCH
}

void ElDistMatrixUpdateImagPart_z
( ElDistMatrix_z* AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCDDM_z(AHandle)->ADM->UpdateImagPart(i,j,alpha); }
    CATCH
}

// void DistMatrix<T,U,V>::MakeReal( Int i, Int j )
// ------------------------------------------------
void ElDistMatrixMakeReal_c( ElDistMatrix_c* AHandle, ElInt i, ElInt j )
{
    try { RCDDM_c(AHandle)->ADM->MakeReal(i,j); }
    CATCH
}

void ElDistMatrixMakeReal_z( ElDistMatrix_z* AHandle, ElInt i, ElInt j )
{
    try { RCDDM_z(AHandle)->ADM->MakeReal(i,j); }
    CATCH
}

// void DistMatrix<T,U,V>::Conjugate( Int i, Int j )
// -------------------------------------------------
void ElDistMatrixConjugate_c( ElDistMatrix_c* AHandle, ElInt i, ElInt j )
{
    try { RCDDM_c(AHandle)->ADM->Conjugate(i,j); }
    CATCH
}

void ElDistMatrixConjugate_z( ElDistMatrix_z* AHandle, ElInt i, ElInt j )
{
    try { RCDDM_z(AHandle)->ADM->Conjugate(i,j); }
    CATCH
}

// DistMatrix<T,UDiag,VDiag> DistMatrix<T,U,V>::GetDiagonal( Int offset ) const
// ----------------------------------------------------------------------------
ElDistMatrix_s* ElDistMatrixGetDiagonal_s
( const ElDistMatrix_s* AHandle, ElInt offset )
{
    ElDistMatrix_s* dHandle = 0;
    try 
    {
        auto ADyn = RCDDM_s_const(AHandle);
        ElDist U = static_cast<ElDist>(DiagColDist( ADyn->U, ADyn->V ));
        ElDist V = static_cast<ElDist>(DiagRowDist( ADyn->U, ADyn->V ));
        const ElGrid* gridHandle = ElDistMatrixGrid_s( AHandle );
        dHandle = ElDistMatrixCreateSpecific_s( U, V, gridHandle );
    }
    CATCH
    return dHandle;
}

ElDistMatrix_d* ElDistMatrixGetDiagonal_d
( const ElDistMatrix_d* AHandle, ElInt offset )
{
    ElDistMatrix_d* dHandle = 0;
    try 
    {
        auto ADyn = RCDDM_d_const(AHandle);
        ElDist U = static_cast<ElDist>(DiagColDist( ADyn->U, ADyn->V ));
        ElDist V = static_cast<ElDist>(DiagRowDist( ADyn->U, ADyn->V ));
        const ElGrid* gridHandle = ElDistMatrixGrid_d( AHandle );
        dHandle = ElDistMatrixCreateSpecific_d( U, V, gridHandle );
    }
    CATCH
    return dHandle;
}

ElDistMatrix_c* ElDistMatrixGetDiagonal_c
( const ElDistMatrix_c* AHandle, ElInt offset )
{
    ElDistMatrix_c* dHandle = 0;
    try 
    {
        auto ADyn = RCDDM_c_const(AHandle);
        ElDist U = static_cast<ElDist>(DiagColDist( ADyn->U, ADyn->V ));
        ElDist V = static_cast<ElDist>(DiagRowDist( ADyn->U, ADyn->V ));
        const ElGrid* gridHandle = ElDistMatrixGrid_c( AHandle );
        dHandle = ElDistMatrixCreateSpecific_c( U, V, gridHandle );
    }
    CATCH
    return dHandle;
}

ElDistMatrix_z* ElDistMatrixGetDiagonal_z
( const ElDistMatrix_z* AHandle, ElInt offset )
{
    ElDistMatrix_z* dHandle = 0;
    try 
    {
        auto ADyn = RCDDM_z_const(AHandle);
        ElDist U = static_cast<ElDist>(DiagColDist( ADyn->U, ADyn->V ));
        ElDist V = static_cast<ElDist>(DiagRowDist( ADyn->U, ADyn->V ));
        const ElGrid* gridHandle = ElDistMatrixGrid_z( AHandle );
        dHandle = ElDistMatrixCreateSpecific_z( U, V, gridHandle );
    }
    CATCH
    return dHandle;
}

// TODO: More diagonal manipulation
// ================================

// DistMatrix<T,STAR,STAR> DistMatrix<T,U,V>::GetSubmatrix
// ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
// --------------------------------------------------------------------------
ElDistMatrix_s* ElDistMatrixGetSubmatrix_s
( const ElDistMatrix_s* AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds )
{
    ElDistMatrix_s* ASubHandle = 0;
    try
    {
        const ElGrid* gridHandle = ElDistMatrixGrid_s( AHandle );
        ASubHandle = ElDistMatrixCreateSpecific_s
                     ( EL_STAR, EL_STAR, gridHandle );
        auto ASubDyn = RCDDM_s(ASubHandle);

        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        RCDDM_s_const(AHandle)->ADM->GetSubmatrix
        ( rowIndVec, colIndVec, 
          *reinterpret_cast<DistMatrix<float,STAR,STAR>*>(ASubDyn->ADM) );
    }
    CATCH
    return ASubHandle;
}

ElDistMatrix_d* ElDistMatrixGetSubmatrix_d
( const ElDistMatrix_d* AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds )
{
    ElDistMatrix_d* ASubHandle = 0;
    try
    {
        const ElGrid* gridHandle = ElDistMatrixGrid_d( AHandle );
        ASubHandle = ElDistMatrixCreateSpecific_d
                     ( EL_STAR, EL_STAR, gridHandle );
        auto ASubDyn = RCDDM_d(ASubHandle);

        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        RCDDM_d_const(AHandle)->ADM->GetSubmatrix
        ( rowIndVec, colIndVec, 
          *reinterpret_cast<DistMatrix<double,STAR,STAR>*>(ASubDyn->ADM) );
    }
    CATCH
    return ASubHandle;
}

ElDistMatrix_c* ElDistMatrixGetSubmatrix_c
( const ElDistMatrix_c* AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds )
{
    ElDistMatrix_c* ASubHandle = 0;
    try
    {
        const ElGrid* gridHandle = ElDistMatrixGrid_c( AHandle );
        ASubHandle = ElDistMatrixCreateSpecific_c
                     ( EL_STAR, EL_STAR, gridHandle );
        auto ASubDyn = RCDDM_c(ASubHandle);

        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        RCDDM_c_const(AHandle)->ADM->GetSubmatrix
        ( rowIndVec, colIndVec, 
          *reinterpret_cast<DistMatrix<Complex<float>,STAR,STAR>*>
          (ASubDyn->ADM) );
    }
    CATCH
    return ASubHandle;
}

ElDistMatrix_z* ElDistMatrixGetSubmatrix_z
( const ElDistMatrix_z* AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds )
{
    ElDistMatrix_z* ASubHandle = 0;
    try
    {
        const ElGrid* gridHandle = ElDistMatrixGrid_z( AHandle );
        ASubHandle = ElDistMatrixCreateSpecific_z
                     ( EL_STAR, EL_STAR, gridHandle );
        auto ASubDyn = RCDDM_z(ASubHandle);

        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        RCDDM_z_const(AHandle)->ADM->GetSubmatrix
        ( rowIndVec, colIndVec, 
          *reinterpret_cast<DistMatrix<Complex<double>,STAR,STAR>*>
          (ASubDyn->ADM) );
    }
    CATCH
    return ASubHandle;
}

/// TODO: Many more member functions

} // extern "C"
