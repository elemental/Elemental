/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.
   Copyright (c) 2014, Jed Brown
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El-C.h"
using namespace El;

template<typename T>
ElError ElDistMatrixCreateSpecific
( ElDist U_C, ElDist V_C, ElConstGrid gridHandle, AbstractDistMatrix<T>** A )
{
    try 
    {
        Dist U = Reinterpret(U_C);
        Dist V = Reinterpret(V_C);
        auto grid = Reinterpret(gridHandle);

        #define GUARD(CDIST,RDIST) U == CDIST && V == RDIST
        #define PAYLOAD(CDIST,RDIST) *A = new DistMatrix<T,CDIST,RDIST>(*grid);
        #include "El/macros/GuardAndPayload.h"
        #undef GUARD
        #undef PAYLOAD
    }
    EL_CATCH
    return EL_SUCCESS;
}

// DistMatrix<T,UDiag,VDiag> DistMatrix<T,U,V>::GetDiagonal( Int offset ) const
// ----------------------------------------------------------------------------
template<typename T>
ElError ElDistMatrixGetDiagonal
( const AbstractDistMatrix<T>* AAbs, Int offset, 
        AbstractDistMatrix<T>** dAbs )
{
    const Grid& g = AAbs->Grid();
    Dist U = AAbs->DistData().colDist;
    Dist V = AAbs->DistData().rowDist;
    try 
    {
        #define GUARD(CDIST,RDIST) U == CDIST && V == RDIST
        #define PAYLOAD(CDIST,RDIST) \
            auto A = dynamic_cast<const DistMatrix<T,CDIST,RDIST>*>(AAbs); \
            DynamicCastCheck(A); \
            auto* d = new DistMatrix<T,DiagColDist<CDIST,RDIST>(), \
                                       DiagRowDist<CDIST,RDIST>()>(g); \
            A->GetDiagonal( *d, offset ); \
            *dAbs = d;
        #include "El/macros/GuardAndPayload.h"
        #undef GUARD
        #undef PAYLOAD
    }
    EL_CATCH
    return EL_SUCCESS;
}

// DistMatrix<Base<T>,UDiag,VDiag> 
// DistMatrix<T,U,V>::GetRealPartOfDiagonal( Int offset ) const
// ------------------------------------------------------------
template<typename T>
ElError ElDistMatrixGetRealPartOfDiagonal
( const AbstractDistMatrix<Complex<T>>* AAbs, Int offset, 
        AbstractDistMatrix<T>** dAbs )
{
    const Grid& g = AAbs->Grid();
    Dist U = AAbs->DistData().colDist;
    Dist V = AAbs->DistData().rowDist;
    try 
    {
        #define GUARD(CDIST,RDIST) U == CDIST && V == RDIST
        #define PAYLOAD(CDIST,RDIST) \
            auto A = \
              dynamic_cast<const DistMatrix<Complex<T>,CDIST,RDIST>*>(AAbs); \
            DynamicCastCheck(A); \
            auto* d = new DistMatrix<T,DiagColDist<CDIST,RDIST>(), \
                                       DiagRowDist<CDIST,RDIST>()>(g); \
            A->GetRealPartOfDiagonal( *d, offset ); \
            *dAbs = d;
        #include "El/macros/GuardAndPayload.h"
        #undef GUARD
        #undef PAYLOAD
    }
    EL_CATCH
    return EL_SUCCESS;
}

// DistMatrix<Base<T>,UDiag,VDiag> 
// DistMatrix<T,U,V>::GetImagPartOfDiagonal( Int offset ) const
// ------------------------------------------------------------
template<typename T>
ElError ElDistMatrixGetImagPartOfDiagonal
( const AbstractDistMatrix<Complex<T>>* AAbs, Int offset, 
        AbstractDistMatrix<T>** dAbs )
{
    const Grid& g = AAbs->Grid();
    Dist U = AAbs->DistData().colDist;
    Dist V = AAbs->DistData().rowDist;
    try 
    {
        #define GUARD(CDIST,RDIST) U == CDIST && V == RDIST
        #define PAYLOAD(CDIST,RDIST) \
            auto A = \
              dynamic_cast<const DistMatrix<Complex<T>,CDIST,RDIST>*>(AAbs); \
            DynamicCastCheck(A); \
            auto* d = new DistMatrix<T,DiagColDist<CDIST,RDIST>(), \
                                       DiagRowDist<CDIST,RDIST>()>(g); \
            A->GetImagPartOfDiagonal( *d, offset ); \
            *dAbs = d;
        #include "El/macros/GuardAndPayload.h"
        #undef GUARD
        #undef PAYLOAD
    }
    EL_CATCH
    return EL_SUCCESS;
}

extern "C" {

#define DISTMATRIX_CREATE(SIG,T) \
  /* DistMatrix<T,MC,MR>::DistMatrix( const Grid& g ) */ \
  ElError ElDistMatrixCreate_ ## SIG \
  ( ElConstGrid gridHandle, ElDistMatrix_ ## SIG *AHandle ) \
  { EL_TRY( auto grid = Reinterpret(gridHandle); \
            *AHandle = Reinterpret( new DistMatrix<T>(*grid) ) ) } \
  /* DistMatrix<T,U,V>::DistMatrix( const Grid& g ) */ \
  ElError ElDistMatrixCreateSpecific_ ## SIG \
  ( ElDist U, ElDist V, ElConstGrid gridHandle, \
    ElDistMatrix_ ## SIG *AHandle ) \
  { \
    AbstractDistMatrix<T>* ADM; \
    ElError error = ElDistMatrixCreateSpecific( U, V, gridHandle, &ADM ); \
    *AHandle = Reinterpret(ADM); \
    return error; \
  } \
  /* DistMatrix<T,U,V>::~DistMatrix() */ \
  ElError ElDistMatrixDestroy_ ## SIG ( ElConstDistMatrix_ ## SIG AHandle ) \
  { EL_TRY( Reinterpret(AHandle) ) } 

#define DISTMATRIX_RECONFIG(SIG,T) \
  /* void DistMatrix<T,U,V>::Empty() */ \
  ElError ElDistMatrixEmpty_ ## SIG ( ElDistMatrix_ ## SIG AHandle ) \
  { EL_TRY( Reinterpret(AHandle)->Empty() ) } \
  /* void DistMatrix<T,U,V>::EmptyData() */ \
  ElError ElDistMatrixEmptyData_ ## SIG ( ElDistMatrix_ ## SIG AHandle ) \
  { EL_TRY( Reinterpret(AHandle)->EmptyData() ) } \
  /* void DistMatrix<T,U,V>::SetGrid( const Grid& g ) */ \
  ElError ElDistMatrixSetGrid_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElConstGrid gridHandle ) \
  { EL_TRY( Reinterpret(AHandle)->SetGrid(*Reinterpret(gridHandle)) ) } \
  /* void DistMatrix<T,U,V>::Resize( Int height, Int width ) */ \
  ElError ElDistMatrixResize_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt height, ElInt width ) \
  { EL_TRY( Reinterpret(AHandle)->Resize(height,width) ) } \
  /* void DistMatrix<T,U,V>::Resize( Int height, Int width, Int ldim ) */ \
  ElError ElDistMatrixResizeWithLDim_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt height, ElInt width, ElInt ldim ) \
  { EL_TRY( Reinterpret(AHandle)->Resize(height,width,ldim) ) } \
  /* void DistMatrix<T,U,V>::MakeConsistent() */ \
  ElError ElDistMatrixMakeConsistent_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, bool includeViewers ) \
  { EL_TRY( Reinterpret(AHandle)->MakeConsistent(includeViewers) ) } \
  /* void DistMatrix<T,U,V>::MakeSizeConsistent() */ \
  ElError ElDistMatrixMakeSizeConsistent_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, bool includeViewers ) \
  { EL_TRY( Reinterpret(AHandle)->MakeSizeConsistent(includeViewers) ) } \
  /* void DistMatrix<T,U,V>::Align \
     ( Int colAlign, Int rowAlign, bool constrain ) */ \
  ElError ElDistMatrixAlign_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt colAlign, ElInt rowAlign, \
    bool constrain ) \
  { EL_TRY( Reinterpret(AHandle)->Align(colAlign,rowAlign,constrain) ) } \
  /* void DistMatrix<T,U,V>::AlignCols( Int colAlign, bool constrain ) */ \
  ElError ElDistMatrixAlignCols_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt colAlign, bool constrain ) \
  { EL_TRY( Reinterpret(AHandle)->AlignCols(colAlign,constrain) ) } \
   /* void DistMatrix<T,U,V>::AlignRows( Int rowAlign, bool constrain ) */ \
  ElError ElDistMatrixAlignRows_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt rowAlign, bool constrain ) \
  { EL_TRY( Reinterpret(AHandle)->AlignRows(rowAlign,constrain) ) } \
  /* void DistMatrix<T,U,V>::FreeAlignments() */ \
  ElError ElDistMatrixFreeAlignments_ ## SIG ( ElDistMatrix_ ## SIG AHandle ) \
  { EL_TRY( Reinterpret(AHandle)->FreeAlignments() ) } \
  /* void DistMatrix<T,U,V>::SetRoot( Int root ) */ \
  ElError ElDistMatrixSetRoot_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt root ) \
  { EL_TRY( Reinterpret(AHandle)->SetRoot(root) ) } \
  /* void DistMatrix<T,U,V>::AlignWith \
     ( const DistData& data, bool constrain ) */ \
  ElError ElDistMatrixAlignWith_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElDistData distData, bool constrain ) \
  { EL_TRY( Reinterpret(AHandle)->AlignWith \
            ( Reinterpret(distData), constrain ) ) } \
  /* void DistMatrix<T,U,V>::AlignColsWith \
     ( const DistData& data, bool constrain ) */ \
  ElError ElDistMatrixAlignColsWith_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElDistData distData, bool constrain ) \
  { EL_TRY( Reinterpret(AHandle)->AlignColsWith \
            ( Reinterpret(distData), constrain ) ) } \
  /* void DistMatrix<T,U,V>::AlignRowsWith \
     ( const DistData& data, bool constrain ) */ \
  ElError ElDistMatrixAlignRowsWith_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElDistData distData, bool constrain ) \
  { EL_TRY( Reinterpret(AHandle)->AlignRowsWith \
            ( Reinterpret(distData), constrain ) ) } \
  /* void DistMatrix<T,U,V>::AlignAndResize \
     ( Int colAlign, Int rowAlign, Int height, Int width, \
       bool force, bool constrain ) */ \
  ElError ElDistMatrixAlignAndResize_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, \
    ElInt colAlign, ElInt rowAlign, ElInt height, ElInt width, \
    bool force, bool constrain ) \
  { EL_TRY( Reinterpret(AHandle)->AlignAndResize \
            (colAlign,rowAlign,height,width,force,constrain) ) } \
  /* void DistMatrix<T,U,V>::AlignColsAndResize
     ( Int colAlign, Int height, Int width, bool force, bool constrain ) */ \
  ElError ElDistMatrixAlignColsAndResize_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, \
    ElInt colAlign, ElInt height, ElInt width, bool force, bool constrain ) \
  { EL_TRY( Reinterpret(AHandle)->AlignColsAndResize \
            (colAlign,height,width,force,constrain) ) } \
  /* void DistMatrix<T,U,V>::AlignRowsAndResize
     ( Int rowAlign, Int height, Int width, bool force, bool constrain ) */ \
  ElError ElDistMatrixAlignRowsAndResize_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, \
    ElInt rowAlign, ElInt height, ElInt width, bool force, bool constrain ) \
  { EL_TRY( Reinterpret(AHandle)->AlignRowsAndResize \
            (rowAlign,height,width,force,constrain) ) } \
  /* void DistMatrix<T,U,V>::Attach
     ( Int height, Int width, const Grid& grid, Int colAlign, Int rowAlign, 
       T* buffer, Int ldim, Int root ) */ \
  ElError ElDistMatrixAttach_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt height, ElInt width, \
    ElConstGrid gridHandle, ElInt colAlign, ElInt rowAlign, \
    CREFLECT(T)* buffer, ElInt ldim, ElInt root ) \
  { EL_TRY( Reinterpret(AHandle)->Attach \
            (height,width,*Reinterpret(gridHandle),colAlign,rowAlign, \
             Reinterpret(buffer),ldim,root) ) } \
  /* void DistMatrix<T,U,V>::LockedAttach
     ( Int height, Int width, const Grid& grid, Int colAlign, Int rowAlign, 
       const T* buffer, Int ldim, Int root ) */ \
  ElError ElDistMatrixLockedAttach_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt height, ElInt width, \
    ElConstGrid gridHandle, ElInt colAlign, ElInt rowAlign, \
    const CREFLECT(T)* buffer, ElInt ldim, ElInt root ) \
  { EL_TRY( Reinterpret(AHandle)->LockedAttach \
            (height,width,*Reinterpret(gridHandle),colAlign,rowAlign, \
             Reinterpret(buffer),ldim,root) ) }

#define DISTMATRIX_BASIC(SIG,T) \
  /* Int DistMatrix<T,U,V>::Height() const */ \
  ElError ElDistMatrixHeight_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* height ) \
  { EL_TRY( *height = Reinterpret(AHandle)->Height() ) } \
  /* Int DistMatrix<T,U,V>::Width() const */ \
  ElError ElDistMatrixWidth_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* width ) \
  { EL_TRY( *width = Reinterpret(AHandle)->Width() ) } \
  /* Int DistMatrix<T,U,V>::DiagonalLength( Int offset ) const */ \
  ElError ElDistMatrixDiagonalLength_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt offset, ElInt* length ) \
  { EL_TRY( *length = Reinterpret(AHandle)->DiagonalLength(offset) ) } \
  /* bool DistMatrix<T,U,V>::Viewing() const */ \
  ElError ElDistMatrixViewing_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, bool* viewing ) \
  { EL_TRY( *viewing = Reinterpret(AHandle)->Viewing() ) } \
  /* bool DistMatrix<T,U,V>::Locked() const */ \
  ElError ElDistMatrixLocked_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, bool* locked ) \
  { EL_TRY( *locked = Reinterpret(AHandle)->Locked() ) } \
  /* Int DistMatrix<T,U,V>::LocalHeight() const */ \
  ElError ElDistMatrixLocalHeight_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* localHeight ) \
  { EL_TRY( *localHeight = Reinterpret(AHandle)->LocalHeight() ) } \
  /* Int DistMatrix<T,U,V>::LocalWidth() const */ \
  ElError ElDistMatrixLocalWidth_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* localWidth ) \
  { EL_TRY( *localWidth = Reinterpret(AHandle)->LocalWidth() ) } \
  /* Int DistMatrix<T,U,V>::LDim() const */ \
  ElError ElDistMatrixLDim_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* ldim ) \
  { EL_TRY( *ldim = Reinterpret(AHandle)->LDim() ) } \
  /* Matrix<T>& DistMatrix<T,U,V>::Matrix() */ \
  ElError ElDistMatrixMatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElMatrix_ ## SIG *ALocHandle ) \
  { EL_TRY( *ALocHandle = Reinterpret(&Reinterpret(AHandle)->Matrix()) ) } \
  /* const Matrix<T>& DistMatrix<T,U,V>::LockedMatrix() const */ \
  ElError ElDistMatrixLockedMatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElConstMatrix_ ## SIG *ALocHandle ) \
  { EL_TRY( *ALocHandle = Reinterpret \
                          (&Reinterpret(AHandle)->LockedMatrix()) ) } \
  /* size_t DistMatrix<T,U,V>::AllocatedMemory() const */ \
  ElError ElDistMatrixAllocatedMemory_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, size_t* mem ) \
  { EL_TRY( *mem = Reinterpret(AHandle)->AllocatedMemory() ) } \
  /* T* DistMatrix<T,U,V>::Buffer() */ \
  ElError ElDistMatrixBuffer_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, CREFLECT(T)** buffer ) \
  { EL_TRY( *buffer = Reinterpret(Reinterpret(AHandle)->Buffer()) ) } \
  /* const T* DistMatrix<T,U,V>::LockedBuffer() const */ \
  ElError ElDistMatrixLockedBuffer_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, const CREFLECT(T)** buffer ) \
  { EL_TRY( *buffer = Reinterpret(Reinterpret(AHandle)->LockedBuffer()) ) } \
  /* const Grid& DistMatrix<T,U,V>::Grid() const */ \
  ElError ElDistMatrixGrid_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElConstGrid* gridHandle ) \
  { EL_TRY( *gridHandle = Reinterpret(&Reinterpret(AHandle)->Grid()) ) } \
  /* bool DistMatrix<T,U,V>::ColConstrained() const */ \
  ElError ElDistMatrixColConstrained_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, bool* colConst ) \
  { EL_TRY( *colConst = Reinterpret(AHandle)->ColConstrained() ) } \
  /* bool DistMatrix<T,U,V>::RowConstrained() const */ \
  ElError ElDistMatrixRowConstrained_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, bool* rowConst ) \
  { EL_TRY( *rowConst = Reinterpret(AHandle)->RowConstrained() ) } \
  /* bool DistMatrix<T,U,V>::RootConstrained() const */ \
  ElError ElDistMatrixRootConstrained_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, bool* rootConst ) \
  { EL_TRY( *rootConst = Reinterpret(AHandle)->RootConstrained() ) } \
  /* Int DistMatrix<T,U,V>::ColAlign() const */ \
  ElError ElDistMatrixColAlign_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* colAlign ) \
  { EL_TRY( *colAlign = Reinterpret(AHandle)->ColAlign() ) } \
  /* Int DistMatrix<T,U,V>::RowAlign() const */ \
  ElError ElDistMatrixRowAlign_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* rowAlign ) \
  { EL_TRY( *rowAlign = Reinterpret(AHandle)->RowAlign() ) } \
  /* Int DistMatrix<T,U,V>::ColShift() const */ \
  ElError ElDistMatrixColShift_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* colShift ) \
  { EL_TRY( *colShift = Reinterpret(AHandle)->ColShift() ) } \
  /* Int DistMatrix<T,U,V>::RowShift() const */ \
  ElError ElDistMatrixRowShift_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* rowShift ) \
  { EL_TRY( *rowShift = Reinterpret(AHandle)->RowShift() ) } \
  /* Int DistMatrix<T,U,V>::ColRank() const */ \
  ElError ElDistMatrixColRank_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* colRank ) \
  { EL_TRY( *colRank = Reinterpret(AHandle)->ColRank() ) } \
  /* Int DistMatrix<T,U,V>::RowRank() const */ \
  ElError ElDistMatrixRowRank_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* rowRank ) \
  { EL_TRY( *rowRank = Reinterpret(AHandle)->RowRank() ) } \
  /* Int DistMatrix<T,U,V>::PartialColRank() const */ \
  ElError ElDistMatrixPartialColRank_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* rank ) \
  { EL_TRY( *rank = Reinterpret(AHandle)->PartialColRank() ) } \
  /* Int DistMatrix<T,U,V>::PartialRowRank() const */ \
  ElError ElDistMatrixPartialRowRank_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* rank ) \
  { EL_TRY( *rank = Reinterpret(AHandle)->PartialRowRank() ) } \
  /* Int DistMatrix<T,U,V>::PartialUnionColRank() const */ \
  ElError ElDistMatrixPartialUnionColRank_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* rank ) \
  { EL_TRY( *rank = Reinterpret(AHandle)->PartialUnionColRank() ) } \
  /* Int DistMatrix<T,U,V>::PartialUnionRowRank() const */ \
  ElError ElDistMatrixPartialUnionRowRank_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* rank ) \
  { EL_TRY( *rank = Reinterpret(AHandle)->PartialUnionRowRank() ) } \
  /* Int DistMatrix<T,U,V>::DistRank() const */ \
  ElError ElDistMatrixDistRank_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* rank ) \
  { EL_TRY( *rank = Reinterpret(AHandle)->DistRank() ) } \
  /* Int DistMatrix<T,U,V>::CrossRank() const */ \
  ElError ElDistMatrixCrossRank_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* rank ) \
  { EL_TRY( *rank = Reinterpret(AHandle)->CrossRank() ) } \
  /* Int DistMatrix<T,U,V>::RedundantRank() const */ \
  ElError ElDistMatrixRedundantRank_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* rank ) \
  { EL_TRY( *rank = Reinterpret(AHandle)->RedundantRank() ) } \
  /* Int DistMatrix<T,U,V>::Root() const */ \
  ElError ElDistMatrixRoot_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* root ) \
  { EL_TRY( *root = Reinterpret(AHandle)->Root() ) } \
  /* bool DistMatrix<T,U,V>::Participating() const */ \
  ElError ElDistMatrixParticipating_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, bool* participating ) \
  { EL_TRY( *participating = Reinterpret(AHandle)->Participating() ) } \
  /* Int DistMatrix<T,U,V>::RowOwner( Int i ) const */ \
  ElError ElDistMatrixRowOwner_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt i, ElInt* rowOwner ) \
  { EL_TRY( *rowOwner = Reinterpret(AHandle)->RowOwner(i) ) } \
  /* Int DistMatrix<T,U,V>::ColOwner( Int i ) const */ \
  ElError ElDistMatrixColOwner_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt j, ElInt* colOwner ) \
  { EL_TRY( *colOwner = Reinterpret(AHandle)->ColOwner(j) ) } \
  /* Int DistMatrix<T,U,V>::Owner( Int i, Int j ) const */ \
  ElError ElDistMatrixOwner_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt i, ElInt j, ElInt* owner ) \
  { EL_TRY( *owner = Reinterpret(AHandle)->Owner(i,j) ) } \
  /* Int DistMatrix<T,U,V>::LocalRow( Int i ) const */ \
  ElError ElDistMatrixLocalRow_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt i, ElInt* iLoc ) \
  { EL_TRY( *iLoc = Reinterpret(AHandle)->LocalRow(i) ) } \
  /* Int DistMatrix<T,U,V>::LocalCol( Int j ) const */ \
  ElError ElDistMatrixLocalCol_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt j, ElInt* jLoc ) \
  { EL_TRY( *jLoc = Reinterpret(AHandle)->LocalCol(j) ) } \
  /* Int DistMatrix<T,U,V>::LocalRowOffset( Int i ) const */ \
  ElError ElDistMatrixLocalRowOffset_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt i, ElInt* iLoc ) \
  { EL_TRY( *iLoc = Reinterpret(AHandle)->LocalRowOffset(i) ) } \
  /* Int DistMatrix<T,U,V>::LocalColOffset( Int j ) const */ \
  ElError ElDistMatrixLocalColOffset_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt j, ElInt* jLoc ) \
  { EL_TRY( *jLoc = Reinterpret(AHandle)->LocalColOffset(j) ) } \
  /* Int DistMatrix<T,U,V>::GlobalRow( Int iLoc ) const */ \
  ElError ElDistMatrixGlobalRow_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt iLoc, ElInt* i ) \
  { EL_TRY( *i = Reinterpret(AHandle)->GlobalRow(iLoc) ) } \
  /* Int DistMatrix<T,U,V>::GlobalCol( Int j ) const */ \
  ElError ElDistMatrixGlobalCol_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt jLoc, ElInt* j ) \
  { EL_TRY( *j = Reinterpret(AHandle)->GlobalCol(jLoc) ) } \
  /* bool DistMatrix<T,U,V>::IsLocalRow( Int i ) const */ \
  ElError ElDistMatrixIsLocalRow_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt i, bool* isLocal ) \
  { EL_TRY( *isLocal = Reinterpret(AHandle)->IsLocalRow(i) ) } \
  /* bool DistMatrix<T,U,V>::IsLocalCol( Int j ) const */ \
  ElError ElDistMatrixIsLocalCol_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt j, bool* isLocal ) \
  { EL_TRY( *isLocal = Reinterpret(AHandle)->IsLocalCol(j) ) } \
  /* bool DistMatrix<T,U,V>::IsLocal( Int i, Int j ) const */ \
  ElError ElDistMatrixIsLocal_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt i, ElInt j, bool* isLocal ) \
  { EL_TRY( *isLocal = Reinterpret(AHandle)->IsLocal(i,j) ) } \
  /* DistData DistMatrix<T,U,V>::DistData() const */ \
  ElError ElDistMatrixDistData_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElDistData* distData ) \
  { EL_TRY( DistData data = Reinterpret(AHandle)->DistData(); \
            *distData = Reinterpret( data ) ) } \
  /* mpi::Comm DistMatrix<T,U,V>::DistComm() const */ \
  ElError ElDistMatrixDistComm_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, MPI_Comm* comm ) \
  { EL_TRY( *comm = Reinterpret(AHandle)->DistComm().comm ) } \
  /* mpi::Comm DistMatrix<T,U,V>::CrossComm() const */ \
  ElError ElDistMatrixCrossComm_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, MPI_Comm* comm ) \
  { EL_TRY( *comm = Reinterpret(AHandle)->CrossComm().comm ) } \
  /* mpi::Comm DistMatrix<T,U,V>::RedundantComm() const */ \
  ElError ElDistMatrixRedundantComm_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, MPI_Comm* comm ) \
  { EL_TRY( *comm = Reinterpret(AHandle)->RedundantComm().comm ) } \
  /* mpi::Comm DistMatrix<T,U,V>::ColComm() const */ \
  ElError ElDistMatrixColComm_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, MPI_Comm* comm ) \
  { EL_TRY( *comm = Reinterpret(AHandle)->ColComm().comm ) } \
  /* mpi::Comm DistMatrix<T,U,V>::RowComm() const */ \
  ElError ElDistMatrixRowComm_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, MPI_Comm* comm ) \
  { EL_TRY( *comm = Reinterpret(AHandle)->RowComm().comm ) } \
  /* mpi::Comm DistMatrix<T,U,V>::PartialColComm() const */ \
  ElError ElDistMatrixPartialColComm_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, MPI_Comm* comm ) \
  { EL_TRY( *comm = Reinterpret(AHandle)->PartialColComm().comm ) } \
  /* mpi::Comm DistMatrix<T,U,V>::PartialRowComm() const */ \
  ElError ElDistMatrixPartialRowComm_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, MPI_Comm* comm ) \
  { EL_TRY( *comm = Reinterpret(AHandle)->PartialRowComm().comm ) } \
  /* mpi::Comm DistMatrix<T,U,V>::PartialUnionColComm() const */ \
  ElError ElDistMatrixPartialUnionColComm_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, MPI_Comm* comm ) \
  { EL_TRY( *comm = Reinterpret(AHandle)->PartialUnionColComm().comm ) } \
  /* mpi::Comm DistMatrix<T,U,V>::PartialUnionRowComm() const */ \
  ElError ElDistMatrixPartialUnionRowComm_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, MPI_Comm* comm ) \
  { EL_TRY( *comm = Reinterpret(AHandle)->PartialUnionRowComm().comm ) } \
  /* Int DistMatrix<T,U,V>::ColStride() const */ \
  ElError ElDistMatrixColStride_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* stride ) \
  { EL_TRY( *stride = Reinterpret(AHandle)->ColStride() ) } \
  /* Int DistMatrix<T,U,V>::RowStride() const */ \
  ElError ElDistMatrixRowStride_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* stride ) \
  { EL_TRY( *stride = Reinterpret(AHandle)->RowStride() ) } \
  /* Int DistMatrix<T,U,V>::PartialColStride() const */ \
  ElError ElDistMatrixPartialColStride_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* stride ) \
  { EL_TRY( *stride = Reinterpret(AHandle)->PartialColStride() ) } \
  /* Int DistMatrix<T,U,V>::PartialRowStride() const */ \
  ElError ElDistMatrixPartialRowStride_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* stride ) \
  { EL_TRY( *stride = Reinterpret(AHandle)->PartialRowStride() ) } \
  /* Int DistMatrix<T,U,V>::PartialUnionColStride() const */ \
  ElError ElDistMatrixPartialUnionColStride_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* stride ) \
  { EL_TRY( *stride = Reinterpret(AHandle)->PartialUnionColStride() ) } \
  /* Int DistMatrix<T,U,V>::PartialUnionRowStride() const */ \
  ElError ElDistMatrixPartialUnionRowStride_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* stride ) \
  { EL_TRY( *stride = Reinterpret(AHandle)->PartialUnionRowStride() ) } \
  /* Int DistMatrix<T,U,V>::DistSize() const */ \
  ElError ElDistMatrixDistSize_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* commSize ) \
  { EL_TRY( *commSize = Reinterpret(AHandle)->DistSize() ) } \
  /* Int DistMatrix<T,U,V>::CrossSize() const */ \
  ElError ElDistMatrixCrossSize_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* commSize ) \
  { EL_TRY( *commSize = Reinterpret(AHandle)->CrossSize() ) } \
  /* Int DistMatrix<T,U,V>::RedundantSize() const */ \
  ElError ElDistMatrixRedundantSize_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt* commSize ) \
  { EL_TRY( *commSize = Reinterpret(AHandle)->RedundantSize() ) }

#define DISTMATRIX_SINGLEENTRY(SIG,T) \
  /* T DistMatrix<T,U,V>::Get( Int i, Int j ) const */ \
  ElError ElDistMatrixGet_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt i, ElInt j, CREFLECT(T)* val ) \
  { EL_TRY( *val = Reinterpret(Reinterpret(AHandle)->Get(i,j)) ) } \
  /* void DistMatrix<T,U,V>::Set( Int i, Int j, T alpha ) */ \
  ElError ElDistMatrixSet_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt i, ElInt j, CREFLECT(T) alpha ) \
  { EL_TRY( Reinterpret(AHandle)->Set(i,j,Reinterpret(alpha)) ) } \
  /* void DistMatrix<T,U,V>::Update( Int i, Int j, T alpha ) */ \
  ElError ElDistMatrixUpdate_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt i, ElInt j, CREFLECT(T) alpha ) \
  { EL_TRY( Reinterpret(AHandle)->Update(i,j,Reinterpret(alpha)) ) } \
  /* T DistMatrix<T,U,V>::GetLocal( Int iLoc, Int jLoc ) const */ \
  ElError ElDistMatrixGetLocal_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt iLoc, ElInt jLoc, \
    CREFLECT(T)* val ) \
  { EL_TRY( *val = Reinterpret(Reinterpret(AHandle)->GetLocal(iLoc,jLoc)) ) } \
  /* void DistMatrix<T,U,V>::SetLocal( Int iLoc, Int jLoc, T alpha ) */ \
  ElError ElDistMatrixSetLocal_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt iLoc, ElInt jLoc, CREFLECT(T) alpha ) \
  { EL_TRY( Reinterpret(AHandle)->SetLocal(iLoc,jLoc,Reinterpret(alpha)) ) } \
  /* void DistMatrix<T,U,V>::UpdateLocal( Int iLoc, Int jLoc, T alpha ) */ \
  ElError ElDistMatrixUpdateLocal_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt iLoc, ElInt jLoc, CREFLECT(T) alpha ) \
  { EL_TRY( Reinterpret(AHandle)->UpdateLocal \
            (iLoc,jLoc,Reinterpret(alpha)) ) }

#define DISTMATRIX_SINGLEENTRY_COMPLEX(SIG,SIGBASE,T) \
  /* Base<T> DistMatrix<T,U,V>::GetRealPart( Int i, Int j ) const */ \
  ElError ElDistMatrixGetRealPart_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt i, ElInt j, Base<T>* val ) \
  { EL_TRY( *val = Reinterpret(AHandle)->GetRealPart(i,j) ) } \
  /* Base<T> DistMatrix<T,U,V>::GetImagPart( Int i, Int j ) const */ \
  ElError ElDistMatrixGetImagPart_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt i, ElInt j, Base<T>* val ) \
  { EL_TRY( *val = Reinterpret(AHandle)->GetImagPart(i,j) ) } \
  /* void DistMatrix<T,U,V>::SetRealPart( Int i, Int j, Base<T> alpha ) */ \
  ElError ElDistMatrixSetRealPart_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt i, ElInt j, Base<T> alpha ) \
  { EL_TRY( Reinterpret(AHandle)->SetRealPart(i,j,alpha) ) } \
  /* void DistMatrix<T,U,V>::SetImagPart( Int i, Int j, Base<T> alpha ) */ \
  ElError ElDistMatrixSetImagPart_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt i, ElInt j, Base<T> alpha ) \
  { EL_TRY( Reinterpret(AHandle)->SetImagPart(i,j,alpha) ) } \
  /* void DistMatrix<T,U,V>::UpdateRealPart( Int i, Int j, Base<T> alpha ) */ \
  ElError ElDistMatrixUpdateRealPart_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt i, ElInt j, Base<T> alpha ) \
  { EL_TRY( Reinterpret(AHandle)->UpdateRealPart(i,j,alpha) ) } \
  /* void DistMatrix<T,U,V>::UpdateImagPart( Int i, Int j, Base<T> alpha ) */ \
  ElError ElDistMatrixUpdateImagPart_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt i, ElInt j, Base<T> alpha ) \
  { EL_TRY( Reinterpret(AHandle)->UpdateImagPart(i,j,alpha) ) } \
  /* void DistMatrix<T,U,V>::MakeReal( Int i, Int j ) */ \
  ElError ElDistMatrixMakeReal_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt i, ElInt j ) \
  { EL_TRY( Reinterpret(AHandle)->MakeReal(i,j) ) } \
  /* void DistMatrix<T,U,V>::Conjugate( Int i, Int j ) */ \
  ElError ElDistMatrixConjugate_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt i, ElInt j ) \
  { EL_TRY( Reinterpret(AHandle)->Conjugate(i,j) ) } \
  /* Base<T> DistMatrix<T,U,V>::GetLocalRealPart \
     ( Int iLoc, Int jLoc ) const */ \
  ElError ElDistMatrixGetLocalRealPart_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt iLoc, ElInt jLoc, Base<T>* val ) \
  { EL_TRY( *val = Reinterpret(AHandle)->GetLocalRealPart(iLoc,jLoc) ) } \
  /* Base<T> DistMatrix<T,U,V>::GetLocalImagPart \
     ( Int iLoc, Int jLoc ) const */ \
  ElError ElDistMatrixGetLocalImagPart_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt iLoc, ElInt jLoc, Base<T>* val ) \
  { EL_TRY( *val = Reinterpret(AHandle)->GetLocalImagPart(iLoc,jLoc) ) } \
  /* void DistMatrix<T,U,V>::SetLocalRealPart
     ( Int iLoc, Int jLoc, Base<T> alpha ) */ \
  ElError ElDistMatrixSetLocalRealPart_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt iLoc, ElInt jLoc, Base<T> alpha ) \
  { EL_TRY( Reinterpret(AHandle)->SetLocalRealPart(iLoc,jLoc,alpha) ) } \
  /* void DistMatrix<T,U,V>::SetLocalImagPart
     ( Int iLoc, Int jLoc, Base<T> alpha ) */ \
  ElError ElDistMatrixSetLocalImagPart_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt iLoc, ElInt jLoc, Base<T> alpha ) \
  { EL_TRY( Reinterpret(AHandle)->SetLocalImagPart(iLoc,jLoc,alpha) ) } \
  /* void DistMatrix<T,U,V>::UpdateLocalRealPart
     ( Int iLoc, Int jLoc, Base<T> alpha ) */ \
  ElError ElDistMatrixUpdateLocalRealPart_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt iLoc, ElInt jLoc, Base<T> alpha ) \
  { EL_TRY( Reinterpret(AHandle)->UpdateLocalRealPart(iLoc,jLoc,alpha) ) } \
  /* void DistMatrix<T,U,V>::UpdateLocalImagPart
     ( Int iLoc, Int jLoc, Base<T> alpha ) */ \
  ElError ElDistMatrixUpdateLocalImagPart_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt iLoc, ElInt jLoc, Base<T> alpha ) \
  { EL_TRY( Reinterpret(AHandle)->UpdateLocalImagPart(iLoc,jLoc,alpha) ) }

#define DISTMATRIX_DIAGONAL(SIG,T) \
  /* bool DistMatrix<T,U,V>::DiagonalAlignedWith
     ( const DistData& data, Int offset ) const */ \
  ElError ElDistMatrixDiagonalAlignedWith_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElDistData distData, ElInt offset, \
    bool* aligned ) \
  { EL_TRY( *aligned = Reinterpret(AHandle)->DiagonalAlignedWith \
                       (Reinterpret(distData),offset) ) } \
  /* Int DistMatrix<T,U,V>::DiagonalRoot( Int offset ) const */ \
  ElError ElDistMatrixDiagonalRoot_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt offset, ElInt* root ) \
  { EL_TRY( *root = Reinterpret(AHandle)->DiagonalRoot(offset) ) } \
  /* Int DistMatrix<T,U,V>::DiagonalAlign( Int offset ) const */ \
  ElError ElDistMatrixDiagonalAlign_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt offset, ElInt* align ) \
  { EL_TRY( *align = Reinterpret(AHandle)->DiagonalAlign(offset) ) } \
  /* DistMatrix<T,UDiag,VDiag> 
     DistMatrix<T,U,V>::GetDiagonal( Int offset ) const */ \
  ElError ElDistMatrixGetDiagonal_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt offset, \
    ElDistMatrix_ ## SIG *dHandle ) \
  { auto AAbs = Reinterpret(AHandle); \
    AbstractDistMatrix<T>* dAbs; \
    ElError error = ElDistMatrixGetDiagonal( AAbs, offset, &dAbs ); \
    *dHandle = Reinterpret(dAbs); \
    return error; }

#define DISTMATRIX_DIAGONAL_COMPLEX(SIG,SIGBASE,T) \
  /* void DistMatrix<T,U,V>::MakeDiagonalReal( Int offset ) */ \
  ElError ElDistMatrixMakeDiagonalReal_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt offset ) \
  { EL_TRY( Reinterpret(AHandle)->MakeDiagonalReal(offset) ) } \
  /* void DistMatrix<T,U,V>::ConjugateDiagonal( Int offset ) */ \
  ElError ElDistMatrixConjugateDiagonal_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt offset ) \
  { EL_TRY( Reinterpret(AHandle)->ConjugateDiagonal(offset) ) } \
  /* DistMatrix<Base<T>,UDiag,VDiag> 
     DistMatrix<T,U,V>::GetRealPartOfDiagonal( Int offset ) const */ \
  ElError ElDistMatrixGetRealPartOfDiagonal_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt offset, \
    ElDistMatrix_ ## SIGBASE *dHandle ) \
  { auto AAbs = Reinterpret(AHandle); \
    AbstractDistMatrix<Base<T>>* dAbs; \
    ElError error = ElDistMatrixGetRealPartOfDiagonal( AAbs, offset, &dAbs ); \
    *dHandle = Reinterpret(dAbs); \
    return error; } \
  /* DistMatrix<Base<T>,UDiag,VDiag> 
     DistMatrix<T,U,V>::GetImagPartOfDiagonal( Int offset ) const */ \
  ElError ElDistMatrixGetImagPartOfDiagonal_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElInt offset, \
    ElDistMatrix_ ## SIGBASE *dHandle ) \
  { auto AAbs = Reinterpret(AHandle); \
    AbstractDistMatrix<Base<T>>* dAbs; \
    ElError error = ElDistMatrixGetImagPartOfDiagonal( AAbs, offset, &dAbs ); \
    *dHandle = Reinterpret(dAbs); \
    return error; }

#define DISTMATRIX_SUBMATRIX(SIG,T) \
  /* DistMatrix<T,STAR,STAR> DistMatrix<T,U,V>::GetSubmatrix
     ( const std::vector<Int>& rowInds, 
       const std::vector<Int>& colInds ) const */ \
  ElError ElDistMatrixGetSubmatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, \
    ElInt numRowInds, const ElInt* rowInds, \
    ElInt numColInds, const ElInt* colInds, ElDistMatrix_ ## SIG *ASubHandle ) \
  { try { auto A = Reinterpret(AHandle); \
          auto ASub = new DistMatrix<T,STAR,STAR>(A->Grid()); \
          std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds); \
          std::vector<Int> colIndVec(colInds,colInds+numColInds); \
          A->GetSubmatrix( rowIndVec, colIndVec, *ASub ); \
          *ASubHandle = Reinterpret(ASub); } EL_CATCH; return EL_SUCCESS; } \
  /* void DistMatrix<T,U,V>::SetSubmatrix
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
      const DistMatrix<T,STAR,STAR>& ASub ); */ \
  ElError ElDistMatrixSetSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, const ElInt* rowInds, const ElInt* colInds, \
    ElConstDistMatrix_ ## SIG ASubHandle ) \
  { try { \
      auto A = Reinterpret(AHandle); \
      auto ASubADM = Reinterpret(ASubHandle); \
      auto ASub = dynamic_cast<const DistMatrix<T,STAR,STAR>*>(ASubADM); \
      DynamicCastCheck(ASub); \
      const Int numRowInds = ASub->Height(); \
      const Int numColInds = ASub->Width(); \
      std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds), \
                       colIndVec(colInds,colInds+numColInds); \
      A->SetSubmatrix( rowIndVec, colIndVec, *ASub ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void DistMatrix<T,U,V>::UpdateSubmatrix
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
      T alpha, const DistMatrix<T,STAR,STAR>& ASub ); */ \
  ElError ElDistMatrixUpdateSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, const ElInt* rowInds, const ElInt* colInds, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG ASubHandle ) \
  { try { \
      auto A = Reinterpret(AHandle); \
      auto ASubADM = Reinterpret(ASubHandle); \
      auto ASub = dynamic_cast<const DistMatrix<T,STAR,STAR>*>(ASubADM); \
      DynamicCastCheck(ASub); \
      const Int numRowInds = ASub->Height(); \
      const Int numColInds = ASub->Width(); \
      std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds), \
                       colIndVec(colInds,colInds+numColInds); \
      A->UpdateSubmatrix( rowIndVec, colIndVec, Reinterpret(alpha), *ASub ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Matrix<T> DistMatrix<T,U,V>::GetLocalSubmatrix
     ( const std::vector<Int>& rowIndsLoc, 
       const std::vector<Int>& colIndsLoc ) const */ \
  ElError ElDistMatrixGetLocalSubmatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, \
    ElInt numRowInds, const ElInt* rowIndsLoc, \
    ElInt numColInds, const ElInt* colIndsLoc, ElMatrix_ ## SIG *ASubHandle ) \
  { try { \
      auto A = Reinterpret(AHandle); \
      auto ASub = new Matrix<T>; \
      std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds), \
                       colIndVec(colIndsLoc,colIndsLoc+numColInds); \
      A->GetLocalSubmatrix( rowIndVec, colIndVec, *ASub ); \
      *ASubHandle = Reinterpret(ASub); } EL_CATCH; return EL_SUCCESS; } \
  /* void DistMatrix<T,U,V>::SetLocalSubmatrix
     ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
       const Matrix<T>& ASub ); */ \
  ElError ElDistMatrixSetLocalSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, \
    const ElInt* rowIndsLoc, const ElInt* colIndsLoc, \
    ElConstMatrix_ ## SIG ASubHandle ) \
  { try { \
        auto A = Reinterpret(AHandle); \
        auto ASub = Reinterpret(ASubHandle); \
        const Int numRowInds = ASub->Height(); \
        const Int numColInds = ASub->Width(); \
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds), \
                         colIndVec(colIndsLoc,colIndsLoc+numColInds); \
        A->SetLocalSubmatrix( rowIndVec, colIndVec, *ASub ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void DistMatrix<T,U,V>::UpdateLocalSubmatrix
     ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
       const Matrix<T>& ASub ); */ \
  ElError ElDistMatrixUpdateLocalSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, \
    const ElInt* rowIndsLoc, const ElInt* colIndsLoc, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG ASubHandle ) \
  { try { \
        auto A = Reinterpret(AHandle); \
        auto ASub = Reinterpret(ASubHandle); \
        const Int numRowInds = ASub->Height(); \
        const Int numColInds = ASub->Width(); \
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds), \
                         colIndVec(colIndsLoc,colIndsLoc+numColInds); \
        A->UpdateLocalSubmatrix \
        ( rowIndVec, colIndVec, Reinterpret(alpha), *ASub ); \
    } EL_CATCH; return EL_SUCCESS; } 

#define DISTMATRIX_SUBMATRIX_COMPLEX(SIG,SIGBASE,T) \
  /* DistMatrix<Base<T>,STAR,STAR> DistMatrix<T,U,V>::GetRealPartOfSubmatrix
     ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
  */ \
  ElError ElDistMatrixGetRealPartOfSubmatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, \
  ElInt numRowInds, const ElInt* rowInds, \
  ElInt numColInds, const ElInt* colInds, \
  ElDistMatrix_ ## SIGBASE *ASubHandle ) \
  { try { auto A = Reinterpret(AHandle); \
          auto ASub = new DistMatrix<Base<T>,STAR,STAR>(A->Grid()); \
          std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds); \
          std::vector<Int> colIndVec(colInds,colInds+numColInds); \
          A->GetRealPartOfSubmatrix( rowIndVec, colIndVec, *ASub ); \
          *ASubHandle = Reinterpret(ASub); } EL_CATCH; return EL_SUCCESS; } \
  /* DistMatrix<Base<T>,STAR,STAR> DistMatrix<T,U,V>::GetImagPartOfSubmatrix
     ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
  */ \
  ElError ElDistMatrixGetImagPartOfSubmatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, \
  ElInt numRowInds, const ElInt* rowInds, \
  ElInt numColInds, const ElInt* colInds, \
  ElDistMatrix_ ## SIGBASE *ASubHandle ) \
  { try { auto A = Reinterpret(AHandle); \
          auto ASub = new DistMatrix<Base<T>,STAR,STAR>(A->Grid()); \
          std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds); \
          std::vector<Int> colIndVec(colInds,colInds+numColInds); \
          A->GetImagPartOfSubmatrix( rowIndVec, colIndVec, *ASub ); \
          *ASubHandle = Reinterpret(ASub); } EL_CATCH; return EL_SUCCESS; } \
  /* void DistMatrix<T,U,V>::SetRealPartOfSubmatrix
     ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
       const DistMatrix<Base<T>,STAR,STAR>& ASub ); */ \
  ElError ElDistMatrixSetRealPartOfSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, const ElInt* rowInds, const ElInt* colInds, \
    ElConstDistMatrix_ ## SIGBASE ASubHandle ) \
  { try { \
      auto A = Reinterpret(AHandle); \
      auto ASubADM = Reinterpret(ASubHandle); \
      auto ASub = dynamic_cast<const DistMatrix<Base<T>,STAR,STAR>*>(ASubADM); \
      DynamicCastCheck(ASub); \
      const Int numRowInds = ASub->Height(); \
      const Int numColInds = ASub->Width(); \
      std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds), \
                       colIndVec(colInds,colInds+numColInds); \
      A->SetRealPartOfSubmatrix( rowIndVec, colIndVec, *ASub ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void DistMatrix<T,U,V>::SetImagPartOfSubmatrix
     ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
       const DistMatrix<Base<T>,STAR,STAR>& ASub ); */ \
  ElError ElDistMatrixSetImagPartOfSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, const ElInt* rowInds, const ElInt* colInds, \
    ElConstDistMatrix_ ## SIGBASE ASubHandle ) \
  { try { \
      auto A = Reinterpret(AHandle); \
      auto ASubADM = Reinterpret(ASubHandle); \
      auto ASub = dynamic_cast<const DistMatrix<Base<T>,STAR,STAR>*>(ASubADM); \
      DynamicCastCheck(ASub); \
      const Int numRowInds = ASub->Height(); \
      const Int numColInds = ASub->Width(); \
      std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds), \
                       colIndVec(colInds,colInds+numColInds); \
      A->SetImagPartOfSubmatrix( rowIndVec, colIndVec, *ASub ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void DistMatrix<T,U,V>::UpdateRealPartOfSubmatrix
     ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
       Base<T> alpha, const DistMatrix<Base<T>,STAR,STAR>& ASub ); */ \
  ElError ElDistMatrixUpdateRealPartOfSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, const ElInt* rowInds, const ElInt* colInds, \
    Base<T> alpha, ElConstDistMatrix_ ## SIGBASE ASubHandle ) \
  { try { \
      auto A = Reinterpret(AHandle); \
      auto ASubADM = Reinterpret(ASubHandle); \
      auto ASub = dynamic_cast<const DistMatrix<Base<T>,STAR,STAR>*>(ASubADM); \
      DynamicCastCheck(ASub); \
      const Int numRowInds = ASub->Height(); \
      const Int numColInds = ASub->Width(); \
      std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds), \
                       colIndVec(colInds,colInds+numColInds); \
      A->UpdateRealPartOfSubmatrix( rowIndVec, colIndVec, alpha, *ASub ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void DistMatrix<T,U,V>::UpdateImagPartOfSubmatrix
     ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
       Base<T> alpha, const DistMatrix<Base<T>,STAR,STAR>& ASub ); */ \
  ElError ElDistMatrixUpdateImagPartOfSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, const ElInt* rowInds, const ElInt* colInds, \
    Base<T> alpha, ElConstDistMatrix_ ## SIGBASE ASubHandle ) \
  { try { \
      auto A = Reinterpret(AHandle); \
      auto ASubADM = Reinterpret(ASubHandle); \
      auto ASub = dynamic_cast<const DistMatrix<Base<T>,STAR,STAR>*>(ASubADM); \
      DynamicCastCheck(ASub); \
      const Int numRowInds = ASub->Height(); \
      const Int numColInds = ASub->Width(); \
      std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds), \
                       colIndVec(colInds,colInds+numColInds); \
      A->UpdateImagPartOfSubmatrix( rowIndVec, colIndVec, alpha, *ASub ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void DistMatrix<T,U,V>::MakeSubmatrixReal
     ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) */ \
  ElError ElDistMatrixMakeSubmatrixReal_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt numRowInds, const ElInt* rowInds, \
                                  ElInt numColInds, const ElInt* colInds ) \
  { try { \
      std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds), \
                       colIndVec(colInds,colInds+numColInds); \
          Reinterpret(AHandle)->MakeSubmatrixReal( rowIndVec, colIndVec ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void DistMatrix<T,U,V>::ConjugateSubmatrix
     ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) */ \
  ElError ElDistMatrixConjugateSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, ElInt numRowInds, const ElInt* rowInds, \
                                  ElInt numColInds, const ElInt* colInds ) \
  { try { \
      std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds), \
                       colIndVec(colInds,colInds+numColInds); \
      Reinterpret(AHandle)->ConjugateSubmatrix( rowIndVec, colIndVec ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Matrix<Base<T>> DistMatrix<T,U,V>::GetRealPartOfLocalSubmatrix
     ( const std::vector<Int>& rowInds, 
       const std::vector<Int>& colInds ) const */ \
  ElError ElDistMatrixGetRealPartOfLocalSubmatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, \
    ElInt numRowInds, const ElInt* rowIndsLoc, \
    ElInt numColInds, const ElInt* colIndsLoc, \
    ElMatrix_ ## SIGBASE *ASubHandle ) \
  { try { \
        auto A = Reinterpret(AHandle); \
        auto ASub = new Matrix<Base<T>>; \
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds), \
                         colIndVec(colIndsLoc,colIndsLoc+numColInds); \
        A->GetRealPartOfLocalSubmatrix( rowIndVec, colIndVec, *ASub ); \
        *ASubHandle = Reinterpret(ASub); } EL_CATCH; return EL_SUCCESS; } \
  /* Matrix<Base<T>> DistMatrix<T,U,V>::GetImagPartOfLocalSubmatrix
     ( const std::vector<Int>& rowInds, 
       const std::vector<Int>& colInds ) const */ \
  ElError ElDistMatrixGetImagPartOfLocalSubmatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, \
    ElInt numRowInds, const ElInt* rowIndsLoc, \
    ElInt numColInds, const ElInt* colIndsLoc, \
    ElMatrix_ ## SIGBASE *ASubHandle ) \
  { try { \
        auto A = Reinterpret(AHandle); \
        auto ASub = new Matrix<Base<T>>; \
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds), \
                         colIndVec(colIndsLoc,colIndsLoc+numColInds); \
        A->GetImagPartOfLocalSubmatrix( rowIndVec, colIndVec, *ASub ); \
        *ASubHandle = Reinterpret(ASub); } EL_CATCH; return EL_SUCCESS; } \
  /* void DistMatrix<T,U,V>::SetRealPartOfLocalSubmatrix
     ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
       const Matrix<Base<T>>& ASub ); */ \
  ElError ElDistMatrixSetRealPartOfLocalSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, \
    const ElInt* rowIndsLoc, const ElInt* colIndsLoc, \
    ElConstMatrix_ ## SIGBASE ASubHandle ) \
  { try { \
        auto A = Reinterpret(AHandle); \
        auto ASub = Reinterpret(ASubHandle); \
        const Int numRowInds = ASub->Height(); \
        const Int numColInds = ASub->Width(); \
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds), \
                         colIndVec(colIndsLoc,colIndsLoc+numColInds); \
        A->SetRealPartOfLocalSubmatrix( rowIndVec, colIndVec, *ASub ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void DistMatrix<T,U,V>::SetImagPartOfLocalSubmatrix
     ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
       const Matrix<Base<T>>& ASub ); */ \
  ElError ElDistMatrixSetImagPartOfLocalSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, \
    const ElInt* rowIndsLoc, const ElInt* colIndsLoc, \
    ElConstMatrix_ ## SIGBASE ASubHandle ) \
  { try { \
        auto A = Reinterpret(AHandle); \
        auto ASub = Reinterpret(ASubHandle); \
        const Int numRowInds = ASub->Height(); \
        const Int numColInds = ASub->Width(); \
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds), \
                         colIndVec(colIndsLoc,colIndsLoc+numColInds); \
        A->SetImagPartOfLocalSubmatrix( rowIndVec, colIndVec, *ASub ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void DistMatrix<T,U,V>::UpdateRealPartOfLocalSubmatrix
     ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
       Base<T> alpha, const Matrix<Base<T>>& ASub ); */ \
  ElError ElDistMatrixUpdateRealPartOfLocalSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, \
    const ElInt* rowIndsLoc, const ElInt* colIndsLoc, \
    Base<T> alpha, ElConstMatrix_ ## SIGBASE ASubHandle ) \
  { try { \
        auto A = Reinterpret(AHandle); \
        auto ASub = Reinterpret(ASubHandle); \
        const Int numRowInds = ASub->Height(); \
        const Int numColInds = ASub->Width(); \
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds), \
                         colIndVec(colIndsLoc,colIndsLoc+numColInds); \
        A->UpdateRealPartOfLocalSubmatrix \
        ( rowIndVec, colIndVec, alpha, *ASub ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void DistMatrix<T,U,V>::UpdateImagPartOfLocalSubmatrix
     ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
       Base<T> alpha, const Matrix<Base<T>>& ASub ); */ \
  ElError ElDistMatrixUpdateImagPartOfLocalSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, \
    const ElInt* rowIndsLoc, const ElInt* colIndsLoc, \
    Base<T> alpha, ElConstMatrix_ ## SIGBASE ASubHandle ) \
  { try { \
        auto A = Reinterpret(AHandle); \
        auto ASub = Reinterpret(ASubHandle); \
        const Int numRowInds = ASub->Height(); \
        const Int numColInds = ASub->Width(); \
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds), \
                         colIndVec(colIndsLoc,colIndsLoc+numColInds); \
        A->UpdateImagPartOfLocalSubmatrix \
        ( rowIndVec, colIndVec, alpha, *ASub ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void DistMatrix<T,U,V>::MakeLocalSubmatrixReal
     ( const std::vector<Int>& rowIndsLoc, 
       const std::vector<Int>& colIndsLoc ) */ \
  ElError ElDistMatrixMakeLocalSubmatrixReal_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, \
    ElInt numRowInds, const ElInt* rowIndsLoc, \
    ElInt numColInds, const ElInt* colIndsLoc ) \
  { try { \
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds), \
                         colIndVec(colIndsLoc,colIndsLoc+numColInds); \
        Reinterpret(AHandle)->MakeLocalSubmatrixReal( rowIndVec, colIndVec ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void DistMatrix<T,U,V>::ConjugateLocalSubmatrix
     ( const std::vector<Int>& rowIndsLoc, 
       const std::vector<Int>& colIndsLoc ) */ \
  ElError ElDistMatrixConjugateLocalSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG AHandle, \
    ElInt numRowInds, const ElInt* rowIndsLoc, \
    ElInt numColInds, const ElInt* colIndsLoc ) \
  { try { \
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds), \
                         colIndVec(colIndsLoc,colIndsLoc+numColInds); \
        Reinterpret(AHandle)->ConjugateLocalSubmatrix( rowIndVec, colIndVec ); \
    } EL_CATCH; return EL_SUCCESS; }

#define C_PROTO(SIG,T) \
  DISTMATRIX_CREATE(SIG,T) \
  DISTMATRIX_RECONFIG(SIG,T) \
  DISTMATRIX_BASIC(SIG,T) \
  DISTMATRIX_SINGLEENTRY(SIG,T) \
  DISTMATRIX_DIAGONAL(SIG,T) \
  DISTMATRIX_SUBMATRIX(SIG,T)

#define C_PROTO_COMPLEX(SIG,SIGBASE,T) \
  C_PROTO(SIG,T) \
  DISTMATRIX_SINGLEENTRY_COMPLEX(SIG,SIGBASE,T) \
  DISTMATRIX_DIAGONAL_COMPLEX(SIG,SIGBASE,T) \
  DISTMATRIX_SUBMATRIX_COMPLEX(SIG,SIGBASE,T)

#include "El/macros/CInstantiate.h"

// TODO: More diagonal manipulation
// ================================

} // extern "C"
