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
#include "El.h"
using namespace El;

template<typename T>
ElError ElDistMatrixCreateSpecific
( ElDist U_C, ElDist V_C, ElConstGrid grid, AbstractDistMatrix<T>** A )
{
    try 
    {
        Dist U = CReflect(U_C);
        Dist V = CReflect(V_C);
        auto gridPtr = CReflect(grid);

        #define GUARD(CDIST,RDIST) U == CDIST && V == RDIST
        #define PAYLOAD(CDIST,RDIST) \
          *A = new DistMatrix<T,CDIST,RDIST>(*gridPtr);
        #include "El/macros/GuardAndPayload.h"
        #undef GUARD
        #undef PAYLOAD
    }
    EL_CATCH
    return EL_SUCCESS;
}

extern "C" {

#define DISTMATRIX_CREATE(SIG,SIGBASE,T) \
  /* DistMatrix<T,MC,MR>::DistMatrix( const Grid& g ) */ \
  ElError ElDistMatrixCreate_ ## SIG \
  ( ElConstGrid grid, ElDistMatrix_ ## SIG *A ) \
  { EL_TRY( auto gridPtr = CReflect(grid); \
            *A = CReflect( new DistMatrix<T>(*gridPtr) ) ) } \
  /* DistMatrix<T,U,V>::DistMatrix( const Grid& g ) */ \
  ElError ElDistMatrixCreateSpecific_ ## SIG \
  ( ElDist U, ElDist V, ElConstGrid grid, \
    ElDistMatrix_ ## SIG *A ) \
  { \
    AbstractDistMatrix<T>* ADM; \
    ElError error = ElDistMatrixCreateSpecific( U, V, grid, &ADM ); \
    *A = CReflect(ADM); \
    return error; \
  } \
  /* DistMatrix<T,U,V>::~DistMatrix() */ \
  ElError ElDistMatrixDestroy_ ## SIG ( ElConstDistMatrix_ ## SIG A ) \
  { EL_TRY( CReflect(A) ) } 

#define DISTMATRIX_RECONFIG(SIG,SIGBASE,T) \
  /* void DistMatrix<T,U,V>::Empty() */ \
  ElError ElDistMatrixEmpty_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( CReflect(A)->Empty() ) } \
  /* void DistMatrix<T,U,V>::EmptyData() */ \
  ElError ElDistMatrixEmptyData_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( CReflect(A)->EmptyData() ) } \
  /* void DistMatrix<T,U,V>::SetGrid( const Grid& g ) */ \
  ElError ElDistMatrixSetGrid_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElConstGrid grid ) \
  { EL_TRY( CReflect(A)->SetGrid(*CReflect(grid)) ) } \
  /* void DistMatrix<T,U,V>::Resize( Int height, Int width ) */ \
  ElError ElDistMatrixResize_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt height, ElInt width ) \
  { EL_TRY( CReflect(A)->Resize(height,width) ) } \
  /* void DistMatrix<T,U,V>::Resize( Int height, Int width, Int ldim ) */ \
  ElError ElDistMatrixResizeWithLDim_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt height, ElInt width, ElInt ldim ) \
  { EL_TRY( CReflect(A)->Resize(height,width,ldim) ) } \
  /* void DistMatrix<T,U,V>::MakeConsistent() */ \
  ElError ElDistMatrixMakeConsistent_ ## SIG \
  ( ElDistMatrix_ ## SIG A, bool includeViewers ) \
  { EL_TRY( CReflect(A)->MakeConsistent(includeViewers) ) } \
  /* void DistMatrix<T,U,V>::MakeSizeConsistent() */ \
  ElError ElDistMatrixMakeSizeConsistent_ ## SIG \
  ( ElDistMatrix_ ## SIG A, bool includeViewers ) \
  { EL_TRY( CReflect(A)->MakeSizeConsistent(includeViewers) ) } \
  /* void DistMatrix<T,U,V>::Align \
     ( Int colAlign, Int rowAlign, bool constrain ) */ \
  ElError ElDistMatrixAlign_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt colAlign, ElInt rowAlign, \
    bool constrain ) \
  { EL_TRY( CReflect(A)->Align(colAlign,rowAlign,constrain) ) } \
  /* void DistMatrix<T,U,V>::AlignCols( Int colAlign, bool constrain ) */ \
  ElError ElDistMatrixAlignCols_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt colAlign, bool constrain ) \
  { EL_TRY( CReflect(A)->AlignCols(colAlign,constrain) ) } \
   /* void DistMatrix<T,U,V>::AlignRows( Int rowAlign, bool constrain ) */ \
  ElError ElDistMatrixAlignRows_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt rowAlign, bool constrain ) \
  { EL_TRY( CReflect(A)->AlignRows(rowAlign,constrain) ) } \
  /* void DistMatrix<T,U,V>::FreeAlignments() */ \
  ElError ElDistMatrixFreeAlignments_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( CReflect(A)->FreeAlignments() ) } \
  /* void DistMatrix<T,U,V>::SetRoot( Int root, bool constrain ) */ \
  ElError ElDistMatrixSetRoot_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt root, bool constrain ) \
  { EL_TRY( CReflect(A)->SetRoot(root,constrain) ) } \
  /* void DistMatrix<T,U,V>::AlignWith \
     ( const DistData& data, bool constrain ) */ \
  ElError ElDistMatrixAlignWith_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistData distData, bool constrain ) \
  { EL_TRY( CReflect(A)->AlignWith \
            ( CReflect(distData), constrain ) ) } \
  /* void DistMatrix<T,U,V>::AlignColsWith \
     ( const DistData& data, bool constrain ) */ \
  ElError ElDistMatrixAlignColsWith_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistData distData, bool constrain ) \
  { EL_TRY( CReflect(A)->AlignColsWith \
            ( CReflect(distData), constrain ) ) } \
  /* void DistMatrix<T,U,V>::AlignRowsWith \
     ( const DistData& data, bool constrain ) */ \
  ElError ElDistMatrixAlignRowsWith_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistData distData, bool constrain ) \
  { EL_TRY( CReflect(A)->AlignRowsWith \
            ( CReflect(distData), constrain ) ) } \
  /* void DistMatrix<T,U,V>::AlignAndResize \
     ( Int colAlign, Int rowAlign, Int height, Int width, \
       bool force, bool constrain ) */ \
  ElError ElDistMatrixAlignAndResize_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElInt colAlign, ElInt rowAlign, ElInt height, ElInt width, \
    bool force, bool constrain ) \
  { EL_TRY( CReflect(A)->AlignAndResize \
            (colAlign,rowAlign,height,width,force,constrain) ) } \
  /* void DistMatrix<T,U,V>::AlignColsAndResize
     ( Int colAlign, Int height, Int width, bool force, bool constrain ) */ \
  ElError ElDistMatrixAlignColsAndResize_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElInt colAlign, ElInt height, ElInt width, bool force, bool constrain ) \
  { EL_TRY( CReflect(A)->AlignColsAndResize \
            (colAlign,height,width,force,constrain) ) } \
  /* void DistMatrix<T,U,V>::AlignRowsAndResize
     ( Int rowAlign, Int height, Int width, bool force, bool constrain ) */ \
  ElError ElDistMatrixAlignRowsAndResize_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElInt rowAlign, ElInt height, ElInt width, bool force, bool constrain ) \
  { EL_TRY( CReflect(A)->AlignRowsAndResize \
            (rowAlign,height,width,force,constrain) ) } \
  /* void DistMatrix<T,U,V>::Attach
     ( Int height, Int width, const Grid& grid, Int colAlign, Int rowAlign, 
       T* buffer, Int ldim, Int root ) */ \
  ElError ElDistMatrixAttach_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt height, ElInt width, \
    ElConstGrid grid, ElInt colAlign, ElInt rowAlign, \
    CREFLECT(T)* buffer, ElInt ldim, ElInt root ) \
  { EL_TRY( CReflect(A)->Attach \
            (height,width,*CReflect(grid),colAlign,rowAlign, \
             CReflect(buffer),ldim,root) ) } \
  /* void DistMatrix<T,U,V>::LockedAttach
     ( Int height, Int width, const Grid& grid, Int colAlign, Int rowAlign, 
       const T* buffer, Int ldim, Int root ) */ \
  ElError ElDistMatrixLockedAttach_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt height, ElInt width, \
    ElConstGrid grid, ElInt colAlign, ElInt rowAlign, \
    const CREFLECT(T)* buffer, ElInt ldim, ElInt root ) \
  { EL_TRY( CReflect(A)->LockedAttach \
            (height,width,*CReflect(grid),colAlign,rowAlign, \
             CReflect(buffer),ldim,root) ) }

#define DISTMATRIX_BASIC(SIG,SIGBASE,T) \
  /* Int DistMatrix<T,U,V>::Height() const */ \
  ElError ElDistMatrixHeight_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* height ) \
  { EL_TRY( *height = CReflect(A)->Height() ) } \
  /* Int DistMatrix<T,U,V>::Width() const */ \
  ElError ElDistMatrixWidth_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* width ) \
  { EL_TRY( *width = CReflect(A)->Width() ) } \
  /* Int DistMatrix<T,U,V>::DiagonalLength( Int offset ) const */ \
  ElError ElDistMatrixDiagonalLength_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt offset, ElInt* length ) \
  { EL_TRY( *length = CReflect(A)->DiagonalLength(offset) ) } \
  /* bool DistMatrix<T,U,V>::Viewing() const */ \
  ElError ElDistMatrixViewing_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, bool* viewing ) \
  { EL_TRY( *viewing = CReflect(A)->Viewing() ) } \
  /* bool DistMatrix<T,U,V>::Locked() const */ \
  ElError ElDistMatrixLocked_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, bool* locked ) \
  { EL_TRY( *locked = CReflect(A)->Locked() ) } \
  /* Int DistMatrix<T,U,V>::LocalHeight() const */ \
  ElError ElDistMatrixLocalHeight_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* localHeight ) \
  { EL_TRY( *localHeight = CReflect(A)->LocalHeight() ) } \
  /* Int DistMatrix<T,U,V>::LocalWidth() const */ \
  ElError ElDistMatrixLocalWidth_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* localWidth ) \
  { EL_TRY( *localWidth = CReflect(A)->LocalWidth() ) } \
  /* Int DistMatrix<T,U,V>::LDim() const */ \
  ElError ElDistMatrixLDim_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* ldim ) \
  { EL_TRY( *ldim = CReflect(A)->LDim() ) } \
  /* Matrix<T>& DistMatrix<T,U,V>::Matrix() */ \
  ElError ElDistMatrixMatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElMatrix_ ## SIG *ALoc ) \
  { EL_TRY( *ALoc = CReflect(&CReflect(A)->Matrix()) ) } \
  /* const Matrix<T>& DistMatrix<T,U,V>::LockedMatrix() const */ \
  ElError ElDistMatrixLockedMatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstMatrix_ ## SIG *ALoc ) \
  { EL_TRY( *ALoc = CReflect(&CReflect(A)->LockedMatrix()) ) } \
  /* size_t DistMatrix<T,U,V>::AllocatedMemory() const */ \
  ElError ElDistMatrixAllocatedMemory_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, size_t* mem ) \
  { EL_TRY( *mem = CReflect(A)->AllocatedMemory() ) } \
  /* T* DistMatrix<T,U,V>::Buffer() */ \
  ElError ElDistMatrixBuffer_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T)** buffer ) \
  { EL_TRY( *buffer = CReflect(CReflect(A)->Buffer()) ) } \
  /* const T* DistMatrix<T,U,V>::LockedBuffer() const */ \
  ElError ElDistMatrixLockedBuffer_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, const CREFLECT(T)** buffer ) \
  { EL_TRY( *buffer = CReflect(CReflect(A)->LockedBuffer()) ) } \
  /* const Grid& DistMatrix<T,U,V>::Grid() const */ \
  ElError ElDistMatrixGrid_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstGrid* grid ) \
  { EL_TRY( *grid = CReflect(&CReflect(A)->Grid()) ) } \
  /* bool DistMatrix<T,U,V>::ColConstrained() const */ \
  ElError ElDistMatrixColConstrained_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, bool* colConst ) \
  { EL_TRY( *colConst = CReflect(A)->ColConstrained() ) } \
  /* bool DistMatrix<T,U,V>::RowConstrained() const */ \
  ElError ElDistMatrixRowConstrained_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, bool* rowConst ) \
  { EL_TRY( *rowConst = CReflect(A)->RowConstrained() ) } \
  /* bool DistMatrix<T,U,V>::RootConstrained() const */ \
  ElError ElDistMatrixRootConstrained_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, bool* rootConst ) \
  { EL_TRY( *rootConst = CReflect(A)->RootConstrained() ) } \
  /* Int DistMatrix<T,U,V>::ColAlign() const */ \
  ElError ElDistMatrixColAlign_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* colAlign ) \
  { EL_TRY( *colAlign = CReflect(A)->ColAlign() ) } \
  /* Int DistMatrix<T,U,V>::RowAlign() const */ \
  ElError ElDistMatrixRowAlign_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* rowAlign ) \
  { EL_TRY( *rowAlign = CReflect(A)->RowAlign() ) } \
  /* Int DistMatrix<T,U,V>::ColShift() const */ \
  ElError ElDistMatrixColShift_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* colShift ) \
  { EL_TRY( *colShift = CReflect(A)->ColShift() ) } \
  /* Int DistMatrix<T,U,V>::RowShift() const */ \
  ElError ElDistMatrixRowShift_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* rowShift ) \
  { EL_TRY( *rowShift = CReflect(A)->RowShift() ) } \
  /* Int DistMatrix<T,U,V>::ColRank() const */ \
  ElError ElDistMatrixColRank_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* colRank ) \
  { EL_TRY( *colRank = CReflect(A)->ColRank() ) } \
  /* Int DistMatrix<T,U,V>::RowRank() const */ \
  ElError ElDistMatrixRowRank_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* rowRank ) \
  { EL_TRY( *rowRank = CReflect(A)->RowRank() ) } \
  /* Int DistMatrix<T,U,V>::PartialColRank() const */ \
  ElError ElDistMatrixPartialColRank_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* rank ) \
  { EL_TRY( *rank = CReflect(A)->PartialColRank() ) } \
  /* Int DistMatrix<T,U,V>::PartialRowRank() const */ \
  ElError ElDistMatrixPartialRowRank_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* rank ) \
  { EL_TRY( *rank = CReflect(A)->PartialRowRank() ) } \
  /* Int DistMatrix<T,U,V>::PartialUnionColRank() const */ \
  ElError ElDistMatrixPartialUnionColRank_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* rank ) \
  { EL_TRY( *rank = CReflect(A)->PartialUnionColRank() ) } \
  /* Int DistMatrix<T,U,V>::PartialUnionRowRank() const */ \
  ElError ElDistMatrixPartialUnionRowRank_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* rank ) \
  { EL_TRY( *rank = CReflect(A)->PartialUnionRowRank() ) } \
  /* Int DistMatrix<T,U,V>::DistRank() const */ \
  ElError ElDistMatrixDistRank_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* rank ) \
  { EL_TRY( *rank = CReflect(A)->DistRank() ) } \
  /* Int DistMatrix<T,U,V>::CrossRank() const */ \
  ElError ElDistMatrixCrossRank_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* rank ) \
  { EL_TRY( *rank = CReflect(A)->CrossRank() ) } \
  /* Int DistMatrix<T,U,V>::RedundantRank() const */ \
  ElError ElDistMatrixRedundantRank_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* rank ) \
  { EL_TRY( *rank = CReflect(A)->RedundantRank() ) } \
  /* Int DistMatrix<T,U,V>::Root() const */ \
  ElError ElDistMatrixRoot_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* root ) \
  { EL_TRY( *root = CReflect(A)->Root() ) } \
  /* bool DistMatrix<T,U,V>::Participating() const */ \
  ElError ElDistMatrixParticipating_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, bool* participating ) \
  { EL_TRY( *participating = CReflect(A)->Participating() ) } \
  /* Int DistMatrix<T,U,V>::RowOwner( Int i ) const */ \
  ElError ElDistMatrixRowOwner_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt i, ElInt* rowOwner ) \
  { EL_TRY( *rowOwner = CReflect(A)->RowOwner(i) ) } \
  /* Int DistMatrix<T,U,V>::ColOwner( Int i ) const */ \
  ElError ElDistMatrixColOwner_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt j, ElInt* colOwner ) \
  { EL_TRY( *colOwner = CReflect(A)->ColOwner(j) ) } \
  /* Int DistMatrix<T,U,V>::Owner( Int i, Int j ) const */ \
  ElError ElDistMatrixOwner_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt i, ElInt j, ElInt* owner ) \
  { EL_TRY( *owner = CReflect(A)->Owner(i,j) ) } \
  /* Int DistMatrix<T,U,V>::LocalRow( Int i ) const */ \
  ElError ElDistMatrixLocalRow_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt i, ElInt* iLoc ) \
  { EL_TRY( *iLoc = CReflect(A)->LocalRow(i) ) } \
  /* Int DistMatrix<T,U,V>::LocalCol( Int j ) const */ \
  ElError ElDistMatrixLocalCol_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt j, ElInt* jLoc ) \
  { EL_TRY( *jLoc = CReflect(A)->LocalCol(j) ) } \
  /* Int DistMatrix<T,U,V>::LocalRowOffset( Int i ) const */ \
  ElError ElDistMatrixLocalRowOffset_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt i, ElInt* iLoc ) \
  { EL_TRY( *iLoc = CReflect(A)->LocalRowOffset(i) ) } \
  /* Int DistMatrix<T,U,V>::LocalColOffset( Int j ) const */ \
  ElError ElDistMatrixLocalColOffset_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt j, ElInt* jLoc ) \
  { EL_TRY( *jLoc = CReflect(A)->LocalColOffset(j) ) } \
  /* Int DistMatrix<T,U,V>::GlobalRow( Int iLoc ) const */ \
  ElError ElDistMatrixGlobalRow_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt iLoc, ElInt* i ) \
  { EL_TRY( *i = CReflect(A)->GlobalRow(iLoc) ) } \
  /* Int DistMatrix<T,U,V>::GlobalCol( Int j ) const */ \
  ElError ElDistMatrixGlobalCol_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt jLoc, ElInt* j ) \
  { EL_TRY( *j = CReflect(A)->GlobalCol(jLoc) ) } \
  /* bool DistMatrix<T,U,V>::IsLocalRow( Int i ) const */ \
  ElError ElDistMatrixIsLocalRow_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt i, bool* isLocal ) \
  { EL_TRY( *isLocal = CReflect(A)->IsLocalRow(i) ) } \
  /* bool DistMatrix<T,U,V>::IsLocalCol( Int j ) const */ \
  ElError ElDistMatrixIsLocalCol_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt j, bool* isLocal ) \
  { EL_TRY( *isLocal = CReflect(A)->IsLocalCol(j) ) } \
  /* bool DistMatrix<T,U,V>::IsLocal( Int i, Int j ) const */ \
  ElError ElDistMatrixIsLocal_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt i, ElInt j, bool* isLocal ) \
  { EL_TRY( *isLocal = CReflect(A)->IsLocal(i,j) ) } \
  /* DistData DistMatrix<T,U,V>::DistData() const */ \
  ElError ElDistMatrixDistData_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistData* distData ) \
  { EL_TRY( DistData data = CReflect(A)->DistData(); \
            *distData = CReflect( data ) ) } \
  /* mpi::Comm DistMatrix<T,U,V>::DistComm() const */ \
  ElError ElDistMatrixDistComm_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, MPI_Comm* comm ) \
  { EL_TRY( *comm = CReflect(A)->DistComm().comm ) } \
  /* mpi::Comm DistMatrix<T,U,V>::CrossComm() const */ \
  ElError ElDistMatrixCrossComm_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, MPI_Comm* comm ) \
  { EL_TRY( *comm = CReflect(A)->CrossComm().comm ) } \
  /* mpi::Comm DistMatrix<T,U,V>::RedundantComm() const */ \
  ElError ElDistMatrixRedundantComm_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, MPI_Comm* comm ) \
  { EL_TRY( *comm = CReflect(A)->RedundantComm().comm ) } \
  /* mpi::Comm DistMatrix<T,U,V>::ColComm() const */ \
  ElError ElDistMatrixColComm_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, MPI_Comm* comm ) \
  { EL_TRY( *comm = CReflect(A)->ColComm().comm ) } \
  /* mpi::Comm DistMatrix<T,U,V>::RowComm() const */ \
  ElError ElDistMatrixRowComm_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, MPI_Comm* comm ) \
  { EL_TRY( *comm = CReflect(A)->RowComm().comm ) } \
  /* mpi::Comm DistMatrix<T,U,V>::PartialColComm() const */ \
  ElError ElDistMatrixPartialColComm_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, MPI_Comm* comm ) \
  { EL_TRY( *comm = CReflect(A)->PartialColComm().comm ) } \
  /* mpi::Comm DistMatrix<T,U,V>::PartialRowComm() const */ \
  ElError ElDistMatrixPartialRowComm_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, MPI_Comm* comm ) \
  { EL_TRY( *comm = CReflect(A)->PartialRowComm().comm ) } \
  /* mpi::Comm DistMatrix<T,U,V>::PartialUnionColComm() const */ \
  ElError ElDistMatrixPartialUnionColComm_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, MPI_Comm* comm ) \
  { EL_TRY( *comm = CReflect(A)->PartialUnionColComm().comm ) } \
  /* mpi::Comm DistMatrix<T,U,V>::PartialUnionRowComm() const */ \
  ElError ElDistMatrixPartialUnionRowComm_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, MPI_Comm* comm ) \
  { EL_TRY( *comm = CReflect(A)->PartialUnionRowComm().comm ) } \
  /* Int DistMatrix<T,U,V>::ColStride() const */ \
  ElError ElDistMatrixColStride_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* stride ) \
  { EL_TRY( *stride = CReflect(A)->ColStride() ) } \
  /* Int DistMatrix<T,U,V>::RowStride() const */ \
  ElError ElDistMatrixRowStride_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* stride ) \
  { EL_TRY( *stride = CReflect(A)->RowStride() ) } \
  /* Int DistMatrix<T,U,V>::PartialColStride() const */ \
  ElError ElDistMatrixPartialColStride_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* stride ) \
  { EL_TRY( *stride = CReflect(A)->PartialColStride() ) } \
  /* Int DistMatrix<T,U,V>::PartialRowStride() const */ \
  ElError ElDistMatrixPartialRowStride_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* stride ) \
  { EL_TRY( *stride = CReflect(A)->PartialRowStride() ) } \
  /* Int DistMatrix<T,U,V>::PartialUnionColStride() const */ \
  ElError ElDistMatrixPartialUnionColStride_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* stride ) \
  { EL_TRY( *stride = CReflect(A)->PartialUnionColStride() ) } \
  /* Int DistMatrix<T,U,V>::PartialUnionRowStride() const */ \
  ElError ElDistMatrixPartialUnionRowStride_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* stride ) \
  { EL_TRY( *stride = CReflect(A)->PartialUnionRowStride() ) } \
  /* Int DistMatrix<T,U,V>::DistSize() const */ \
  ElError ElDistMatrixDistSize_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* commSize ) \
  { EL_TRY( *commSize = CReflect(A)->DistSize() ) } \
  /* Int DistMatrix<T,U,V>::CrossSize() const */ \
  ElError ElDistMatrixCrossSize_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* commSize ) \
  { EL_TRY( *commSize = CReflect(A)->CrossSize() ) } \
  /* Int DistMatrix<T,U,V>::RedundantSize() const */ \
  ElError ElDistMatrixRedundantSize_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt* commSize ) \
  { EL_TRY( *commSize = CReflect(A)->RedundantSize() ) }

#define DISTMATRIX_SINGLEENTRY(SIG,SIGBASE,T) \
  /* T DistMatrix<T,U,V>::Get( Int i, Int j ) const */ \
  ElError ElDistMatrixGet_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt i, ElInt j, CREFLECT(T)* val ) \
  { EL_TRY( *val = CReflect(CReflect(A)->Get(i,j)) ) } \
  /* void DistMatrix<T,U,V>::Set( Int i, Int j, T alpha ) */ \
  ElError ElDistMatrixSet_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt i, ElInt j, CREFLECT(T) alpha ) \
  { EL_TRY( CReflect(A)->Set(i,j,CReflect(alpha)) ) } \
  /* void DistMatrix<T,U,V>::Update( Int i, Int j, T alpha ) */ \
  ElError ElDistMatrixUpdate_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt i, ElInt j, CREFLECT(T) alpha ) \
  { EL_TRY( CReflect(A)->Update(i,j,CReflect(alpha)) ) } \
  /* T DistMatrix<T,U,V>::GetLocal( Int iLoc, Int jLoc ) const */ \
  ElError ElDistMatrixGetLocal_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt iLoc, ElInt jLoc, \
    CREFLECT(T)* val ) \
  { EL_TRY( *val = CReflect(CReflect(A)->GetLocal(iLoc,jLoc)) ) } \
  /* void DistMatrix<T,U,V>::SetLocal( Int iLoc, Int jLoc, T alpha ) */ \
  ElError ElDistMatrixSetLocal_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt iLoc, ElInt jLoc, CREFLECT(T) alpha ) \
  { EL_TRY( CReflect(A)->SetLocal(iLoc,jLoc,CReflect(alpha)) ) } \
  /* void DistMatrix<T,U,V>::UpdateLocal( Int iLoc, Int jLoc, T alpha ) */ \
  ElError ElDistMatrixUpdateLocal_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt iLoc, ElInt jLoc, CREFLECT(T) alpha ) \
  { EL_TRY( CReflect(A)->UpdateLocal \
            (iLoc,jLoc,CReflect(alpha)) ) }

#define DISTMATRIX_SINGLEENTRY_COMPLEX(SIG,SIGBASE,T) \
  /* Base<T> DistMatrix<T,U,V>::GetRealPart( Int i, Int j ) const */ \
  ElError ElDistMatrixGetRealPart_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt i, ElInt j, Base<T>* val ) \
  { EL_TRY( *val = CReflect(A)->GetRealPart(i,j) ) } \
  /* Base<T> DistMatrix<T,U,V>::GetImagPart( Int i, Int j ) const */ \
  ElError ElDistMatrixGetImagPart_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt i, ElInt j, Base<T>* val ) \
  { EL_TRY( *val = CReflect(A)->GetImagPart(i,j) ) } \
  /* void DistMatrix<T,U,V>::SetRealPart( Int i, Int j, Base<T> alpha ) */ \
  ElError ElDistMatrixSetRealPart_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt i, ElInt j, Base<T> alpha ) \
  { EL_TRY( CReflect(A)->SetRealPart(i,j,alpha) ) } \
  /* void DistMatrix<T,U,V>::SetImagPart( Int i, Int j, Base<T> alpha ) */ \
  ElError ElDistMatrixSetImagPart_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt i, ElInt j, Base<T> alpha ) \
  { EL_TRY( CReflect(A)->SetImagPart(i,j,alpha) ) } \
  /* void DistMatrix<T,U,V>::UpdateRealPart( Int i, Int j, Base<T> alpha ) */ \
  ElError ElDistMatrixUpdateRealPart_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt i, ElInt j, Base<T> alpha ) \
  { EL_TRY( CReflect(A)->UpdateRealPart(i,j,alpha) ) } \
  /* void DistMatrix<T,U,V>::UpdateImagPart( Int i, Int j, Base<T> alpha ) */ \
  ElError ElDistMatrixUpdateImagPart_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt i, ElInt j, Base<T> alpha ) \
  { EL_TRY( CReflect(A)->UpdateImagPart(i,j,alpha) ) } \
  /* void DistMatrix<T,U,V>::MakeReal( Int i, Int j ) */ \
  ElError ElDistMatrixMakeReal_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt i, ElInt j ) \
  { EL_TRY( CReflect(A)->MakeReal(i,j) ) } \
  /* void DistMatrix<T,U,V>::Conjugate( Int i, Int j ) */ \
  ElError ElDistMatrixConjugate_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt i, ElInt j ) \
  { EL_TRY( CReflect(A)->Conjugate(i,j) ) } \
  /* Base<T> DistMatrix<T,U,V>::GetLocalRealPart \
     ( Int iLoc, Int jLoc ) const */ \
  ElError ElDistMatrixGetLocalRealPart_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt iLoc, ElInt jLoc, Base<T>* val ) \
  { EL_TRY( *val = CReflect(A)->GetLocalRealPart(iLoc,jLoc) ) } \
  /* Base<T> DistMatrix<T,U,V>::GetLocalImagPart \
     ( Int iLoc, Int jLoc ) const */ \
  ElError ElDistMatrixGetLocalImagPart_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt iLoc, ElInt jLoc, Base<T>* val ) \
  { EL_TRY( *val = CReflect(A)->GetLocalImagPart(iLoc,jLoc) ) } \
  /* void DistMatrix<T,U,V>::SetLocalRealPart
     ( Int iLoc, Int jLoc, Base<T> alpha ) */ \
  ElError ElDistMatrixSetLocalRealPart_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt iLoc, ElInt jLoc, Base<T> alpha ) \
  { EL_TRY( CReflect(A)->SetLocalRealPart(iLoc,jLoc,alpha) ) } \
  /* void DistMatrix<T,U,V>::SetLocalImagPart
     ( Int iLoc, Int jLoc, Base<T> alpha ) */ \
  ElError ElDistMatrixSetLocalImagPart_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt iLoc, ElInt jLoc, Base<T> alpha ) \
  { EL_TRY( CReflect(A)->SetLocalImagPart(iLoc,jLoc,alpha) ) } \
  /* void DistMatrix<T,U,V>::UpdateLocalRealPart
     ( Int iLoc, Int jLoc, Base<T> alpha ) */ \
  ElError ElDistMatrixUpdateLocalRealPart_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt iLoc, ElInt jLoc, Base<T> alpha ) \
  { EL_TRY( CReflect(A)->UpdateLocalRealPart(iLoc,jLoc,alpha) ) } \
  /* void DistMatrix<T,U,V>::UpdateLocalImagPart
     ( Int iLoc, Int jLoc, Base<T> alpha ) */ \
  ElError ElDistMatrixUpdateLocalImagPart_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt iLoc, ElInt jLoc, Base<T> alpha ) \
  { EL_TRY( CReflect(A)->UpdateLocalImagPart(iLoc,jLoc,alpha) ) }

#define DISTMATRIX_DIAGONAL(SIG,SIGBASE,T) \
  /* bool AbstractDistMatrix<T>::DiagonalAlignedWith
     ( const DistData& data, Int offset ) const */ \
  ElError ElDistMatrixDiagonalAlignedWith_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistData distData, ElInt offset, \
    bool* aligned ) \
  { EL_TRY( *aligned = CReflect(A)->DiagonalAlignedWith \
                       (CReflect(distData),offset) ) } \
  /* Int AbstractDistMatrix<T>::DiagonalRoot( Int offset ) const */ \
  ElError ElDistMatrixDiagonalRoot_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt offset, ElInt* root ) \
  { EL_TRY( *root = CReflect(A)->DiagonalRoot(offset) ) } \
  /* Int AbstractDistMatrix<T>::DiagonalAlign( Int offset ) const */ \
  ElError ElDistMatrixDiagonalAlign_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt offset, ElInt* align ) \
  { EL_TRY( *align = CReflect(A)->DiagonalAlign(offset) ) }

#define DISTMATRIX_SUBMATRIX(SIG,SIGBASE,T) \
  /* void AbstractDistMatrix<T>::GetSubmatrix
     ( const std::vector<Int>& I, const std::vector<Int>& J,
       AbstractDistMatrix<T>& ASub ) const */ \
  ElError ElDistMatrixGetSubmatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElInt numRowInds, const ElInt* I, \
    ElInt numColInds, const ElInt* J, ElDistMatrix_ ## SIG ASub ) \
  { EL_TRY( std::vector<Int> IVec(I,I+numRowInds); \
            std::vector<Int> JVec(J,J+numColInds); \
            CReflect(A)->GetSubmatrix(IVec,JVec,*CReflect(ASub)) ) } \
  /* void AbstractDistMatrix<T>::GetImagPartOfSubmatrix
     ( const std::vector<Int>& I, const std::vector<Int>& J,
       AbstractDistMatrix<Base<T>>& ASub ) const */ \
  ElError ElDistMatrixGetImagPartOfSubmatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElInt numRowInds, const ElInt* I, \
    ElInt numColInds, const ElInt* J, ElDistMatrix_ ## SIGBASE ASub ) \
  { EL_TRY( std::vector<Int> IVec(I,I+numRowInds); \
            std::vector<Int> JVec(J,J+numColInds); \
            CReflect(A)->GetImagPartOfSubmatrix(IVec,JVec,*CReflect(ASub)) ) } \
  /* void AbstractDistMatrix<T>::SetSubmatrix
    ( const std::vector<Int>& I, const std::vector<Int>& J,
      const AbstractDistMatrix<T>& ASub ); */ \
  ElError ElDistMatrixSetSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, const ElInt* I, const ElInt* J, \
    ElConstDistMatrix_ ## SIG ASub ) \
  { try { \
      const Int numRowInds = CReflect(ASub)->Height(); \
      const Int numColInds = CReflect(ASub)->Width(); \
      std::vector<Int> IVec(I,I+numRowInds), \
                       JVec(J,J+numColInds); \
      CReflect(A)->SetSubmatrix(IVec,JVec,*CReflect(ASub)); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void AbstractDistMatrix<T>::UpdateSubmatrix
    ( const std::vector<Int>& I, const std::vector<Int>& J,
      T alpha, const AbstractDistMatrix<T>& ASub ); */ \
  ElError ElDistMatrixUpdateSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, const ElInt* I, const ElInt* J, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG ASub ) \
  { try { \
      const Int numRowInds = CReflect(ASub)->Height(); \
      const Int numColInds = CReflect(ASub)->Width(); \
      std::vector<Int> IVec(I,I+numRowInds), \
                       JVec(J,J+numColInds); \
      CReflect(A)->UpdateSubmatrix(IVec,JVec,CReflect(alpha),*CReflect(ASub)); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void AbstractDistMatrix<T>::GetLocalSubmatrix
     ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc,
       Matrix<T>& ASub ) const */ \
  ElError ElDistMatrixGetLocalSubmatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElInt numRowInds, const ElInt* ILoc, \
    ElInt numColInds, const ElInt* JLoc, ElMatrix_ ## SIG ASub ) \
  { try { \
      std::vector<Int> IVec(ILoc,ILoc+numRowInds), \
                       JVec(JLoc,JLoc+numColInds); \
      CReflect(A)->GetLocalSubmatrix( IVec, JVec, *CReflect(ASub) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void AbstractDistMatrix<T>::GetImagPartOfLocalSubmatrix
     ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc,
       Matrix<Base<T>>& ASub ) const */ \
  ElError ElDistMatrixGetImagPartOfLocalSubmatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElInt numRowInds, const ElInt* ILoc, \
    ElInt numColInds, const ElInt* JLoc, ElMatrix_ ## SIGBASE ASub ) \
  { try { \
      std::vector<Int> IVec(ILoc,ILoc+numRowInds), \
                       JVec(JLoc,JLoc+numColInds); \
      CReflect(A)->GetImagPartOfLocalSubmatrix( IVec, JVec, *CReflect(ASub) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void AbstractDistMatrix<T>::SetLocalSubmatrix
     ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc,
       const Matrix<T>& ASub ); */ \
  ElError ElDistMatrixSetLocalSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    const ElInt* ILoc, const ElInt* JLoc, ElConstMatrix_ ## SIG ASub ) \
  { try { \
        const Int numRowInds = CReflect(ASub)->Height(); \
        const Int numColInds = CReflect(ASub)->Width(); \
        std::vector<Int> IVec(ILoc,ILoc+numRowInds), \
                         JVec(JLoc,JLoc+numColInds); \
        CReflect(A)->SetLocalSubmatrix(IVec,JVec,*CReflect(ASub)); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void AbstractDistMatrix<T>::UpdateLocalSubmatrix
     ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc,
       T alpha, const Matrix<T>& ASub ); */ \
  ElError ElDistMatrixUpdateLocalSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    const ElInt* ILoc, const ElInt* JLoc, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG ASub ) \
  { try { \
        const Int numRowInds = CReflect(ASub)->Height(); \
        const Int numColInds = CReflect(ASub)->Width(); \
        std::vector<Int> IVec(ILoc,ILoc+numRowInds), \
                         JVec(JLoc,JLoc+numColInds); \
        CReflect(A)->UpdateLocalSubmatrix \
        (IVec,JVec,CReflect(alpha),*CReflect(ASub)); \
    } EL_CATCH; return EL_SUCCESS; }

#define DISTMATRIX_SUBMATRIX_COMPLEX(SIG,SIGBASE,T) \
  /* void AbstractDistMatrix<T>::GetRealPartOfSubmatrix
     ( const std::vector<Int>& I, const std::vector<Int>& J,
       AbstractDistMatrix<Base<T>>& ASub ) const */ \
  ElError ElDistMatrixGetRealPartOfSubmatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElInt numRowInds, const ElInt* I, \
    ElInt numColInds, const ElInt* J, ElDistMatrix_ ## SIGBASE ASub ) \
  { try { std::vector<Int> IVec(I,I+numRowInds); \
          std::vector<Int> JVec(J,J+numColInds); \
          CReflect(A)->GetRealPartOfSubmatrix(IVec,JVec,*CReflect(ASub)); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void AbstractDistMatrix<T>::SetRealPartOfSubmatrix
    ( const std::vector<Int>& I, const std::vector<Int>& J,
      const AbstractDistMatrix<Base<T>>& ASub ); */ \
  ElError ElDistMatrixSetRealPartOfSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, const ElInt* I, const ElInt* J, \
    ElConstDistMatrix_ ## SIGBASE ASub ) \
  { try { \
      const Int numRowInds = CReflect(ASub)->Height(); \
      const Int numColInds = CReflect(ASub)->Width(); \
      std::vector<Int> IVec(I,I+numRowInds), \
                       JVec(J,J+numColInds); \
      CReflect(A)->SetRealPartOfSubmatrix(IVec,JVec,*CReflect(ASub)); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void AbstractDistMatrix<T>::SetImagPartOfSubmatrix
    ( const std::vector<Int>& I, const std::vector<Int>& J,
      const AbstractDistMatrix<Base<T>>& ASub ); */ \
  ElError ElDistMatrixSetImagPartOfSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, const ElInt* I, const ElInt* J, \
    ElConstDistMatrix_ ## SIGBASE ASub ) \
  { try { \
      const Int numRowInds = CReflect(ASub)->Height(); \
      const Int numColInds = CReflect(ASub)->Width(); \
      std::vector<Int> IVec(I,I+numRowInds), \
                       JVec(J,J+numColInds); \
      CReflect(A)->SetImagPartOfSubmatrix(IVec,JVec,*CReflect(ASub)); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void AbstractDistMatrix<T>::UpdateRealPartOfSubmatrix
    ( const std::vector<Int>& I, const std::vector<Int>& J,
      Base<T> alpha, const AbstractDistMatrix<Base<T>>& ASub ); */ \
  ElError ElDistMatrixUpdateRealPartOfSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, const ElInt* I, const ElInt* J, \
    Base<T> alpha, ElConstDistMatrix_ ## SIGBASE ASub ) \
  { try { \
      const Int numRowInds = CReflect(ASub)->Height(); \
      const Int numColInds = CReflect(ASub)->Width(); \
      std::vector<Int> IVec(I,I+numRowInds), \
                       JVec(J,J+numColInds); \
      CReflect(A)->UpdateRealPartOfSubmatrix(IVec,JVec,alpha,*CReflect(ASub)); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void AbstractDistMatrix<T>::UpdateImagPartOfSubmatrix
    ( const std::vector<Int>& I, const std::vector<Int>& J,
      Base<T> alpha, const AbstractDistMatrix<Base<T>>& ASub ); */ \
  ElError ElDistMatrixUpdateImagPartOfSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, const ElInt* I, const ElInt* J, \
    Base<T> alpha, ElConstDistMatrix_ ## SIGBASE ASub ) \
  { try { \
      const Int numRowInds = CReflect(ASub)->Height(); \
      const Int numColInds = CReflect(ASub)->Width(); \
      std::vector<Int> IVec(I,I+numRowInds), \
                       JVec(J,J+numColInds); \
      CReflect(A)->UpdateImagPartOfSubmatrix(IVec,JVec,alpha,*CReflect(ASub)); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void AbstractDistMatrix<T>::MakeSubmatrixReal
     ( const std::vector<Int>& I, const std::vector<Int>& J ) */ \
  ElError ElDistMatrixMakeSubmatrixReal_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt numRowInds, const ElInt* I, \
                                  ElInt numColInds, const ElInt* J ) \
  { try { \
      std::vector<Int> IVec(I,I+numRowInds), \
                       JVec(J,J+numColInds); \
          CReflect(A)->MakeSubmatrixReal( IVec, JVec ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void AbstractDistMatrix<T>::ConjugateSubmatrix
     ( const std::vector<Int>& I, const std::vector<Int>& J ) */ \
  ElError ElDistMatrixConjugateSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt numRowInds, const ElInt* I, \
                                  ElInt numColInds, const ElInt* J ) \
  { try { \
      std::vector<Int> IVec(I,I+numRowInds), \
                       JVec(J,J+numColInds); \
      CReflect(A)->ConjugateSubmatrix( IVec, JVec ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void AbstractDistMatrix<T>::GetRealPartOfLocalSubmatrix
     ( const std::vector<Int>& I, const std::vector<Int>& J,
       Matrix<Base<T>>& ASub ) const */ \
  ElError ElDistMatrixGetRealPartOfLocalSubmatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElInt numRowInds, const ElInt* ILoc, \
    ElInt numColInds, const ElInt* JLoc, \
    ElMatrix_ ## SIGBASE ASub ) \
  { try { \
        std::vector<Int> IVec(ILoc,ILoc+numRowInds), \
                         JVec(JLoc,JLoc+numColInds); \
        CReflect(A)->GetRealPartOfLocalSubmatrix(IVec,JVec,*CReflect(ASub)); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void AbstractDistMatrix<T>::SetRealPartOfLocalSubmatrix
     ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc,
       const Matrix<Base<T>>& ASub ); */ \
  ElError ElDistMatrixSetRealPartOfLocalSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    const ElInt* ILoc, const ElInt* JLoc, \
    ElConstMatrix_ ## SIGBASE ASub ) \
  { try { \
        const Int numRowInds = CReflect(ASub)->Height(); \
        const Int numColInds = CReflect(ASub)->Width(); \
        std::vector<Int> IVec(ILoc,ILoc+numRowInds), \
                         JVec(JLoc,JLoc+numColInds); \
        CReflect(A)->SetRealPartOfLocalSubmatrix(IVec,JVec,*CReflect(ASub)); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void AbstractDistMatrix<T>::SetImagPartOfLocalSubmatrix
     ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc,
       const Matrix<Base<T>>& ASub ); */ \
  ElError ElDistMatrixSetImagPartOfLocalSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    const ElInt* ILoc, const ElInt* JLoc, \
    ElConstMatrix_ ## SIGBASE ASub ) \
  { try { \
        const Int numRowInds = CReflect(ASub)->Height(); \
        const Int numColInds = CReflect(ASub)->Width(); \
        std::vector<Int> IVec(ILoc,ILoc+numRowInds), \
                         JVec(JLoc,JLoc+numColInds); \
        CReflect(A)->SetImagPartOfLocalSubmatrix(IVec,JVec,*CReflect(ASub)); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void AbstractDistMatrix<T>::UpdateRealPartOfLocalSubmatrix
     ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc,
       Base<T> alpha, const Matrix<Base<T>>& ASub ); */ \
  ElError ElDistMatrixUpdateRealPartOfLocalSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    const ElInt* ILoc, const ElInt* JLoc, \
    Base<T> alpha, ElConstMatrix_ ## SIGBASE ASub ) \
  { try { \
        const Int numRowInds = CReflect(ASub)->Height(); \
        const Int numColInds = CReflect(ASub)->Width(); \
        std::vector<Int> IVec(ILoc,ILoc+numRowInds), \
                         JVec(JLoc,JLoc+numColInds); \
        CReflect(A)->UpdateRealPartOfLocalSubmatrix \
        (IVec,JVec,alpha,*CReflect(ASub)); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void AbstractDistMatrix<T>::UpdateImagPartOfLocalSubmatrix
     ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc,
       Base<T> alpha, const Matrix<Base<T>>& ASub ); */ \
  ElError ElDistMatrixUpdateImagPartOfLocalSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    const ElInt* ILoc, const ElInt* JLoc, \
    Base<T> alpha, ElConstMatrix_ ## SIGBASE ASub ) \
  { try { \
        const Int numRowInds = CReflect(ASub)->Height(); \
        const Int numColInds = CReflect(ASub)->Width(); \
        std::vector<Int> IVec(ILoc,ILoc+numRowInds), \
                         JVec(JLoc,JLoc+numColInds); \
        CReflect(A)->UpdateImagPartOfLocalSubmatrix \
        (IVec,JVec,alpha,*CReflect(ASub)); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void AbstractDistMatrix<T>::MakeLocalSubmatrixReal
     ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc ) */ \
  ElError ElDistMatrixMakeLocalSubmatrixReal_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElInt numRowInds, const ElInt* ILoc, \
    ElInt numColInds, const ElInt* JLoc ) \
  { try { \
        std::vector<Int> IVec(ILoc,ILoc+numRowInds), \
                         JVec(JLoc,JLoc+numColInds); \
        CReflect(A)->MakeLocalSubmatrixReal( IVec, JVec ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* void AbstractDistMatrix<T>::ConjugateLocalSubmatrix
     ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc ) */ \
  ElError ElDistMatrixConjugateLocalSubmatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElInt numRowInds, const ElInt* ILoc, \
    ElInt numColInds, const ElInt* JLoc ) \
  { try { \
        std::vector<Int> IVec(ILoc,ILoc+numRowInds), \
                         JVec(JLoc,JLoc+numColInds); \
        CReflect(A)->ConjugateLocalSubmatrix( IVec, JVec ); \
    } EL_CATCH; return EL_SUCCESS; }

#define DISTMATRIX_SUM(SIG,SIGBASE,T) \
  ElError ElDistMatrixSumOver_ ## SIG \
  ( ElDistMatrix_ ## SIG A, MPI_Comm comm ) \
  { EL_TRY( CReflect(A)->SumOver( comm ) ) }

#define C_PROTO(SIG,SIGBASE,T) \
  DISTMATRIX_CREATE(SIG,SIGBASE,T) \
  DISTMATRIX_RECONFIG(SIG,SIGBASE,T) \
  DISTMATRIX_BASIC(SIG,SIGBASE,T) \
  DISTMATRIX_SINGLEENTRY(SIG,SIGBASE,T) \
  DISTMATRIX_DIAGONAL(SIG,SIGBASE,T) \
  DISTMATRIX_SUBMATRIX(SIG,SIGBASE,T) \
  DISTMATRIX_SUM(SIG,SIGBASE,T)

#define C_PROTO_COMPLEX(SIG,SIGBASE,T) \
  C_PROTO(SIG,SIGBASE,T) \
  DISTMATRIX_SINGLEENTRY_COMPLEX(SIG,SIGBASE,T) \
  DISTMATRIX_SUBMATRIX_COMPLEX(SIG,SIGBASE,T)

#include "El/macros/CInstantiate.h"

// TODO: More diagonal manipulation
// ================================

} // extern "C"
