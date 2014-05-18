/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLOCKDISTMATRIX_GENERAL_DECL_HPP
#define EL_BLOCKDISTMATRIX_GENERAL_DECL_HPP

namespace El {

template<typename T,Dist U,Dist V> 
class GeneralBlockDistMatrix : public AbstractBlockDistMatrix<T>
{
public:
    // Typedefs
    // ========
    typedef AbstractBlockDistMatrix<T> absType;
    typedef GeneralBlockDistMatrix<T,U,V> type;
    static constexpr Dist UDiag = DiagColDist<U,V>();
    static constexpr Dist VDiag = DiagRowDist<U,V>();
    static constexpr Dist UGath = GatheredDist<U>();
    static constexpr Dist VGath = GatheredDist<V>();
    static constexpr Dist UPart = PartialDist<U>();
    static constexpr Dist VPart = PartialDist<V>();
    static constexpr Dist UScat = ScatteredColDist<U,V>();
    static constexpr Dist VScat = ScatteredRowDist<U,V>();

    // Constructors and destructors
    // ============================
    // Move constructor
    GeneralBlockDistMatrix( type&& A ) EL_NOEXCEPT;

    // Assignment and reconfiguration
    // ==============================
    // Move assignment
    type& operator=( type&& A );

    void AlignColsWith
    ( const El::BlockDistData& data, bool constrain=true ) override;
    void AlignRowsWith
    ( const El::BlockDistData& data, bool constrain=true ) override;

    void Translate( BlockDistMatrix<T,U,V>& A ) const;

    void AllGather( BlockDistMatrix<T,UGath,VGath>& A ) const;
    void ColAllGather( BlockDistMatrix<T,UGath,V>& A ) const;
    void RowAllGather( BlockDistMatrix<T,U,VGath>& A ) const;
    void PartialColAllGather( BlockDistMatrix<T,UPart,V>& A ) const;
    void PartialRowAllGather( BlockDistMatrix<T,U,VPart>& A ) const;

    void FilterFrom( const BlockDistMatrix<T,UGath,VGath>& A );
    void ColFilterFrom( const BlockDistMatrix<T,UGath,V>& A );
    void RowFilterFrom( const BlockDistMatrix<T,U,VGath>& A );
    void PartialColFilterFrom( const BlockDistMatrix<T,UPart,V>& A );
    void PartialRowFilterFrom( const BlockDistMatrix<T,U,VPart>& A );

    void PartialColAllToAllFrom( const BlockDistMatrix<T,UPart,VScat>& A );
    void PartialRowAllToAllFrom( const BlockDistMatrix<T,UScat,VPart>& A );
    void PartialColAllToAll( BlockDistMatrix<T,UPart,VScat>& A ) const;
    void PartialRowAllToAll( BlockDistMatrix<T,UScat,VPart>& A ) const;

    void SumScatterFrom( const BlockDistMatrix<T,UGath,VGath>& A );
    void RowSumScatterFrom( const BlockDistMatrix<T,U,VGath>& A );
    void ColSumScatterFrom( const BlockDistMatrix<T,UGath,V>& A );
    void PartialRowSumScatterFrom( const BlockDistMatrix<T,U,VPart>& A );
    void PartialColSumScatterFrom( const BlockDistMatrix<T,UPart,V>& A );

    void SumScatterUpdate( T alpha, const BlockDistMatrix<T,UGath,VGath>& A );
    void RowSumScatterUpdate( T alpha, const BlockDistMatrix<T,U,VGath>& A );
    void ColSumScatterUpdate( T alpha, const BlockDistMatrix<T,UGath,V>& A );
    void PartialRowSumScatterUpdate
    ( T alpha, const BlockDistMatrix<T,U,VPart>& A );
    void PartialColSumScatterUpdate
    ( T alpha, const BlockDistMatrix<T,UPart,V>& A );

    // Transposed variants of some of the above routines which avoid 
    // large amounts of non-uniform data access
    // -------------------------------------------------------------
    void TransposeColAllGather
    ( BlockDistMatrix<T,V,UGath>& A, bool conjugate=false ) const;
    void TransposePartialColAllGather
    ( BlockDistMatrix<T,V,UPart>& A, bool conjugate=false ) const;
    void AdjointColAllGather( BlockDistMatrix<T,V,UGath>& A ) const;
    void AdjointPartialColAllGather( BlockDistMatrix<T,V,UPart>& A ) const;

    void TransposeColFilterFrom
    ( const BlockDistMatrix<T,V,UGath>& A, bool conjugate=false );
    void TransposeRowFilterFrom
    ( const BlockDistMatrix<T,VGath,U>& A, bool conjugate=false );
    void TransposePartialColFilterFrom
    ( const BlockDistMatrix<T,V,UPart>& A, bool conjugate=false );
    void TransposePartialRowFilterFrom
    ( const BlockDistMatrix<T,VPart,U>& A, bool conjugate=false );
    void AdjointColFilterFrom( const BlockDistMatrix<T,V,UGath>& A );
    void AdjointRowFilterFrom( const BlockDistMatrix<T,VGath,U>& A );
    void AdjointPartialColFilterFrom( const BlockDistMatrix<T,V,UPart>& A );
    void AdjointPartialRowFilterFrom( const BlockDistMatrix<T,VPart,U>& A );

    void TransposeColSumScatterFrom
    ( const BlockDistMatrix<T,V,UGath>& A, bool conjugate=false );
    void TransposePartialColSumScatterFrom
    ( const BlockDistMatrix<T,V,UPart>& A, bool conjugate=false );
    void AdjointColSumScatterFrom( const BlockDistMatrix<T,V,UGath>& A );
    void AdjointPartialColSumScatterFrom( const BlockDistMatrix<T,V,UPart>& A );

    void TransposeColSumScatterUpdate
    ( T alpha, const BlockDistMatrix<T,V,UGath>& A, bool conjugate=false );
    void TransposePartialColSumScatterUpdate
    ( T alpha, const BlockDistMatrix<T,V,UPart>& A, bool conjugate=false );
    void AdjointColSumScatterUpdate
    ( T alpha, const BlockDistMatrix<T,V,UGath>& A );
    void AdjointPartialColSumScatterUpdate
    ( T alpha, const BlockDistMatrix<T,V,UPart>& A );

    // Diagonal manipulation
    // =====================
    bool DiagonalAlignedWith
    ( const El::BlockDistData& d, Int offset=0 ) const override;
    Int DiagonalRoot( Int offset=0 ) const override;
    Int DiagonalAlign( Int offset=0 ) const override;

    void GetDiagonal( BlockDistMatrix<T,UDiag,VDiag>& d, Int offset=0 ) const;
    void GetRealPartOfDiagonal
    ( BlockDistMatrix<Base<T>,UDiag,VDiag>& d, Int offset=0 ) const;
    void GetImagPartOfDiagonal
    ( BlockDistMatrix<Base<T>,UDiag,VDiag>& d, Int offset=0 ) const;

    BlockDistMatrix<T,UDiag,VDiag> GetDiagonal( Int offset=0 ) const;
    BlockDistMatrix<Base<T>,UDiag,VDiag> 
    GetRealPartOfDiagonal( Int offset=0 ) const;
    BlockDistMatrix<Base<T>,UDiag,VDiag> 
    GetImagPartOfDiagonal( Int offset=0 ) const;

    void SetDiagonal( const BlockDistMatrix<T,UDiag,VDiag>& d, Int offset=0 );
    void SetRealPartOfDiagonal
    ( const BlockDistMatrix<Base<T>,UDiag,VDiag>& d, Int offset=0 );
    void SetImagPartOfDiagonal
    ( const BlockDistMatrix<Base<T>,UDiag,VDiag>& d, Int offset=0 );

    void UpdateDiagonal
    ( T alpha, const BlockDistMatrix<T,UDiag,VDiag>& d, Int offset=0 );
    void UpdateRealPartOfDiagonal
    ( Base<T> alpha, const BlockDistMatrix<Base<T>,UDiag,VDiag>& d, 
      Int offset=0 );
    void UpdateImagPartOfDiagonal
    ( Base<T> alpha, const BlockDistMatrix<Base<T>,UDiag,VDiag>& d, 
      Int offset=0 );

protected:

    // Private constructors
    // ====================

    // Inherited constructors are part of C++11 but not yet widely supported.
    //using AbstractDistMatrix<T>::AbstractDistMatrix;

    GeneralBlockDistMatrix( const El::Grid& g=DefaultGrid(), Int root=0 );
    GeneralBlockDistMatrix
    ( const El::Grid& g, Int blockHeight, Int blockWidth, Int root=0 );

    // Diagonal helper routines
    // ========================
    template<typename S,class Function>
    void GetDiagonalHelper
    ( BlockDistMatrix<S,UDiag,VDiag>& d, Int offset, Function func ) const;
    template<typename S,class Function>
    void SetDiagonalHelper
    ( const BlockDistMatrix<S,UDiag,VDiag>& d, Int offset, Function func );

    // Friend declarations
    // ===================
    template<typename S,Dist J,Dist K> friend class DistMatrix;
    template<typename S,Dist J,Dist K> friend class BlockDistMatrix;
};

} // namespace El

#endif // ifndef EL_BLOCKDISTMATRIX_GENERAL_DECL_HPP
