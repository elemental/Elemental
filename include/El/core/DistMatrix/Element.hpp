/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DISTMATRIX_ELEMENTAL_HPP
#define EL_DISTMATRIX_ELEMENTAL_HPP

namespace El {

template<typename T> 
class ElementalMatrix : public AbstractDistMatrix<T>
{
public:
    // Typedefs
    // ========
    typedef ElementalMatrix<T> type;
    typedef AbstractDistMatrix<T> absType;

    // Constructors and destructors
    // ============================
    // Move constructor
    ElementalMatrix( type&& A ) EL_NO_EXCEPT;

    virtual ~ElementalMatrix();

    virtual ElementalMatrix<T>* Construct
    ( const El::Grid& g, int root ) const = 0;
    virtual ElementalMatrix<T>* ConstructTranspose
    ( const El::Grid& g, int root ) const = 0;
    virtual ElementalMatrix<T>* ConstructDiagonal
    ( const El::Grid& g, int root ) const = 0;
    // TODO: ConstructPartialCol and friends?

    // Assignment and reconfiguration
    // ==============================
    void Resize( Int height, Int width ) override;
    void Resize( Int height, Int width, Int ldim ) override;

    void MakeConsistent( bool includingViewers=false );

    // Realignment
    // -----------
    void Align( int colAlign, int rowAlign, bool constrain=true );
    void AlignCols( int colAlign, bool constrain=true );
    void AlignRows( int rowAlign, bool constrain=true );

    void AlignWith
    ( const El::DistData& data,
      bool constrain=true, bool allowMismatch=false ) override;
    void AlignColsWith
    ( const El::DistData& data,
      bool constrain=true, bool allowMismatch=false ) override;
    void AlignRowsWith
    ( const El::DistData& data,
      bool constrain=true, bool allowMismatch=false ) override;

    void AlignWith
    ( const El::ElementalData& data,
      bool constrain=true, bool allowMismatch=false );
    void AlignColsWith
    ( const El::ElementalData& data,
      bool constrain=true, bool allowMismatch=false );
    void AlignRowsWith
    ( const El::ElementalData& data,
      bool constrain=true, bool allowMismatch=false );

    void AlignAndResize
    ( int colAlign, int rowAlign, Int height, Int width, 
      bool force=false, bool constrain=true );
    void AlignColsAndResize
    ( int colAlign, Int height, Int width, 
      bool force=false, bool constrain=true );
    void AlignRowsAndResize
    ( int rowAlign, Int height, Int width, 
      bool force=false, bool constrain=true );

    // Buffer attachment
    // -----------------
    // (Immutable) view of a distributed matrix's buffer
    void Attach
    ( Int height, Int width, const El::Grid& grid, 
      int colAlign, int rowAlign, T* buffer, Int ldim, int root=0 );
    void LockedAttach
    ( Int height, Int width, const El::Grid& grid,
      int colAlign, int rowAlign, const T* buffer, Int ldim, int root=0 );
    void Attach
    ( Int height, Int width, const El::Grid& grid,
      int colAlign, int rowAlign, El::Matrix<T>& A, int root=0 );
    void LockedAttach
    ( Int height, Int width, const El::Grid& grid,
      int colAlign, int rowAlign, const El::Matrix<T>& A, int root=0 );
    // (Immutable) view of a local matrix's buffer
    void Attach( const El::Grid& grid, El::Matrix<T>& A );
    void LockedAttach( const El::Grid& grid, const El::Matrix<T>& A );

    // Operator overloading
    // ====================

    // Copy
    // ----
    const type& operator=( const type& A );
    // TODO: Eliminate this routine
    const type& operator=( const DistMultiVec<T>& A );

    // Addition/subtraction
    // --------------------
    const type& operator+=( const type& A );
    const type& operator-=( const type& A );

    // Move assignment
    // ---------------
    type& operator=( type&& A );

    // Basic queries
    // =============
    DistWrap Wrap() const override EL_NO_EXCEPT { return ELEMENT; }
    virtual El::ElementalData DistData() const = 0;

    Int BlockHeight() const override EL_NO_EXCEPT { return 1; }
    Int BlockWidth()  const override EL_NO_EXCEPT { return 1; }
    Int ColCut()      const override EL_NO_EXCEPT { return 0; }
    Int RowCut()      const override EL_NO_EXCEPT { return 0; }

    int  RowOwner( Int i )       const override EL_NO_EXCEPT;
    int  ColOwner( Int j )       const override EL_NO_EXCEPT;
    Int  LocalRowOffset( Int i ) const override EL_NO_EXCEPT;
    Int  LocalColOffset( Int j ) const override EL_NO_EXCEPT;
    Int  GlobalRow( Int iLoc )   const override EL_NO_EXCEPT;
    Int  GlobalCol( Int jLoc )   const override EL_NO_EXCEPT;

    // Diagonal manipulation
    // =====================
    bool DiagonalAlignedWith
    ( const El::DistData& d, Int offset=0 ) const override EL_NO_EXCEPT;
    bool DiagonalAlignedWith
    ( const El::ElementalData& d, Int offset=0 ) const EL_NO_EXCEPT;
    int DiagonalRoot( Int offset=0 ) const override EL_NO_EXCEPT;
    int DiagonalAlign( Int offset=0 ) const override EL_NO_EXCEPT;

protected:
    // Protected constructors
    // ======================
    // Create a 0 x 0 distributed matrix
    ElementalMatrix( const El::Grid& g=DefaultGrid(), int root=0 );

private:
    // Exchange metadata with another matrix
    // =====================================
    void ShallowSwap( type& A );

    template<typename S,Dist J,Dist K,DistWrap wrap> friend class DistMatrix;
};

template<typename T>
void AssertConforming1x2
( const ElementalMatrix<T>& AL, const ElementalMatrix<T>& AR );

template<typename T>
void AssertConforming2x1
( const ElementalMatrix<T>& AT, const ElementalMatrix<T>& AB );

template<typename T>
void AssertConforming2x2
( const ElementalMatrix<T>& ATL, const ElementalMatrix<T>& ATR,
  const ElementalMatrix<T>& ABL, const ElementalMatrix<T>& ABR );

} // namespace El

#endif // ifndef EL_DISTMATRIX_ELEMENTAL_HPP
