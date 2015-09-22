/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DISTMATRIX_BLOCK_HPP
#define EL_DISTMATRIX_BLOCK_HPP

namespace El {

template<typename T> 
class BlockMatrix : public AbstractDistMatrix<T>
{
public:
    // Typedefs
    // ========
    typedef BlockMatrix<T> type;
    typedef AbstractDistMatrix<T> absType;

    // Constructors and destructors
    // ============================
    // Move constructor
    BlockMatrix( type&& A ) EL_NO_EXCEPT;

    virtual ~BlockMatrix();

    virtual BlockMatrix<T>* Construct
    ( const El::Grid& g, int root ) const = 0;
    virtual BlockMatrix<T>* ConstructTranspose
    ( const El::Grid& g, int root ) const = 0;
    virtual BlockMatrix<T>* ConstructDiagonal
    ( const El::Grid& g, int root ) const = 0;

    // Assignment and reconfiguration
    // ==============================
    void Empty() override;
    void Resize( Int height, Int width ) override;
    void Resize( Int height, Int width, Int ldim ) override;

    void MakeConsistent( bool includingViewers=false );

    // Realignment
    // -----------
    void Align
    ( Int blockHeight, Int blockWidth, 
      int colAlign, int rowAlign, Int colCut=0, Int rowCut=0, 
      bool constrain=true );
    void AlignCols
    ( Int blockHeight, int colAlign, Int colCut=0, bool constrain=true );
    void AlignRows
    ( Int blockWidth, int rowAlign, Int rowCut=0, bool constrain=true );

    void AlignWith
    ( const El::DistData& data,
      bool constrain=true, bool allowMismatch=false ) override;
    void AlignColsWith
    ( const El::DistData& data,
      bool constrain=true, bool allowMismatch=false ) override;
    void AlignRowsWith
    ( const El::DistData& data,
      bool constrain=true, bool allowMismatch=false ) override;

    // TODO: The interface for these routines could be improved
    void AlignAndResize
    ( Int blockHeight, Int blockWidth, 
      int colAlign, int rowAlign, Int colCut, Int rowCut, 
      Int height, Int width, bool force=false, bool constrain=true );
    void AlignColsAndResize
    ( Int blockHeight, int colAlign, Int colCut, Int height, Int width, 
      bool force=false, bool constrain=true );
    void AlignRowsAndResize
    ( Int blockWidth, int rowAlign, Int rowCut, Int height, Int width, 
      bool force=false, bool constrain=true );

    // Buffer attachment
    // -----------------
    // (Immutable) view of a distributed matrix's buffer
    void Attach
    ( Int height, Int width, const El::Grid& g, 
      Int blockHeight, Int blockWidth,
      int colAlign, int rowAlign, Int colCut, Int rowCut,
      T* buffer, Int ldim, int root=0 );
    void LockedAttach
    ( Int height, Int width, const El::Grid& g,
      Int blockHeight, Int blockWidth,
      int colAlign, int rowAlign, Int colCut, Int rowCut,
      const T* buffer, Int ldim, int root=0 );
    void Attach
    ( Int height, Int width, const El::Grid& g,
      Int blockHeight, Int blockWidth,
      int colAlign, int rowAlign, Int colCut, Int rowCut,
      El::Matrix<T>& A, int root=0 );
    void LockedAttach
    ( Int height, Int width, const El::Grid& g,
      Int blockHeight, Int blockWidth,
      int colAlign, int rowAlign, Int colCut, Int rowCut,
      const El::Matrix<T>& A, int root=0 );

    // Operator overloading
    // ====================

    // Copy
    // ----
    const type& operator=( const type& A );

    // Addition/subtraction
    // --------------------
    const type& operator+=( const type& A );
    const type& operator-=( const type& A );

    // Move assignment
    // ---------------
    type& operator=( type&& A );

    // Basic queries
    // =============
    DistWrap Wrap() const override EL_NO_EXCEPT { return BLOCK; }
    virtual El::DistData DistData() const = 0;

    // Distribution information
    // ------------------------
    Int BlockHeight() const override EL_NO_EXCEPT;
    Int BlockWidth()  const override EL_NO_EXCEPT;
    Int ColCut()      const override EL_NO_EXCEPT;
    Int RowCut()      const override EL_NO_EXCEPT;

    int RowOwner( Int i )       const override EL_NO_EXCEPT;
    int ColOwner( Int j )       const override EL_NO_EXCEPT;
    Int LocalRowOffset( Int i ) const override EL_NO_EXCEPT;
    Int LocalColOffset( Int j ) const override EL_NO_EXCEPT;
    Int GlobalRow( Int iLoc )   const override EL_NO_EXCEPT;
    Int GlobalCol( Int jLoc )   const override EL_NO_EXCEPT;

    // Diagonal manipulation
    // =====================
    bool DiagonalAlignedWith
    ( const El::DistData& d, Int offset=0 ) const override;
    int DiagonalRoot( Int offset=0 ) const override;
    int DiagonalAlign( Int offset=0 ) const override;

protected:
    // Member variables
    // ================

    // Process grid and distribution metadata
    // --------------------------------------
    Int blockHeight_, blockWidth_;
    Int colCut_,   rowCut_;

    // Private constructors
    // ====================
    BlockMatrix( const El::Grid& g=DefaultGrid(),  int root=0 );
    BlockMatrix
    ( const El::Grid& g, Int blockHeight, Int blockWidth, int root=0 );

private:
    Int NewLocalHeight( Int height ) const;
    Int NewLocalWidth( Int width ) const;

    Int NewLocalHeight_( Int height ) const EL_NO_EXCEPT;
    Int NewLocalWidth_( Int width ) const EL_NO_EXCEPT;

    // Exchange metadata with another matrix
    // =====================================
    void ShallowSwap( type& A );

    template<typename S,Dist J,Dist K,DistWrap wrap> friend class DistMatrix;
};

template<typename T>
void AssertConforming1x2
( const BlockMatrix<T>& AL,
  const BlockMatrix<T>& AR );

template<typename T>
void AssertConforming2x1
( const BlockMatrix<T>& AT,
  const BlockMatrix<T>& AB );

template<typename T>
void AssertConforming2x2
( const BlockMatrix<T>& ATL, 
  const BlockMatrix<T>& ATR,
  const BlockMatrix<T>& ABL, 
  const BlockMatrix<T>& ABR );

} // namespace El

#endif // ifndef EL_DISTMATRIX_BLOCK_HPP
