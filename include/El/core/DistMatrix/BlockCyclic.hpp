/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DISTMATRIX_BLOCKCYCLIC_HPP
#define EL_DISTMATRIX_BLOCKCYCLIC_HPP

namespace El {

struct BlockCyclicData
{
    Dist colDist, rowDist;
    Int blockHeight, blockWidth;
    int colAlign, rowAlign;
    Int colCut, rowCut;
    int root;  // relevant for [o ,o ]/[MD,* ]/[* ,MD]
    const Grid* grid;

    BlockCyclicData() { }

    template<typename T>
    BlockCyclicData( const BlockCyclicMatrix<T>& A )
    : colDist(A.ColDist()), rowDist(A.RowDist()),
      blockHeight(A.BlockHeight()), blockWidth(A.BlockWidth()),
      colAlign(A.ColAlign()), rowAlign(A.RowAlign()),
      colCut(A.ColCut()), rowCut(A.RowCut()),
      root(A.Root()), grid(&A.Grid())
    { }
};
inline bool operator==( const BlockCyclicData& A, const BlockCyclicData& B )
{ return A.colDist     == B.colDist &&
         A.rowDist     == B.rowDist &&
         A.blockHeight == B.blockHeight &&
         A.blockWidth  == B.blockWidth &&
         A.colAlign    == B.colAlign &&
         A.rowAlign    == B.rowAlign &&
         A.root        == B.root &&
         A.grid        == B.grid; }

template<typename T> 
class BlockCyclicMatrix : public AbstractDistMatrix<T>
{
public:
    // Typedefs
    // ========
    typedef BlockCyclicMatrix<T> type;
    typedef AbstractDistMatrix<T> absType;

    // Constructors and destructors
    // ============================
    // Move constructor
    BlockCyclicMatrix( type&& A ) EL_NO_EXCEPT;

    virtual ~BlockCyclicMatrix();

    virtual BlockCyclicMatrix<T>* Construct
    ( const El::Grid& g, int root ) const = 0;
    virtual BlockCyclicMatrix<T>* ConstructTranspose
    ( const El::Grid& g, int root ) const = 0;
    virtual BlockCyclicMatrix<T>* ConstructDiagonal
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
    void FreeAlignments();
    void AlignWith( const El::BlockCyclicData& data, bool constrain=true );
    void AlignColsWith( const El::BlockCyclicData& data, bool constrain=true );
    void AlignRowsWith( const El::BlockCyclicData& data, bool constrain=true );
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
    DistWrap Wrap() const override EL_NO_EXCEPT { return BLOCK_CYCLIC; }
    virtual El::BlockCyclicData DistData() const = 0;

    // Distribution information
    // ------------------------
    Int BlockHeight() const EL_NO_EXCEPT;
    Int BlockWidth()  const EL_NO_EXCEPT;
    Int ColCut()      const EL_NO_EXCEPT;
    Int RowCut()      const EL_NO_EXCEPT;

    int RowOwner( Int i )       const override EL_NO_EXCEPT;
    int ColOwner( Int j )       const override EL_NO_EXCEPT;
    Int LocalRowOffset( Int i ) const override EL_NO_EXCEPT;
    Int LocalColOffset( Int j ) const override EL_NO_EXCEPT;
    Int GlobalRow( Int iLoc )   const override EL_NO_EXCEPT;
    Int GlobalCol( Int jLoc )   const override EL_NO_EXCEPT;

    // Diagonal manipulation
    // =====================
    bool DiagonalAlignedWith
    ( const El::BlockCyclicData& d, Int offset=0 ) const;
    int DiagonalRoot( Int offset=0 ) const;
    int DiagonalAlign( Int offset=0 ) const;

protected:
    // Member variables
    // ================

    // Process grid and distribution metadata
    // --------------------------------------
    Int blockHeight_, blockWidth_;
    Int colCut_,   rowCut_;

    // Private constructors
    // ====================
    BlockCyclicMatrix( const El::Grid& g=DefaultGrid(),  int root=0 );
    BlockCyclicMatrix
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
( const BlockCyclicMatrix<T>& AL,
  const BlockCyclicMatrix<T>& AR );

template<typename T>
void AssertConforming2x1
( const BlockCyclicMatrix<T>& AT,
  const BlockCyclicMatrix<T>& AB );

template<typename T>
void AssertConforming2x2
( const BlockCyclicMatrix<T>& ATL, 
  const BlockCyclicMatrix<T>& ATR,
  const BlockCyclicMatrix<T>& ABL, 
  const BlockCyclicMatrix<T>& ABR );

} // namespace El

#endif // ifndef EL_DISTMATRIX_BLOCKCYCLIC_HPP
