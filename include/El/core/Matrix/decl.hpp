/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_MATRIX_DECL_HPP
#define EL_MATRIX_DECL_HPP

#include <El/core/Grid.hpp>

namespace El {

// Matrix base for arbitrary rings
template<typename Ring>
class Matrix
{
public:    
    // Constructors and destructors
    // ============================

    // Create a 0x0 matrix
    Matrix();

    // Create a matrix with the specified dimensions
    Matrix( Int height, Int width );

    // Create a matrix with the specified dimensions and leading dimension
    Matrix( Int height, Int width, Int leadingDimension );

    // Construct a matrix around an existing (possibly immutable) buffer
    Matrix
    ( Int height,
      Int width,
      const Ring* buffer,
      Int leadingDimension );
    Matrix
    ( Int height,
      Int width,
      Ring* buffer,
      Int leadingDimension );

    // Create a copy of a matrix
    Matrix( const Matrix<Ring>& A );

    // Move the metadata from a given matrix
    Matrix( Matrix<Ring>&& A ) EL_NO_EXCEPT;

    // Destructor
    ~Matrix();

    // Assignment and reconfiguration
    // ==============================

    void Empty( bool freeMemory=true );
    void Resize( Int height, Int width );
    void Resize( Int height, Int width, Int leadingDimension );

    // Reconfigure around the given buffer, but do not assume ownership
    void Attach
    ( Int height, Int width, Ring* buffer, Int leadingDimension );
    void LockedAttach
    ( Int height, Int width, const Ring* buffer, Int leadingDimension );

    // Reconfigure around the given buffer and assume ownership
    void Control
    ( Int height, Int width, Ring* buffer, Int leadingDimension );

    // Force the size to remain constant (but allow the entries to be modified).
    void FixSize() EL_NO_EXCEPT;

    // Operator overloading
    // ====================

    // Return a view
    // -------------
          Matrix<Ring> operator()( Range<Int> I, Range<Int> J );
    const Matrix<Ring> operator()( Range<Int> I, Range<Int> J ) const;

    // Return a copy of (potentially non-contiguous) subset of indices
    // ---------------------------------------------------------------
    Matrix<Ring>
    operator()( Range<Int> I, const vector<Int>& J ) const;
    Matrix<Ring>
    operator()( const vector<Int>& I, Range<Int> J ) const;
    Matrix<Ring>
    operator()( const vector<Int>& I, const vector<Int>& J ) const;

    // Make a copy
    // -----------
    const Matrix<Ring>& operator=( const Matrix<Ring>& A );

    // Move assignment
    // ---------------
    Matrix<Ring>& operator=( Matrix<Ring>&& A );

    // Rescaling
    // ---------
    const Matrix<Ring>& operator*=( const Ring& alpha );

    // Addition/substraction
    // ---------------------
    const Matrix<Ring>& operator+=( const Matrix<Ring>& A );
    const Matrix<Ring>& operator-=( const Matrix<Ring>& A );

    // Basic queries
    // =============
    Int Height() const EL_NO_EXCEPT;
    Int Width() const EL_NO_EXCEPT;
    Int LDim() const EL_NO_EXCEPT;
    Int MemorySize() const EL_NO_EXCEPT;
    Int DiagonalLength( Int offset=0 ) const EL_NO_EXCEPT;

    Ring* Buffer() EL_NO_RELEASE_EXCEPT;
    Ring* Buffer( Int i, Int j ) EL_NO_RELEASE_EXCEPT;
    const Ring* LockedBuffer() const EL_NO_EXCEPT;
    const Ring* LockedBuffer( Int i, Int j ) const EL_NO_EXCEPT;

    bool Viewing()   const EL_NO_EXCEPT;
    bool FixedSize() const EL_NO_EXCEPT;
    bool Locked()    const EL_NO_EXCEPT;
    // Advanced
    // --------
    void SetViewType( El::ViewType viewType ) EL_NO_EXCEPT;
    El::ViewType ViewType() const EL_NO_EXCEPT;

    // Single-entry manipulation
    // =========================
    Ring Get( Int i, Int j=0 ) const EL_NO_RELEASE_EXCEPT;
    Base<Ring> GetRealPart( Int i, Int j=0 ) const EL_NO_RELEASE_EXCEPT;
    Base<Ring> GetImagPart( Int i, Int j=0 ) const EL_NO_RELEASE_EXCEPT;

    void Set( Int i, Int j, const Ring& alpha ) EL_NO_RELEASE_EXCEPT;
    void Set( const Entry<Ring>& entry ) EL_NO_RELEASE_EXCEPT;

    void SetRealPart
    ( Int i, Int j, const Base<Ring>& alpha ) EL_NO_RELEASE_EXCEPT;
    void SetImagPart
    ( Int i, Int j, const Base<Ring>& alpha ) EL_NO_RELEASE_EXCEPT;

    void SetRealPart
    ( const Entry<Base<Ring>>& entry ) EL_NO_RELEASE_EXCEPT;
    void SetImagPart
    ( const Entry<Base<Ring>>& entry ) EL_NO_RELEASE_EXCEPT;

    void Update( Int i, Int j, const Ring& alpha ) EL_NO_RELEASE_EXCEPT;
    void Update( const Entry<Ring>& entry ) EL_NO_RELEASE_EXCEPT;

    void UpdateRealPart
    ( Int i, Int j, const Base<Ring>& alpha ) EL_NO_RELEASE_EXCEPT;
    void UpdateImagPart
    ( Int i, Int j, const Base<Ring>& alpha ) EL_NO_RELEASE_EXCEPT;

    void UpdateRealPart
    ( const Entry<Base<Ring>>& entry ) EL_NO_RELEASE_EXCEPT;
    void UpdateImagPart
    ( const Entry<Base<Ring>>& entry ) EL_NO_RELEASE_EXCEPT;

    void MakeReal( Int i, Int j ) EL_NO_RELEASE_EXCEPT;
    void Conjugate( Int i, Int j ) EL_NO_RELEASE_EXCEPT;

    // Return a reference to a single entry without error-checking
    // -----------------------------------------------------------
    inline const Ring& CRef( Int i, Int j=0 ) const EL_NO_RELEASE_EXCEPT;
    inline const Ring& operator()( Int i, Int j=0 ) const EL_NO_RELEASE_EXCEPT;

    inline Ring& Ref( Int i, Int j=0 ) EL_NO_RELEASE_EXCEPT;
    inline Ring& operator()( Int i, Int j=0 ) EL_NO_RELEASE_EXCEPT;

private:
    // Member variables
    // ================
    El::ViewType viewType_=OWNER;
    Int height_=0, width_=0, leadingDimension_=1;

    Memory<Ring> memory_;
    // Const-correctness is internally managed to avoid the need for storing
    // two separate pointers with different 'const' attributes
    Ring* data_=nullptr;

    // Exchange metadata with another matrix
    // =====================================
    void ShallowSwap( Matrix<Ring>& A );

    // Reconfigure without error-checking
    // ==================================
    void Empty_( bool freeMemory=true );
    void Resize_( Int height, Int width );
    void Resize_( Int height, Int width, Int leadingDimension );

    void Control_
    ( Int height, Int width, Ring* buffer, Int leadingDimension );
    void Attach_
    ( Int height, Int width, Ring* buffer, Int leadingDimension );
    void LockedAttach_
    ( Int height, Int width, const Ring* buffer, Int leadingDimension );

    // Assertions
    // ==========
    void AssertValidDimensions( Int height, Int width ) const;
    void AssertValidDimensions
    ( Int height, Int width, Int leadingDimension ) const;
    void AssertValidEntry( Int i, Int j ) const;
   
    // Friend declarations
    // ===================
    template<typename S> friend class Matrix;
    template<typename S> friend class AbstractDistMatrix;
    template<typename S> friend class ElementalMatrix;
    template<typename S> friend class BlockMatrix;

    // For supporting duck typing
    // ==========================
    // The following are provided in order to aid duck-typing over
    // {Matrix, DistMatrix, DistMultiVec, etc.}.

    // This is equivalent to the trivial constructor in functionality
    // (though an error is thrown if 'grid' is not equal to 'Grid::Trivial()').
    explicit Matrix( const El::Grid& grid );

    // This is a no-op
    // (though an error is thrown if 'grid' is not equal to 'Grid::Trivial()').
    void SetGrid( const El::Grid& grid );

    // This always returns 'Grid::Trivial()'.
    const El::Grid& Grid() const;

    // This is a no-op
    // (though an error is thrown if 'colAlign' or 'rowAlign' is not zero).
    void Align( Int colAlign, Int rowAlign, bool constrain=true );

    // These always return 0.
    int ColAlign() const EL_NO_EXCEPT;
    int RowAlign() const EL_NO_EXCEPT;
};

} // namespace El

#endif // ifndef EL_MATRIX_DECL_HPP
