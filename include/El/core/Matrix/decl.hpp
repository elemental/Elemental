/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_MATRIX_DECL_HPP
#define EL_MATRIX_DECL_HPP

namespace El {

// Matrix base for arbitrary rings
template<typename scalarType>
class Matrix
{
public:    
    // Constructors and destructors
    // ============================

    // Create a 0x0 matrix
    Matrix( bool fixed=false );

    // Create a matrix with the specified dimensions
    Matrix( Int height, Int width, bool fixed=false );

    // Create a matrix with the specified dimensions and leading dimension
    Matrix( Int height, Int width, Int ldim, bool fixed=false );

    // Construct a matrix around an existing (possibly immutable) buffer
    Matrix
    ( Int height,
      Int width,
      const scalarType* buffer,
      Int ldim,
      bool fixed=false );
    Matrix
    ( Int height,
      Int width,
      scalarType* buffer,
      Int ldim,
      bool fixed=false );

    // Create a copy of a matrix
    Matrix( const Matrix<scalarType>& A );

    // Move the metadata from a given matrix
    Matrix( Matrix<scalarType>&& A ) EL_NO_EXCEPT;

    // Destructor
    ~Matrix();

    // Assignment and reconfiguration
    // ==============================

    void Empty( bool freeMemory=true );
    void Resize( Int height, Int width );
    void Resize( Int height, Int width, Int ldim );

    // Reconfigure around the given buffer, but do not assume ownership
    void Attach
    ( Int height, Int width, scalarType* buffer, Int ldim );
    void LockedAttach
    ( Int height, Int width, const scalarType* buffer, Int ldim );

    // Reconfigure around the given buffer and assume ownership
    void Control( Int height, Int width, scalarType* buffer, Int ldim );

    // Operator overloading
    // ====================

    // Return a view
    // -------------
          Matrix<scalarType> operator()( Range<Int> I, Range<Int> J );
    const Matrix<scalarType> operator()( Range<Int> I, Range<Int> J ) const;

    // Return a copy of (potentially non-contiguous) subset of indices
    // ---------------------------------------------------------------
    Matrix<scalarType>
    operator()( Range<Int> I, const vector<Int>& J ) const;
    Matrix<scalarType>
    operator()( const vector<Int>& I, Range<Int> J ) const;
    Matrix<scalarType>
    operator()( const vector<Int>& I, const vector<Int>& J ) const;

    // Make a copy
    // -----------
    const Matrix<scalarType>& operator=( const Matrix<scalarType>& A );

    // Move assignment
    // ---------------
    Matrix<scalarType>& operator=( Matrix<scalarType>&& A );

    // Rescaling
    // ---------
    const Matrix<scalarType>& operator*=( scalarType alpha );

    // Addition/substraction
    // ---------------------
    const Matrix<scalarType>& operator+=( const Matrix<scalarType>& A );
    const Matrix<scalarType>& operator-=( const Matrix<scalarType>& A );

    // Basic queries
    // =============
    Int Height() const EL_NO_EXCEPT;
    Int Width() const EL_NO_EXCEPT;
    Int LDim() const EL_NO_EXCEPT;
    Int MemorySize() const EL_NO_EXCEPT;
    Int DiagonalLength( Int offset=0 ) const EL_NO_EXCEPT;

    scalarType* Buffer() EL_NO_RELEASE_EXCEPT;
    scalarType* Buffer( Int i, Int j ) EL_NO_RELEASE_EXCEPT;
    const scalarType* LockedBuffer() const EL_NO_EXCEPT;
    const scalarType* LockedBuffer( Int i, Int j ) const EL_NO_EXCEPT;

    bool Viewing()   const EL_NO_EXCEPT;
    bool FixedSize() const EL_NO_EXCEPT;
    bool Locked()    const EL_NO_EXCEPT;
    // Advanced
    // --------
    void SetViewType( El::ViewType viewType ) EL_NO_EXCEPT;
    El::ViewType ViewType() const EL_NO_EXCEPT;

    // Single-entry manipulation
    // =========================
    scalarType Get( Int i, Int j ) const EL_NO_RELEASE_EXCEPT;
    Base<scalarType> GetRealPart( Int i, Int j ) const EL_NO_RELEASE_EXCEPT;
    Base<scalarType> GetImagPart( Int i, Int j ) const EL_NO_RELEASE_EXCEPT;

    void Set( Int i, Int j, scalarType alpha ) EL_NO_RELEASE_EXCEPT;
    void Set( const Entry<scalarType>& entry ) EL_NO_RELEASE_EXCEPT;

    void SetRealPart
    ( Int i, Int j, Base<scalarType> alpha ) EL_NO_RELEASE_EXCEPT;
    void SetImagPart
    ( Int i, Int j, Base<scalarType> alpha ) EL_NO_RELEASE_EXCEPT;

    void SetRealPart
    ( const Entry<Base<scalarType>>& entry ) EL_NO_RELEASE_EXCEPT;
    void SetImagPart
    ( const Entry<Base<scalarType>>& entry ) EL_NO_RELEASE_EXCEPT;

    void Update( Int i, Int j, scalarType alpha ) EL_NO_RELEASE_EXCEPT;
    void Update( const Entry<scalarType>& entry ) EL_NO_RELEASE_EXCEPT;

    void UpdateRealPart
    ( Int i, Int j, Base<scalarType> alpha ) EL_NO_RELEASE_EXCEPT;
    void UpdateImagPart
    ( Int i, Int j, Base<scalarType> alpha ) EL_NO_RELEASE_EXCEPT;

    void UpdateRealPart
    ( const Entry<Base<scalarType>>& entry ) EL_NO_RELEASE_EXCEPT;
    void UpdateImagPart
    ( const Entry<Base<scalarType>>& entry ) EL_NO_RELEASE_EXCEPT;

    void MakeReal( Int i, Int j ) EL_NO_RELEASE_EXCEPT;
    void Conjugate( Int i, Int j ) EL_NO_RELEASE_EXCEPT;

private:
    // Member variables
    // ================
    El::ViewType viewType_;
    Int height_, width_, ldim_;

    Memory<scalarType> memory_;
    // Const-correctness is internally managed to avoid the need for storing
    // two separate pointers with different 'const' attributes
    scalarType* data_;

    // Exchange metadata with another matrix
    // =====================================
    void ShallowSwap( Matrix<scalarType>& A );

    // Reconfigure without error-checking
    // ==================================
    void Empty_( bool freeMemory=true );
    void Resize_( Int height, Int width );
    void Resize_( Int height, Int width, Int ldim );

    void Control_
    ( Int height, Int width, scalarType* buffer, Int ldim );
    void Attach_
    ( Int height, Int width, scalarType* buffer, Int ldim );
    void LockedAttach_
    ( Int height, Int width, const scalarType* buffer, Int ldim );

    // Return a reference to a single entry without error-checking
    // ===========================================================
    const scalarType& Get_( Int i, Int j ) const EL_NO_RELEASE_EXCEPT;
    scalarType& Set_( Int i, Int j ) EL_NO_RELEASE_EXCEPT;

    // Assertions
    // ==========
    void AssertValidDimensions( Int height, Int width ) const;
    void AssertValidDimensions( Int height, Int width, Int ldim ) const;
    void AssertValidEntry( Int i, Int j ) const;
   
    // Friend declarations
    // ===================
    template<typename S> friend class Matrix;
    template<typename S> friend class AbstractDistMatrix;
    template<typename S> friend class ElementalMatrix;
    template<typename S> friend class BlockMatrix;
};

} // namespace El

#endif // ifndef EL_MATRIX_DECL_HPP
