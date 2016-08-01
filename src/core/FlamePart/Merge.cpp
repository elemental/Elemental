/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>

namespace El {

// Horizontally merge two contiguous matrices
// ==========================================

template<typename T>
void Merge1x2( Matrix<T>& A, Matrix<T>& BL, Matrix<T>& BR )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( BL.Locked() || BR.Locked() )
          LogicError("Cannot grab an unlocked view of a locked matrix");
      if( BL.Height() != BR.Height() )
          LogicError("1x2 must have consistent height to combine");
      if( BL.LDim() != BR.LDim() )
          LogicError("1x2 must have consistent ldims to combine");
      if( BR.Buffer() != (BL.Buffer()+BL.LDim()*BL.Width()) )
          LogicError("1x2 must have contiguous memory");
    )
    A.Attach( BL.Height(), BL.Width()+BR.Width(), BL.Buffer(), BL.LDim() );
}

template<typename T>
void Merge1x2
( ElementalMatrix<T>& A, 
  ElementalMatrix<T>& BL,
  ElementalMatrix<T>& BR )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( A, BL, BR );
      AssertSameDists( A, BL, BR );
      AssertConforming1x2( BL, BR );
    )
    A.Attach
    ( BL.Height(), BL.Width()+BR.Width(), BL.Grid(), 
      BL.ColAlign(), BL.RowAlign(), BL.Buffer(), BL.LDim(), BL.Root() );
}

template<typename T>
void LockedMerge1x2( Matrix<T>& A, const Matrix<T>& BL, const Matrix<T>& BR )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( BL.Height() != BR.Height() )
          LogicError("1x2 must have consistent height to combine");
      if( BL.LDim() != BR.LDim() )
          LogicError("1x2 must have consistent ldims to combine");
      if( BR.LockedBuffer() != (BL.LockedBuffer()+BL.LDim()*BL.Width()) )
          LogicError("1x2 must have contiguous memory");
    )
    A.LockedAttach
    ( BL.Height(), BL.Width()+BR.Width(), BL.LockedBuffer(), BL.LDim() );
}

template<typename T>
void LockedMerge1x2
(       ElementalMatrix<T>& A,
  const ElementalMatrix<T>& BL, const ElementalMatrix<T>& BR )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( A, BL, BR );
      AssertSameDists( A, BL, BR );
      AssertConforming1x2( BL, BR );
    )
    A.LockedAttach
    ( BL.Height(), BL.Width()+BR.Width(), BL.Grid(), 
      BL.ColAlign(), BL.RowAlign(), BL.LockedBuffer(), BL.LDim(), BL.Root() );
}

// Return by value
// ^^^^^^^^^^^^^^^

template<typename T>
Matrix<T> Merge1x2( Matrix<T>& BL, Matrix<T>& BR )
{
    Matrix<T> A;
    Merge1x2( A, BL, BR );
    return A;
}

template<typename T,Dist U,Dist V>
DistMatrix<T,U,V> Merge1x2( DistMatrix<T,U,V>& BL, DistMatrix<T,U,V>& BR )
{
    DistMatrix<T,U,V> A(BL.Grid());
    Merge1x2( A, BL, BR );
    return A;
}

template<typename T>
Matrix<T> LockedMerge1x2( const Matrix<T>& BL, const Matrix<T>& BR )
{
    Matrix<T> A;
    LockedMerge1x2( A, BL, BR );
    return A;
}

template<typename T,Dist U,Dist V>
DistMatrix<T,U,V> LockedMerge1x2
( const DistMatrix<T,U,V>& BL, const DistMatrix<T,U,V>& BR )
{
    DistMatrix<T,U,V> A(BL.Grid());
    LockedMerge1x2( A, BL, BR );
    return A;
}

// Vertically merge two contiguous matrices
// ========================================

template<typename T>
void Merge2x1( Matrix<T>& A, Matrix<T>& BT, Matrix<T>& BB )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( BT.Locked() || BB.Locked() )
          LogicError("Cannot grab an unlocked view of a locked matrix");
      if( BT.Width() != BB.Width() )
          LogicError("2x1 must have consistent width to combine");
      if( BT.LDim() != BB.LDim() )
          LogicError("2x1 must have consistent ldim to combine");
      if( BB.Buffer() != (BT.Buffer() + BT.Height()) )
          LogicError("2x1 must have contiguous memory");
    )
    A.Attach( BT.Height()+BB.Height(), BT.Width(), BT.Buffer(), BT.LDim() );
}

template<typename T>
void Merge2x1
( ElementalMatrix<T>& A, 
  ElementalMatrix<T>& BT,
  ElementalMatrix<T>& BB )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( A, BT, BB );
      AssertSameDists( A, BT, BB );
      AssertConforming2x1( BT, BB );
    )
    A.Attach
    ( BT.Height()+BB.Height(), BT.Width(), BT.Grid(), 
      BT.ColAlign(), BT.RowAlign(), BT.Buffer(), BT.LDim(), BT.Root() );
}

template<typename T>
void LockedMerge2x1( Matrix<T>& A, const Matrix<T>& BT, const Matrix<T>& BB )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( BT.Width() != BB.Width() )
          LogicError("2x1 must have consistent width to combine");
      if( BT.LDim() != BB.LDim() )
          LogicError("2x1 must have consistent ldim to combine");
      if( BB.LockedBuffer() != (BT.LockedBuffer() + BT.Height()) )
          LogicError("2x1 must have contiguous memory");
    )
    A.LockedAttach
    ( BT.Height()+BB.Height(), BT.Width(), BT.LockedBuffer(), BT.LDim() );
}

template<typename T>
void LockedMerge2x1
(       ElementalMatrix<T>& A,
  const ElementalMatrix<T>& BT,
  const ElementalMatrix<T>& BB )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( A, BT, BB );
      AssertSameDists( A, BT, BB );
      AssertConforming2x1( BT, BB );
    )
    A.LockedAttach
    ( BT.Height()+BB.Height(), BT.Width(), BT.Grid(), 
      BT.ColAlign(), BT.RowAlign(), BT.LockedBuffer(), BT.LDim(), BT.Root() );
}

// Return by value
// ^^^^^^^^^^^^^^^

template<typename T>
Matrix<T> Merge2x1( Matrix<T>& BT, Matrix<T>& BB )
{
    Matrix<T> A;
    Merge2x1( A, BT, BB );
    return A;
}

template<typename T,Dist U,Dist V>
DistMatrix<T,U,V> Merge2x1( DistMatrix<T,U,V>& BT, DistMatrix<T,U,V>& BB )
{
    DistMatrix<T,U,V> A(BT.Grid());
    Merge2x1( A, BT, BB );
    return A;
}

template<typename T>
Matrix<T> LockedMerge2x1( const Matrix<T>& BT, const Matrix<T>& BB )
{
    Matrix<T> A;
    LockedMerge2x1( A, BT, BB );
    return A;
}

template<typename T,Dist U,Dist V>
DistMatrix<T,U,V> LockedMerge2x1
( const DistMatrix<T,U,V>& BT, const DistMatrix<T,U,V>& BB )
{
    DistMatrix<T,U,V> A(BT.Grid());
    LockedMerge2x1( A, BT, BB );
    return A;
}

// Merge a contiguous 2x2 block of matrices
// ========================================

template<typename T>
void Merge2x2
( Matrix<T>& A,
  Matrix<T>& BTL,
  Matrix<T>& BTR,
  Matrix<T>& BBL,
  Matrix<T>& BBR )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( BTL.Locked() || BTR.Locked() || BBL.Locked() || BBR.Locked() )
          LogicError("Cannot grab an unlocked view of a locked matrix");
      if( BTL.Width() != BBL.Width()   ||
          BTR.Width() != BBR.Width()   ||
          BTL.Height() != BTR.Height() ||
          BBL.Height() != BBR.Height()   )
          LogicError("2x2 must conform to combine");
      if( BTL.LDim() != BTR.LDim() ||
          BTR.LDim() != BBL.LDim() ||
          BBL.LDim() != BBR.LDim()   )
          LogicError("2x2 must have consistent ldims to combine");
      if( BBL.Buffer() != (BTL.Buffer() + BTL.Height()) ||
          BBR.Buffer() != (BTR.Buffer() + BTR.Height()) ||
          BTR.Buffer() != (BTL.Buffer() + BTL.LDim()*BTL.Width()) )
          LogicError("2x2 must have contiguous memory");
    )
    A.Attach
    ( BTL.Height()+BBL.Height(), BTL.Width()+BTR.Width(), 
      BTL.Buffer(), BTL.LDim() );
}

template<typename T>
void Merge2x2
( ElementalMatrix<T>& A,
  ElementalMatrix<T>& BTL,
  ElementalMatrix<T>& BTR,
  ElementalMatrix<T>& BBL,
  ElementalMatrix<T>& BBR )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( A, BTL, BTR, BBL, BBR );
      AssertSameDists( A, BTL, BTR, BBL, BBR );
      AssertConforming2x2( BTL, BTR, BBL, BBR );
    )
    A.Attach
    ( BTL.Height()+BBL.Height(), BTL.Width()+BTR.Width(), BTL.Grid(),
      BTL.ColAlign(), BTL.RowAlign(), BTL.Buffer(), BTL.LDim(), BTL.Root() );
}

template<typename T>
void LockedMerge2x2
(       Matrix<T>& A,
  const Matrix<T>& BTL,
  const Matrix<T>& BTR,
  const Matrix<T>& BBL,
  const Matrix<T>& BBR )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( BTL.Width() != BBL.Width()   ||
          BTR.Width() != BBR.Width()   ||
          BTL.Height() != BTR.Height() ||
          BBL.Height() != BBR.Height()   )
          LogicError("2x2 must conform to combine");
      if( BTL.LDim() != BTR.LDim() ||
          BTR.LDim() != BBL.LDim() ||
          BBL.LDim() != BBR.LDim()   )
          LogicError("2x2 must have consistent ldims to combine");
      if( BBL.LockedBuffer() != (BTL.LockedBuffer()+BTL.Height()) ||
          BBR.LockedBuffer() != (BTR.LockedBuffer()+BTR.Height()) ||
          BTR.LockedBuffer() != (BTL.LockedBuffer()+BTL.LDim()*BTL.Width()) )
          LogicError("2x2 must have contiguous memory");
    )
    A.LockedAttach
    ( BTL.Height()+BBL.Height(), BTL.Width()+BTR.Width(), 
      BTL.LockedBuffer(), BTL.LDim() );
}

template<typename T>
void LockedMerge2x2
(       ElementalMatrix<T>& A,
  const ElementalMatrix<T>& BTL,
  const ElementalMatrix<T>& BTR,
  const ElementalMatrix<T>& BBL,
  const ElementalMatrix<T>& BBR )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( A, BTL, BTR, BBL, BBR );
      AssertSameDists( A, BTL, BTR, BBL, BBR );
      AssertConforming2x2( BTL, BTR, BBL, BBR );
    )
    A.LockedAttach
    ( BTL.Height()+BBL.Height(), BTL.Width()+BTR.Width(), BTL.Grid(),
      BTL.ColAlign(), BTL.RowAlign(), BTL.LockedBuffer(), BTL.LDim(),
      BTL.Root() );
}

// Return by value
// ^^^^^^^^^^^^^^^

template<typename T>
Matrix<T> Merge2x2
( Matrix<T>& BTL, Matrix<T>& BTR,
  Matrix<T>& BBL, Matrix<T>& BBR )
{
    Matrix<T> A;
    Merge2x2( A, BTL, BTR, BBL, BBR );
    return A;
}

template<typename T,Dist U,Dist V>
DistMatrix<T,U,V> Merge2x2
( DistMatrix<T,U,V>& BTL, DistMatrix<T,U,V>& BTR,
  DistMatrix<T,U,V>& BBL, DistMatrix<T,U,V>& BBR )
{
    DistMatrix<T,U,V> A(BTL.Grid());
    Merge2x2( A, BTL, BTR, BBL, BBR );
    return A;
}

template<typename T>
Matrix<T> LockedMerge2x2
( const Matrix<T>& BTL, const Matrix<T>& BTR,
  const Matrix<T>& BBL, const Matrix<T>& BBR )
{
    Matrix<T> A;
    LockedMerge2x2( A, BTL, BTR, BBL, BBR );
    return A;
}

template<typename T,Dist U,Dist V>
DistMatrix<T,U,V> LockedMerge2x2
( const DistMatrix<T,U,V>& BTL, const DistMatrix<T,U,V>& BTR,
  const DistMatrix<T,U,V>& BBL, const DistMatrix<T,U,V>& BBR )
{
    DistMatrix<T,U,V> A(BTL.Grid());
    LockedMerge2x2( A, BTL, BTR, BBL, BBR );
    return A;
}

#define PROTO_DIST(T,U,V) \
  template DistMatrix<T,U,V> Merge1x2 \
  ( DistMatrix<T,U,V>& BL, DistMatrix<T,U,V>& BR ); \
  template DistMatrix<T,U,V> LockedMerge1x2 \
  ( const DistMatrix<T,U,V>& BL, const DistMatrix<T,U,V>& BR ); \
  template DistMatrix<T,U,V> Merge2x1 \
  ( DistMatrix<T,U,V>& BT, DistMatrix<T,U,V>& BB ); \
  template DistMatrix<T,U,V> LockedMerge2x1 \
  ( const DistMatrix<T,U,V>& BT, const DistMatrix<T,U,V>& BB ); \
  template DistMatrix<T,U,V> Merge2x2 \
  ( DistMatrix<T,U,V>& BTL, DistMatrix<T,U,V>& BTR, \
    DistMatrix<T,U,V>& BBL, DistMatrix<T,U,V>& BBR ); \
  template DistMatrix<T,U,V> LockedMerge2x2 \
  ( const DistMatrix<T,U,V>& BTL, const DistMatrix<T,U,V>& BTR, \
    const DistMatrix<T,U,V>& BBL, const DistMatrix<T,U,V>& BBR );

#define PROTO(T) \
  template void Merge1x2( Matrix<T>& A, Matrix<T>& BL, Matrix<T>& BR ); \
  template void Merge1x2 \
  ( ElementalMatrix<T>& A, \
    ElementalMatrix<T>& BL, ElementalMatrix<T>& BR ); \
  template void LockedMerge1x2 \
  ( Matrix<T>& A, const Matrix<T>& BL, const Matrix<T>& BR ); \
  template void LockedMerge1x2 \
  ( ElementalMatrix<T>& A, \
    const ElementalMatrix<T>& BL, const ElementalMatrix<T>& BR ); \
  template void Merge2x1( Matrix<T>& A, Matrix<T>& BT, Matrix<T>& BB ); \
  template void Merge2x1 \
  ( ElementalMatrix<T>& A, \
    ElementalMatrix<T>& BT, ElementalMatrix<T>& BB ); \
  template void LockedMerge2x1 \
  ( Matrix<T>& A, const Matrix<T>& BT, const Matrix<T>& BB ); \
  template void LockedMerge2x1 \
  ( ElementalMatrix<T>& A, \
    const ElementalMatrix<T>& BT, const ElementalMatrix<T>& BB ); \
  template void Merge2x2 \
  ( Matrix<T>& A, \
    Matrix<T>& BTL, Matrix<T>& BTR, \
    Matrix<T>& BBL, Matrix<T>& BBR ); \
  template void Merge2x2 \
  ( ElementalMatrix<T>& A, \
    ElementalMatrix<T>& BTL, ElementalMatrix<T>& BTR, \
    ElementalMatrix<T>& BBL, ElementalMatrix<T>& BBR ); \
  template void LockedMerge2x2 \
  ( Matrix<T>& A, \
    const Matrix<T>& BTL, const Matrix<T>& BTR, \
    const Matrix<T>& BBL, const Matrix<T>& BBR ); \
  template void LockedMerge2x2 \
  ( ElementalMatrix<T>& A, \
    const ElementalMatrix<T>& BTL, const ElementalMatrix<T>& BTR, \
    const ElementalMatrix<T>& BBL, const ElementalMatrix<T>& BBR ); \
  template Matrix<T> Merge1x2( Matrix<T>& BL, Matrix<T>& BR );  \
  template Matrix<T> LockedMerge1x2 \
  ( const Matrix<T>& BL, const Matrix<T>& BR ); \
  template Matrix<T> Merge2x1( Matrix<T>& BT, Matrix<T>& BB ); \
  template Matrix<T> LockedMerge2x1 \
  ( const Matrix<T>& BT, const Matrix<T>& BB ); \
  template Matrix<T> Merge2x2 \
  ( Matrix<T>& BTL, Matrix<T>& BTR, Matrix<T>& BBL, Matrix<T>& BBR ); \
  template Matrix<T> LockedMerge2x2 \
  ( const Matrix<T>& BTL, const Matrix<T>& BTR, \
    const Matrix<T>& BBL, const Matrix<T>& BBR ); \
  PROTO_DIST(T,CIRC,CIRC) \
  PROTO_DIST(T,MC,  MR  ) \
  PROTO_DIST(T,MC,  STAR) \
  PROTO_DIST(T,MD,  STAR) \
  PROTO_DIST(T,MR,  MC  ) \
  PROTO_DIST(T,MR,  STAR) \
  PROTO_DIST(T,STAR,MC  ) \
  PROTO_DIST(T,STAR,MD  ) \
  PROTO_DIST(T,STAR,MR  ) \
  PROTO_DIST(T,STAR,STAR) \
  PROTO_DIST(T,STAR,VC  ) \
  PROTO_DIST(T,STAR,VR  ) \
  PROTO_DIST(T,VC,  STAR) \
  PROTO_DIST(T,VR,  STAR)

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
