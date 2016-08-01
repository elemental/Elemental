/*
   Copyright 2009-2011, Jack Poulson.
   All rights reserved.

   Copyright 2011-2012, Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin.
   All rights reserved.

   Copyright 2013, Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright 2013-2014, Jack Poulson and The Georgia Institute of Technology.
   All rights reserved.

   Copyright 2014-2015, Jack Poulson and Stanford University.
   All rights reserved.
   
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace ldl {

template<typename T>
MatrixNode<T>::MatrixNode( MatrixNode<T>* parentNode )
: parent(parentNode), duplicateMat(nullptr), duplicateMV(nullptr)
{ }

template<typename T>
MatrixNode<T>::MatrixNode( DistMatrixNode<T>* dupNode )
: parent(nullptr), duplicateMat(dupNode), duplicateMV(nullptr)
{ }

template<typename T>
MatrixNode<T>::MatrixNode( DistMultiVecNode<T>* dupNode )
: parent(nullptr), duplicateMat(nullptr), duplicateMV(dupNode)
{ }

template<typename T>
MatrixNode<T>::MatrixNode
( const vector<Int>& invMap, const NodeInfo& info, const Matrix<T>& X )
: parent(nullptr), duplicateMat(nullptr), duplicateMV(nullptr)
{ 
    DEBUG_CSE
    Pull( invMap, info, X ); 
}

template<typename T>
MatrixNode<T>::~MatrixNode()
{
    for( auto* child : children )
        delete child;
}

template<typename T>
const MatrixNode<T>& MatrixNode<T>::operator=( const MatrixNode<T>& X )
{
    DEBUG_CSE
    matrix = X.matrix; 

    // Clean up any pre-existing children if not the right amount
    const Int numChildren = X.children.size();
    if( children.size() != X.children.size() )
    {   
        for( auto* child : children )
            delete child;
        children.resize( numChildren );
        for( Int c=0; c<numChildren; ++c )
            children[c] = new MatrixNode<T>(this);
    }
 
    for( Int c=0; c<numChildren; ++c )
        *children[c] = *X.children[c];

    return *this;
}

template<typename T>
void MatrixNode<T>::Pull
( const vector<Int>& invMap, const NodeInfo& info, const Matrix<T>& X )
{
    DEBUG_CSE
 
    const Int width = X.Width();
    matrix.Resize( info.size, width );
    for( Int t=0; t<info.size; ++t )
    {
        const Int i = invMap[info.off+t];
        for( Int j=0; j<width; ++j )
            matrix(t,j) = X(i,j);
    }

    // Clean up any pre-existing children if not the right amount
    const Int numChildren = info.children.size();
    if( children.size() != info.children.size() )
    {
        for( auto* child : children )
            delete child;
        children.resize( numChildren );
        for( Int c=0; c<numChildren; ++c )
            children[c] = new MatrixNode<T>(this);
    }

    for( Int c=0; c<numChildren; ++c )
        children[c]->Pull( invMap, *info.children[c], X );
}

template<typename T>
void MatrixNode<T>::Push
( const vector<Int>& invMap, const NodeInfo& info, Matrix<T>& X ) const
{
    DEBUG_CSE

    const Int width = matrix.Width();
    X.Resize( info.off+info.size, width );

    function<void(const MatrixNode<T>&,const NodeInfo&)> push = 
      [&]( const MatrixNode<T>& matNode, const NodeInfo& infoNode ) 
      {
          const Int numChildren = infoNode.children.size();
          for( Int c=0; c<numChildren; ++c )
              push( *matNode.children[c], *infoNode.children[c] );

          for( Int t=0; t<infoNode.size; ++t )
          {
              const Int i = invMap[infoNode.off+t];
              for( Int j=0; j<width; ++j )
                  X(i,j) = matNode.matrix(t,j);
          }
      };
    push( *this, info ); 
}

template<typename T>
Int MatrixNode<T>::Height() const
{
    DEBUG_CSE
    Int height = 0;
    function<void(const MatrixNode<T>&)> count = 
      [&]( const MatrixNode<T>& node )
      {
          for( const MatrixNode<T>* child : node.children )
              count( *child );
          height += node.matrix.Height();
      };
    count( *this );
    return height;
}

#define PROTO(T) template struct MatrixNode<T>;
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace ldl
} // namespace El
