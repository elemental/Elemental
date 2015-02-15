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
#include "El.hpp"

namespace El {

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
( const vector<Int>& invMap, const SymmNodeInfo& info, const Matrix<T>& X )
: parent(nullptr), duplicateMat(nullptr), duplicateMV(nullptr)
{ 
    DEBUG_ONLY(CallStackEntry cse("MatrixNode::MatrixNode"))
    Pull( invMap, info, X ); 
}

template<typename T>
const MatrixNode<T>& MatrixNode<T>::operator=( const MatrixNode<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("MatrixNode::operator="))
    matrix = X.matrix; 
 
    const Int numChildren = X.children.size();
    children.resize( numChildren );
    for( Int c=0; c<numChildren; ++c )
    {
        children[c] = new MatrixNode<T>(this);
        *children[c] = *X.children[c];
    }
    return *this;
}

template<typename T>
void MatrixNode<T>::Pull
( const vector<Int>& invMap, const SymmNodeInfo& info, const Matrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("MatrixNode::Pull"))
 
    // Clean up any pre-existing children
    for( MatrixNode<T>* child : children )
        delete child;

    const Int numChildren = info.children.size();
    children.resize( numChildren );
    for( Int c=0; c<numChildren; ++c )
    {
        children[c] = new MatrixNode<T>(this);
        children[c]->Pull( invMap, *info.children[c], X );
    }

    const Int width = X.Width();
    matrix.Resize( info.size, width );
    for( Int t=0; t<info.size; ++t )
    {
        const Int i = invMap[info.off+t];
        for( Int j=0; j<width; ++j )
            matrix.Set( t, j, X.Get(i,j) );
    }
}

template<typename T>
void MatrixNode<T>::Push
( const vector<Int>& invMap, const SymmNodeInfo& info, Matrix<T>& X ) const
{
    DEBUG_ONLY(CallStackEntry cse("MatrixNode::Push"))

    const Int width = matrix.Width();
    X.Resize( info.off+info.size, width );

    function<void(const MatrixNode<T>&,const SymmNodeInfo&)> push = 
      [&]( const MatrixNode<T>& matNode, const SymmNodeInfo& infoNode ) 
      {
          const Int numChildren = infoNode.children.size();
          for( Int c=0; c<numChildren; ++c )
              push( *matNode.children[c], *infoNode.children[c] );

          for( Int t=0; t<infoNode.size; ++t )
          {
              const Int i = invMap[infoNode.off+t];
              for( Int j=0; j<width; ++j )
                  X.Set( i, j, matNode.matrix.Get(t,j) );
          }
      };
    push( *this, info ); 
}

template<typename T>
Int MatrixNode<T>::Height() const
{
    DEBUG_ONLY(CallStackEntry cse("MatrixNode::LocalHeight"))    
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

#define PROTO(T) template class MatrixNode<T>;
#include "El/macros/Instantiate.h"

} // namespace El
