/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CORE_PROXY_HPP
#define EL_CORE_PROXY_HPP

namespace El {

struct ProxyCtrl 
{
    bool colConstrain, rowConstrain, rootConstrain;
    Int colAlign, rowAlign, root;
    Int blockHeight, blockWidth;
    Int colCut, rowCut;

    ProxyCtrl() 
    : colConstrain(false), rowConstrain(false), rootConstrain(false),
      colAlign(0), rowAlign(0), root(0), 
      blockHeight(DefaultBlockHeight()), blockWidth(DefaultBlockWidth()),
      colCut(0), rowCut(0)
    { }
};

struct ElementalProxyCtrl
{
    bool colConstrain, rowConstrain, rootConstrain;
    Int colAlign, rowAlign, root;

    ElementalProxyCtrl() 
    : colConstrain(false), rowConstrain(false), rootConstrain(false),
      colAlign(0), rowAlign(0), root(0) 
    { }
};

// TODO: Detailed description of the allowed (S,T) pairings

// Read proxy
// ==========

// Sequential
// ----------
template<typename T,typename S>
shared_ptr<const Matrix<T>> ReadProxy( const Matrix<S>* A );
template<typename T,typename S>
shared_ptr<Matrix<T>> ReadProxy( Matrix<S>* A );

// Distributed
// -----------
template<typename T,Dist U,Dist V,typename S>
shared_ptr<const DistMatrix<T,U,V>>
ReadProxy
( const AbstractDistMatrix<S>* A,
  const ElementalProxyCtrl& ctrl=ElementalProxyCtrl() );
template<typename T,Dist U,Dist V,typename S>
shared_ptr<DistMatrix<T,U,V>>
ReadProxy
( AbstractDistMatrix<S>* A,
  const ElementalProxyCtrl& ctrl=ElementalProxyCtrl() );

// TODO: General DistMatrix version

// Read-write proxy
// ================

// Sequential
// ----------
template<typename T,typename S>
shared_ptr<Matrix<T>> ReadWriteProxy( Matrix<S>* A );

// Distributed
// -----------
template<typename T,Dist U,Dist V,typename S>
shared_ptr<DistMatrix<T,U,V>> ReadWriteProxy
( AbstractDistMatrix<S>* A,
  const ElementalProxyCtrl& ctrl=ElementalProxyCtrl() );

// Write proxy
// ===========

// Sequential
// ----------
template<typename T,typename S>
shared_ptr<Matrix<T>> WriteProxy( Matrix<S>* A );

// Distributed
// -----------
template<typename T,Dist U,Dist V,typename S>
shared_ptr<DistMatrix<T,U,V>> WriteProxy
( AbstractDistMatrix<S>* A,
  const ElementalProxyCtrl& ctrl=ElementalProxyCtrl() );

// TODO: General implementation

} // namespace El

#endif // ifndef EL_CORE_PROXY_HPP
