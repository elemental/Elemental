/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CORE_PROXY_HPP
#define EL_CORE_PROXY_HPP

// TODO: Split into several files

namespace El {

template<typename S,typename T,typename=EnableIf<CanCast<S,T>>>
class MatrixReadProxy 
{
private:
    Matrix<T>* prox_;

public:
    MatrixReadProxy( const Matrix<S>& A )
    { 
        prox_ = new Matrix<T>;
        Copy( A, *prox_ );
    }

    ~MatrixReadProxy() { delete prox_; }

    const Matrix<T>& GetLocked() const { return *prox_; }
          Matrix<T>& Get()             { return *prox_; }
};

template<typename T>
class MatrixReadProxy<T,T,void>
{
private:
    bool locked_;
    Matrix<T>& orig_;

public:
    MatrixReadProxy( const Matrix<T>& A )
    : locked_(true),
      orig_(const_cast<Matrix<T>&>(A))
    { }

    MatrixReadProxy( Matrix<T>& A )
    : locked_(false),
      orig_(A)
    { }

    ~MatrixReadProxy() { }

    const Matrix<T>& GetLocked() const { return orig_; }

    Matrix<T>& Get()
    {
        if( locked_ )
            LogicError("Attempted to extract mutable from immutable");
        return orig_;
    }
};

template<typename S,typename T,typename=EnableIf<CanCast<T,S>>>
class MatrixWriteProxy 
{
private:
    Matrix<S>& orig_;
    Matrix<T>* prox_;

public:
    MatrixWriteProxy( Matrix<S>& A )
    : orig_(A)
    { prox_ = new Matrix<T>( A.Height(), A.Width() ); }

    ~MatrixWriteProxy() 
    { 
        if( !uncaught_exception() )
            Copy( *prox_, orig_ );
        delete prox_;
    }

    const Matrix<T>& GetLocked() const { return *prox_; }
          Matrix<T>& Get()             { return *prox_; }
};

template<typename T>
class MatrixWriteProxy<T,T,void>
{
private:
    Matrix<T>& orig_;

public:
    MatrixWriteProxy( Matrix<T>& A )
    : orig_(A)
    { }

    ~MatrixWriteProxy() { }

    const Matrix<T>& GetLocked() const { return orig_; }
          Matrix<T>& Get()             { return orig_; }
};

template<typename S,typename T,typename=EnableIf<CanBidirectionalCast<S,T>>>
class MatrixReadWriteProxy 
{
private:
    Matrix<S>& orig_;
    Matrix<T>* prox_;

public:
    MatrixReadWriteProxy( Matrix<S>& A )
    : orig_(A)
    { prox_ = new Matrix<T>(A); }

    ~MatrixReadWriteProxy() 
    { 
        if( !uncaught_exception() )
            Copy( *prox_, orig_ );
        delete prox_;
    }

    const Matrix<T>& GetLocked() const { return *prox_; }
          Matrix<T>& Get()             { return *prox_; }
};

template<typename T>
class MatrixReadWriteProxy<T,T,void>
{
private:
    Matrix<T>& orig_;

public:
    MatrixReadWriteProxy( Matrix<T>& A )
    : orig_(A)
    { }

    ~MatrixReadWriteProxy() { }

    const Matrix<T>& GetLocked() const { return orig_; }
          Matrix<T>& Get()             { return orig_; }
};

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

template<typename S,typename T,Dist U=MC,Dist V=MR,DistWrap wrap=ELEMENT,
         typename=EnableIf<CanCast<S,T>>>
class DistMatrixReadProxy; 

template<typename S,typename T,Dist U,Dist V>
class DistMatrixReadProxy<S,T,U,V,ELEMENT,void>
{
private:
    typedef DistMatrix<T,U,V,ELEMENT> proxType;

    bool locked_;
    proxType* prox_;

public:
    // TODO: Add try-catch statements into constructors?

    DistMatrixReadProxy
    ( const AbstractDistMatrix<S>& A,
      const ElementalProxyCtrl& ctrl=ElementalProxyCtrl() )
    { 
        prox_ = new proxType(A.Grid());
        if( ctrl.rootConstrain )
            prox_->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            prox_->AlignCols( ctrl.colAlign );    
        if( ctrl.rowConstrain )
            prox_->AlignRows( ctrl.rowAlign );
        Copy( A, *prox_ );
    }

    DistMatrixReadProxy
    ( AbstractDistMatrix<S>& A,
      const ElementalProxyCtrl& ctrl=ElementalProxyCtrl() )
    { 
        prox_ = new proxType(A.Grid());
        if( ctrl.rootConstrain )
            prox_->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            prox_->AlignCols( ctrl.colAlign );    
        if( ctrl.rowConstrain )
            prox_->AlignRows( ctrl.rowAlign );
        Copy( A, *prox_ );
    }

    ~DistMatrixReadProxy() { delete prox_; }

    const proxType& GetLocked() const { return *prox_; }
          proxType& Get()             { return *prox_; }
};

template<typename T,Dist U,Dist V>
class DistMatrixReadProxy<T,T,U,V,ELEMENT,void>
{
private:
    typedef DistMatrix<T,U,V,ELEMENT> proxType;

    bool locked_;
    bool madeCopy_;
    proxType* prox_;

public:
    // TODO: Add try-catch statements into constructors?

    DistMatrixReadProxy
    ( const AbstractDistMatrix<T>& A,
      const ElementalProxyCtrl& ctrl=ElementalProxyCtrl() )
    { 
        if( A.ColDist() == U && A.RowDist() == V && A.Wrap() == ELEMENT )
        {
            const bool colMisalign = 
              ( ctrl.colConstrain && A.ColAlign() != ctrl.colAlign );
            const bool rowMisalign = 
              ( ctrl.rowConstrain && A.RowAlign() != ctrl.rowAlign );
            const bool rootMisalign = 
              ( ctrl.rootConstrain && A.Root() != ctrl.root );
            if( !colMisalign && !rowMisalign && !rootMisalign )
            {
                locked_ = true;
                madeCopy_ = false;
                auto ACast = static_cast<const proxType*>(&A);
                prox_ = const_cast<proxType*>(ACast);
                return;
            }
        }
        locked_ = false;
        madeCopy_ = true;
        prox_ = new proxType(A.Grid());
        if( ctrl.rootConstrain )
            prox_->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            prox_->AlignCols( ctrl.colAlign );    
        if( ctrl.rowConstrain )
            prox_->AlignRows( ctrl.rowAlign );
        Copy( A, *prox_ );
    }

    DistMatrixReadProxy
    ( AbstractDistMatrix<T>& A,
      const ElementalProxyCtrl& ctrl=ElementalProxyCtrl() )
    : locked_(false)
    { 
        if( A.ColDist() == U && A.RowDist() == V && A.Wrap() == ELEMENT )
        {
            const bool colMisalign = 
              ( ctrl.colConstrain && A.ColAlign() != ctrl.colAlign );
            const bool rowMisalign = 
              ( ctrl.rowConstrain && A.RowAlign() != ctrl.rowAlign );
            const bool rootMisalign = 
              ( ctrl.rootConstrain && A.Root() != ctrl.root );
            if( !colMisalign && !rowMisalign && !rootMisalign )
            {
                madeCopy_ = false;
                prox_ = static_cast<proxType*>(&A);
                if( ctrl.rootConstrain )
                    prox_->SetRoot( ctrl.root );
                if( ctrl.colConstrain )
                    prox_->AlignCols( ctrl.colAlign );    
                if( ctrl.rowConstrain )
                    prox_->AlignRows( ctrl.rowAlign );
                return;
            }
        }
        madeCopy_ = true;
        prox_ = new proxType(A.Grid());
        if( ctrl.rootConstrain )
            prox_->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            prox_->AlignCols( ctrl.colAlign );    
        if( ctrl.rowConstrain )
            prox_->AlignRows( ctrl.rowAlign );
        Copy( A, *prox_ );
    }

    ~DistMatrixReadProxy() 
    { 
        if( madeCopy_ )
            delete prox_;
    }

    const proxType& GetLocked() const { return *prox_; }

    proxType& Get()
    {
        if( locked_ )
            LogicError("Attempted to extract mutable from immutable");
        return *prox_;
    }
};

template<typename S,typename T,Dist U,Dist V>
class DistMatrixReadProxy<S,T,U,V,BLOCK,void>
{
private:
    typedef DistMatrix<T,U,V,BLOCK> proxType;
    proxType* prox_;

public:
    // TODO: Add try-catch statements into constructors?

    DistMatrixReadProxy
    ( const AbstractDistMatrix<S>& A,
      const ProxyCtrl& ctrl=ProxyCtrl() )
    { 
        prox_ = new proxType(A.Grid());
        if( ctrl.rootConstrain )
            prox_->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            prox_->AlignCols( ctrl.blockHeight, ctrl.colAlign, ctrl.colCut );
        if( ctrl.rowConstrain )
            prox_->AlignRows( ctrl.blockWidth, ctrl.rowAlign, ctrl.rowCut );
        Copy( A, *prox_ );
    }

    DistMatrixReadProxy
    ( AbstractDistMatrix<S>& A,
      const ProxyCtrl& ctrl=ProxyCtrl() )
    { 
        prox_ = new proxType(A.Grid());
        if( ctrl.rootConstrain )
            prox_->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            prox_->AlignCols( ctrl.blockHeight, ctrl.colAlign, ctrl.colCut );
        if( ctrl.rowConstrain )
            prox_->AlignRows( ctrl.blockWidth, ctrl.rowAlign, ctrl.rowCut );
        Copy( A, *prox_ );
    }

    ~DistMatrixReadProxy() { delete prox_; }

    const proxType& GetLocked() const { return *prox_; }
          proxType& Get()             { return *prox_; }
};

template<typename T,Dist U,Dist V>
class DistMatrixReadProxy<T,T,U,V,BLOCK,void>
{
private:
    typedef DistMatrix<T,U,V,BLOCK> proxType;

    bool locked_;
    bool madeCopy_;
    proxType* prox_;

public:
    // TODO: Add try-catch statements into constructors?

    DistMatrixReadProxy
    ( const AbstractDistMatrix<T>& A,
      const ProxyCtrl& ctrl=ProxyCtrl() )
    { 
        if( A.ColDist() == U && A.RowDist() == V && A.Wrap() == BLOCK )
        {
            const bool colMisalign = 
              ( ctrl.colConstrain && 
                (A.ColAlign() != ctrl.colAlign ||
                 A.BlockHeight() != ctrl.blockHeight ||
                 A.ColCut() != ctrl.colCut) );
            const bool rowMisalign = 
              ( ctrl.rowConstrain && 
                (A.RowAlign() != ctrl.rowAlign ||
                 A.BlockWidth() != ctrl.blockWidth ||
                 A.RowCut() != ctrl.rowCut) );
            const bool rootMisalign = 
              ( ctrl.rootConstrain && A.Root() != ctrl.root );
            if( !colMisalign && !rowMisalign && !rootMisalign )
            {
                locked_ = true;
                madeCopy_ = false;
                auto ACast = static_cast<const proxType*>(&A);
                prox_ = const_cast<proxType*>(ACast);
                return;
            }
        }
        locked_ = false;
        madeCopy_ = true;
        prox_ = new proxType(A.Grid());
        if( ctrl.rootConstrain )
            prox_->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            prox_->AlignCols( ctrl.blockHeight, ctrl.colAlign, ctrl.colCut );
        if( ctrl.rowConstrain )
            prox_->AlignRows( ctrl.blockWidth, ctrl.rowAlign, ctrl.rowCut );
        Copy( A, *prox_ );
    }

    DistMatrixReadProxy
    ( AbstractDistMatrix<T>& A,
      const ProxyCtrl& ctrl=ProxyCtrl() )
    : locked_(false)
    { 
        if( A.ColDist() == U && A.RowDist() == V && A.Wrap() == BLOCK )
        {
            const bool colMisalign = 
              ( ctrl.colConstrain && 
                (A.ColAlign() != ctrl.colAlign ||
                 A.BlockHeight() != ctrl.blockHeight ||
                 A.ColCut() != ctrl.colCut) );
            const bool rowMisalign = 
              ( ctrl.rowConstrain && 
                (A.RowAlign() != ctrl.rowAlign ||
                 A.BlockWidth() != ctrl.blockWidth ||
                 A.RowCut() != ctrl.rowCut) );
            const bool rootMisalign = 
              ( ctrl.rootConstrain && A.Root() != ctrl.root );
            if( !colMisalign && !rowMisalign && !rootMisalign )
            {
                madeCopy_ = false;
                prox_ = static_cast<proxType*>(&A);
                if( ctrl.rootConstrain )
                    prox_->SetRoot( ctrl.root );
                if( ctrl.colConstrain )
                    prox_->AlignCols
                    ( ctrl.blockHeight, ctrl.colAlign, ctrl.colCut ); 
                if( ctrl.rowConstrain )
                    prox_->AlignRows
                    ( ctrl.blockWidth, ctrl.rowAlign, ctrl.rowCut );
                return;
            }
        }
        madeCopy_ = true;
        prox_ = new proxType(A.Grid());
        if( ctrl.rootConstrain )
            prox_->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            prox_->AlignCols( ctrl.blockHeight, ctrl.colAlign, ctrl.colCut );
        if( ctrl.rowConstrain )
            prox_->AlignRows( ctrl.blockWidth, ctrl.rowAlign, ctrl.rowCut );
        Copy( A, *prox_ );
    }

    ~DistMatrixReadProxy() 
    { 
        if( madeCopy_ )
            delete prox_;
    }

    const proxType& GetLocked() const { return *prox_; }

    proxType& Get()
    {
        if( locked_ )
            LogicError("Attempted to extract mutable from immutable");
        return *prox_;
    }
};

template<typename S,typename T,Dist U=MC,Dist V=MR,DistWrap wrap=ELEMENT,
         typename=EnableIf<CanCast<T,S>>>
class DistMatrixWriteProxy;

template<typename S,typename T,Dist U,Dist V>
class DistMatrixWriteProxy<S,T,U,V,ELEMENT,void>
{
private:
    typedef DistMatrix<T,U,V,ELEMENT> proxType;

    AbstractDistMatrix<S>& orig_;
    proxType* prox_;

public:
    // TODO: Add try-catch statements into constructors?
    DistMatrixWriteProxy
    ( AbstractDistMatrix<S>& A,
      const ElementalProxyCtrl& ctrl=ElementalProxyCtrl() )
    : orig_(A)
    { 
        prox_ = new proxType(A.Grid());
        if( ctrl.rootConstrain )
            prox_->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            prox_->AlignCols( ctrl.colAlign );    
        if( ctrl.rowConstrain )
            prox_->AlignRows( ctrl.rowAlign );
        prox_->Resize( A.Height(), A.Width() );
    }

    ~DistMatrixWriteProxy() 
    { 
        if( !uncaught_exception() ) 
            Copy( *prox_, orig_ );
        delete prox_;
    }

    const proxType& GetLocked() const { return *prox_; }
          proxType& Get()             { return *prox_; }
};

template<typename T,Dist U,Dist V>
class DistMatrixWriteProxy<T,T,U,V,ELEMENT,void>
{
private:
    typedef DistMatrix<T,U,V,ELEMENT> proxType;

    bool madeCopy_;
    AbstractDistMatrix<T>& orig_;
    proxType* prox_;

public:
    // TODO: Add try-catch statements into constructors?
    DistMatrixWriteProxy
    ( AbstractDistMatrix<T>& A,
      const ElementalProxyCtrl& ctrl=ElementalProxyCtrl() )
    : orig_(A)
    { 
        if( A.ColDist() == U && A.RowDist() == V && A.Wrap() == ELEMENT )
        {
            const bool colMisalign = 
              ( ctrl.colConstrain && A.ColAlign() != ctrl.colAlign );
            const bool rowMisalign = 
              ( ctrl.rowConstrain && A.RowAlign() != ctrl.rowAlign );
            const bool rootMisalign = 
              ( ctrl.rootConstrain && A.Root() != ctrl.root );
            if( !colMisalign && !rowMisalign && !rootMisalign )
            {
                madeCopy_ = false;
                prox_ = static_cast<proxType*>(&A);
                if( ctrl.rootConstrain )
                    prox_->SetRoot( ctrl.root );
                if( ctrl.colConstrain )
                    prox_->AlignCols( ctrl.colAlign );    
                if( ctrl.rowConstrain )
                    prox_->AlignRows( ctrl.rowAlign );
                return;
            }
        }
        madeCopy_ = true;
        prox_ = new proxType(A.Grid());
        if( ctrl.rootConstrain )
            prox_->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            prox_->AlignCols( ctrl.colAlign );    
        if( ctrl.rowConstrain )
            prox_->AlignRows( ctrl.rowAlign );
        prox_->Resize( A.Height(), A.Width() );
    }

    ~DistMatrixWriteProxy() 
    { 
        if( madeCopy_ )
        {
            if( !uncaught_exception() )
                Copy( *prox_, orig_ );
            delete prox_;
        }
    }

    const proxType& GetLocked() const { return *prox_; }
          proxType& Get()             { return *prox_; }
};

template<typename S,typename T,Dist U,Dist V>
class DistMatrixWriteProxy<S,T,U,V,BLOCK,void>
{
private:
    typedef DistMatrix<T,U,V,BLOCK> proxType;

    AbstractDistMatrix<S>& orig_;
    proxType* prox_;

public:
    // TODO: Add try-catch statements into constructors?
    DistMatrixWriteProxy
    ( AbstractDistMatrix<S>& A,
      const ProxyCtrl& ctrl=ProxyCtrl() )
    : orig_(A)
    { 
        prox_ = new proxType(A.Grid());
        if( ctrl.rootConstrain )
            prox_->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            prox_->AlignCols( ctrl.blockHeight, ctrl.colAlign, ctrl.colCut );
        if( ctrl.rowConstrain )
            prox_->AlignRows( ctrl.blockWidth, ctrl.rowAlign, ctrl.rowCut );
        prox_->Resize( A.Height(), A.Width() );
    }

    ~DistMatrixWriteProxy() 
    { 
        if( !uncaught_exception() ) 
            Copy( *prox_, orig_ );
        delete prox_;
    }

    const proxType& GetLocked() const { return *prox_; }
          proxType& Get()             { return *prox_; }
};

template<typename T,Dist U,Dist V>
class DistMatrixWriteProxy<T,T,U,V,BLOCK,void>
{
private:
    typedef DistMatrix<T,U,V,BLOCK> proxType;

    bool madeCopy_;
    AbstractDistMatrix<T>& orig_;
    proxType* prox_;

public:
    // TODO: Add try-catch statements into constructors?
    DistMatrixWriteProxy
    ( AbstractDistMatrix<T>& A,
      const ProxyCtrl& ctrl=ProxyCtrl() )
    : orig_(A)
    { 
        if( A.ColDist() == U && A.RowDist() == V && A.Wrap() == BLOCK )
        {
            const bool colMisalign = 
              ( ctrl.colConstrain && 
                (A.ColAlign() != ctrl.colAlign ||
                 A.BlockHeight() != ctrl.blockHeight ||
                 A.ColCut() != ctrl.colCut) );
            const bool rowMisalign = 
              ( ctrl.rowConstrain && 
                (A.RowAlign() != ctrl.rowAlign ||
                 A.BlockWidth() != ctrl.blockWidth ||
                 A.RowCut() != ctrl.rowCut) );
            const bool rootMisalign = 
              ( ctrl.rootConstrain && A.Root() != ctrl.root );
            if( !colMisalign && !rowMisalign && !rootMisalign )
            {
                madeCopy_ = false;
                prox_ = static_cast<proxType*>(&A);
                if( ctrl.rootConstrain )
                    prox_->SetRoot( ctrl.root );
                if( ctrl.colConstrain )
                    prox_->AlignCols
                    ( ctrl.blockHeight, ctrl.colAlign, ctrl.colCut ); 
                if( ctrl.rowConstrain )
                    prox_->AlignRows
                    ( ctrl.blockWidth, ctrl.rowAlign, ctrl.rowCut );
                return;
            }
        }
        madeCopy_ = true;
        prox_ = new proxType(A.Grid());
        if( ctrl.rootConstrain )
            prox_->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            prox_->AlignCols( ctrl.blockHeight, ctrl.colAlign, ctrl.colCut );
        if( ctrl.rowConstrain )
            prox_->AlignRows( ctrl.blockWidth, ctrl.rowAlign, ctrl.rowCut );
        prox_->Resize( A.Height(), A.Width() );
    }

    ~DistMatrixWriteProxy() 
    { 
        if( madeCopy_ )
        {
            if( !uncaught_exception() )
                Copy( *prox_, orig_ );
            delete prox_;
        }
    }

    const proxType& GetLocked() const { return *prox_; }
          proxType& Get()             { return *prox_; }
};

template<typename S,typename T,Dist U=MC,Dist V=MR,DistWrap wrap=ELEMENT,
         typename=EnableIf<CanBidirectionalCast<S,T>>>
class DistMatrixReadWriteProxy; 

template<typename S,typename T,Dist U,Dist V>
class DistMatrixReadWriteProxy<S,T,U,V,ELEMENT,void>
{
private:
    typedef DistMatrix<T,U,V,ELEMENT> proxType;

    AbstractDistMatrix<S>& orig_;
    proxType* prox_;

public:
    // TODO: Add try-catch statements into constructors?
    DistMatrixReadWriteProxy
    ( AbstractDistMatrix<S>& A,
      const ElementalProxyCtrl& ctrl=ElementalProxyCtrl() )
    : orig_(A)
    { 
        prox_ = new proxType(A.Grid());
        if( ctrl.rootConstrain )
            prox_->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            prox_->AlignCols( ctrl.colAlign );    
        if( ctrl.rowConstrain )
            prox_->AlignRows( ctrl.rowAlign );
        Copy( A, *prox_ );
    }

    ~DistMatrixReadWriteProxy() 
    { 
        if( !uncaught_exception() )
            Copy( *prox_, orig_ );
        delete prox_;
    }

    const proxType& GetLocked() const { return *prox_; }
          proxType& Get()             { return *prox_; }
};

template<typename T,Dist U,Dist V>
class DistMatrixReadWriteProxy<T,T,U,V,ELEMENT,void>
{
private:
    typedef DistMatrix<T,U,V,ELEMENT> proxType;

    bool madeCopy_;
    AbstractDistMatrix<T>& orig_;
    proxType* prox_;

public:
    // TODO: Add try-catch statements into constructors?
    DistMatrixReadWriteProxy
    ( AbstractDistMatrix<T>& A,
      const ElementalProxyCtrl& ctrl=ElementalProxyCtrl() )
    : orig_(A)
    { 
        if( A.ColDist() == U && A.RowDist() == V && A.Wrap() == ELEMENT )
        {
            const bool colMisalign = 
              ( ctrl.colConstrain && A.ColAlign() != ctrl.colAlign );
            const bool rowMisalign = 
              ( ctrl.rowConstrain && A.RowAlign() != ctrl.rowAlign );
            const bool rootMisalign = 
              ( ctrl.rootConstrain && A.Root() != ctrl.root );
            if( !colMisalign && !rowMisalign && !rootMisalign )
            {
                madeCopy_ = false;
                prox_ = static_cast<proxType*>(&A);
                if( ctrl.rootConstrain )
                    prox_->SetRoot( ctrl.root );
                if( ctrl.colConstrain )
                    prox_->AlignCols( ctrl.colAlign );    
                if( ctrl.rowConstrain )
                    prox_->AlignRows( ctrl.rowAlign );
                return;
            }
        }
        madeCopy_ = true;
        prox_ = new proxType(A.Grid());
        if( ctrl.rootConstrain )
            prox_->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            prox_->AlignCols( ctrl.colAlign );    
        if( ctrl.rowConstrain )
            prox_->AlignRows( ctrl.rowAlign );
        Copy( A, *prox_ );
    }

    ~DistMatrixReadWriteProxy() 
    { 
        if( madeCopy_ )
        {
            if( !uncaught_exception() )
                Copy( *prox_, orig_ );
            delete prox_;
        }
    }

    const proxType& GetLocked() const { return *prox_; }
          proxType& Get()             { return *prox_; }
};

template<typename S,typename T,Dist U,Dist V>
class DistMatrixReadWriteProxy<S,T,U,V,BLOCK,void>
{
private:
    typedef DistMatrix<T,U,V,BLOCK> proxType;

    AbstractDistMatrix<S>& orig_;
    proxType* prox_;

public:
    // TODO: Add try-catch statements into constructors?
    DistMatrixReadWriteProxy
    ( AbstractDistMatrix<S>& A,
      const ProxyCtrl& ctrl=ProxyCtrl() )
    : orig_(A)
    { 
        prox_ = new proxType(A.Grid());
        if( ctrl.rootConstrain )
            prox_->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            prox_->AlignCols( ctrl.blockHeight, ctrl.colAlign, ctrl.colCut );
        if( ctrl.rowConstrain )
            prox_->AlignRows( ctrl.blockWidth, ctrl.rowAlign, ctrl.rowCut );
        Copy( A, *prox_ );
    }

    ~DistMatrixReadWriteProxy() 
    { 
        if( !uncaught_exception() )
            Copy( *prox_, orig_ );
        delete prox_;
    }

    const proxType& GetLocked() const { return *prox_; }
          proxType& Get()             { return *prox_; }
};

template<typename T,Dist U,Dist V>
class DistMatrixReadWriteProxy<T,T,U,V,BLOCK,void>
{
private:
    typedef DistMatrix<T,U,V,BLOCK> proxType;

    bool madeCopy_;
    AbstractDistMatrix<T>& orig_;
    proxType* prox_;

public:
    // TODO: Add try-catch statements into constructors?
    DistMatrixReadWriteProxy
    ( AbstractDistMatrix<T>& A,
      const ProxyCtrl& ctrl=ProxyCtrl() )
    : orig_(A)
    { 
        if( A.ColDist() == U && A.RowDist() == V && A.Wrap() == BLOCK )
        {
            const bool colMisalign = 
              ( ctrl.colConstrain && 
                (A.ColAlign() != ctrl.colAlign ||
                 A.BlockHeight() != ctrl.blockHeight ||
                 A.ColCut() != ctrl.colCut) );
            const bool rowMisalign = 
              ( ctrl.rowConstrain && 
                (A.RowAlign() != ctrl.rowAlign ||
                 A.BlockWidth() != ctrl.blockWidth ||
                 A.RowCut() != ctrl.rowCut) );
            const bool rootMisalign = 
              ( ctrl.rootConstrain && A.Root() != ctrl.root );
            if( !colMisalign && !rowMisalign && !rootMisalign )
            {
                madeCopy_ = false;
                prox_ = static_cast<proxType*>(&A);
                if( ctrl.rootConstrain )
                    prox_->SetRoot( ctrl.root );
                if( ctrl.colConstrain )
                    prox_->AlignCols
                    ( ctrl.blockHeight, ctrl.colAlign, ctrl.colCut ); 
                if( ctrl.rowConstrain )
                    prox_->AlignRows
                    ( ctrl.blockWidth, ctrl.rowAlign, ctrl.rowCut );
                return;
            }
        }
        madeCopy_ = true;
        prox_ = new proxType(A.Grid());
        if( ctrl.rootConstrain )
            prox_->SetRoot( ctrl.root );
        if( ctrl.colConstrain )
            prox_->AlignCols( ctrl.blockHeight, ctrl.colAlign, ctrl.colCut );
        if( ctrl.rowConstrain )
            prox_->AlignRows( ctrl.blockWidth, ctrl.rowAlign, ctrl.rowCut );
        Copy( A, *prox_ );
    }

    ~DistMatrixReadWriteProxy() 
    { 
        if( madeCopy_ )
        {
            if( !uncaught_exception() )
                Copy( *prox_, orig_ );
            delete prox_;
        }
    }

    const proxType& GetLocked() const { return *prox_; }
          proxType& Get()             { return *prox_; }
};

} // namespace El

#endif // ifndef EL_CORE_PROXY_HPP
