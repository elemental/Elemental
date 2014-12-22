/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T,Dist U,Dist V>
GeneralDistMatrix<T,U,V>::GeneralDistMatrix( const El::Grid& grid, Int root )
: AbstractDistMatrix<T>(grid,root)
{ }

template<typename T,Dist U,Dist V>
GeneralDistMatrix<T,U,V>::GeneralDistMatrix( GeneralDistMatrix<T,U,V>&& A ) 
EL_NOEXCEPT
: AbstractDistMatrix<T>(std::move(A))
{ }

// Assignment and reconfiguration
// ==============================

template<typename T,Dist U,Dist V>
GeneralDistMatrix<T,U,V>& 
GeneralDistMatrix<T,U,V>::operator=( GeneralDistMatrix<T,U,V>&& A )
{
    AbstractDistMatrix<T>::operator=( std::move(A) );
    return *this;
}

// Basic queries
// =============
// Distribution information
// ------------------------
template<typename T,Dist U,Dist V>
Dist GeneralDistMatrix<T,U,V>::ColDist() const { return U; }
template<typename T,Dist U,Dist V>
Dist GeneralDistMatrix<T,U,V>::RowDist() const { return V; }
template<typename T,Dist U,Dist V>
Dist GeneralDistMatrix<T,U,V>::PartialColDist() const 
{ return Partial<U>(); }
template<typename T,Dist U,Dist V>
Dist GeneralDistMatrix<T,U,V>::PartialRowDist() const 
{ return Partial<V>(); }
template<typename T,Dist U,Dist V>
Dist GeneralDistMatrix<T,U,V>::PartialUnionColDist() const 
{ return PartialUnionCol<U,V>(); }
template<typename T,Dist U,Dist V>
Dist GeneralDistMatrix<T,U,V>::PartialUnionRowDist() const 
{ return PartialUnionRow<U,V>(); }
template<typename T,Dist U,Dist V>
Dist GeneralDistMatrix<T,U,V>::CollectedColDist() const
{ return Collect<U>(); }
template<typename T,Dist U,Dist V>
Dist GeneralDistMatrix<T,U,V>::CollectedRowDist() const
{ return Collect<V>(); }

// Instantiations for {Int,Real,Complex<Real>} for each Real in {float,double}
// ###########################################################################

#define DISTPROTO(T,U,V) template class GeneralDistMatrix<T,U,V>
  
#define PROTO(T)\
  DISTPROTO(T,CIRC,CIRC);\
  DISTPROTO(T,MC,  MR  );\
  DISTPROTO(T,MC,  STAR);\
  DISTPROTO(T,MD,  STAR);\
  DISTPROTO(T,MR,  MC  );\
  DISTPROTO(T,MR,  STAR);\
  DISTPROTO(T,STAR,MC  );\
  DISTPROTO(T,STAR,MD  );\
  DISTPROTO(T,STAR,MR  );\
  DISTPROTO(T,STAR,STAR);\
  DISTPROTO(T,STAR,VC  );\
  DISTPROTO(T,STAR,VR  );\
  DISTPROTO(T,VC,  STAR);\
  DISTPROTO(T,VR,  STAR);

#include "El/macros/Instantiate.h"

} // namespace El
