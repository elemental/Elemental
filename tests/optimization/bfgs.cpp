/*
   Copyright (c) 2009-2016, Ryan H. Lewis
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <utility>
#include <El.hpp>
#include <El/optimization/bfgs.hpp>

using namespace El;

template< typename T>
std::pair< DistMatrix<T>, T>
LogisticRegression( const DistMatrix<T>& X, const DistMatrix<T>& y, const T lambda=1){
  const std::function< T(const DistMatrix<T>&)>
  logistic_function = [&](const DistMatrix<T>& theta){
     auto nrm = El::Norm(theta);
     El::DistMatrix<T> y(theta);
     El::DistMatrix<T> r(X.Height(), 1);
     El::Gemv(El::NORMAL, T(-1), X, theta, T(0), r);
     for( auto i = 0; i < y.Height(); ++i){ r.Set(i, 0, y.Get(i,0));  }

     T v = lambda*nrm*nrm;
     for(auto i = 0; i < y.Height(); ++i){ v += El::Log(T(1) + El::Exp(-1*y.Get(i,0)*r.Get(i,0))); }
     return v;
  };

  const std::function< DistMatrix<T>(const DistMatrix<T>&, DistMatrix<T>&)>
  logistic_gradient = [&](const DistMatrix<T>& theta, El::DistMatrix<T>& y){
      El::DistMatrix<T> r(X.Height(), 1);
      El::Gemv(El::NORMAL, T(-1), X, theta, T(0), r);
      for(auto i = 0; i < y.Height(); ++i){
         //TODO: Finish this.
         y.Set(i, 0, T(1)/(T(1) + El::Exp(-1*y.Get(i,0)*r.Get(i,0))));
      }
      return y;
  };
  DistMatrix<T> x0( X.Width(), 1);
  Gaussian( x0, X.Width(), 1);
  auto val = El::BFGS( x0, logistic_function, logistic_gradient);
  return std::make_pair( x0, val);
}

template< typename T>
std::pair< DistMatrix<T>, T>
QuadraticFunction( const Int & N){
  const std::function< T(const DistMatrix<T>&)>
  quadratic_function = [&](const DistMatrix<T>& theta){
        return .5*Dot(theta,theta);
  };

  const std::function< DistMatrix<T>(const DistMatrix<T>&, DistMatrix<T>&)>
  gradient = [&](const DistMatrix<T>& theta, El::DistMatrix<T>& y){
      y = theta;
      return y;
  };
  DistMatrix<T> x0( N, 1);
  x0.Set(0,0,5);
  auto val = El::BFGS( x0, quadratic_function,  gradient);
  return std::make_pair( x0, val);
}



template< typename T>
void TestBFGS(){
   //DistMatrix<T> X(1000, 30);
   //DistMatrix<T> y(X.Height(), 1);
   //for(int i = 0; i < y.Height(); ++i){
   //    y.Set(i, 0, 2*(i % 2)-1);
   //}
   //auto opt = LogisticRegression( X, y, T(2));
   //std::cout << opt.second << std::endl;
    for(Int i = 1; i < 2; ++i){
        std::cout << QuadraticFunction<T>(i).second << std::endl;
    }

}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );

    try
    {
        if( commRank == 0 )
            Output("Testing with doubles:");
        TestBFGS<double>();
    }
    catch( exception& e ) { ReportException(e); }
    return 0;
}