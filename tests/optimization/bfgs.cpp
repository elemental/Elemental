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
#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

using namespace El;

template< typename T>
std::pair< DistMatrix<T>, T>
SimpleQuadraticBFGSTest( const Int & N){
  const std::function< T(const DistMatrix<T>&)>
  quadratic_function = [&](const DistMatrix<T>& theta){
        return Dot(theta,theta)/T(2);
  };

  const std::function< DistMatrix<T>(const DistMatrix<T>&, DistMatrix<T>&)>
  gradient = [&](const DistMatrix<T>& theta, El::DistMatrix<T>& y){
      y = theta;
      return y;
  };
  DistMatrix<T> x0( N, 1);
  Ones(x0, N, 1);
  auto val = El::BFGS( x0, quadratic_function,  gradient);
  return std::make_pair( x0, val);
}
/***
 *  f(x) = .5*(x'Ax - x'b)
 *  grad(f) = Ax - b
 * @param A
 * @param x
 * @return
 */
template< typename T>
std::pair< DistMatrix<T>, T>
QuadraticFunction( const DistMatrix<T> & A, const DistMatrix<T>& b)
{
  const std::function< T(const DistMatrix<T>&)>
  quadratic_function = [&](const DistMatrix<T>& theta)
  {
        DistMatrix<T> y( theta);
        Gemv(El::NORMAL, T(1), A, theta, T(0), y);
        return (Dot(y,theta) - Dot(b, theta))/T(2);
  };

  const std::function< DistMatrix<T>(const DistMatrix<T>&, DistMatrix<T>&)>
  gradient = [&](const DistMatrix<T>& theta, El::DistMatrix<T>& y)
  {
      Gemv(El::NORMAL, T(1), A, theta, T(0), y);
      Axpy(T(-1), b, y);
      return y;
  };
  DistMatrix<T> x0( A.Height(), 1);
  Gaussian( x0, A.Height(), 1);
  auto val = El::BFGS( x0, quadratic_function,  gradient);
  return std::make_pair( x0, val);
}

template< typename T>
std::pair< DistMatrix<T>, T>
QuadraticBFGSTest(Int N){
 DistMatrix<T> A,B;
 Laplacian(A, N);
 for(int i = 0; i < N; ++i){ A.Set(i,i,Abs(A.Get(i,i))); }
 Print(A);
 auto nrm = FrobeniusNorm(A);
 Gemm(El::NORMAL, El::ADJOINT, T(1), A, A, B);
 A = B;
 A *= T(1)/(nrm*nrm);
 DistMatrix<T> b( N, 1);
 Gaussian(b, N, 1);
 b *= (T(1)/FrobeniusNorm(b));
 return QuadraticFunction<T>(A, b);
}

template<typename T>
std::function< T(const DistMatrix<T>&)>
rosenbrock_function(){
std::function< T(const DistMatrix<T>&)>
          rosenbrock = [&](const DistMatrix<T>& theta)
{
    auto x1 = theta.Get(0,0);
    auto x2 = theta.Get(1,0);
    auto t1 = (x2 -x1*x1);
    auto t2 = (T(1)-x1);
    //f(x1,x2) = 100*(x2-x1^2)^2 + (1-x1)^2
    return T(100)*t1*t1 + t2*t2;
};
return rosenbrock;
}


template< typename T>
const std::function< DistMatrix<T>(const DistMatrix<T>&, DistMatrix<T>&)>
rosenbrock_gradient(){
const std::function< DistMatrix<T>(const DistMatrix<T>&, DistMatrix<T>&)>
      gradient = [&](const DistMatrix<T>& theta, El::DistMatrix<T>& y)
{
    auto x1 = theta.Get(0,0);
    auto x2 = theta.Get(1,0);

    T g1 = T(2)*(T(200)*x1*x1*x1 - T(200)*x1*x2  + x1 - T(1));
    T g2 = T(2)*(x2-x1*x1);

    y.Set(0, 0, g1);
    y.Set(1, 0, g2);
    return y;
};
return gradient;
}

template< typename T>
std::pair< DistMatrix<T>, T>
RosenbrockTestEasyStart(){
    DistMatrix<T> x0( 2, 1);
    x0.Set(0, 0, 1.2);
    x0.Set(1, 0, 1.2);
    auto val = El::BFGS( x0, rosenbrock_function<T>(),  rosenbrock_gradient<T>());
    return std::make_pair(x0,val);
}

template< typename T>
std::pair< DistMatrix<T>, T>
RosenbrockTestHardStart(){
    DistMatrix<T> x0( 2, 1);
    x0.Set(0, 0, -1.2);
    x0.Set(1, 0, 1.0);
    auto val = El::BFGS( x0, rosenbrock_function<T>(),  rosenbrock_gradient<T>());
    return std::make_pair(x0,val);
}

template< typename T>
std::pair< DistMatrix<T>, T>
RosenbrockTest(){
    DistMatrix<T> x0( 2, 1);
    Gaussian( x0, 2, 1);
    auto val = El::BFGS( x0, rosenbrock_function<T>(),  rosenbrock_gradient<T>());
    return std::make_pair(x0,val);

}

//TEST_CASE( "Can Minimize Simple Quadratic", "[BFGS]" )
//{
//    for(int i = 0; i < 10; ++i) {
//        auto pair = SimpleQuadraticBFGSTest<double>(10);
//        auto x0 = pair.first;
//        auto val = pair.second;
//        REQUIRE(El::Norm(x0) < 1e-15);
//        REQUIRE(val < 1e-15);
//    }
//}

TEST_CASE( "Can Minimize Laplacian Quadratic", "[BFGS]")
{
    for(int i = 0; i < 10; ++i) {
        std::cout << i << std::endl;
        auto pair = QuadraticBFGSTest<double>(10);
        auto x = pair.first;
        auto val = pair.second;
        REQUIRE(El::MaxNorm(x) < 941.5901860155);
        REQUIRE(val < 1e-11);
    }
}

//TEST_CASE( "Can Minimize rosenbrock", "[BFGS]" )
//{
//    for(int i = 0; i < 10; ++i){
//        auto pair = RosenbrockTest<double>();
//        auto x0 = pair.first;
//        auto val = pair.second;
//        REQUIRE( El::Abs(x0.Get(0,0) - 1) < 1e-14);
//        REQUIRE( El::Abs(x0.Get(1,0) - 1) < 1e-14);
//        REQUIRE( val < 1e-16);
//    }
//    {
//        auto pair = RosenbrockTestEasyStart<double>();
//        auto x0 = pair.first;
//        auto val = pair.second;
//        REQUIRE( El::Abs(x0.Get(0,0) - 1) < 1e-14);
//        REQUIRE( El::Abs(x0.Get(1,0) - 1) < 1e-14);
//        REQUIRE( val < 1e-16);
//
//    }
//    {
//        auto pair = RosenbrockTestHardStart<double>();
//        auto x0 = pair.first;
//        auto val = pair.second;
//        REQUIRE( El::Abs(x0.Get(0,0) - 1) < 1e-14);
//        REQUIRE( El::Abs(x0.Get(1,0) - 1) < 1e-14);
//        REQUIRE( val < 1e-16);
//
    }
}


int main( int argc, char* argv[])  {
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );
    return Catch::Session().run( 0, (const char**)(argv));
}
