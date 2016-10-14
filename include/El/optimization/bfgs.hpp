/*
   Copyright (c) 2009-2016, Ryan H. Lewis
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_OPTIMIZATION_BFGS_HPP
#define EL_OPTIMIZATION_BFGS_HPP
#include <utility>
#include <El.hpp>

namespace El {
/***
 * This function approximately minimizes the function min_\alpha f(x0 + alpha*z);
 * @param x0
 * @param z
 * @param f
 * @param grad
 * @return
 */
template< typename T, typename Function, typename Gradient>
T lineSearch( const DistMatrix<T>& x0, const DistMatrix<T>& z, 
              const Function& f, const Gradient& grad){
    
    DistMatrix<T> x(x0);
    T f0 = f(x0);
    DistMatrix<T> g;
    g = grad(x, g);
    T f0_dash = Dot(z,g);
    T alpha =1;
    bool decrease_direction = true;
    DistMatrix<T> x_candidate;
    DistMatrix<T> grad2;
    for(std::size_t iter = 0; iter < 100; ++iter){
        x_candidate = x;
        Axpy(alpha, z, x_candidate);
        auto fval = f(x_candidate);
        if ( fval > f0 + 0.0001*alpha*f0_dash){
            alpha *= 0.5;
            decrease_direction = false;
        } else {
            grad2 = grad( x_candidate, g);
            auto fval_dash = Dot(z, grad2);
            if ( fval_dash < 0.9 < f0_dash || !decrease_direction){
                alpha *= 4.0;
            } else {
                return alpha;
            }
        }
    }
    return alpha;
}

namespace detail {
/***
 * This class stores the updates to the inverse hessian of F, and provides the ability to apply inverse hessian.
 */
template< typename T>
struct HessianInverseOperator {
    DistMatrix<T> operator*( const DistMatrix<T>& x){    
        //Initially this is just the identity matrix;
        if (_hessian_data.size() == 0) { return x; }
        DistMatrix<T> z( x); 
        for( auto& rank_two_update: _hessian_data){
            //we maintain that z = B_{k-1}*x
            auto coeff = std::get<0>(rank_two_update);
            const auto& s = std::get<1>(rank_two_update);
            const auto& y = std::get<2>(rank_two_update);
            const auto& By = std::get<3>(rank_two_update);
            const auto& sy = std::get<4>(rank_two_update);
            //s'x
            auto sx = Dot(s,x);
            //y_k'z = y_kB_{k-1}*x
            auto yz = Dot(y,z);
            auto coeff2 = -sx/sy;
            auto coeff3 = -yz/sy;
            // Rank-1 update:
            //[(sy + yBy)/(sy*sy)](s'x)*s
            Axpy(coeff*sx, s, z);

            //Rank-1 update: [(s'x)/(s'y)]*(By)
            Axpy(coeff2, By, z);
            //Rank-1 update: (y_kB_{k-1}*x)/(s'y)
            Axpy(coeff3, s, z);
        }
        return z;
    }
    void update( DistMatrix<T>& s, DistMatrix<T>& y){
        auto sy = Dot(s,y);
        auto By = (*this)*y;
        auto yBy = Dot(By,y);
        T coeff = (sy + yBy)/(sy*sy);
        _hessian_data.emplace_back( coeff, s, y, By, sy);
    }

private:
    std::vector< std::tuple< T, DistMatrix<T>, DistMatrix<T>, DistMatrix<T>, T> > _hessian_data;
};
} //end namespace detail

/**
 * Performs BFGS to minimize a function f: R^n --> R with gradient g: R^n --> R^n
 * @param x
 * @param f
 * @param gradient
 */
template< typename Vector, typename T>
T BFGS( Vector& x, const std::function< T(const Vector&)>& f, const std::function< Vector(const Vector&, Vector&)>& gradient){
    const Int DIM = x.Height();
    Vector g(x);
    Vector g_old(x);
    Vector s(x);
    Vector x_old(x);
    Vector y(g);
    detail::HessianInverseOperator<T> H;
    for( std::size_t iter=0; (InfinityNorm(g) > T(1e-16) && iter < DIM); ++iter){
        gradient(x, g);
        //Construct the quasi-newton step
        auto p = (H*g); p*=(-1);
        //Evaluate the wolf conditions..
        const T stepSize = lineSearch(x, p, f, gradient);
        //This is the step
        //Update x with the step
        //x = x + stepSize*p
        Axpy(stepSize, x,  p);
        x = std::move( s); //swap x and s
        //Hold on to old gradient
        g_old = g;
        //Re-evaluate
        gradient(x, g);

        //Evaluate change in gradient
        y = g;
        //y = g - g_old
        Axpy(T(-1), y, g_old);
        H.update(s,y);
        x_old = x;
        ++iter;
    }
    return f(x);
}

} // namespace El

#endif // ifndef EL_OPTIMIZATION_PROX_HPP
