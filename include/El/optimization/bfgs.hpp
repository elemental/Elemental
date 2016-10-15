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
 * This function approximately minimizes the function min_\alpha f(x0 + alpha*p)
 * In some sense minimizing this function is "ideal" however, this may take a lot of work, and since
 * we are using this inside of a larger algorithm on an approximate problem, there is no reaosn to believe
 * that this is the "best" step. Instead we just hope for a step with "sufficient" decrease, and then we
 * hope to prove that if {x_i} approaches the minimum that this function converges.
 * This algorithm attempts to satsify the weak wolf conditions.
 */
template< typename T, typename Function, typename Gradient>
T lineSearch( const Function& f, const Gradient& gradient,
              const DistMatrix<T>& g, Int D,
              const DistMatrix<T>& x0, const DistMatrix<T>& p,
              Int maxIter=100, T c1=1e-4, T c2=0.9){

        if( c1 > c2) { std::swap(c1, c2); }

        T f0 = f(x0);
        DistMatrix<T> g2(D, 1);
        T f0_dash = Dot(p,g);
        T  alpha = 0;
        T  t = 1;
        T  beta = El::limits::Infinity<T>();
        DistMatrix<T> x_candidate(x0);
        for(std::size_t iter = 0; iter < maxIter; ++iter){
            x_candidate = x0;
            Axpy(alpha, p, x_candidate);
            auto fval = f(x_candidate);
            //The quantity on the LHS is less than f0
            //So we have found a direction where the fval
            //has not gone down sufficiently.
            if ( fval > f0 + c1*alpha*f0_dash){
                beta = t;
                alpha = T(0.5)*(alpha + beta);
            } else { //The function has gone down sufficiently
            //We now seek to satisfy curvature constraints.
            gradient( x_candidate, g2);
            auto fval_dash = Dot(p, g2);
            if ( fval_dash < c2*f0_dash){
                alpha = t;
                if (El::limits::IsFinite(beta)){
                    t = (alpha+beta)*T(.5);
                } else {
                    t = 2*alpha;
                }
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

    typedef El::DistMatrix<T> Matrix;
    typedef std::tuple< T, Matrix, Matrix, Matrix, T> Update;
    typedef std::vector< Update > Hessian_updates;
    /**
     * This method applies y = H_k*x, eg. solves B_ky = x;
     * (H_k)x = z + alpha*(s'x)s - beta*H_{k-1}y- rho*s
     * where alpha = (s'y + y'H_{k-1}y)/(s'y)^2
     * where beta = (s'x)/(s'y)
     * where rho = (y'z)/(s'y)
     *
     * @param x
     * @return H_k*x
     */
    Matrix operator*( const Matrix& x){
        //Initially this is just the identity matrix;
        //H_0 = I, so H_0*x = x;
        if (_hessian_data.size() == 0) { return x; }
        DistMatrix<T> z( x);

        //we maintain that z = H_{k-1}*x
        // H_kx = z + alpha*(s'x)s - (s'x)/(s'y)H_{k-1}y  (y'z)/(s'y)*s]/(s'y)
        for( auto& rank_two_update_data: _hessian_data){
            auto alpha = std::get<0>(rank_two_update_data);
            const auto& s = std::get<1>(rank_two_update_data);
            const auto& y = std::get<2>(rank_two_update_data);
            const auto& Hy = std::get<3>(rank_two_update_data);
            const auto& sy = std::get<4>(rank_two_update_data);
            //s'x
            auto sx = Dot(s,x);
            //y'z
            auto yz = Dot(y,z);
            auto beta = -sx/sy;
            auto rho = -yz/sy;
            // (H_k)x = z + alpha*(s'x)s + beta*H_{k-1}y + rho*s
            //We now view:
            //  (H_k)x = z + alpha*s + beta*Hy + rho*s
            //         = z + (alpha+rho)*s - beta*Hy
            // Rank-1 update: z <- z+(alpha+rho)*s
            Axpy((alpha*sx+rho), s, z);
            //Rank-1 update: z <- z + beta*Hy
            Axpy( beta, Hy, z);
        }
        return z;
    }
    /**
     * The forward update is:
     * B_{k+1} = B_k + yy'/(y's) + B_ksks'B_k/(s'Bs)
     * So the backwards update is:
     * H_{k+1} = H_k + alpha*(ss') - [H_kys' + sy'H_k]/(s'y)
     * where alpha = [(s'y + y'H_ky)/(s'y)^2];
     * This method advances the operator to represent H_{k+1}
     * @param s
     * @param y
     */
    void update( DistMatrix<T>& s, DistMatrix<T>& y){
        // We store s,y, alpha, H*y, and s'y.
        auto sy = Dot(s,y);
        auto Hy = (*this)*y;
        auto yHy = Dot(Hy,y);
        T alpha = (sy + yHy)/(sy*sy);
        _hessian_data.emplace_back( alpha, s, y, Hy, sy);
    }

private:
    Hessian_updates _hessian_data;
};
} //end namespace detail

/**
 * Performs BFGS to minimize a function f: R^n --> R with gradient g: R^m --> R^n
 * @param x
 * @param f
 * @param gradient
 */
template< typename Vector, typename T>
T BFGS( Vector& x, const std::function< T(const Vector&)>& f,
        const std::function< Vector(const Vector&, Vector&)>& gradient){
    const Int D = x.Height();
    Vector g(x);
    Vector g_old(x);
    Vector y(g);
    detail::HessianInverseOperator<T> Hinv;
    for( std::size_t iter=0; (InfinityNorm(g) > limits::Epsilon<T>()); ++iter){
        gradient(x, g);
        //El::Display(x, "Iterate");
        //El::Display(g, "Gradient");
        //Construct the quasi-newton step
        auto p = Hinv*g; p *= -1;
        //El::Display(p," Descent direction");
        //Evaluate the wolf conditions..
        const T stepSize = lineSearch(f, gradient, g, D, x, p);
        //std::cout << "Step size: " << stepSize << std::endl;
        //Update x with the step
        //s = stepSize*p;
        //x = x + s
        p *= stepSize;
        Axpy(T(1), p, x);
        //Hold on to old gradient
        g_old = g;
        //Re-evaluate
        gradient(x, g);

        //Evaluate change in gradient
        y = g;
        //y = g - g_old
        Axpy(T(-1), g_old, y);
        Hinv.update(p, y);
    }
    return f(x);
}

} // namespace El

#endif // ifndef EL_OPTIMIZATION_PROX_HPP
