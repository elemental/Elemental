/*
   Copyright (c) 2016, Ryan H. Lewis
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_OPTIMIZATION_BFGS_HPP
#define EL_OPTIMIZATION_BFGS_HPP
#include <utility>
#include <El.hpp>
#include <El.h>

namespace El {

/***
 * Source:
 * NONSMOOTH OPTIMIZATION VIA BFGS
 * ADRIAN S. LEWIS AND MICHAEL L. OVERTON
 * http://www.cs.nyu.edu/~overton/papers/pdffiles/bfgs_inexactLS.pdf
 *
 * This function approximately minimizes the function min_\alpha f(x0 + alpha*p)
 * The routine lineSearch is gaurunteed to find a step length \alpha satisfying the Weak Wolfe Conditions:
 * 1. h(\alpha) = f(x + \alpha*p)-f(x) <=  c1*s*alpha (The sufficient decrease condition)
 *      where s = lim sup {t -> 0 from above}  h(t)/t
 *      We assume that s < 0
 * 2. h is differentiable at alpha with h'(\alpha) > c2*s
 * For c1 and c2 fixed: 0 < c1 < c2 < 1
 * We require that f is locally lipschitz or semi-algebraic
 */
template< typename T, typename Vector, typename Gradient>
T lineSearch( const std::function< T(const Vector&)>& f, const Gradient& gradient,
              const DistMatrix<T>& g,
              const DistMatrix<T>& x0, const DistMatrix<T>& p,
              Int maxIter=100,
              T armijoParameter=Pow(limits::Epsilon< Base<T> >(), Base<T>(1)/Base<T>(4)),
              T wolfeParameter=T(9)/T(10)) {
        T c1 = armijoParameter;
        T c2 = wolfeParameter;
        if (c1 > c2) { throw ArgException("c1 must be less than c2"); }
        if (c2 > 1 || c1 < 0) { throw ArgException("0 < c1 < c2 < 1"); }
        T alpha = 0;
        T beta = limits::Infinity<T>();
        T t = 1;
        DistMatrix<T> g2(g);
        DistMatrix<T> x_candidate(x0);
        const T s = Dot(p,g);
        const T f_x0 = f(x0);
        auto h = [&](const T &stepSize) {
            x_candidate = x0;
            Axpy(stepSize, p, x_candidate);
            return f(x_candidate) - f_x0;
        };
        auto phiPrime = [&]() {
            gradient(x_candidate, g2);
            return Dot(p, g2);
        };
        Int numIter=0;
        do {
            if (!(h(t) < c1 * s * t)) { beta = t; }
            else if (!(phiPrime() > c2 * s)) { alpha = t; }
            else { return t; }

            if (limits::IsFinite(beta)) { t = (alpha + beta) / T(2); }
            else { t = 2 * alpha; }
            ++numIter;
        }while( numIter < maxIter);
        return t;
    }

namespace detail {
/***
 * This class stores the updates to the inverse hessian of F, and provides the ability to apply inverse hessian.
 */
template< typename T>
struct HessianInverseOperator {

    typedef DistMatrix<T> DMatrix;
    typedef std::tuple<T, DMatrix, DMatrix, DMatrix, T> UpdateTerm;
    typedef std::vector<UpdateTerm> Hessian_updates;

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
    DMatrix operator*( const DMatrix& x)
    {
        //Initially this is just the identity matrix;
        //H_0 = I, so H_0*x = x;
        if (hessian_data.size() == 0) { return x; }
        DistMatrix<T> z( x);

        // We maintain that z = H_{k-1}*x
        // H_kx = z + alpha*(s'x)s - (s'x)/(s'y)H_{k-1}y  (y'z)/(s'y)*s]/(s'y)
        for( auto& rank_two_update_data: hessian_data)
        {
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
    void Update( DMatrix& s, DMatrix& y)
    {
        // We store s,y, alpha, H*y, and s'y.
        auto sy = Dot(s,y);
        auto Hy = (*this)*y;
        auto yHy = Dot(Hy,y);
        T alpha = (sy + yHy)/(sy*sy);
        hessian_data.emplace_back( alpha, s, y, Hy, sy);
    }

private:
    Hessian_updates hessian_data;
};
} //end namespace detail

/**
 * Performs BFGS to minimize a function f: R^n --> R with gradient g: R^m --> R^n
 * @param x initial position
 * @param F T valued function
 * @param gradient of funtion
 */
template< typename Vector, typename T>
T BFGS( Vector& x, const std::function< T(const Vector&)>& F,
        const std::function< Vector(const Vector&, Vector&)>& gradient)
    {
    const Int D = x.Height();
    Vector g(x); //g stores the gradient at x_n
    Vector g_old(x); //this one the old gradient x_{n-1}
    Vector y(g); //y stores g - g_old
    gradient(x, g);
    detail::HessianInverseOperator<T> Hinv;
    auto norm_g = InfinityNorm(g);
    for( std::size_t iter=0; (norm_g > T(100)*limits::Epsilon<Base<T>>()); ++iter)
    {
        //std::cout << "iter: " << iter << std::endl;
        //Display(x, "Iterate");
        //Display(g, "Gradient");
        //Construct the quasi-newton step
        auto p = Hinv*g; p *= T(-1);
        // Display(p," Descent direction");
        // Evaluate the wolf conditions..
        const T stepSize = lineSearch(F, gradient, g, x, p);
        // std::cout << "Step size: " << stepSize << std::endl;
        // Update x with the step
        // s = stepSize*p;
        // x = x + s
        p *= stepSize;
        Axpy(T(1), p, x);
        // Hold on to old gradient
        g_old = g;
        // Re-evaluate
        gradient(x, g);
        norm_g = InfinityNorm(g);
        if( norm_g < T(100)*limits::Epsilon<Base<T>>()){ return F(x); }
        // Evaluate change in gradient
        y = g;
        // y = g - g_old
        Axpy(T(-1), g_old, y);
        Hinv.Update(p, y);
         std::cout << iter << " ||g||_inf = " << norm_g << std::endl;
    }
    return F(x);
}

} //end namespace El

#endif // ifndef EL_OPTIMIZATION_PROX_HPP
