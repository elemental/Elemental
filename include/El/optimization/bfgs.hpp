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

namespace El {
/***
 * Finds the min of the quadratic interpolation
 * phi_q(alpha) = a*alpha^2 + b*alpha + c
 * where
 * a =[fval_low -f0 - alpha_low*f0_dash]/alpha_low^2
 * b = f0_dash
 * c = f0;
 * then the minimizer is:
 */
template< typename T>
T MinQuadraticInterpolate( T alpha_low, T  f0_dash, T f0, T fval_low) {
        return (f0_dash * alpha_low * alpha_low) / (T(-2) * (fval_low - f0 - f0_dash * alpha_low));
}
 /***
 * Finds the min of the cubic interpolation
 * phi_q(alpha) = a*alpha^3 + b*alpha^2 + alpha*f0_dash + f0
 * a and b are too complicated to encode here, please follow the source code.
 *
 */
template< typename T>
T MinCubicInterpolate(T f0, T f0_dash, T alpha_0, T f_alpha0, T alpha_1, T f_alpha1){
    T alpha_0_sq = alpha_0*alpha_0; //alpha_0^2
    T alpha_0_cb = alpha_0*alpha_0*alpha_0; //alpha_0^3
    T alpha_1_sq = alpha_1*alpha_1; //alpha_1^2
    T alpha_1_cb = alpha_1*alpha_1*alpha_1; //alpha_1^3
    T x = f_alpha1 - f0 - f0_dash*alpha_1;
    T y = f_alpha0 - f0 - f0_dash*alpha_0;
    T a = alpha_0_sq*x - alpha_1_sq*y;
    T b = -alpha_0_cb*x + alpha_1_cb*y;
    T denom = alpha_0_sq*alpha_1_sq*(alpha_1-alpha_0);
    a /= denom;
    b /= denom;
    return (-b + El::Sqrt(b*b - 3*a*f0_dash))/(T(3)*a);
}
/**
 * Preconditions:
 * -------------
 * 1. The interval bounded by alpha_low and alpha_high contains a value alpha satisfying the strong wolfe conditions.
 * 2. alpha_low satisfies the sufficient decrease condition.
 * 3. f(x + alpha_low*p)*(alpha_high - alpha_low) < 0
 *
 * Note: it is _not_ required that alpha_low < alpha_high.
 *
 * @param f
 * @param gradient
 * @param f0
 * @param x0
 * @param p
 * @param alpha_low
 * @param alpha_high
 * @param c1
 * @param c2
 * @return
 */
template< typename T, typename Function, typename Gradient>
T Zoom( const Function& f, const Gradient& gradient,
        T f0, T f0_dash, T fval_low, T flow_dash, const DistMatrix<T>& x0, const DistMatrix<T>& p,
        T alpha_low, T alpha_high, T c1, T c2)
{
    DistMatrix<T> x_j(x0);
    DistMatrix<T> g2(p.Height(), 1, x0.Grid());
    Int iter = 0;
    T alpha;
    while (alpha_high - alpha_low > T(10)*limits::Epsilon<Base<T>>())
    {
        // BEGIN Interpolation
        /**
         * The treatment of the interpolation step of zoom is somewhat vague. We have chosen to build a quadratic interpolant,
         * and if sufficient decrease is not met, we then follow up by a cubic interpolant involving the minimizer of the quadratic.
         */
        alpha = MinQuadraticInterpolate(alpha_low, f0_dash, f0, fval_low);
        x_j = x0;
        Axpy(alpha, p, x_j);
        T fval = f(x_j);
        if ( fval > f0 + c1*alpha*f0_dash)
        {
            alpha = MinCubicInterpolate(f0, f0_dash, alpha_low, fval_low, alpha, fval);
            x_j = x0;
            Axpy(alpha, p, x_j);
            fval = f(x_j);
        }
        // End Interpolation
        //x_j = x0;
        //Axpy(alpha, p, x_j);
        //fval = f(x_j);
        if ((fval > f0 + c1*alpha*f0_dash) || (fval >= fval_low)){
            alpha_high = alpha;
        }else { //alpha certainly satisfies sufficient decrease
            auto fval_dash = Dot(p, g2);
            gradient( x_j, g2); //evaluate fval_dash = phi'(alpha_j)
            //If this value of alpha also satsifies the curvature conditions
            if(Abs(fval_dash)<= -c2*f0) //then return.
            {
                return alpha;
            }
            if(  fval_dash*(alpha_high - alpha_low) >= 0) //maintain invariant (3)
            {
                alpha_high = alpha_low;
            }
            /**
             * alpha_low is maintained as of all steps satisfying suff. dec.
             * it is the one with smallest function value.
             */
            fval_low = fval;
            alpha_low = alpha;
        }
        ++iter;
    }
    return alpha;
}
/***
 * This function approximately minimizes the function min_\alpha f(x0 + alpha*p)
 * The routine lineSearch is gaurunteed to find a step length \alpha satisfying the Strong Wolfe Conditions:
 * 1. f(x + \alpha*p) <= f(x) + c1*\alpha*p'gradient(x) (The sufficient decrease condition)
 * 2. c2|p'gradient(x)| <= |p'gradient(x + \alpha*p)| (The curvature condition)
 * For c1 and c2 fixed: 0 < c1 < c2 < 1
 * The preconditions for the method are that:
 *  1. p is a descent direction (p'g(x) < 0)
 *  2. f is bounded below along the direction p
 *
 */
template< typename T, typename Function, typename Gradient>
T lineSearch( const Function& f, const Gradient& gradient,
              const DistMatrix<T>& g, Int D,
              const DistMatrix<T>& x0, const DistMatrix<T>& p,
              Int maxIter=100,
              T c1=Pow(limits::Epsilon< Base<T> >(), Base<T>(1)/Base<T>(4)), T c2=T(9)/T(10))
    {

        if(c1 > c2) {
            throw ArgException("c1 must be less than c2");
        }

        T f0 = f(x0);
        DistMatrix<T> g2(D, 1, x0.Grid());
        T f0_dash = Dot(p,g);
        T  alpha(1);
        T  alpha_prev(0);
        T  alphaMax(1e3);
        T fvalPrev(0);
        T fvalPrev_dash = f0_dash;
        T fval(0);
        DistMatrix<T> x_candidate(x0);
        // We is maintain that alpha_prev < alpha
        for(std::size_t iter = 0; iter < maxIter; ++iter)
        {
            //We work out f(x + alpha*p)
            x_candidate = x0;
            Axpy(alpha, p, x_candidate);
            fval = f(x_candidate);
            // alpha violates the sufficient descrease condition
            // Since p is a descent direction, f(x+\epsilon*p) < f(x) < f(x+alpha*p)
            // Since f is bounded below we bracket an interval satisfying the Strong wolfe conditions.
            // We use the algorithm Zoom to find a step length satisfying condition 2.
            if ( fval > f0 + c1*alpha*f0_dash || (iter > 0 && fval > fvalPrev) )
            {
                return Zoom(f, gradient, f0, f0_dash, fvalPrev_dash, fvalPrev, x0, p, alpha_prev, alpha, c1, c2);
            }
            // x_candidate satsifies condition (1)
            gradient( x_candidate, g2);
            auto fval_dash = Dot(p, g2);
            // If we satisfy condition (2) then terminate.
            if(  Abs(fval_dash) <= -c2*f0_dash)
            {
                return alpha;
            }
            // Otherwise the curvature condition is violated.
            // If p'g2 > 0 then we are no longer in a descent direction
            // So again the interval (alpha_prev,alpha) must contain a value
            // Which satisfies the curvature conditions.
            if( fval_dash > 0)
            {
                return Zoom(f, gradient, f0, f0_dash, fval, fval_dash, x0, p, alpha, alpha_prev, c1, c2);
            }
            alpha_prev = alpha;
            fvalPrev = fval;
            fvalPrev_dash = fval_dash;
            alpha = Min(2*alpha, alphaMax);
        }
    return alpha;
}

namespace detail {
/***
 * This class stores the updates to the inverse hessian of F, and provides the ability to apply inverse hessian.
 */
template< typename T>
struct HessianInverseOperator {

    typedef DistMatrix<T> DMatrix;
    typedef std::tuple<T, DMatrix, DMatrix, DMatrix, T> Update;
    typedef std::vector<Update> Hessian_updates;
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
    void update( DMatrix& s, DMatrix& y)
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
    Vector g(x);
    Vector g_old(x);
    Vector y(g);
    gradient(x, g);
    detail::HessianInverseOperator<T> Hinv;
    auto norm_g = InfinityNorm(g);
    for( std::size_t iter=0; (norm_g > T(100)*limits::Epsilon<Base<T>>()); ++iter)
    {
        // std::cout << "iter: " << iter << std::endl;
        // Display(x, "Iterate");
        // Display(g, "Gradient");
        //Construct the quasi-newton step
        auto p = Hinv*g; p *= T(-1);
        // Display(p," Descent direction");
        // Evaluate the wolf conditions..
        const T stepSize = lineSearch(F, gradient, g, D, x, p);
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
        Hinv.update(p, y);
        // std::cout << "||g||_inf = " << norm_g << std::endl;
    }
    return F(x);
}

} //end namespace El

#endif // ifndef EL_OPTIMIZATION_PROX_HPP
