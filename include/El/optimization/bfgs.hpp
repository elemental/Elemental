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
 * Finds the minimizer of the quadratic that interpolates (x,fx) and (y, fy) with derivative (x, gx)
 *  f(x) = ax^2 + bx + c
 *  the minimizer = -b / 2*a
 *  So we need only to determine b and a.
 *  it turns out that b = [x*(fx-fy)+gx*(x^2 - y2)]/(x*y-y^2)
 *  and a = [(fx-fy) + gx*(y-x)]/(x*y-y^2)
 *  -b/2*a = -[x*(fx-fy)+gx*(x^2 - y2)]/(2 [(fx-fy) + gx*(y-x)])
 */
template< typename T>
inline T MinQuadraticInterpolate( T x, T fX, T y, T fY, T gX) { return -(x*(fX-fY) + gX*(x*x - y*y))/(T(2)*(fX-fY) + gX*(y-x)); }
/**
 * Finds the minimizer of the quadratic that interpolates (l,fl) with derivative (l, gL) and (t,gT)
 *  f(x) = ax^2 + bx + c
 *  the minimizer = -b / 2*a
 *  It turns out b = [1/(l-t)]*[l*gT - t*gL];
 *  It turns out a = [1/(l-t)]*[-1*gT - gL];
 *  -b/2*a = [l*gT - t*gL]/(T(2)(gT+gL));
 * @param l
 * @param fL
 * @param gL
 * @param t
 * @param gT
 * @return
 */
template< typename T>
inline T MinQuadraticInterpolateTwoDerivatives( T l, T gL, T t, T gT) { return(l*gT - t*gL)/(T(2)*(gT + gL)); }



 /**
 * Finds the min of the cubic interpolation
 * f(x) = a*x^3 + b*x^2 + c*x + d
 * this thing has solutions:
 * x = -b + sqrt(b^2 - 3*a*c) / 3*a
 * and
 * x = -b - sqrt(b^2 - 3*a*c) / 3*a
 * One can find the solutions for a,b,c
 * (a,b,c,d) = (invert {{x^3, x^2, x, 1},{t^3, t^2, t, 1}, {3x^2, 2x, 1, 0}, {3t^2, 2t, 1, 0}})*{{fX},{fT},{gX},{gT}}
 */
template< typename T, typename Function>
inline T MinCubicInterpolate(T x, T fX, T t, T fT, T gX, T gT, const Function& phi){
     T xSquared = x*x;
     T tSquared = t*t;
     T invDet = (T(1)/((t-x)*(t-x)));
     T a = ((T(2)*fX - T(2)*fT)/(t-x) + gX + gT);
     T b = (((t+x)/(t-x))*(T(-3)*fX + 3*gT) + (T(-2)*t-x)*gX  + (-t-T(2)*x)*gT);
     T c = ((T(6)*t*x/(t-x))*(fX-fT) + (tSquared + 2*x*t)*gX  + (xSquared + 2*t*x)*gT);
     T discriminant = sqrt(b*b - 3*a*c)*(1/(t-x));
     b*=invDet;
     a*=invDet;
     T root1 = (-b + discriminant)/(T(3)*a);
     T root2 = (-b - discriminant)/(T(3)*a);
     T f1 = phi(root1);
     T f2 = phi(root2);
     if( f1 <= f2){ return root1; }
     return root2;
 }

/**
 * This technique is borrowed from the paper:
 * Line Search Algorithms with Guaranteed Sufficient Decrease
 * By JORGE J. MORE and DAVID J. THUENTE
 * Section 4, Trial Value Selection
 * https://pdfs.semanticscholar.org/d258/1519560dd59e5b92599f9711aa1ab249ff86.pdf

 * @param alphaLow
 * @param alphaHigh
 * @param alpha
 * @param phiLow
 * @param phiHigh
 * @param phiAlpha
 * @param dPhiLow
 * @param dphiAlpha
 * @return
 */
template< typename Function, typename T>
inline T Interpolate(T alphaLow, T alphaHigh, T alpha, T phiLow, T phiHigh, T phiAlpha,
              T dPhiLow, T dPhiHigh, T dPhiAlpha, const Function& phi){

            T alphaC = MinCubicInterpolate(alphaLow, phiLow, alpha, phiAlpha, dPhiLow, dPhiAlpha, phi);
            if( phiAlpha > phiLow){
                T alphaQ = MinQuadraticInterpolate(alphaLow, phiLow, alpha, phiAlpha, dPhiLow);
                if( Abs(alphaC - alphaLow) < Abs(alphaQ - alphaLow)){ return alphaC; }
                return (alphaQ+alphaC)/T(2);
            }
            T alphaS = MinQuadraticInterpolateTwoDerivatives(alphaLow, dPhiLow, alpha, dPhiAlpha);
            if(  phiAlpha < phiLow and dPhiAlpha*dPhiLow < 0){
                if( Abs(alphaC - alpha) >= Abs(alphaS - alpha)){ return alphaC; }
                return alphaS;
            }
            if(  phiAlpha < phiLow and dPhiAlpha*dPhiLow > 0 and Abs(dPhiAlpha) <= Abs(dPhiLow)){
                if( Abs(alphaC - alpha) < Abs(alphaS - alpha)){ return alphaC; }
                return alphaS;
            }
            //Else: phiAlpha <= phiLow and dPhiAlpha*dPhiLow > 0 and Abs(dPhiAlpha) > Abs(dPhiLow))
            return  MinCubicInterpolate(alphaHigh, phiHigh, alpha, phiAlpha, dPhiHigh, dPhiAlpha, phi);
}

/**
 * Preconditions:
 * -------------
 * 1. The interval bounded by alphaLow and alphaHigh contains a value alpha satisfying the strong wolfe conditions.
 * 2. alphaLow satisfies the sufficient decrease condition.
 * 3. f(x + alphaLow*p)*(alphaHigh - alphaLow) < 0
 *
 * Note: it is _not_ required that alphaLow < alphaHigh.
 *
 * @param f
 * @param gradient
 * @param phi
 * @param phi0
 * @param phiPrime0
 * @param alphaLow
 * @param alphaHigh
 * @param phiAlphaLow
 * @param x0
 * @param p
 * @param c1
 * @param c2
 * @return
 */
template< typename T, typename Function, typename Function1, typename Gradient>
T Zoom( const Function& f, const Gradient& gradient,
        const Function1& phi,
        T phi0, T phiPrime0,
        T alphaLow, T alphaHigh,
        const DistMatrix<T>& x0, const DistMatrix<T>& p,
        T c1, T c2)
{
    DistMatrix<T> x_j(x0);
    DistMatrix<T> g2(p);
    Int iter = 0;

    auto phiPrime = [&] (const T& stepSize){
        x_j = x0; //TODO: This axpy can be optimized away later if it is an issue later.
        Axpy(stepSize, p, x_j);
        gradient(x_j, g2);
        return Dot(p, g2);
    };

    //Initial Guess
    T alpha = (alphaHigh+alphaLow)/T(2);
    T phiLow = phi(alphaLow);
    T phiHigh = phi(alphaHigh);
    T phiAlpha = phi(alpha);
    T prevAlpha;
    T dPhiLow = phiPrime(alphaLow);
    T dPhiHigh = phiPrime(alphaHigh);
    T dPhiAlpha = phiPrime(alpha);
    while (alphaHigh - alphaLow > T(100)*limits::Epsilon<Base<T>>() ) {
        alpha = Interpolate(alphaLow, alphaHigh, alpha, phiLow, phiHigh, phiAlpha,
                      dPhiLow, dPhiHigh, dPhiAlpha, phi);
        if( std::isnan(alpha)) { return prevAlpha; }
        if( Abs(alpha-prevAlpha)<c1){ return alpha; }
        //hand inlined phi(alpha);
        x_j = x0;
        Axpy(alpha, p, x_j);
        //Hand inlined phiPrime(alpha);
        phiAlpha = f(x_j);
        gradient(x_j, g2);
        dPhiAlpha = Dot(p, g2);
        // End Interpolation
        if ((phiAlpha > phi0 + c1*alpha*phiPrime0) || (phiAlpha >= phiLow)){
            alphaHigh = alpha;
            phiHigh = phiAlpha;
            dPhiHigh = dPhiAlpha;
        }else { //alpha certainly satisfies sufficient decrease
            //If this value of alpha also satsifies the curvature conditions
            if(Abs(dPhiAlpha)<= -c2*phi0) //then return.
            {
                return alpha;
            }
            if(  dPhiAlpha*(alphaHigh - alphaLow) >= 0) //maintain invariant (3)
            {
                alphaHigh = alphaLow;
                phiHigh = phiLow;
                dPhiHigh = dPhiLow;
            }
            /**
             * alphaLow is maintained as of all steps satisfying suff. dec.
             * it is the one with smallest function value.
             */
            phiLow = phiAlpha;
            alphaLow = alpha;
            dPhiLow = dPhiAlpha;
        }
        prevAlpha = alpha;
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

        if(c1 > c2) { throw ArgException("c1 must be less than c2"); }

        T phi0 = f(x0);
        DistMatrix<T> g2(D, 1, x0.Grid());
        T  phiPrime0 = Dot(p,g);
        T  alpha(1);
        T  alpha_prev(0);
        T  alphaMax(1e3);
        T fvalPrev(0);
        T phiPrimePrev = phiPrime0;
        T fval(0);
        DistMatrix<T> x_candidate(x0);

        auto phi = [&](const T& stepSize){
            x_candidate = x0;
            Axpy(stepSize, p, x_candidate);
            return f(x_candidate);
        };


        auto phiPrime = [&](){
            gradient(x_candidate, g2);
            return Dot(p, g2);
        };


        // We is maintain that alpha_prev < alpha
        for(std::size_t iter = 0; iter < maxIter; ++iter)
        {
            //We work out f(x + alpha*p)
            fval = phi(alpha);
            // alpha violates the sufficient descrease condition
            // Since p is a descent direction, f(x+\epsilon*p) < f(x) < f(x+alpha*p)
            // Since f is bounded below we bracket an interval satisfying the Strong wolfe conditions.
            // We use the algorithm Zoom to find a step length satisfying condition 2.
            if ( fval > phi0 + c1*alpha*phiPrime0 || (iter > 0 && fval > fvalPrev) )
            {
                return Zoom(f, gradient, phi, phi0, phiPrime0, alpha_prev, alpha, x0, p, c1, c2);
            }
            // x_candidate satsifies condition (1)
            auto phiPrimeAlpha = phiPrime();

            // If we satisfy condition (2) then terminate.
            if(  Abs(phiPrimeAlpha) <= -c2*phiPrime0)
            {
                return alpha;
            }
            // Otherwise the curvature condition is violated.
            // If p'g2 > 0 then we are no longer in a descent direction
            // So again the interval (alpha_prev,alpha) must contain a value
            // Which satisfies the curvature conditions.
            if( phiPrimeAlpha > 0)
            {
                return Zoom(f, gradient, phi, phi0, phiPrime0, alpha, alpha_prev, x0, p, c1, c2);
            }
            alpha_prev = alpha;
            fvalPrev = fval;
            phiPrimePrev = phiPrimeAlpha;
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
        Hinv.Update(p, y);
         std::cout << iter << " ||g||_inf = " << norm_g << std::endl;
    }
    return F(x);
}

} //end namespace El

#endif // ifndef EL_OPTIMIZATION_PROX_HPP
