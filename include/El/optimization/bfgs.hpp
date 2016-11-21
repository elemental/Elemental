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
namespace detail {
/***
 * This class stores the updates to the inverse hessian of F, and provides the ability to apply inverse hessian.
 */
template< typename T>
struct HessianInverseOperator {

    typedef DistMatrix<T> DMatrix;
    typedef std::tuple<T, DMatrix, DMatrix, T> UpdateTerm;
    typedef std::vector<UpdateTerm> Hessian_updates;

    /**
    * This method applies y = H_k*x, eg. solves B_ky = x;
    * This is implemented using the procedure from page 779
    * of Updating Quasi-Netwon Matrices With Limited Storage
    * http://www.ii.uib.no/~lennart/drgrad/Nocedal1980.pdf
    * @param x
    * @return H_k*x
    */
    DMatrix operator*( const DMatrix& x)
    {
        std::size_t k = hessian_data.size();
        Int bound = (k <= M) ? k : M;
        Int incr = (k <= M) ? 0 : k-M;
        if( k == 0){ return x; } //H_0 is initially identity
        DMatrix q(x);
        alphaList.resize(bound);
        for(Int i = bound-1; i >= 0; --i){
            Int j = i+incr;
            const UpdateTerm& updateTerm = hessian_data[j];
            const T& stepSize = std::get<0>(updateTerm);
            const auto& p = std::get<1>(updateTerm.get);
            const auto& y = std::get<2>(updateTerm.get);
            const auto& rho = std::get<3>(updateTerm);
            T alpha = stepSize*Dot(p,q)*rho;
            alphaList[i] = alpha;
            Axpy(-alpha, y, q);
        }
        //r_0 = H_0*q_0, current we assume H_0 = I
        //Later we may scale it.
        for(int i = 0; i < bound-1; ++i){
            Int j = i+incr;
            const UpdateTerm& updateTerm = hessian_data[j];
            const T& stepSize = std::get<0>(updateTerm);
            const auto& p = std::get<1>(updateTerm.get);
            const auto& y = std::get<2>(updateTerm.get);
            const auto& rho = std::get<3>(updateTerm);
            const auto& alpha = alphaList[i];
            T beta = rho*Dot(y,q);
            Axpy(stepSize*(alpha-beta),p,q);
        }
        return q;
    }
    /**
    * The forward update is:
    * B_{k+1} = B_k + yy'/(y's) + B_ksks'B_k/(s'Bs)
    * This method advances the operator H_{k+1} = B_{k+1}^{-1}
    * @param s = x_{k+1} - x_{k} = stepSize * p_k
    * @param y = g_{k+1} - g_{k}
    */
    void Update( const T& stepSize, DMatrix& p, DMatrix& y)
    {
        hessian_data.emplace_back( stepSize, p, y, T(1)/stepSize*Dot(y,p));
    }

    private:
        Hessian_updates hessian_data;
        Int M = hessian_data.max_size(); //for future LBFGS
        std::vector< T> alphaList;
}; //end class HessianInverseOperator
} //end namespace detail

template< typename T>
inline bool IsNaN(const T& t){ return (!limits::IsFinite<T>(t) && t != limits::Infinity<T>()); }

/***
 * Source:
 * NONSMOOTH OPTIMIZATION VIA BFGS
 * ADRIAN S. LEWIS AND MICHAEL L. OVERTON
 * http://www.cs.nyu.edu/~overton/papers/pdffiles/bfgs_inexactLS.pdf
 *
 * This function approximately minimizes the function min_\alpha f(x0 + alpha*p)
 * The routine lineSearch is gaurunteed to find a step length \alpha satisfying the Weak Wolfe Conditions:
 * 1. h(t) = f(x + t*p)-f(x) <=  c1*s*t (The sufficient decrease condition)
 *      where s = lim sup {t -> 0 from above}  h(t)/t
 *      We assume that s < 0
 * 2. h is differentiable at alpha with h'(t) > c2*s
 * For c1 and c2 fixed: 0 < c1 < c2 < 1
 * We require that f is locally lipschitz or semi-algebraic
 */
template< typename T, typename Vector, typename Gradient>
T lineSearch( const std::function< T(const Vector&)>& f, const Gradient& gradient,
              const DistMatrix<T>& g,
              const DistMatrix<T>& x0, const DistMatrix<T>& p,
              T armijoParameter=Pow(limits::Epsilon< Base<T> >(), Base<T>(1)/Base<T>(4)),
              T wolfeParameter=T(9)/T(10)) {
        T c1 = armijoParameter;
        T c2 = wolfeParameter;
        if (c1 > c2) { throw ArgException("c1 must be less than c2"); }
        if (c2 > 1 || c1 < 0) { throw ArgException("0 < c1 < c2 < 1"); }
        T pnorm = Norm(p);
        if (pnorm <= limits::Epsilon<T>()){ throw ArgException("p is too close to zero"); }
        T nBisectMax = Max(T(30), Round(Log2(1e5*pnorm))); // allows more if ||p|| big
        T nExpandMax = Max(T(10), Round(Log2(1e5/pnorm))); // allows more if ||p|| small

        T alpha = 0;
        T beta = limits::Infinity<T>();
        T t = 1;
        DistMatrix<T> g2(g);
        DistMatrix<T> x_candidate(x0);
        const T f_x0 = f(x0);
        const T s = Dot(p,g);
        Int nBisect=0;
        Int nExpand=0;
        std::cout << "t \t f((x+alpha*t) \t f'(x+alpha*t)" << std::endl;
        bool done = false;
        do {
            x_candidate = x0;
            Axpy(t, p, x_candidate);
            T ft = f(x_candidate);
            T gtp = 0;
            gradient(x_candidate, g2);
            gtp = Dot(p, g2);
            std::cout << t << "\t" << ft  << "\t" << gtp << std::endl;
            //Armijo rule is violated.
            if (ft >= f_x0 + c1*s*t || IsNaN(ft) ) { beta = t; }
            else {
                //wolfe condition violated
                if (gtp <= c2*s || IsNaN(gtp)) { alpha = t; }
                else { return t; }
            }
            //Adjust t for the next function evaluation
            if (beta < limits::Infinity<T>()) {
                if( nBisect < nBisectMax){
                    t = (alpha + beta) / T(2);
                    ++nBisect;
                }else { done = true; }
            }
            else {
                if( nExpand < nExpandMax){
                    t = T(2) * alpha;
                    ++nExpand;
                } else {
                    done = true;
                }
            }
        }while( !done);
        if( !limits::IsFinite(beta)){ RuntimeError("Line search failed to brack point satisfying weak wolfe conditions. Function may be unbounded below"); }
        RuntimeError("[",alpha,",",beta,"] brackets an interval containing a point satisfying WWC");
        return t;
    }


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
    auto norm_g = Norm(g);
    auto norm_x = Norm(x);
    T convergenceEpsilon = T(1e-5);
    for( std::size_t iter=0; (norm_g > Max(norm_x,T(1))*convergenceEpsilon); ++iter)
    {
        //std::cout << "iter: " << iter << std::endl;
        //Display(x, "Iterate");
        //Display(g, "Gradient");
        //Construct the quasi-newton step
        auto p = Hinv*g; p *= T(-1);
        //Display(p," Descent direction");
        auto s = Dot(p,g);
        if ( s >= 0) { RuntimeError("p is not a descent direction"); };
        // Evaluate the wolf conditions..
        T stepSize = lineSearch(F, gradient, g, x, p);
        stepSize = Min(Max(stepSize, T(1e-20)),T(1e20)); //These are the defaults used in LBFGS
        // std::cout << "Step size: " << stepSize << std::endl;
        // Update x with the step
        // x = x + stepSize*p;
        Axpy(stepSize, p, x);
        // Hold on to old gradient
        g_old = g;
        // Re-evaluate
        gradient(x, g);
        norm_g = Norm(g);
        norm_x = Norm(x);
        if( norm_g < Max(norm_x,T(1))*convergenceEpsilon){ return F(x); }
        // Evaluate change in gradient
        y = g;
        // y = g - g_old
        Axpy(T(-1), g_old, y);
        Hinv.Update(stepSize, p, y);
        std::cout << iter << " ||g||_inf = " << norm_g << std::endl;
    }
    return F(x);
}

} //end namespace El

#endif // ifndef EL_OPTIMIZATION_PROX_HPP
