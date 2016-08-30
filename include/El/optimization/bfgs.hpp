/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_OPTIMIZATION_BFGS_HPP
#define EL_OPTIMIZATION_BFGS_HPP

namespace El {

template< typename T, typename Function, typename Gradient>
T lineSearch( const DistMatrix<T>& x0, const DistMatrix<T>& z, Function& f, gradient& grad){
    
    DistMatrix<T> x(x0);
    T f0 = f(x0);
    auto g = grad(x);
    T f0_dash = Dot(z,grad);
    T alpha =1;
    bool decrease_direction = true;
    DistMatrix<T> x_candidate;
    DistMatrix<T> grad2;
    for(std::size_t iter = 0; iter < 100; ++iter){
        x_candidate = x;
        Axpy(x_candidate,alpha,z);
        auto fval = f(x_candidate);
        if ( fval > f0 + 0.0001*alpha*f0_dash){
            alpha *= 0.5;
            decrease_direction = false;
        } else {
            grad = grad( x_candidate);
            auto fval_dash = Dot(z,grad2);
            if ( fval_dash < 0.9 < phi0_dash || !decrease_direction){
                alpha *= 4.0;
            } else {
                return alpha;
            }
        }
    }
    return alpha;
}
namespace detail {
template< typename T>
class HessianInverseOperator {
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
            Axpy(z,coeff*sx,s);
            //Rank-1 update: [(s'x)/(s'y)]*(By)
            Axpy(z,coeff2,By);
            //Rank-1 update: (y_kB_{k-1}*x)/(s'y)
            Axpy(z,coeff3,s);
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

    std::vector< std::tuple< T, DistMatrix<T>, DistMatrix<T>, DistMatrix<T>, T> > _hessian_data;
};
} //end namespace detail

template< typename Vector, typename T>
void BFGS( Vector& x, 
           const std::function< const Vector&, T>& f,
           const std::function< const Vector&, Vector>& gradient){
    const Int DIM = x.Height();
    Vector g(x);
    Vector g_old(x);
    Vector s(x);
    Vector x_old(x);
    Vector y(g);
    detail::HessianInverseOperator<T> H;
    do
    {
        gradient(x, g);
        //Construct the quasi-newton step
        auto p = (H*g); p*=(-1);
        //Evaluate the wolf conditions..
        const T stepSize = lineSearch(x,p,f,gradient);
        //This is the step
        s = (stepSize)*p;
        //Update x with the step
        //x = x + stepSize*p
        Axpy(x,T(1),s);
        //Hold on to old gradient
        g_old = g;
        //Re-evaluate
        gradient(x, g);

        //Evaluate change in gradient
        y = g;
        //y = g - g_old
        Axpy(y,T(-1),g_old);
        H.update(s,y);
        x_old = x;
        ++iter;
    }while( MaxNorm(g) > ctrl.tolerance && iter < ctrl.maxIter);

} // namespace El

#endif // ifndef EL_OPTIMIZATION_PROX_HPP
