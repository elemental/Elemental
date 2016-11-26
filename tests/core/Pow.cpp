/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int alphaMin = Input("--alphaMin","minimum base",-5);
        const Int alphaMax = Input("--alphaMax","maximum base",5);
        const Int betaMax = Input("--betaMax","maximum exponent",13);
        ProcessInput();

        for( Int alpha=alphaMin; alpha<=alphaMax; ++alpha )
        {
            for( Int beta=0; beta<=betaMax; ++beta )
            {
                const Int result = Pow( alpha, beta );
                Int resultMan = 1;
                for( Int exponent=0; exponent<beta; ++exponent )
                    resultMan *= alpha;
                if( result != resultMan )
                    LogicError(alpha,"^",beta,"=",resultMan,"!=",result);
            }
        }
        Output
        ("Integer exponentiation tests passed for bases in [",alphaMin,",",
         alphaMax,"] and exponents in [0,",betaMax,"]");
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
