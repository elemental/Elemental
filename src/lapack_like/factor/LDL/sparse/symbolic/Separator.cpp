/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2012 Jack Poulson, Lexing Ying, and
   The University of Texas at Austin.
   All rights reserved.

   Copyright (c) 2013 Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright (c) 2014 Jack Poulson and The Georgia Institute of Technology.
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace ldl {

Separator::Separator( Separator* parentNode )
: parent(parentNode)
{ }

Separator::Separator( DistSeparator* dupNode )
: duplicate(dupNode)
{
    EL_DEBUG_CSE
    off = duplicate->off;
    inds = duplicate->inds;
}

Separator::~Separator()
{
    EL_DEBUG_CSE
    if( uncaught_exception() )
    {
        cerr << "Uncaught exception" << endl;
        EL_DEBUG_ONLY(DumpCallStack())
        return;
    }
    for( const Separator* child : children )
        delete child;
}

DistSeparator::DistSeparator( DistSeparator* parentNode )
: parent(parentNode)
{ }

DistSeparator::~DistSeparator()
{
    EL_DEBUG_CSE
    if( uncaught_exception() )
    {
        cerr << "Uncaught exception" << endl;
        EL_DEBUG_ONLY(DumpCallStack())
        return;
    }
    delete child;
    delete duplicate;
}

} // namespace ldl
} // namespace El
