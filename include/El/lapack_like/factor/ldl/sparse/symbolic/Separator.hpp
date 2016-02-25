/*
   Copyright 2009-2011, Jack Poulson.
   All rights reserved.

   Copyright 2011-2012, Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin.
   All rights reserved.

   Copyright 2013, Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright 2013-2014, Jack Poulson and The Georgia Institute of Technology.
   All rights reserved.

   Copyright 2014-2015, Jack Poulson and Stanford University.
   All rights reserved.
   
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_UTIL_SEPARATOR_HPP
#define EL_UTIL_SEPARATOR_HPP

namespace El {
namespace ldl {

struct DistSeparator;

struct Separator
{
    Int off;
    vector<Int> inds;

    Separator* parent; 
    vector<Separator*> children;
    DistSeparator* duplicate; 

    Separator( Separator* parentNode=nullptr )
    : parent(parentNode), duplicate(nullptr)
    { }

    Separator( DistSeparator* dupNode );

    ~Separator()
    {
        if( uncaught_exception() )
        {
            cerr << "Uncaught exception" << endl;
            DEBUG_ONLY(DumpCallStack())
            return;
        }

        for( const Separator* child : children )
            delete child;
    }
};

struct DistSeparator
{
    mpi::Comm comm;
    Int off;
    vector<Int> inds;

    DistSeparator* parent; 
    DistSeparator* child;  
    Separator* duplicate;  

    DistSeparator( DistSeparator* parentNode=nullptr )
    : comm(mpi::COMM_WORLD),
      parent(parentNode), child(nullptr), duplicate(nullptr)
    { }

    ~DistSeparator()
    {
        if( uncaught_exception() )
        {
            cerr << "Uncaught exception" << endl;
            DEBUG_ONLY(DumpCallStack())
            return;
        }

        delete child;
        delete duplicate;

        if( comm != mpi::COMM_WORLD )
            mpi::Free( comm );
    }
};

inline Separator::Separator( DistSeparator* dupNode )
: parent(nullptr), duplicate(dupNode)
{
    off = duplicate->off;
    inds = duplicate->inds;
}

} // namespace ldl
} // namespace El

#endif // ifndef EL_UTIL_SEPARATOR_HPP
