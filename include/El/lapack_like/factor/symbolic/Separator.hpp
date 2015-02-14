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
#pragma once
#ifndef EL_SYMBOLIC_SEPARATOR_HPP
#define EL_SYMBOLIC_SEPARATOR_HPP

namespace El {

struct DistSeparator;

struct Separator
{
    Int off;
    vector<Int> inds;

    Separator* parent; 
    vector<Separator*> children;
    DistSeparator* duplicate; 

    Separator( Separator* parentNode=nullptr );
    Separator( DistSeparator* dupNode );

    ~Separator();
};

struct DistSeparator
{
    mpi::Comm comm;
    Int off;
    vector<Int> inds;

    DistSeparator* parent; 
    DistSeparator* child;  
    Separator* duplicate;  

    DistSeparator( DistSeparator* parentNode=nullptr );
    ~DistSeparator();
};

} // namespace El

#endif // ifndef EL_SYMBOLIC_SEPARATOR_HPP
