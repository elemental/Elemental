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
#ifndef EL_FACTOR_LDL_SPARSE_SYMBOLIC_SEPARATOR_HPP
#define EL_FACTOR_LDL_SPARSE_SYMBOLIC_SEPARATOR_HPP

namespace El {
namespace ldl {

struct DistSeparator;

struct Separator
{
    Int off;
    vector<Int> inds;

    Separator* parent=nullptr;
    vector<Separator*> children;
    DistSeparator* duplicate=nullptr;

    Separator( Separator* parentNode=nullptr );
    Separator( DistSeparator* dupNode );
    ~Separator();

    void BuildMap( vector<Int>& map ) const;
};

struct DistNodeInfo;

struct DistSeparator
{
    Int off;
    vector<Int> inds;

    DistSeparator* parent=nullptr;
    DistSeparator* child=nullptr;
    Separator* duplicate=nullptr;

    DistSeparator( DistSeparator* parentNode=nullptr );
    ~DistSeparator();

    void BuildMap( const DistNodeInfo& info, DistMap& map ) const;
};

} // namespace ldl
} // namespace El

#endif // ifndef EL_FACTOR_LDL_SPARSE_SYMBOLIC_SEPARATOR_HPP
