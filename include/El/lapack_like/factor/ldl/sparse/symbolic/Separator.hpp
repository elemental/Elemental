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

    // An observing pointer to the parent (should it exist).
    Separator* parent=nullptr;

    // An observing pointer to the distributed duplicate (should it exist).
    DistSeparator* duplicate=nullptr;

    // Unique pointers to the childen (should they exist).
    vector<unique_ptr<Separator>> children;

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

    // An observing pointer to the parent (should it exist).
    DistSeparator* parent=nullptr;

    // A unique pointer to the distributed child node shared by this process
    // (should it exist).
    unique_ptr<DistSeparator> child;

    // A unique pointer to the sequential duplicate node (should it exist).
    unique_ptr<Separator> duplicate;

    DistSeparator( DistSeparator* parentNode=nullptr );
    ~DistSeparator();

    void BuildMap( const DistNodeInfo& info, DistMap& map ) const;
};

} // namespace ldl
} // namespace El

#endif // ifndef EL_FACTOR_LDL_SPARSE_SYMBOLIC_SEPARATOR_HPP
