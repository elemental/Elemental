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
#ifndef EL_CORE_DISTMAP_DECL_HPP
#define EL_CORE_DISTMAP_DECL_HPP

namespace El {

// Use a simple 1d distribution where each process owns a fixed number of
// indices,
//     if last process,  height - (commSize-1)*floor(height/commSize)
//     otherwise,        floor(height/commSize)
class DistMap
{
public:
    // Constructors and destructors
    DistMap( const El::Grid& grid=El::Grid::Default() );
    DistMap( Int numSources, const El::Grid& grid=El::Grid::Default() );
    // TODO(poulson): Constructor for building from a DistMap
    ~DistMap();

    // Map manipulation
    // Collectively map each process's local set of indices
    void Translate( vector<Int>& localInds ) const;
    void Translate
    ( vector<Int>& localInds, const vector<int>& origOwners ) const;

    // composite(i) := second(first(i))
    void Extend( DistMap& firstMap ) const;
    void Extend( const DistMap& firstMap, DistMap& compositeMap ) const;

    // High-level information
    Int NumSources() const;

    // Distribution information
    void SetGrid( const El::Grid& grid );
    const El::Grid& Grid() const;
    Int Blocksize() const;
    Int FirstLocalSource() const;
    Int NumLocalSources() const;
    int RowOwner( Int i ) const;

    // Local data
    Int GetLocal( Int localSource ) const;
    void SetLocal( Int localSource, Int target );
    Int* Buffer();
    const Int* Buffer() const;
    vector<Int>& Map();
    const vector<Int>& Map() const;

    // For modifying the size of the map
    void Empty();
    void Resize( Int numSources );

    // Assignment
    const DistMap& operator=( const DistMap& map );

private:
    Int numSources_;

    // An observing pointer to a pre-existing Grid.
    const El::Grid* grid_=nullptr;

    Int blocksize_;

    vector<Int> map_;

    void InitializeLocalData();
};

void InvertMap( const vector<Int>& map, vector<Int>& inverseMap );
void InvertMap( const DistMap& map, DistMap& inverseMap );

} // namespace El

#endif // ifndef EL_CORE_DISTMAP_DECL_HPP
