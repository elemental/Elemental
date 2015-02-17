/*
   Copyright (c) 2009-2015, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
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
    DistMap();
    DistMap( mpi::Comm comm );
    DistMap( Int numSources, mpi::Comm comm );
    // TODO: Constructor for building from a DistMap
    ~DistMap();

    // If the total number of sources is partitioned among the processes, 
    // calling this routine will have the DistMap map the indices to the 
    // owning process
    void StoreOwners( Int numSource, vector<Int>& localInds, mpi::Comm comm );

    // Map manipulation
    // Collectively map each process's local set of indices
    void Translate( vector<Int>& localInds ) const;
   // composite(i) := second(first(i))
    void Extend( DistMap& firstMap ) const;
    void Extend( const DistMap& firstMap, DistMap& compositeMap ) const;

    // High-level information
    Int NumSources() const;

    // Communicator management
    void SetComm( mpi::Comm comm );
    mpi::Comm Comm() const;

    // Distribution information
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

    mpi::Comm comm_;

    Int blocksize_;
    Int firstLocalSource_;

    vector<Int> map_;
};

void InvertMap( const vector<Int>& map, vector<Int>& inverseMap );
void InvertMap( const DistMap& map, DistMap& inverseMap );
  
} // namespace El

#endif // ifndef EL_CORE_DISTMAP_DECL_HPP
