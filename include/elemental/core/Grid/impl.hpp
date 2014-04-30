/*
   Copyright (c) 2009-2014, Jack Poulson
                      2013, Jed Brown 
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_GRID_IMPL_HPP
#define ELEM_GRID_IMPL_HPP

namespace elem {

inline int
Grid::FindFactor( int p )
{
    int factor = int(sqrt(double(p)));
    while( p % factor != 0 )
        ++factor;
    return factor;
}

inline 
Grid::Grid( mpi::Comm comm, GridOrder order )
: haveViewers_(false), order_(order)
{
    DEBUG_ONLY(CallStackEntry cse("Grid::Grid"))

    // Extract our rank, the underlying group, and the number of processes
    mpi::Dup( comm, viewingComm_ );
    mpi::CommGroup( viewingComm_, viewingGroup_ );
    size_ = mpi::Size( viewingComm_ );

    // All processes own the grid, so we have to trivially split viewingGroup_
    owningGroup_ = viewingGroup_;

    // Factor p
    height_ = FindFactor( size_ );
    SetUpGrid();
}

inline 
Grid::Grid( mpi::Comm comm, int height, GridOrder order )
: haveViewers_(false), order_(order)
{
    DEBUG_ONLY(CallStackEntry cse("Grid::Grid"))

    // Extract our rank, the underlying group, and the number of processes
    mpi::Dup( comm, viewingComm_ );
    mpi::CommGroup( viewingComm_, viewingGroup_ );
    size_ = mpi::Size( viewingComm_ );

    // All processes own the grid, so we have to trivially split viewingGroup_
    owningGroup_ = viewingGroup_;

    height_ = height;
    if( height_ < 0 )
        LogicError("Process grid dimensions must be non-negative");

    SetUpGrid();
}

inline void 
Grid::SetUpGrid()
{
    DEBUG_ONLY(CallStackEntry cse("Grid::SetUpGrid"))
    if( size_ % height_ != 0 )
        LogicError
        ("Grid height, ",height_,", does not evenly divide grid size, ",size_);

    const int width = size_ / height_;
    gcd_ = elem::GCD( height_, width );
    int lcm = size_ / gcd_;

    // Create the communicator for the owning group (mpi::COMM_NULL otherwise)
    mpi::Create( viewingComm_, owningGroup_, owningComm_ );

    vectorColToViewingMap_.resize(size_);
    diagPathsAndRanks_.resize(2*size_);
    MemZero( diagPathsAndRanks_.data(), 2*size_ );
    const bool colMajor = (order_==COLUMN_MAJOR);
    if( InGrid() )
    {
        // Create a cartesian communicator
        int dimensions[2];
        if( colMajor )
        {
            dimensions[0] = width;
            dimensions[1] = height_;
        }
        else
        {
            dimensions[0] = height_;
            dimensions[1] = width;
        }
        int periods[2] = { true, true };
        bool reorder = false;
        mpi::CartCreate
        ( owningComm_, 2, dimensions, periods, reorder, cartComm_ );

        // Set up the MatrixCol and MatrixRow communicators
        int remainingDimensions[2];
        remainingDimensions[0] = ( colMajor ? false : true  );
        remainingDimensions[1] = ( colMajor ? true  : false );
        mpi::CartSub( cartComm_, remainingDimensions, matrixColComm_ );
        remainingDimensions[0] = ( colMajor ? true  : false );
        remainingDimensions[1] = ( colMajor ? false : true  );
        mpi::CartSub( cartComm_, remainingDimensions, matrixRowComm_ );
        const int matrixColRank = mpi::Rank( matrixColComm_ );
        const int matrixRowRank = mpi::Rank( matrixRowComm_ );

        // Set up the VectorCol and VectorRow communicators
        const int vectorColRank = matrixColRank + height_*matrixRowRank;
        const int vectorRowRank = matrixRowRank + width*matrixColRank;
        mpi::Split( cartComm_, 0, vectorColRank, vectorColComm_ );
        mpi::Split( cartComm_, 0, vectorRowRank, vectorRowComm_ );

        // Set up the map from the VC group to the viewingGroup_ ranks.
        mpi::Group vectorColGroup;
        mpi::CommGroup( vectorColComm_, vectorColGroup ); 
        std::vector<int> ranks(size_);
        for( int i=0; i<size_; ++i )
            ranks[i] = i;
        mpi::Translate
        ( vectorColGroup, size_, ranks.data(), 
          viewingGroup_,         vectorColToViewingMap_.data() );
        mpi::Free( vectorColGroup );

        // Compute which diagonal 'path' we're in, and what our rank is, then
        // perform AllGather world to store everyone's info
        std::vector<int> myDiagPathAndRank(2);
        myDiagPathAndRank[0] = Mod(matrixRowRank-matrixColRank,gcd_);
        int diagPathRank = 0;
        int row = 0;
        int col = myDiagPathAndRank[0];
        for( int j=0; j<lcm; ++j )
        {
            if( row == matrixColRank && col == matrixRowRank )
            {
                myDiagPathAndRank[1] = diagPathRank;
                break;
            }
            else
            {
                row = (row + 1) % height_;
                col = (col + 1) % width;
                ++diagPathRank;
            }
        }
        mpi::AllGather
        ( myDiagPathAndRank.data(),  2, 
          diagPathsAndRanks_.data(), 2, vectorColComm_ );

        mpi::Split( cartComm_, DiagPath(), DiagPathRank(), matrixDiagComm_ );
        mpi::Split
        ( cartComm_, DiagPathRank(), DiagPath(), matrixDiagPerpComm_ );

        DEBUG_ONLY(
            mpi::ErrorHandlerSet( matrixColComm_,      mpi::ERRORS_RETURN );
            mpi::ErrorHandlerSet( matrixRowComm_,      mpi::ERRORS_RETURN );
            mpi::ErrorHandlerSet( vectorColComm_,      mpi::ERRORS_RETURN );
            mpi::ErrorHandlerSet( vectorRowComm_,      mpi::ERRORS_RETURN );
            mpi::ErrorHandlerSet( matrixDiagComm_,     mpi::ERRORS_RETURN );
            mpi::ErrorHandlerSet( matrixDiagPerpComm_, mpi::ERRORS_RETURN );
        )
    }
    else
    {
        matrixColComm_      = mpi::COMM_NULL;
        matrixRowComm_      = mpi::COMM_NULL;
        vectorColComm_      = mpi::COMM_NULL;
        vectorRowComm_      = mpi::COMM_NULL;
        matrixDiagComm_     = mpi::COMM_NULL; 
        matrixDiagPerpComm_ = mpi::COMM_NULL;
        // diag paths and ranks are implicitly set to undefined
    }
    // Translate the rank of the root process of the owningGroup so that we can
    // broadcast data
    int owningRoot = mpi::Translate( owningGroup_, 0, viewingGroup_ );
    mpi::Broadcast
    ( vectorColToViewingMap_.data(), size_, owningRoot, viewingComm_ );
    mpi::Broadcast
    ( diagPathsAndRanks_.data(), 2*size_, owningRoot, viewingComm_ );
}

inline 
Grid::~Grid()
{
    if( !mpi::Finalized() )
    {
        if( InGrid() )
        {
            mpi::Free( matrixDiagComm_ );
            mpi::Free( matrixDiagPerpComm_ );
            mpi::Free( matrixColComm_ );
            mpi::Free( matrixRowComm_ );
            mpi::Free( vectorColComm_ );
            mpi::Free( vectorRowComm_ );
            mpi::Free( cartComm_ );
            mpi::Free( owningComm_ );
        }
        mpi::Free( viewingComm_ );
        if( HaveViewers() )
            mpi::Free( owningGroup_ );
        mpi::Free( viewingGroup_ );
    }
}

inline int Grid::MCRank() const { return mpi::Rank(matrixColComm_); }
inline int Grid::MRRank() const { return mpi::Rank(matrixRowComm_); }
inline int Grid::VCRank() const { return mpi::Rank(vectorColComm_); }
inline int Grid::VRRank() const { return mpi::Rank(vectorRowComm_); }

inline int Grid::MCSize() const { return height_;       }
inline int Grid::MRSize() const { return size_/height_; }
inline int Grid::VCSize() const { return size_;         }
inline int Grid::VRSize() const { return size_;         }

inline mpi::Comm Grid::MCComm()     const { return matrixColComm_;      }
inline mpi::Comm Grid::MRComm()     const { return matrixRowComm_;      }
inline mpi::Comm Grid::VCComm()     const { return vectorColComm_;      }
inline mpi::Comm Grid::VRComm()     const { return vectorRowComm_;      }
inline mpi::Comm Grid::MDComm()     const { return matrixDiagComm_;     }
inline mpi::Comm Grid::MDPerpComm() const { return matrixDiagPerpComm_; }

// Provided for simplicity, but redundant
// ======================================

inline int 
Grid::Rank() const { return ( order_==COLUMN_MAJOR ? VCRank() : VRRank() ); }

inline int Grid::Height() const { return MCSize(); }
inline int Grid::Width()  const { return MRSize(); }
inline int Grid::Size()   const { return VCSize(); }

inline GridOrder Grid::Order() const { return order_; }

inline int Grid::Row() const { return MCRank(); }
inline int Grid::Col() const { return MRRank(); }
inline mpi::Comm Grid::ColComm() const { return MCComm(); }
inline mpi::Comm Grid::RowComm() const { return MRComm(); }
inline mpi::Comm Grid::Comm() const
{ return ( order_==COLUMN_MAJOR ? VCComm() : VRComm() ); }

// Advanced routines
// =================

// Currently forces a columnMajor absolute rank on the grid
inline 
Grid::Grid( mpi::Comm viewers, mpi::Group owners, int height, GridOrder order )
: haveViewers_(true), order_(order)
{
    DEBUG_ONLY(CallStackEntry cse("Grid::Grid"))

    // Extract our rank and the underlying group from the viewing comm
    mpi::Dup( viewers, viewingComm_ );
    mpi::CommGroup( viewingComm_, viewingGroup_ );

    // Extract our rank and the number of processes from the owning group
    mpi::Dup( owners, owningGroup_ );
    size_ = mpi::Size( owningGroup_ );

    height_ = height;
    if( height_ < 0 )
        LogicError("Process grid dimensions must be non-negative");

    SetUpGrid();
}

inline int Grid::GCD() const { return gcd_; }
inline int Grid::LCM() const { return size_/gcd_; }

inline bool Grid::HaveViewers() const { return haveViewers_; }
inline bool Grid::InGrid() const 
{ return mpi::Rank(owningGroup_) != mpi::UNDEFINED; }

inline int Grid::OwningRank() const { return mpi::Rank(owningGroup_); }
inline int Grid::ViewingRank() const { return mpi::Rank(viewingComm_); }

inline int Grid::VCToViewingMap( int VCRank ) const
{ return vectorColToViewingMap_[VCRank]; }

inline mpi::Group Grid::OwningGroup() const { return owningGroup_; }
inline mpi::Comm Grid::OwningComm()  const { return owningComm_; }
inline mpi::Comm Grid::ViewingComm() const { return viewingComm_; }

inline int 
Grid::DiagPath() const
{ 
    const int vcRank = VCRank();
    if( vcRank != mpi::UNDEFINED )
        return diagPathsAndRanks_[2*vcRank];
    else
        return mpi::UNDEFINED;
}

inline int 
Grid::DiagPath( int vectorColRank ) const
{ 
    if( vectorColRank != mpi::UNDEFINED )
        return diagPathsAndRanks_[2*vectorColRank]; 
    else
        return mpi::UNDEFINED;
}

inline int 
Grid::DiagPathRank() const
{ 
    const int vcRank = VCRank();
    if( vcRank != mpi::UNDEFINED )
        return diagPathsAndRanks_[2*vcRank+1];
    else
        return mpi::UNDEFINED;
}

inline int 
Grid::DiagPathRank( int vectorColRank ) const
{ 
    if( vectorColRank != mpi::UNDEFINED )
        return diagPathsAndRanks_[2*vectorColRank+1]; 
    else
        return mpi::UNDEFINED;
}

inline int Grid::FirstVCRank( int diagPath ) const{ return diagPath*height_; }

// Comparison functions
// ====================

inline bool operator==( const Grid& A, const Grid& B ) { return &A == &B; }
inline bool operator!=( const Grid& A, const Grid& B ) { return &A != &B; }

} // namespace elem

#endif // ifndef ELEM_GRID_IMPL_HPP
