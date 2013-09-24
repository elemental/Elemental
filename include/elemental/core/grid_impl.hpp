/*
   Copyright (c) 2009-2013, Jack Poulson
                      2013, Jed Brown 
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_GRID_IMPL_HPP
#define ELEM_CORE_GRID_IMPL_HPP

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
Grid::Grid( mpi::Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("Grid::Grid");
#endif
    inGrid_ = true; // this is true by assumption for this constructor

    // Extract our rank, the underlying group, and the number of processes
    mpi::CommDup( comm, viewingComm_ );
    mpi::CommGroup( viewingComm_, viewingGroup_ );
    viewingRank_ = mpi::CommRank( viewingComm_ );
    size_ = mpi::CommSize( viewingComm_ );

    // All processes own the grid, so we have to trivially split viewingGroup_
    owningGroup_ = viewingGroup_;
    notOwningGroup_ = mpi::GROUP_EMPTY;
    owningRank_ = viewingRank_;

    // Factor p
    height_ = FindFactor( size_ );
    width_ = size_ / height_;

    SetUpGrid();
}

inline 
Grid::Grid( mpi::Comm comm, int height )
{
#ifndef RELEASE
    CallStackEntry entry("Grid::Grid");
#endif
    inGrid_ = true; // this is true by assumption for this constructor

    // Extract our rank, the underlying group, and the number of processes
    mpi::CommDup( comm, viewingComm_ );
    mpi::CommGroup( viewingComm_, viewingGroup_ );
    viewingRank_ = mpi::CommRank( viewingComm_ );
    size_ = mpi::CommSize( viewingComm_ );

    // All processes own the grid, so we have to trivially split viewingGroup_
    owningGroup_ = viewingGroup_;
    notOwningGroup_ = mpi::GROUP_EMPTY;
    owningRank_ = viewingRank_;

    height_ = height;
    width_ = size_ /  height_;
    if( height_ < 0 )
        LogicError("Process grid dimensions must be non-negative");

    SetUpGrid();
}

inline void 
Grid::SetUpGrid()
{
#ifndef RELEASE
    CallStackEntry entry("Grid::SetUpGrid");
#endif
    if( size_ != height_*width_ )
    {
        std::ostringstream msg;
        msg << "Number of processes must match grid size:\n"
            << "  size=" << size_ << ", (height,width)=(" 
            << height_ << "," << width_ << ")";
        LogicError( msg.str() );
    }

    gcd_ = elem::GCD( height_, width_ );
    int lcm = size_ / gcd_;

    // Split the viewing comm into the owning and not owning subsets
    mpi::CommSplit( viewingComm_, inGrid_, owningRank_, owningComm_ );

    // Set up the map from the VC group to the viewingGroup_ ranks.
    // Since the VC communicator preserves the ordering of the owningGroup_
    // ranks, we can simply translate from owningGroup_.
    std::vector<int> ranks(size_);
    for( int i=0; i<size_; ++i )
        ranks[i] = i;
    vectorColToViewingMap_.resize(size_);
    mpi::GroupTranslateRanks
    ( owningGroup_, size_, ranks.data(), viewingGroup_, 
      vectorColToViewingMap_.data() );

    diagPathsAndRanks_.resize(2*size_);
    MemZero( diagPathsAndRanks_.data(), 2*size_ );
    if( inGrid_ )
    {
        // Create a cartesian communicator
        int dimensions[2] = { width_, height_ };
        int periods[2] = { true, true };
        bool reorder = false;
        mpi::CartCreate
        ( owningComm_, 2, dimensions, periods, reorder, cartComm_ );

        // Set up the MatrixCol and MatrixRow communicators
        int remainingDimensions[2];
        remainingDimensions[0] = false;
        remainingDimensions[1] = true;
        mpi::CartSub( cartComm_, remainingDimensions, matrixColComm_ );
        remainingDimensions[0] = true;
        remainingDimensions[1] = false;
        mpi::CartSub( cartComm_, remainingDimensions, matrixRowComm_ );
        matrixColRank_ = mpi::CommRank( matrixColComm_ );
        matrixRowRank_ = mpi::CommRank( matrixRowComm_ );

        // Set up the VectorCol and VectorRow communicators
        vectorColRank_ = matrixColRank_ + height_*matrixRowRank_;
        vectorRowRank_ = matrixRowRank_ + width_*matrixColRank_;
        mpi::CommSplit( cartComm_, 0, vectorColRank_, vectorColComm_ );
        mpi::CommSplit( cartComm_, 0, vectorRowRank_, vectorRowComm_ );

        // Compute which diagonal 'path' we're in, and what our rank is, then
        // perform AllGather world to store everyone's info
        std::vector<int> myDiagPathAndRank(2);
        myDiagPathAndRank[0] = (matrixRowRank_+height_-matrixColRank_) % gcd_;
        int diagPathRank = 0;
        int row = 0;
        int col = myDiagPathAndRank[0];
        for( int j=0; j<lcm; ++j )
        {
            if( row == matrixColRank_ && col == matrixRowRank_ )
            {
                myDiagPathAndRank[1] = diagPathRank;
                break;
            }
            else
            {
                row = (row + 1) % height_;
                col = (col + 1) % width_;
                ++diagPathRank;
            }
        }
        mpi::AllGather
        ( myDiagPathAndRank.data(), 2, 
          diagPathsAndRanks_.data(), 2, vectorColComm_ );

#ifndef RELEASE
        mpi::ErrorHandlerSet
        ( matrixColComm_, mpi::ERRORS_RETURN );
        mpi::ErrorHandlerSet
        ( matrixRowComm_, mpi::ERRORS_RETURN );
        mpi::ErrorHandlerSet
        ( vectorColComm_, mpi::ERRORS_RETURN );
        mpi::ErrorHandlerSet
        ( vectorRowComm_, mpi::ERRORS_RETURN );
#endif
    }
    else
    {
        // diag paths and ranks are implicitly set to undefined
        matrixColRank_ = mpi::UNDEFINED;
        matrixRowRank_ = mpi::UNDEFINED;
        vectorColRank_ = mpi::UNDEFINED;
        vectorRowRank_ = mpi::UNDEFINED;
    }
    mpi::Broadcast
    ( diagPathsAndRanks_.data(), 2*size_, 
      vectorColToViewingMap_[0], viewingComm_ );
}

inline 
Grid::~Grid()
{
    if( !mpi::Finalized() )
    {
        if( inGrid_ )
        {
            mpi::CommFree( matrixColComm_ );
            mpi::CommFree( matrixRowComm_ );
            mpi::CommFree( vectorColComm_ );
            mpi::CommFree( vectorRowComm_ );
            mpi::CommFree( cartComm_ );
        }

        mpi::CommFree( owningComm_ );
        if( notOwningGroup_ != mpi::GROUP_EMPTY )
            mpi::GroupFree( notOwningGroup_ );

        mpi::CommFree( viewingComm_ );
        mpi::GroupFree( viewingGroup_ );
    }
}

inline int 
Grid::MCRank() const
{ return matrixColRank_; }

inline int 
Grid::MRRank() const
{ return matrixRowRank_; }

inline int 
Grid::VCRank() const
{ return vectorColRank_; }

inline int 
Grid::VRRank() const
{ return vectorRowRank_; }

inline int 
Grid::MCSize() const
{ return height_; }

inline int 
Grid::MRSize() const
{ return width_; }

inline int 
Grid::VCSize() const
{ return size_; }

inline int 
Grid::VRSize() const
{ return size_; }

inline mpi::Comm 
Grid::MCComm() const
{ return matrixColComm_; }

inline mpi::Comm 
Grid::MRComm() const
{ return matrixRowComm_; }

inline mpi::Comm 
Grid::VCComm() const
{ return vectorColComm_; }

inline mpi::Comm 
Grid::VRComm() const
{ return vectorRowComm_; }

//
// Provided for simplicity, but redundant
//

inline int 
Grid::Row() const
{ return matrixColRank_; }

inline int 
Grid::Col() const
{ return matrixRowRank_; }

inline int 
Grid::Rank() const
{ return vectorColRank_; }

inline int 
Grid::Height() const
{ return height_; }

inline int 
Grid::Width() const
{ return width_; }

inline int 
Grid::Size() const
{ return size_; }

inline mpi::Comm 
Grid::ColComm() const
{ return matrixColComm_; }

inline mpi::Comm 
Grid::RowComm() const
{ return matrixRowComm_; }

inline mpi::Comm 
Grid::Comm() const
{ return vectorColComm_; }

//
// Advanced routines
//

// Currently forces a columnMajor absolute rank on the grid
inline 
Grid::Grid( mpi::Comm viewers, mpi::Group owners, int height )
{
#ifndef RELEASE
    CallStackEntry entry("Grid::Grid");
#endif

    // Extract our rank and the underlying group from the viewing comm
    mpi::CommDup( viewers, viewingComm_ );
    mpi::CommGroup( viewingComm_, viewingGroup_ );
    viewingRank_ = mpi::CommRank( viewingComm_ );

    // Extract our rank and the number of processes from the owning group
    owningGroup_ = owners;
    size_ = mpi::GroupSize( owningGroup_ );
    owningRank_ = mpi::GroupRank( owningGroup_ );
    inGrid_ = ( owningRank_ != mpi::UNDEFINED );

    // Create the complement of the owning group
    mpi::GroupDifference( viewingGroup_, owningGroup_, notOwningGroup_ );

    height_ = height;
    width_ = size_ / height;

    if( height_ < 0 )
        LogicError("Process grid dimensions must be non-negative");

    SetUpGrid();
}

inline int 
Grid::GCD() const
{ return gcd_; }

inline int 
Grid::LCM() const
{ return size_/gcd_; }

inline bool 
Grid::InGrid() const
{ return inGrid_; }

inline int 
Grid::OwningRank() const
{ return owningRank_; }

inline int 
Grid::ViewingRank() const
{ return viewingRank_; }

inline int 
Grid::VCToViewingMap( int VCRank ) const
{ return vectorColToViewingMap_[VCRank]; }

inline mpi::Group 
Grid::OwningGroup() const
{ return owningGroup_; }

inline mpi::Comm 
Grid::OwningComm() const
{ return owningComm_; }

inline mpi::Comm 
Grid::ViewingComm() const
{ return viewingComm_; }

inline int 
Grid::DiagPath() const
{ 
    if( inGrid_ )
        return diagPathsAndRanks_[2*vectorColRank_];
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
    if( inGrid_ )
        return diagPathsAndRanks_[2*vectorColRank_+1];
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

inline int
Grid::FirstVCRank( int diagPath ) const
{ return diagPath*height_; }

//
// Comparison functions
//

inline bool 
operator==( const Grid& A, const Grid& B )
{ return &A == &B; }

inline bool 
operator!=( const Grid& A, const Grid& B )
{ return &A != &B; }

} // namespace elem

#endif // ifndef ELEM_CORE_GRID_IMPL_HPP
