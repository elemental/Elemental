/*
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef ELEMENTAL_PROCESSGRID_HPP
#define ELEMENTAL_PROCESSGRID_HPP 1

namespace elemental {

class Grid
{
public:
    // Basic routines
    Grid( mpi::Comm comm=mpi::COMM_WORLD );
    Grid( mpi::Comm comm, int height, int width );
    ~Grid();
    bool InGrid() const;
    int Size() const;
    int Height() const;
    int Width() const;
    int GCD() const;
    int LCM() const;
    int Rank() const; // same as VCRank(), but provided for simplicity
    int MCRank() const;
    int MRRank() const;
    int VCRank() const;
    int VRRank() const;
    mpi::Comm Comm() const; // same as VCComm(), but provided for simplicity
    mpi::Comm MCComm() const;
    mpi::Comm MRComm() const;
    mpi::Comm VCComm() const;
    mpi::Comm VRComm() const;

    // Advanced routines
    Grid( mpi::Comm viewingComm, mpi::Group owningGroup );
    Grid( mpi::Comm viewingComm, mpi::Group owningGroup, 
          int height, int width );
    int OwningRank() const;
    int ViewingRank() const;
    int VCToViewingMap( int VCRank ) const;
    mpi::Group OwningGroup() const;
    mpi::Comm OwningComm() const;
    mpi::Comm ViewingComm() const;
    int DiagPath() const;
    int DiagPath( int vectorColRank ) const;
    int DiagPathRank() const;
    int DiagPathRank( int vectorColRank ) const;

private:
    int height_, width_, size_, gcd_;
    int matrixColRank_;
    int matrixRowRank_;
    int vectorColRank_;
    int vectorRowRank_;
    std::vector<int> diagPathsAndRanks_;

    mpi::Comm viewingComm_; // all processes that create the grid
    mpi::Group viewingGroup_;
    int viewingRank_; // our rank in the viewing communicator

    mpi::Group owningGroup_; // the processes that can own data
    mpi::Group notOwningGroup_; // contains the remaining processes

    std::vector<int> vectorColToViewingMap_;

    // Keep track of whether or not our process is in the grid. This is 
    // necessary to avoid calls like MPI_Comm_size when we're not in the
    // communicator's group. Note that we can__ call MPI_Group_rank when not 
    // in the group and that the result is MPI_UNDEFINED.
    bool inGrid_;

    // Create a communicator for the processes that are in the process grid
    mpi::Comm owningComm_;
    mpi::Comm notOwningComm_; // necessary complimentary communicator
    int owningRank_;

    // These will only be valid if we are in the grid
    mpi::Comm cartComm_;  // the processes that are in the grid
    mpi::Comm matrixColComm_;
    mpi::Comm matrixRowComm_;
    mpi::Comm vectorColComm_;
    mpi::Comm vectorRowComm_;

    void SetUpGrid();

    // Disable copying this class due to MPI_Comm/MPI_Group ownership issues
    // and potential performance loss from duplicating MPI communicators, e.g.,
    // on Blue Gene/P there is supposedly a performance loss
    const Grid& operator=( Grid& );
    Grid( const Grid& );
};

bool operator== ( const Grid& A, const Grid& B );
bool operator!= ( const Grid& A, const Grid& B );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline bool Grid::InGrid() const
{ return inGrid_; }

inline int Grid::Size() const
{ return size_; }

inline int Grid::Height() const
{ return height_; }

inline int Grid::Width() const
{ return width_; }

inline int Grid::GCD() const
{ return gcd_; }

inline int Grid::LCM() const
{ return size_/gcd_; }

inline int Grid::DiagPath() const
{ 
    if( inGrid_ )
        return diagPathsAndRanks_[2*vectorColRank_];
    else
        return mpi::UNDEFINED;
}

inline int Grid::DiagPath( int vectorColRank ) const
{ 
    if( vectorColRank != mpi::UNDEFINED )
        return diagPathsAndRanks_[2*vectorColRank]; 
    else
        return mpi::UNDEFINED;
}

inline int Grid::DiagPathRank() const
{ 
    if( inGrid_ )
        return diagPathsAndRanks_[2*vectorColRank_+1];
    else
        return mpi::UNDEFINED;
}

inline int Grid::DiagPathRank( int vectorColRank ) const
{ 
    if( vectorColRank != mpi::UNDEFINED )
        return diagPathsAndRanks_[2*vectorColRank+1]; 
    else
        return mpi::UNDEFINED;
}

inline int Grid::Rank() const
{ return vectorColRank_; }

inline int Grid::MCRank() const
{ return matrixColRank_; }

inline int Grid::MRRank() const
{ return matrixRowRank_; }

inline int Grid::VCRank() const
{ return vectorColRank_; }

inline int Grid::VRRank() const
{ return vectorRowRank_; }

inline int Grid::OwningRank() const
{ return owningRank_; }

inline int Grid::ViewingRank() const
{ return viewingRank_; }

inline int Grid::VCToViewingMap( int VCRank ) const
{ return vectorColToViewingMap_[VCRank]; }

inline mpi::Group Grid::OwningGroup() const
{ return owningGroup_; }

inline mpi::Comm Grid::OwningComm() const
{ return owningComm_; }

inline mpi::Comm Grid::ViewingComm() const
{ return viewingComm_; }

inline mpi::Comm Grid::Comm() const
{ return vectorColComm_; }

inline mpi::Comm Grid::MCComm() const
{ return matrixColComm_; }

inline mpi::Comm Grid::MRComm() const
{ return matrixRowComm_; }

inline mpi::Comm Grid::VCComm() const
{ return vectorColComm_; }

inline mpi::Comm Grid::VRComm() const
{ return vectorRowComm_; }

inline bool operator== ( const Grid& A, const Grid& B )
{ return &A == &B; }

inline bool operator!= ( const Grid& A, const Grid& B )
{ return &A != &B; }

inline Grid::Grid( mpi::Comm comm )
{
#ifndef RELEASE
    PushCallStack("Grid::Grid");
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
    height_ = static_cast<int>(sqrt(static_cast<double>(size_)));
    while( size_ % height_ != 0 )
        ++height_;
    width_ = size_ / height_;

    SetUpGrid();

#ifndef RELEASE
    PopCallStack();
#endif
}

inline Grid::Grid( mpi::Comm comm, int height, int width )
{
#ifndef RELEASE
    PushCallStack("Grid::Grid");
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
    width_ = width;
    if( height_ < 0 || width_ < 0 )
        throw std::logic_error
        ("Process grid dimensions must be non-negative");

    SetUpGrid();

#ifndef RELEASE
    PopCallStack();
#endif
}

inline Grid::Grid( mpi::Comm viewingComm, mpi::Group owningGroup )
{
#ifndef RELEASE
    PushCallStack("Grid::Grid");
#endif

    // Extract our rank and the underlying group from the viewing comm
    mpi::CommDup( viewingComm, viewingComm_ );
    mpi::CommGroup( viewingComm_, viewingGroup_ );
    viewingRank_ = mpi::CommRank( viewingComm_ );

    // Extract our rank and the number of processes from the owning group
    owningGroup_ = owningGroup;
    size_ = mpi::GroupSize( owningGroup_ );
    owningRank_ = mpi::GroupRank( owningGroup_ );
    inGrid_ = ( owningRank_ != mpi::UNDEFINED );

    // Create the complement of the owning group
    mpi::GroupDifference( viewingGroup_, owningGroup_, notOwningGroup_ );

    // Factor the grid size
    height_ = static_cast<int>(sqrt(static_cast<double>(size_)));
    while( size_ % height_ != 0 )
        ++height_;
    width_ = size_ / height_;

    SetUpGrid();

#ifndef RELEASE
    PopCallStack();
#endif
}

// Currently forces a columnMajor absolute rank on the grid
inline Grid::Grid
( mpi::Comm viewingComm, mpi::Group owningGroup, int height, int width )
{
#ifndef RELEASE
    PushCallStack("Grid::Grid");
#endif

    // Extract our rank and the underlying group from the viewing comm
    mpi::CommDup( viewingComm, viewingComm_ );
    mpi::CommGroup( viewingComm_, viewingGroup_ );
    viewingRank_ = mpi::CommRank( viewingComm_ );

    // Extract our rank and the number of processes from the owning group
    owningGroup_ = owningGroup;
    size_ = mpi::GroupSize( owningGroup_ );
    owningRank_ = mpi::GroupRank( owningGroup_ );
    inGrid_ = ( owningRank_ != mpi::UNDEFINED );

    // Create the complement of the owning group
    mpi::GroupDifference( viewingGroup_, owningGroup_, notOwningGroup_ );

    height_ = height;
    width_ = width;

    if( height_ < 0 || width_ < 0 )
        throw std::logic_error
        ("Process grid dimensions must be non-negative");

    SetUpGrid();

#ifndef RELEASE
    PopCallStack();
#endif
}

inline void Grid::SetUpGrid()
{
#ifndef RELEASE
    PushCallStack("Grid::SetUpGrid");
#endif
    if( size_ != height_*width_ )
    {
        std::ostringstream msg;
        msg << "Number of processes must match grid size:\n"
            << "  size=" << size_ << ", (height,width)=(" 
            << height_ << "," << width_ << ")";
        throw std::logic_error( msg.str().c_str() );
    }

    gcd_ = elemental::GCD( height_, width_ );
    int lcm = size_ / gcd_;

    // Split the viewing comm into the owning and not owning subsets
    if( inGrid_ )
        mpi::CommSplit( viewingComm_, true, owningRank_, owningComm_ );
    else
        mpi::CommSplit( viewingComm_, false, 0, notOwningComm_ );

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
        diagPathsAndRanks_.resize(2*size_);
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
        ( &myDiagPathAndRank[0], 2, &diagPathsAndRanks_[0], 2, vectorColComm_ );

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

    // Set up the map from the VC group to the viewingGroup_ ranks.
    // Since the VC communicator preserves the ordering of the owningGroup_
    // ranks, we can simply translate from owningGroup_.
    std::vector<int> ranks(size_);
    for( int i=0; i<size_; ++i )
        ranks[i] = i;
    vectorColToViewingMap_.resize(size_);
    mpi::GroupTranslateRanks
    ( owningGroup_, size_, &ranks[0], viewingGroup_, 
      &vectorColToViewingMap_[0] );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline Grid::~Grid()
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

        if( inGrid_ )
            mpi::CommFree( owningComm_ );
        else
            mpi::CommFree( notOwningComm_ );

        if( notOwningGroup_ != mpi::GROUP_EMPTY )
            mpi::GroupFree( notOwningGroup_ );

        mpi::CommFree( viewingComm_ );
    }
}

} // namespace elemental

#endif /* ELEMENTAL_PROCESSGRID_HPP */

