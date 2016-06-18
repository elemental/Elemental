/*
   Copyright (c) 2009-2016, Jack Poulson
                      2013, Jed Brown 
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>

namespace El {

Grid* Grid::defaultGrid = 0;

void Grid::InitializeDefault()
{
    defaultGrid = new Grid( mpi::COMM_WORLD );    
}

void Grid::FinalizeDefault()
{
    delete defaultGrid;
    defaultGrid = 0;
}

const Grid& Grid::Default() EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( defaultGrid == 0 )
          LogicError
          ("Attempted to return a non-existant default grid. Please ensure "
           "that Elemental is initialized before creating a DistMatrix.");
    )
    return *defaultGrid;
}

int Grid::FindFactor( int p ) EL_NO_EXCEPT
{
    int factor = int(sqrt(double(p)));
    while( p % factor != 0 )
        ++factor;
    return factor;
}

Grid::Grid( mpi::Comm comm, GridOrder order )
: haveViewers_(false), order_(order)
{
    DEBUG_CSE

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

Grid::Grid( mpi::Comm comm, int height, GridOrder order )
: haveViewers_(false), order_(order)
{
    DEBUG_CSE

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

void Grid::SetUpGrid()
{
    DEBUG_CSE
    if( size_ % height_ != 0 )
        LogicError
        ("Grid height, ",height_,", does not evenly divide grid size, ",size_);
    owningRank_ = mpi::Rank( owningGroup_ );
    viewingRank_ = mpi::Rank( viewingComm_ );
    inGrid_ = ( owningRank_ != mpi::UNDEFINED );

    const int width = size_ / height_;
    gcd_ = El::GCD( height_, width );
    int lcm = size_ / gcd_;

    // Create the communicator for the owning group (mpi::COMM_NULL otherwise)
    mpi::Create( viewingComm_, owningGroup_, owningComm_ );

    vcToViewing_.resize(size_);
    diagsAndRanks_.resize(2*size_);
    MemZero( diagsAndRanks_.data(), 2*size_ );
    const bool colMajor = (order_==COLUMN_MAJOR);
    if( InGrid() )
    {
        // Create a cartesian communicator
        // TODO: Allow for reordering and non-periodicity
        int dims[2];
        if( colMajor )
        {
            dims[0] = width;
            dims[1] = height_;
        }
        else
        {
            dims[0] = height_;
            dims[1] = width;
        }
        int periods[2] = { true, true };
        bool reorder = false;
        mpi::CartCreate( owningComm_, 2, dims, periods, reorder, cartComm_ );

        // Set up the MatrixCol and MatrixRow communicators
        int remainingDims[2];
        remainingDims[0] = ( colMajor ? false : true  );
        remainingDims[1] = ( colMajor ? true  : false );
        mpi::CartSub( cartComm_, remainingDims, mcComm_ );
        remainingDims[0] = ( colMajor ? true  : false );
        remainingDims[1] = ( colMajor ? false : true  );
        mpi::CartSub( cartComm_, remainingDims, mrComm_ );
        mcRank_ = mpi::Rank( mcComm_ );
        mrRank_ = mpi::Rank( mrComm_ );

        // Set up the VectorCol and VectorRow communicators
        vcRank_ = mcRank_ + height_*mrRank_;
        vrRank_ = mrRank_ + width*mcRank_;
        mpi::Split( cartComm_, 0, vcRank_, vcComm_ );
        mpi::Split( cartComm_, 0, vrRank_, vrComm_ );

        // Set up the map from the VC group to the viewingGroup_ ranks.
        vector<int> ranks(size_);
        for( int i=0; i<size_; ++i )
            ranks[i] = i;
        mpi::Translate
        ( vcComm_, size_, ranks.data(), viewingGroup_,  vcToViewing_.data() );

        // Compute which diagonal we're in, and what our rank is, then
        // perform AllGather world to store everyone's info
        vector<int> myDiagAndRank(2);
        myDiagAndRank[0] = mdPerpRank_ = Mod(mrRank_-mcRank_,gcd_);
        int diagRank = 0;
        int row = 0;
        int col = myDiagAndRank[0];
        for( int j=0; j<lcm; ++j )
        {
            if( row == mcRank_ && col == mrRank_ )
            {
                myDiagAndRank[1] = mdRank_ = diagRank;
                break;
            }
            else
            {
                row = (row + 1) % height_;
                col = (col + 1) % width;
                ++diagRank;
            }
        }
        mpi::AllGather
        ( myDiagAndRank.data(),  2, 
          diagsAndRanks_.data(), 2, vcComm_ );

        mpi::Split( cartComm_, mdPerpRank_, mdRank_,     mdComm_     );
        mpi::Split( cartComm_, mdRank_,     mdPerpRank_, mdPerpComm_ );

        DEBUG_ONLY(
          mpi::ErrorHandlerSet( mcComm_,     mpi::ERRORS_RETURN );
          mpi::ErrorHandlerSet( mrComm_,     mpi::ERRORS_RETURN );
          mpi::ErrorHandlerSet( vcComm_,     mpi::ERRORS_RETURN );
          mpi::ErrorHandlerSet( vrComm_,     mpi::ERRORS_RETURN );
          mpi::ErrorHandlerSet( mdComm_,     mpi::ERRORS_RETURN );
          mpi::ErrorHandlerSet( mdPerpComm_, mpi::ERRORS_RETURN );
        )
    }
    else
    {
        mcComm_     = mpi::COMM_NULL;
        mrComm_     = mpi::COMM_NULL;
        mdComm_     = mpi::COMM_NULL; 
        mdPerpComm_ = mpi::COMM_NULL;
        vcComm_     = mpi::COMM_NULL;
        vrComm_     = mpi::COMM_NULL;

        mcRank_     = mpi::UNDEFINED;
        mrRank_     = mpi::UNDEFINED;
        mdRank_     = mpi::UNDEFINED;
        mdPerpRank_ = mpi::UNDEFINED;
        vcRank_     = mpi::UNDEFINED;
        vrRank_     = mpi::UNDEFINED;

        // diags and ranks are implicitly set to undefined
    }
    // Translate the rank of the root process of the owningGroup so that we can
    // broadcast data
    int owningRoot = mpi::Translate( owningGroup_, 0, viewingGroup_ );
    mpi::Broadcast( vcToViewing_.data(), size_, owningRoot, viewingComm_ );
    mpi::Broadcast( diagsAndRanks_.data(), 2*size_, owningRoot, viewingComm_ );
}

Grid::~Grid()
{
    if( !mpi::Finalized() )
    {
        if( InGrid() )
        {
            mpi::Free( mdComm_ );
            mpi::Free( mdPerpComm_ );
            mpi::Free( mcComm_ );
            mpi::Free( mrComm_ );
            mpi::Free( vcComm_ );
            mpi::Free( vrComm_ );
            mpi::Free( cartComm_ );
            mpi::Free( owningComm_ );
        }
        mpi::Free( viewingComm_ );
        if( HaveViewers() )
            mpi::Free( owningGroup_ );
        mpi::Free( viewingGroup_ );
    }
}

int Grid::MCRank()     const EL_NO_RELEASE_EXCEPT { return mcRank_; }
int Grid::MRRank()     const EL_NO_RELEASE_EXCEPT { return mrRank_; }
int Grid::MDRank()     const EL_NO_RELEASE_EXCEPT { return mdRank_; }
int Grid::MDPerpRank() const EL_NO_RELEASE_EXCEPT { return mdPerpRank_; }
int Grid::VCRank()     const EL_NO_RELEASE_EXCEPT { return vcRank_; }
int Grid::VRRank()     const EL_NO_RELEASE_EXCEPT { return vrRank_; }

int Grid::MCSize()     const EL_NO_EXCEPT { return height_;       }
int Grid::MRSize()     const EL_NO_EXCEPT { return size_/height_; }
int Grid::MDSize()     const EL_NO_EXCEPT { return size_/gcd_;    }
int Grid::MDPerpSize() const EL_NO_EXCEPT { return gcd_;          }
int Grid::VCSize()     const EL_NO_EXCEPT { return size_;         }
int Grid::VRSize()     const EL_NO_EXCEPT { return size_;         }

mpi::Comm Grid::MCComm()     const EL_NO_EXCEPT { return mcComm_;     }
mpi::Comm Grid::MRComm()     const EL_NO_EXCEPT { return mrComm_;     }
mpi::Comm Grid::MDComm()     const EL_NO_EXCEPT { return mdComm_;     }
mpi::Comm Grid::MDPerpComm() const EL_NO_EXCEPT { return mdPerpComm_; }
mpi::Comm Grid::VCComm()     const EL_NO_EXCEPT { return vcComm_;     }
mpi::Comm Grid::VRComm()     const EL_NO_EXCEPT { return vrComm_;     }

// Provided for simplicity, but redundant
// ======================================
int Grid::Height() const EL_NO_EXCEPT { return MCSize(); }
int Grid::Width()  const EL_NO_EXCEPT { return MRSize(); }
int Grid::Size()   const EL_NO_EXCEPT { return VCSize(); }
int Grid::Rank()   const EL_NO_RELEASE_EXCEPT { return OwningRank(); }

GridOrder Grid::Order() const EL_NO_EXCEPT { return order_; }

int Grid::Row() const EL_NO_RELEASE_EXCEPT { return MCRank(); }
int Grid::Col() const EL_NO_RELEASE_EXCEPT { return MRRank(); }
mpi::Comm Grid::ColComm() const EL_NO_EXCEPT { return MCComm(); }
mpi::Comm Grid::RowComm() const EL_NO_EXCEPT { return MRComm(); }
mpi::Comm Grid::Comm() const EL_NO_EXCEPT
{ return ( order_==COLUMN_MAJOR ? VCComm() : VRComm() ); }

// Advanced routines
// =================

// Currently forces a columnMajor absolute rank on the grid
Grid::Grid( mpi::Comm viewers, mpi::Group owners, int height, GridOrder order )
: haveViewers_(true), order_(order)
{
    DEBUG_CSE

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

int Grid::GCD() const EL_NO_EXCEPT { return gcd_; }
int Grid::LCM() const EL_NO_EXCEPT { return size_/gcd_; }

bool Grid::HaveViewers() const EL_NO_EXCEPT { return haveViewers_; }
bool Grid::InGrid() const EL_NO_RELEASE_EXCEPT { return inGrid_; }

int Grid::OwningRank() const EL_NO_RELEASE_EXCEPT { return owningRank_; }
int Grid::ViewingRank() const EL_NO_RELEASE_EXCEPT { return viewingRank_; }

int Grid::VCToVR( int vcRank ) const EL_NO_EXCEPT
{
    const int height = Height();
    const int width = Width();
    const int mcRank = vcRank % height;
    const int mrRank = vcRank / height;
    return mrRank + mcRank*width;
}

int Grid::VRToVC( int vrRank ) const EL_NO_EXCEPT
{
    const int height = Height();
    const int width = Width();
    const int mcRank = vrRank / width;
    const int mrRank = vrRank % width;
    return mcRank + mrRank*height;
}

int Grid::CoordsToVC
( Dist colDist, Dist rowDist, 
  int distRank, int crossRank, int redundantRank ) const EL_NO_RELEASE_EXCEPT
{
    if( colDist == CIRC && rowDist == CIRC )
    {
        return crossRank;
    }
    else if( colDist == MC && rowDist == MR )
    {
        return distRank;
    }
    else if( (colDist == MC && rowDist == STAR) || 
             (rowDist == MC && colDist == STAR) )
    {
        return distRank + redundantRank*Height();
    }
    else if( (colDist == MD && rowDist == STAR) ||
             (rowDist == MD && colDist == STAR) )
    {
        const int row = distRank % Height();
        const int col = (crossRank+distRank) % Width();
        return row + col*Height();
    }
    else if( colDist == MR && rowDist == MC )
    {
        return VRToVC(distRank);
    }
    else if( (colDist == MR && rowDist == STAR) ||
             (rowDist == MR && colDist == STAR) )
    {
        return redundantRank + distRank*Height();
    }
    else if( colDist == STAR && rowDist == STAR )
    {
        return redundantRank;
    }
    else if( (colDist == STAR && rowDist == VC) ||
             (rowDist == STAR && colDist == VC) )
    {
        return distRank;
    }
    else if( (colDist == STAR && rowDist == VR) ||
             (rowDist == STAR && colDist == VR) )
    {
        return VRToVC(distRank);
    }
    DEBUG_ONLY(else LogicError("Invalid data distribution"))
    return -1;
}

int Grid::VCToViewing( int vcRank ) const EL_NO_EXCEPT
{ return vcToViewing_[vcRank]; }

mpi::Group Grid::OwningGroup() const EL_NO_EXCEPT { return owningGroup_; }
mpi::Comm Grid::OwningComm()  const EL_NO_EXCEPT { return owningComm_; }
mpi::Comm Grid::ViewingComm() const EL_NO_EXCEPT { return viewingComm_; }

int Grid::Diag() const EL_NO_RELEASE_EXCEPT
{ 
    const int vcRank = VCRank();
    if( vcRank != mpi::UNDEFINED )
        return diagsAndRanks_[2*vcRank];
    else
        return mpi::UNDEFINED;
}

int Grid::Diag( int vcRank ) const EL_NO_EXCEPT
{ 
    if( vcRank != mpi::UNDEFINED )
        return diagsAndRanks_[2*vcRank]; 
    else
        return mpi::UNDEFINED;
}

int Grid::DiagRank() const EL_NO_RELEASE_EXCEPT
{ 
    const int vcRank = VCRank();
    if( vcRank != mpi::UNDEFINED )
        return diagsAndRanks_[2*vcRank+1];
    else
        return mpi::UNDEFINED;
}

int Grid::DiagRank( int vcRank ) const EL_NO_EXCEPT
{ 
    if( vcRank != mpi::UNDEFINED )
        return diagsAndRanks_[2*vcRank+1]; 
    else
        return mpi::UNDEFINED;
}

// Comparison functions
// ====================

bool operator==( const Grid& A, const Grid& B ) EL_NO_EXCEPT
{ return &A == &B; }
bool operator!=( const Grid& A, const Grid& B ) EL_NO_EXCEPT
{ return &A != &B; }

} // namespace El
