/*
   Copyright (c) 2009-2015, Jack Poulson
                      2013, Jed Brown 
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

int Grid::FindFactor( int p )
{
    int factor = int(sqrt(double(p)));
    while( p % factor != 0 )
        ++factor;
    return factor;
}

Grid::Grid( mpi::Comm comm, GridOrder order )
: haveViewers_(false), order_(order)
{
    DEBUG_ONLY(CSE cse("Grid::Grid"))

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
    DEBUG_ONLY(CSE cse("Grid::Grid"))

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
    DEBUG_ONLY(CSE cse("Grid::SetUpGrid"))
    if( size_ % height_ != 0 )
        LogicError
        ("Grid height, ",height_,", does not evenly divide grid size, ",size_);

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
        const int mcRank = mpi::Rank( mcComm_ );
        const int mrRank = mpi::Rank( mrComm_ );

        // Set up the VectorCol and VectorRow communicators
        const int vcRank = mcRank + height_*mrRank;
        const int vrRank = mrRank + width*mcRank;
        mpi::Split( cartComm_, 0, vcRank, vcComm_ );
        mpi::Split( cartComm_, 0, vrRank, vrComm_ );

        // Set up the map from the VC group to the viewingGroup_ ranks.
        mpi::Group vcGroup;
        mpi::CommGroup( vcComm_, vcGroup ); 
        vector<int> ranks(size_);
        for( int i=0; i<size_; ++i )
            ranks[i] = i;
        mpi::Translate
        ( vcGroup, size_, ranks.data(), viewingGroup_,  vcToViewing_.data() );
        mpi::Free( vcGroup );

        // Compute which diagonal we're in, and what our rank is, then
        // perform AllGather world to store everyone's info
        vector<int> myDiagAndRank(2);
        myDiagAndRank[0] = Mod(mrRank-mcRank,gcd_);
        int diagRank = 0;
        int row = 0;
        int col = myDiagAndRank[0];
        for( int j=0; j<lcm; ++j )
        {
            if( row == mcRank && col == mrRank )
            {
                myDiagAndRank[1] = diagRank;
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

        mpi::Split( cartComm_, Diag(),     DiagRank(), mdComm_     );
        mpi::Split( cartComm_, DiagRank(), Diag(),     mdPerpComm_ );

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
        vcComm_     = mpi::COMM_NULL;
        vrComm_     = mpi::COMM_NULL;
        mdComm_     = mpi::COMM_NULL; 
        mdPerpComm_ = mpi::COMM_NULL;
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

int Grid::MCRank() const { return mpi::Rank(mcComm_); }
int Grid::MRRank() const { return mpi::Rank(mrComm_); }
int Grid::VCRank() const { return mpi::Rank(vcComm_); }
int Grid::VRRank() const { return mpi::Rank(vrComm_); }

int Grid::MCSize() const { return height_;       }
int Grid::MRSize() const { return size_/height_; }
int Grid::VCSize() const { return size_;         }
int Grid::VRSize() const { return size_;         }

mpi::Comm Grid::MCComm()     const { return mcComm_;     }
mpi::Comm Grid::MRComm()     const { return mrComm_;     }
mpi::Comm Grid::VCComm()     const { return vcComm_;     }
mpi::Comm Grid::VRComm()     const { return vrComm_;     }
mpi::Comm Grid::MDComm()     const { return mdComm_;     }
mpi::Comm Grid::MDPerpComm() const { return mdPerpComm_; }

// Provided for simplicity, but redundant
// ======================================
int Grid::Height() const { return MCSize(); }
int Grid::Width()  const { return MRSize(); }
int Grid::Size()   const { return VCSize(); }
int Grid::Rank() const { return OwningRank(); }

GridOrder Grid::Order() const { return order_; }

int Grid::Row() const { return MCRank(); }
int Grid::Col() const { return MRRank(); }
mpi::Comm Grid::ColComm() const { return MCComm(); }
mpi::Comm Grid::RowComm() const { return MRComm(); }
mpi::Comm Grid::Comm() const
{ return ( order_==COLUMN_MAJOR ? VCComm() : VRComm() ); }

// Advanced routines
// =================

// Currently forces a columnMajor absolute rank on the grid
Grid::Grid( mpi::Comm viewers, mpi::Group owners, int height, GridOrder order )
: haveViewers_(true), order_(order)
{
    DEBUG_ONLY(CSE cse("Grid::Grid"))

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

int Grid::GCD() const { return gcd_; }
int Grid::LCM() const { return size_/gcd_; }

bool Grid::HaveViewers() const { return haveViewers_; }
bool Grid::InGrid() const 
{ return mpi::Rank(owningGroup_) != mpi::UNDEFINED; }

int Grid::OwningRank() const { return mpi::Rank(owningGroup_); }
int Grid::ViewingRank() const { return mpi::Rank(viewingComm_); }

int Grid::VCToVR( int vcRank ) const
{
    const int height = Height();
    const int width = Width();
    const int mcRank = vcRank % height;
    const int mrRank = vcRank / height;
    return mrRank + mcRank*width;
}

int Grid::VRToVC( int vrRank ) const
{
    const int height = Height();
    const int width = Width();
    const int mcRank = vrRank / width;
    const int mrRank = vrRank % width;
    return mcRank + mrRank*height;
}

int Grid::CoordsToVC
( Dist colDist, Dist rowDist, 
  int distRank, int crossRank, int redundantRank ) const
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
        return distRank + redundantRank*Width();
    }
    else if( (colDist == MD && rowDist == STAR) ||
             (rowDist == MD && colDist == STAR) )
    {
        const int row =         DiagRank()  % Height();
        const int col = (Diag()+DiagRank()) % Width();
        return row + col*Height();
    }
    else if( colDist == MR && rowDist == MC )
    {
        return VRToVC(distRank);
    }
    else if( (colDist == MR && rowDist == STAR) ||
             (rowDist == MR && colDist == STAR) )
    {
        return redundantRank + distRank*Width();
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
    else
        LogicError("Invalid data distribution");
    return -1;
}

int Grid::VCToViewing( int vcRank ) const
{ return vcToViewing_[vcRank]; }

mpi::Group Grid::OwningGroup() const { return owningGroup_; }
mpi::Comm Grid::OwningComm()  const { return owningComm_; }
mpi::Comm Grid::ViewingComm() const { return viewingComm_; }

int Grid::Diag() const
{ 
    const int vcRank = VCRank();
    if( vcRank != mpi::UNDEFINED )
        return diagsAndRanks_[2*vcRank];
    else
        return mpi::UNDEFINED;
}

int Grid::Diag( int vcRank ) const
{ 
    if( vcRank != mpi::UNDEFINED )
        return diagsAndRanks_[2*vcRank]; 
    else
        return mpi::UNDEFINED;
}

int Grid::DiagRank() const
{ 
    const int vcRank = VCRank();
    if( vcRank != mpi::UNDEFINED )
        return diagsAndRanks_[2*vcRank+1];
    else
        return mpi::UNDEFINED;
}

int Grid::DiagRank( int vcRank ) const
{ 
    if( vcRank != mpi::UNDEFINED )
        return diagsAndRanks_[2*vcRank+1]; 
    else
        return mpi::UNDEFINED;
}

// Comparison functions
// ====================

bool operator==( const Grid& A, const Grid& B ) { return &A == &B; }
bool operator!=( const Grid& A, const Grid& B ) { return &A != &B; }

} // namespace El
