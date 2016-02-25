/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_PSEUDOSPECTRA_UTIL_SNAPSHOT_HPP
#define EL_PSEUDOSPECTRA_UTIL_SNAPSHOT_HPP

namespace El {

namespace pspec {

template<typename Real>
inline void
Snapshot
( const Matrix<Int>& preimage,
  const Matrix<Real>& estimates, 
  const Matrix<Int>& itCounts,
        Int numIts,
        bool deflate,
        SnapshotCtrl& snapCtrl )
{
    DEBUG_ONLY(CSE cse("pspec::Snapshot"));
    auto logMap = []( Real alpha ) { return Log(alpha); };
    if( snapCtrl.realSize != 0 && snapCtrl.imagSize != 0 )
    {
        const bool numSave = 
            ( snapCtrl.numSaveFreq > 0 && 
              snapCtrl.numSaveCount >= snapCtrl.numSaveFreq );
        const bool imgSave = 
            ( snapCtrl.imgSaveFreq > 0 && 
              snapCtrl.imgSaveCount >= snapCtrl.imgSaveFreq );
        const bool imgDisp = 
            ( snapCtrl.imgDispFreq > 0 &&
              snapCtrl.imgDispCount >= snapCtrl.imgDispFreq );
        Matrix<Real> invNorms, estMap;
        Matrix<Int> itCountsReord, itCountMap;
        if( numSave || imgSave || imgDisp )
        {
            invNorms = estimates;
            if( deflate )
                RestoreOrdering( preimage, invNorms );
            ReshapeIntoGrid
            ( snapCtrl.realSize, snapCtrl.imagSize, invNorms, estMap );
            if( snapCtrl.itCounts )
            {
                itCountsReord = itCounts;
                if( deflate )
                    RestoreOrdering( preimage, itCountsReord );
                ReshapeIntoGrid
                ( snapCtrl.realSize, snapCtrl.imagSize, itCountsReord, 
                  itCountMap );
            }
        }
        if( numSave )
        {
            auto title = BuildString( snapCtrl.numBase, "_", numIts );
            Write( estMap, title, snapCtrl.numFormat );
            if( snapCtrl.itCounts )
                Write( itCountMap, title+"_counts", snapCtrl.numFormat );
            snapCtrl.numSaveCount = 0;
        }
        if( imgSave || imgDisp )
            EntrywiseMap( estMap, function<Real(Real)>(logMap) );
        if( imgSave )
        {
            auto title = BuildString( snapCtrl.imgBase, "_", numIts );
            Write( estMap, title, snapCtrl.imgFormat );
            if( snapCtrl.itCounts )
                Write( itCountMap, title+"_counts", snapCtrl.imgFormat );
            auto colorMap = GetColorMap();
            SetColorMap( GRAYSCALE_DISCRETE );
            Write( estMap, title+"_discrete", snapCtrl.imgFormat );
            SetColorMap( colorMap );
            snapCtrl.imgSaveCount = 0;
        }
        if( imgDisp )
        {
            auto title = BuildString( snapCtrl.imgBase, "_", numIts );
            Display( estMap, title );       
            if( snapCtrl.itCounts )
                Display( itCountMap, title+"_counts" );
            auto colorMap = GetColorMap();
            SetColorMap( GRAYSCALE_DISCRETE );
            Display( estMap, title+"_discrete" );
            SetColorMap( colorMap );
            snapCtrl.imgDispCount = 0;
        }
    }
}

template<typename Real>
inline void
FinalSnapshot
( const Matrix<Real>& estimates,
  const Matrix<Int>& itCounts, 
        SnapshotCtrl& snapCtrl )
{
    DEBUG_ONLY(CSE cse("pspec::FinalSnapshot"));
    auto logMap = []( Real alpha ) { return Log(alpha); };
    if( snapCtrl.realSize != 0 && snapCtrl.imagSize != 0 )
    {
        const bool numSave = ( snapCtrl.numSaveFreq >= 0 );
        const bool imgSave = ( snapCtrl.imgSaveFreq >= 0 );
        const bool imgDisp = ( snapCtrl.imgDispFreq >= 0 );
        Matrix<Real> estMap;
        Matrix<Int> itCountMap;
        if( numSave || imgSave || imgDisp )
        {
            ReshapeIntoGrid
            ( snapCtrl.realSize, snapCtrl.imagSize, estimates, estMap );
            if( snapCtrl.itCounts )
                ReshapeIntoGrid
                ( snapCtrl.realSize, snapCtrl.imagSize, itCounts, itCountMap );
        }
        if( numSave )
        {
            string base = snapCtrl.numBase;
            Write( estMap, base, snapCtrl.numFormat );
            if( snapCtrl.itCounts )
                Write( itCountMap, base+"_counts", snapCtrl.numFormat );
        }
        if( imgSave || imgDisp )
            EntrywiseMap( estMap, function<Real(Real)>(logMap) );
        if( imgSave )
        {
            string base = snapCtrl.imgBase;
            Write( estMap, base, snapCtrl.imgFormat );
            if( snapCtrl.itCounts )
                Write( itCountMap, base+"_counts", snapCtrl.imgFormat );
            auto colorMap = GetColorMap();
            SetColorMap( GRAYSCALE_DISCRETE );
            Write( estMap, base+"_discrete", snapCtrl.imgFormat );
            SetColorMap( colorMap );
        }
        if( imgDisp )
        {
            string base = snapCtrl.imgBase;
            Display( estMap, base );       
            if( snapCtrl.itCounts )
                Display( itCountMap, base+"_counts" );
            auto colorMap = GetColorMap();
            SetColorMap( GRAYSCALE_DISCRETE );
            Display( estMap, base+"_discrete" );
            SetColorMap( colorMap );
        }
    }
}

template<typename Real>
inline void
Snapshot
( const DistMatrix<Int,    VR,STAR>& preimage, 
  const DistMatrix<Real,MR,STAR>& estimates, 
  const DistMatrix<Int, VR,STAR>& itCounts,
        Int numIts,
        bool deflate,
        SnapshotCtrl& snapCtrl )
{
    DEBUG_ONLY(CSE cse("pspec::Snapshot"));
    auto logMap = []( Real alpha ) { return Log(alpha); };
    if( snapCtrl.realSize != 0 && snapCtrl.imagSize != 0 )
    {
        const bool numSave = 
            ( snapCtrl.numSaveFreq > 0 && 
              snapCtrl.numSaveCount >= snapCtrl.numSaveFreq );
        const bool imgSave = 
            ( snapCtrl.imgSaveFreq > 0 && 
              snapCtrl.imgSaveCount >= snapCtrl.imgSaveFreq );
        const bool imgDisp =
            ( snapCtrl.imgDispFreq > 0 &&
              snapCtrl.imgDispCount >= snapCtrl.imgDispFreq );
        DistMatrix<Real,VR,STAR> invNorms(estimates.Grid());
        DistMatrix<Real> estMap(estimates.Grid()); 
        DistMatrix<Int, VR,STAR> itCountsReord(itCounts.Grid());
        DistMatrix<Int> itCountMap(itCounts.Grid());
        if( numSave || imgSave || imgDisp )
        {
            invNorms = estimates;
            if( deflate )
                RestoreOrdering( preimage, invNorms );
            ReshapeIntoGrid
            ( snapCtrl.realSize, snapCtrl.imagSize, invNorms, estMap );
            if( snapCtrl.itCounts )
            {
                itCountsReord = itCounts;
                if( deflate )
                    RestoreOrdering( preimage, itCountsReord );
                ReshapeIntoGrid
                ( snapCtrl.realSize, snapCtrl.imagSize, itCountsReord, 
                  itCountMap );
            }
        }
        if( numSave )
        {
            auto title = BuildString( snapCtrl.numBase, "_", numIts );
            Write( estMap, title, snapCtrl.numFormat );
            if( snapCtrl.itCounts )
                Write( itCountMap, title+"_counts", snapCtrl.numFormat );
            snapCtrl.numSaveCount = 0;
        }
        if( imgSave || imgDisp )
            EntrywiseMap( estMap, function<Real(Real)>(logMap) );
        if( imgSave )
        {
            auto title = BuildString( snapCtrl.imgBase, "_", numIts );
            Write( estMap, title, snapCtrl.imgFormat );
            if( snapCtrl.itCounts )
                Write( itCountMap, title+"_counts", snapCtrl.imgFormat );
            auto colorMap = GetColorMap();
            SetColorMap( GRAYSCALE_DISCRETE );
            Write( estMap, title+"_discrete", snapCtrl.imgFormat );
            SetColorMap( colorMap );
            snapCtrl.imgSaveCount = 0;
        }
        if( imgDisp )
        {
            auto title = BuildString( snapCtrl.imgBase, "_", numIts );
            Display( estMap, title );
            if( snapCtrl.itCounts )
                Display( itCountMap, title+"_counts" );
            auto colorMap = GetColorMap();
            SetColorMap( GRAYSCALE_DISCRETE );
            Display( estMap, title+"_discrete" );
            SetColorMap( colorMap );
            snapCtrl.imgDispCount = 0;
        }
    }
}

template<typename Real>
inline void
FinalSnapshot
( const DistMatrix<Real,VR,STAR>& estimates, 
  const DistMatrix<Int, VR,STAR>& itCounts,
        SnapshotCtrl& snapCtrl )
{
    DEBUG_ONLY(CSE cse("pspec::FinalSnapshot"));
    auto logMap = []( Real alpha ) { return Log(alpha); };
    if( snapCtrl.realSize != 0 && snapCtrl.imagSize != 0 )
    {
        const bool numSave = ( snapCtrl.numSaveFreq >= 0 );
        const bool imgSave = ( snapCtrl.imgSaveFreq >= 0 );
        const bool imgDisp = ( snapCtrl.imgDispFreq >= 0 );
        DistMatrix<Real> estMap(estimates.Grid()); 
        DistMatrix<Int> itCountMap(itCounts.Grid());
        if( numSave || imgSave || imgDisp )
        {
            ReshapeIntoGrid
            ( snapCtrl.realSize, snapCtrl.imagSize, estimates, estMap );
            if( snapCtrl.itCounts )
                ReshapeIntoGrid
                ( snapCtrl.realSize, snapCtrl.imagSize, itCounts, itCountMap );
        }
        if( numSave )
        {
            string base = snapCtrl.numBase;
            Write( estMap, base, snapCtrl.numFormat );
            if( snapCtrl.itCounts )
                Write( itCountMap, base+"_counts", snapCtrl.numFormat );
        }
        if( imgSave || imgDisp )
            EntrywiseMap( estMap, function<Real(Real)>(logMap) );
        if( imgSave )
        {
            string base = snapCtrl.imgBase;
            Write( estMap, base, snapCtrl.imgFormat );
            if( snapCtrl.itCounts )
                Write( itCountMap, base+"_counts", snapCtrl.imgFormat );
            auto colorMap = GetColorMap();
            SetColorMap( GRAYSCALE_DISCRETE );
            Write( estMap, base+"_discrete", snapCtrl.imgFormat );
            SetColorMap( colorMap );
        }
        if( imgDisp )
        {
            string base = snapCtrl.imgBase;
            Display( estMap, base );           
            if( snapCtrl.itCounts )
                Display( itCountMap, base+"_counts" );
            auto colorMap = GetColorMap();
            SetColorMap( GRAYSCALE_DISCRETE );
            Display( estMap, base+"_discrete" );
            SetColorMap( colorMap );
        }
    }
}

} // namespace pspec
} // namespace El

#endif // ifndef EL_PSEUDOSPECTRA_UTIL_SNAPSHOT_HPP
