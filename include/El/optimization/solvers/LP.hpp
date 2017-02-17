/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_OPTIMIZATION_SOLVERS_LP_HPP
#define EL_OPTIMIZATION_SOLVERS_LP_HPP

#include <El/optimization/solvers/util.hpp>

namespace El {

namespace LPApproachNS {
enum LPApproach {
  LP_ADMM,
  LP_MEHROTRA
};
} // namespace LPApproachNS
using namespace LPApproachNS;

namespace lp {

namespace direct {

// Attempt to solve a pair of Linear Programs in "direct" conic form:
//
//   min c^T x,
//   s.t. A x = b, x >= 0
//
//   max -b^T y
//   s.t. A^T y -z + c = 0, z >= 0
//

// Control structure for the high-level "direct" conic-form LP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    LPApproach approach=LP_MEHROTRA;
    ADMMCtrl<Real> admmCtrl;
    MehrotraCtrl<Real> mehrotraCtrl;

    Ctrl( bool isSparse )
    { mehrotraCtrl.system = ( isSparse ? AUGMENTED_KKT : NORMAL_KKT ); }
};

} // namespace direct

namespace affine {

// Attempt to solve a pair of Linear Programs in "affine" conic form:
//
//   min c^T x,
//   s.t. A x = b, G x + s = h, s >= 0
//
//   max -b^T y - h^T z
//   s.t. A^T y + G^T z + c = 0, z >= 0
//

// Control structure for the high-level "affine" conic-form LP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    LPApproach approach=LP_MEHROTRA;
    MehrotraCtrl<Real> mehrotraCtrl;
};

} // namespace affine

} // namespace lp

// Direct conic form
// -----------------
template<typename MatrixType,typename VectorType>
struct DirectLPProblem
{
    // The objective is 'c^T x'.
    VectorType c;

    // The primal equality constraint is 'A x = b; x >= 0'.
    MatrixType A;
    VectorType b;
};

// Quack...
template<typename Real>
void ForceSimpleAlignments
( DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const Grid& grid )
{ }
template<typename Real>
void ForceSimpleAlignments
( DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
  const Grid& grid )
{
    problem.c.SetGrid( grid );
    problem.A.SetGrid( grid );
    problem.b.SetGrid( grid );
    problem.c.Align(0,0);
    problem.A.Align(0,0);
    problem.b.Align(0,0);
}
template<typename Real>
void ForceSimpleAlignments
( DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const Grid& grid )
{ }
template<typename Real>
void ForceSimpleAlignments
( DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
  const Grid& grid )
{
    problem.c.SetGrid( grid );
    problem.A.SetGrid( grid );
    problem.b.SetGrid( grid );
}

template<typename Real>
bool SimpleAlignments
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem )
{ return true; }
template<typename Real>
bool SimpleAlignments
( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem )
{
    return problem.c.ColAlign() == 0 && problem.c.RowAlign() == 0 &&
           problem.A.ColAlign() == 0 && problem.A.RowAlign() == 0 &&
           problem.b.ColAlign() == 0 && problem.b.RowAlign() == 0;
}
template<typename Real>
bool SimpleAlignments
( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem )
{ return true; }
template<typename Real>
bool SimpleAlignments
( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem )
{ return true; }

template<typename VectorType>
struct DirectLPSolution
{
    // The primal solution is 'x', where 'A x = b, x >= 0'.
    VectorType x;

    // The dual solution is '(y,z)', where 'A^T y - z + c = 0, z >= 0'.
    VectorType y;
    VectorType z;
};

// Quack...
template<typename Real>
void ForceSimpleAlignments
( DirectLPSolution<Matrix<Real>>& solution, const Grid& grid )
{ }
template<typename Real>
void ForceSimpleAlignments
( DirectLPSolution<DistMatrix<Real>>& solution, const Grid& grid )
{
    solution.x.SetGrid( grid );
    solution.y.SetGrid( grid );
    solution.z.SetGrid( grid );
    solution.x.Align(0,0);
    solution.y.Align(0,0);
    solution.z.Align(0,0);
}
template<typename Real>
void ForceSimpleAlignments
( DirectLPSolution<DistMultiVec<Real>>& solution, const Grid& grid )
{
    solution.x.SetGrid( grid );
    solution.y.SetGrid( grid );
    solution.z.SetGrid( grid );
}

// Quack...
template<typename Real>
bool SimpleAlignments
( const DirectLPSolution<Matrix<Real>>& solution )
{ return true; }
template<typename Real>
bool SimpleAlignments
( const DirectLPSolution<DistMatrix<Real>>& solution )
{
    return solution.x.ColAlign() == 0 && solution.x.RowAlign() == 0 &&
           solution.y.ColAlign() == 0 && solution.y.RowAlign() == 0 &&
           solution.z.ColAlign() == 0 && solution.z.RowAlign() == 0;
}
template<typename Real>
bool SimpleAlignments
( const DirectLPSolution<DistSparseMatrix<Real>>& solution )
{ return true; }

template<typename VectorType>
struct DirectLPResidual
{
    VectorType primalEquality; // This residual is 'A x - b'.

    VectorType dualEquality; // This residual is 'A^T y - z + c'.
    VectorType dualConic; // This residual is 'x o z'.
};

// Quack...
template<typename Real>
void ForceSimpleAlignments
( DirectLPResidual<Matrix<Real>>& residual, const Grid& grid )
{ }
template<typename Real>
void ForceSimpleAlignments
( DirectLPResidual<DistMatrix<Real>>& residual, const Grid& grid )
{
    residual.primalEquality.SetGrid( grid );
    residual.dualEquality.SetGrid( grid );
    residual.dualConic.SetGrid( grid );
    residual.primalEquality.Align(0,0);
    residual.dualEquality.Align(0,0);
    residual.dualConic.Align(0,0);
}
template<typename Real>
void ForceSimpleAlignments
( DirectLPResidual<DistMultiVec<Real>>& residual, const Grid& grid )
{
    residual.primalEquality.SetGrid( grid );
    residual.dualEquality.SetGrid( grid );
    residual.dualConic.SetGrid( grid );
}

template<typename Real>
bool SimpleAlignments
( const DirectLPResidual<DistMatrix<Real>>& residual )
{
    return residual.primalEquality.ColAlign() == 0 &&
           residual.primalEquality.RowAlign() == 0 &&
           residual.dualEquality.ColAlign() == 0 &&
           residual.dualEquality.RowAlign() == 0 &&
           residual.dualConic.ColAlign() == 0 &&
           residual.dualConic.RowAlign() == 0;
}

template<typename Real>
void LP
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        DirectLPSolution<Matrix<Real>>& solution,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(false) );
template<typename Real>
void LP
( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        DirectLPSolution<DistMatrix<Real>>& solution,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(false) );
template<typename Real>
void LP
( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        DirectLPSolution<Matrix<Real>>& solution,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(true) );
template<typename Real>
void LP
( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        DirectLPSolution<DistMultiVec<Real>>& solution,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(true) );

// These interfaces are now deprecated in favor of the above.
template<typename Real>
[[deprecated]]
void LP
( const Matrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(false) );
template<typename Real>
[[deprecated]]
void LP
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,
        AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(false) );
template<typename Real>
[[deprecated]]
void LP
( const SparseMatrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(true) );
template<typename Real>
[[deprecated]]
void LP
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(true) );

// Affine conic form
// -----------------
template<typename MatrixType,typename VectorType>
struct AffineLPProblem
{
    // The objective is 'c^T x'.
    VectorType c;

    // The primal equality constraint is 'A x = b'.
    MatrixType A;
    VectorType b;

    // The primal cone constraint is 'G x + s = h, s >= 0'.
    MatrixType G;
    VectorType h;
};

// Quack...
template<typename Real>
void ForceSimpleAlignments
( AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const Grid& grid )
{ }
template<typename Real>
void ForceSimpleAlignments
( AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
  const Grid& grid )
{
    problem.c.SetGrid( grid );
    problem.A.SetGrid( grid );
    problem.b.SetGrid( grid );
    problem.G.SetGrid( grid );
    problem.h.SetGrid( grid );
    problem.c.Align(0,0);
    problem.A.Align(0,0);
    problem.b.Align(0,0);
    problem.G.Align(0,0);
    problem.h.Align(0,0);
}
template<typename Real>
void ForceSimpleAlignments
( AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const Grid& grid )
{ }
template<typename Real>
void ForceSimpleAlignments
( AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
  const Grid& grid )
{
    problem.c.SetGrid( grid );
    problem.A.SetGrid( grid );
    problem.b.SetGrid( grid );
    problem.G.SetGrid( grid );
    problem.h.SetGrid( grid );
}

// Quack...
template<typename Real>
bool SimpleAlignments
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem )
{ return true; }
template<typename Real>
bool SimpleAlignments
( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem )
{
    return problem.c.ColAlign() == 0 && problem.c.RowAlign() == 0 &&
           problem.A.ColAlign() == 0 && problem.A.RowAlign() == 0 &&
           problem.b.ColAlign() == 0 && problem.b.RowAlign() == 0 &&
           problem.G.ColAlign() == 0 && problem.G.RowAlign() == 0 &&
           problem.h.ColAlign() == 0 && problem.h.RowAlign() == 0;
}
template<typename Real>
bool SimpleAlignments
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem )
{ return true; }
template<typename Real>
bool SimpleAlignments
( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem )
{ return true; }

template<typename VectorType>
struct AffineLPSolution
{
    // The primal solution is '(x,s)', where 'A x = b, G x + s = h, s >= 0'.
    VectorType x;
    VectorType s;

    // The dual solution is '(y,z)', where 'A^T y + G^T z + c = 0, z >= 0'.
    VectorType y;
    VectorType z;
};

// Quack...
template<typename Real>
void ForceSimpleAlignments
( AffineLPSolution<Matrix<Real>>& solution, const Grid& grid )
{ }
template<typename Real>
void ForceSimpleAlignments
( AffineLPSolution<DistMatrix<Real>>& solution, const Grid& grid )
{
    solution.x.SetGrid( grid );
    solution.s.SetGrid( grid );
    solution.y.SetGrid( grid );
    solution.z.SetGrid( grid );
    solution.x.Align(0,0);
    solution.s.Align(0,0);
    solution.y.Align(0,0);
    solution.z.Align(0,0);
}
template<typename Real>
void ForceSimpleAlignments
( AffineLPSolution<DistMultiVec<Real>>& solution, const Grid& grid )
{
    solution.x.SetGrid( grid );
    solution.s.SetGrid( grid );
    solution.y.SetGrid( grid );
    solution.z.SetGrid( grid );
}

// Quack...
template<typename Real>
bool SimpleAlignments
( const AffineLPSolution<Matrix<Real>>& solution )
{ return true; }
template<typename Real>
bool SimpleAlignments
( const AffineLPSolution<DistMatrix<Real>>& solution )
{
    return solution.x.ColAlign() == 0 && solution.x.RowAlign() == 0 &&
           solution.y.ColAlign() == 0 && solution.y.RowAlign() == 0 &&
           solution.z.ColAlign() == 0 && solution.z.RowAlign() == 0 &&
           solution.s.ColAlign() == 0 && solution.s.RowAlign() == 0;
}
template<typename Real>
bool SimpleAlignments
( const AffineLPSolution<DistMultiVec<Real>>& solution )
{ return true; }

template<typename VectorType>
struct AffineLPResidual
{
    VectorType primalEquality; // This residual is 'A x - b'.
    VectorType primalConic; // This residual is 'G x + s - h'.

    VectorType dualEquality; // This residual is 'A^T y + G^T z + c'.
    VectorType dualConic; // This residual is 's o z'.
};

// Quack...
template<typename Real>
void ForceSimpleAlignments
( AffineLPResidual<Matrix<Real>>& residual, const Grid& grid )
{ }
template<typename Real>
void ForceSimpleAlignments
( AffineLPResidual<DistMatrix<Real>>& residual, const Grid& grid )
{
    residual.primalEquality.SetGrid( grid );
    residual.primalConic.SetGrid( grid );
    residual.dualEquality.SetGrid( grid );
    residual.dualConic.SetGrid( grid );
    residual.primalEquality.Align(0,0);
    residual.primalConic.Align(0,0);
    residual.dualEquality.Align(0,0);
    residual.dualConic.Align(0,0);
}
template<typename Real>
void ForceSimpleAlignments
( AffineLPResidual<DistMultiVec<Real>>& residual, const Grid& grid )
{
    residual.primalEquality.SetGrid( grid );
    residual.primalConic.SetGrid( grid );
    residual.dualEquality.SetGrid( grid );
    residual.dualConic.SetGrid( grid );
}

// Quack...
template<typename Real>
bool SimpleAlignments
( const AffineLPResidual<Matrix<Real>>& residual )
{ return true; }
template<typename Real>
bool SimpleAlignments
( const AffineLPResidual<DistMatrix<Real>>& residual )
{
    return residual.primalEquality.ColAlign() == 0 &&
           residual.primalEquality.RowAlign() == 0 &&
           residual.primalConic.ColAlign() == 0 &&
           residual.primalConic.RowAlign() == 0 &&
           residual.dualEquality.ColAlign() == 0 &&
           residual.dualEquality.RowAlign() == 0 &&
           residual.dualConic.ColAlign() == 0 &&
           residual.dualConic.RowAlign() == 0;
}
template<typename Real>
bool SimpleAlignments
( const AffineLPResidual<DistMultiVec<Real>>& residual )
{ return true; }

template<typename Real>
void LP
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
void LP
( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        AffineLPSolution<DistMatrix<Real>>& solution,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
void LP
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
void LP
( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        AffineLPSolution<DistMultiVec<Real>>& solution,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );

// These interfaces are now deprecated in favor of the above.
template<typename Real>
[[deprecated]]
void LP
( const Matrix<Real>& A,
  const Matrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
[[deprecated]]
void LP
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
        AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
[[deprecated]]
void LP
( const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
[[deprecated]]
void LP
( const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
        DistMultiVec<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );

// Mathematical Programming System
// -------------------------------

// Load an affine conic form LP stored in the Mathematical Programming System
// (MPS) format.
template<class MatrixType,class VectorType>
void ReadMPS
( AffineLPProblem<MatrixType,VectorType>& problem,
  const string& filename,
  bool compressed=false,
  bool minimize=true,
  bool keepNonnegativeWithZeroUpperBound=true,
  bool metadataSummary=false );

// Write out an LP in the MPS format.
template<class MatrixType,class VectorType>
void WriteMPS
( const DirectLPProblem<MatrixType,VectorType>& problem,
  const string& filename,
  bool compressed=false );
template<class MatrixType,class VectorType>
void WriteMPS
( const AffineLPProblem<MatrixType,VectorType>& problem,
  const string& filename,
  bool compressed=false );

void CompressMPS
( const string& filename, const string& compressedFilename );
void DecompressMPS
( const string& filename, const string& decompressedFilename );

} // namespace El

#endif // ifndef EL_OPTIMIZATION_SOLVERS_LP_HPP
