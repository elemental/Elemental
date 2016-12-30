/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

enum MPSSection {
  MPS_NAME,
  MPS_ROWS,
  MPS_COLUMNS,
  MPS_RHS,
  MPS_BOUNDS,
  MPS_RANGES,
  MPS_OBJSENSE,
  MPS_END,
  MPS_NONE 
};

enum MPSRowType {
  MPS_LESSER_ROW,
  MPS_GREATER_ROW,
  MPS_EQUALITY_ROW,
  MPS_NONCONSTRAINING_ROW
};

struct MPSRowData
{
  MPSRowType type;
  Int typeIndex;
};

struct MPSVariableData
{
  Int index;
  Int numNonzeros=0;
  bool lowerBounded=false;
  bool upperBounded=false;
  bool fixed=false;
  bool free=false;
  bool nonpositive=false;
  bool nonnegative=false;
};

struct MPSMeta
{
  // From the NAME section
  string name="";

  // From the ROWS section
  string costName="";
  Int numLesserRows=0;
  Int numGreaterRows=0;
  Int numEqualityRows=0;
  Int numNonconstrainingRows=0;
  std::map<string,MPSRowData> rowDict;

  // From the COLUMNS section
  std::map<string,MPSVariableData> variableDict;
  Int numColumnEntries=0;

  // From the BOUNDS section
  string boundName="";
  Int numUpperBounds=0;
  Int numLowerBounds=0;
  Int numFixedBounds=0;
  Int numFreeBounds=0;
  Int numNonpositiveBounds=0;
  Int numNonnegativeBounds=0;

  // From the RANGES section
  // TODO(poulson)

  // From the RHS section
  string rhsName="";
  Int numRHS=0;
};

namespace read_mps {

// NOTE: It is assumed that MPSMeta is in its original state on entry.
void GetMetadata
( const string& filename,
  MPSMeta& meta,
  bool upperBoundImplicitlyNonnegative )
{
    EL_DEBUG_CSE
    std::ifstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    // Rather than assuming that std::map<string,Int>::size() is constant-time,
    // we can maintain counters for the sizes of the variable dictionary.
    Int variableCounter=0;

    // Temporaries for the metadata extraction process.
    string token, rowType, rowName, variableName, boundMark;
    string rhsNameCandidate, boundNameCandidate;
    double value;

    // TODO(poulson): Convert each token to upper-case letters before each
    // comparison. While capital letters are used by convention, they are
    // not required.
    MPSSection section = MPS_NONE;
    string line;
    while( std::getline( file, line ) )
    {
        if( section == MPS_END )
        {
            Output("WARNING: Manually stopping after 'ENDATA'");
            break;
        }
        std::stringstream lineStream( line );

        // The first token on each line should be a string. We will check it for
        // equivalence with each section string.
        if( !(lineStream >> token) )
        {
            // This line only consists of whitespace.
            continue;
        }
        if( token == "NAME" )
        {
            if( meta.name != "" )
                Output("WARNING: Multiple 'NAME' sections");
            if( !(lineStream >> meta.name) )
                LogicError("Missing 'NAME' string");
            section = MPS_NAME;
            continue;
        }
        else if( token == "ROWS" )
        {
            if( meta.numLesserRows > 0 ||
                meta.numGreaterRows > 0 ||
                meta.numEqualityRows > 0 ||
                meta.numNonconstrainingRows > 0 )
                Output("WARNING: Multiple ROWS sections");
            section = MPS_ROWS;
            continue;
        }
        else if( token == "COLUMNS" )
        {
            if( meta.numColumnEntries > 0 )
                Output("WARNING: Multiple 'COLUMNS' sections");
            section = MPS_COLUMNS;
            continue;
        }
        else if( token == "RHS" )
        {
            section = MPS_RHS;
            continue;
        }
        else if( token == "BOUNDS" )
        {
            if( meta.numUpperBounds > 0 ||
                meta.numLowerBounds > 0 ||
                meta.numFixedBounds > 0 ||
                meta.numFreeBounds > 0 ||
                meta.numNonpositiveBounds > 0 ||
                meta.numNonnegativeBounds > 0 )
                Output("WARNING: Multiple 'BOUNDS' sections");
            section = MPS_BOUNDS;
            continue;
        }
        else if( token == "RANGES" )
        {
            section = MPS_RANGES;
            LogicError("MPS 'RANGES' section is not yet supported");
            continue;
        }
        else if( token == "ENDATA" )
        {
            section = MPS_END;
            continue;
        }
        else if( token == "MARKER" )
        {
            LogicError("MPS 'MARKER' section is not yet supported");
        }
        else if( token == "SOS" )
        {
            LogicError("MPS 'SOS' section is not yet supported");
        }

        // No section marker was found, so handle this data line.
        if( section == MPS_ROWS ) 
        {
            rowType = token;
            if( !(lineStream >> rowName) )
                LogicError("Invalid 'ROWS' section");
            MPSRowData rowData;
            if( rowType == "L" )
            {
                rowData.type = MPS_LESSER_ROW;
                rowData.typeIndex = meta.numLesserRows++;
                meta.rowDict[rowName] = rowData;
            }
            else if( rowType == "G" )
            {
                rowData.type = MPS_GREATER_ROW;
                rowData.typeIndex = meta.numGreaterRows++;
                meta.rowDict[rowName] = rowData;
            }
            else if( rowType == "E" )
            {
                rowData.type = MPS_EQUALITY_ROW;
                rowData.typeIndex = meta.numEqualityRows++;
                meta.rowDict[rowName] = rowData;
            }
            else if( rowType == "N" )
            {
                rowData.type = MPS_NONCONSTRAINING_ROW;
                rowData.typeIndex = meta.numNonconstrainingRows++;
                meta.rowDict[rowName] = rowData;
                if( meta.numNonconstrainingRows == 1 )
                    meta.costName = rowName;
            }
            else
                LogicError("Invalid 'ROWS' section");
        }
        else if( section == MPS_COLUMNS )
        {
            variableName = token;
            auto variableIter = meta.variableDict.find( variableName );
            if( variableIter == meta.variableDict.end() )
            {
                MPSVariableData variableData;
                variableData.index = variableCounter++;    
                meta.variableDict[variableName] = variableData;
                variableIter = meta.variableDict.find( variableName );
            }
            MPSVariableData& variableData = variableIter->second;

            // There should be either one or two pairs of entries left to read
            // from this line. We will now read in the first such pair.
            if( !(lineStream >> rowName) )
                LogicError("Invalid 'COLUMNS' section");
            if( !(lineStream >> value) )
                LogicError("Invalid 'COLUMNS' section");
            ++variableData.numNonzeros;
            ++meta.numColumnEntries;

            // We will attempt to read in a second pair.
            if( !(lineStream >> rowName) )
            {
                // There was not a second pair.
                continue;
            }
            if( !(lineStream >> value) )
                LogicError("Invalid 'COLUMNS' section");
            ++variableData.numNonzeros;
            ++meta.numColumnEntries;
        }
        else if( section == MPS_RHS )
        {
            rhsNameCandidate = token;
            if( meta.numRHS == 0 )
            {
                // We should currently have that rhsName == "".
                meta.rhsName = rhsNameCandidate;
                meta.numRHS = 1;
            }
            else if( rhsNameCandidate != meta.rhsName )
                LogicError
                ("Only single problem instances are currently supported "
                 "(multiple right-hand side names were encountered)");
            continue;
        }
        else if( section == MPS_BOUNDS )
        {
            // We already have the first token of a bounding row, which should
            // be of the same general form as
            //
            //   FX BOUNDROW VARIABLENAME 1734.
            //
            // in the case of 'VARIABLENAME' being fixed ('FX') at the value
            // 1734 (with this problem's bound name being 'BOUNDROW').
            boundMark = token;
            if( !(lineStream >> boundNameCandidate) )
                LogicError("Invalid 'BOUNDS' section");
            if( meta.boundName == "" )
                meta.boundName = boundNameCandidate;
            else if( meta.boundName != boundNameCandidate )
                LogicError
                ("Only single problem instances are currently supported "
                 "(multiple bound names were encountered)");
            if( !(lineStream >> variableName) )
                LogicError("Invalid 'BOUNDS' section");
            auto variableIter = meta.variableDict.find( variableName );
            if( variableIter == meta.variableDict.end() )
                LogicError
                ("Invalid 'BOUNDS' section (variable name not found)");
            MPSVariableData& variableData = variableIter->second;
            if( boundMark == "UP" )
                variableData.upperBounded = true;
            else if( boundMark == "LO" )
                variableData.lowerBounded = true;
            else if( boundMark == "FX" )
                variableData.fixed = true;
            else if( boundMark == "FR" )
                variableData.free = true;
            else if( boundMark == "MI" )
                variableData.nonpositive = true;
            else if( boundMark == "PL" )
                variableData.nonnegative = true;
            else
                LogicError("Invalid 'BOUNDS' section (unknown bound mark)");
        }
        else if( section == MPS_RANGES )
        {
            LogicError("The 'RANGES' section is not yet supported");
        }
        else
        {
            LogicError("Invalid MPS file");
        }
    }
    if( meta.name == "" )
        LogicError("No nontrivial 'NAME' was found");
    if( meta.numRHS == 0 )
    {
        // Any unmentioned values are assumed to be zero.
        meta.numRHS = 1;
    }

    // Now iterate through the variable map and make use of the requested
    // conventions for counting the number of bounds of each type.
    // Also warn if there are possibly conflicting bound types.
    for( auto& entry : meta.variableDict )
    {
        // Handle explicit upper and lower bounds.
        if( upperBoundImplicitlyNonnegative )
        {
            if( entry.second.upperBounded && entry.second.lowerBounded )
            {
                ++meta.numUpperBounds;
                ++meta.numLowerBounds;
            }
            else if( entry.second.upperBounded )
            {
                ++meta.numUpperBounds;
                entry.second.nonnegative = true;
            }
            else if( entry.second.lowerBounded )
            {
                ++meta.numLowerBounds;
            }
        }
        else
        {
            if( entry.second.upperBounded )
            {
                ++meta.numUpperBounds;
            }
            if( entry.second.lowerBounded )
            {
                ++meta.numLowerBounds;
            }
        }
        
        // Handle fixed values.
        if( entry.second.fixed )
        {
            if( entry.second.upperBounded || entry.second.lowerBounded )
                LogicError("Invalid bound combination");
            ++meta.numFixedBounds;
        }

        // Handle free values.
        if( entry.second.free )
        {
            if( entry.second.upperBounded || entry.second.lowerBounded )
                LogicError("Invalid bound combination");
            ++meta.numFreeBounds;
        }

        // Handle non-positive values.
        if( entry.second.nonpositive )
        {
            if( entry.second.upperBounded )
                Output
                ("WARNING: Combined nonpositive constraint with upper bound");
            ++meta.numNonpositiveBounds;
        }

        // Handle non-negative values. 
        if( entry.second.nonnegative )
        {
            if( entry.second.lowerBounded )
                Output
                ("WARNING: Combined nonnegative constraint with lower bound");
            ++meta.numNonnegativeBounds;
        }
    }
}

template<typename Real>
void GetMetadata
( AffileLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const string& filename,
  const MPSMeta& meta )
{
    EL_DEBUG_CSE
    const Int n = meta.variableDict.size();

    // We form the primal problem
    //
    //   arginf_{x,s} { c^T x | A x = b, G x + s = h, s >= 0 },
    //
    // where 'A x = b' consists of the 'equality' MPS rows and the 'fixed'
    // rows, e.g.,
    //
    //   | A0 | x = | b0 |,
    //   | A1 |     | b1 |
    //
    // where 'A0 x = b0' the set of 'equality' rows and 'A1 x = b1' is the 
    // set of 'fixed' rows (so that each row of 'A1' is all zeros except for
    // a single one).
    //
    const Int equalityOffset = 0;
    const Int fixedOffset = meta.numEqualityRows;
    const Int m = fixedOffset + meta.numFixedBounds;

    // Furthermore, 'G x + s = h, s >= 0' takes the form
    //
    //   | G0 | x <= | h0 |,
    //   | G1 |      | h1 |
    //   | G2 |      | h2 |
    //   | G3 |      | h3 |
    //
    // where 'G0 x <= h0' consists of the 'lesser' rows, 'G1 x <= h1' is the 
    // negation of the 'greater' rows, 'G2 x <= h2' consists of the upper
    // bounds (each row of 'G0' is all zeros except for a single one),
    // and 'G3 x <= h3' is the negation of the lower bounds (each row of 'G3'
    // is all zeros except for a single negative one).
    //
    const Int lesserOffset = 0;
    const Int greaterOffset = meta.numLesserRows;
    const Int upperBoundOffset = greaterOffset + meta.numGreaterRows;
    const Int lowerBoundOffset = upperBoundOffset + meta.numUpperBounds;
    const Int k = lowerBoundOffset + meta.numLowerBounds;

    // TODO(poulson): Support for reducing the system by eliminating the
    // 'fixed' variables. 

    Zeros( problem.c, n, 1 );
    Zeros( problem.A, m, n );
    Zeros( problem.b, m, 1 );
    Zeros( problem.G, k, n );
    Zeros( problem.h, k, 1 );

    std::ifstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    // Temporaries for the metadata extraction process.
    string token, rowType, rowName, variableName, boundMark;
    double value;

    // TODO(poulson): Convert each token to upper-case letters before each
    // comparison. While capital letters are used by convention, they are
    // not required.
    MPSSection section = MPS_NONE;
    string line;
    while( std::getline( file, line ) )
    {
        if( section == MPS_END )
        {
            Output("WARNING: Manually stopping after 'ENDATA'");
            break;
        }
        std::stringstream lineStream( line );

        // The first token on each line should be a string. We will check it for
        // equivalence with each section string.
        if( !(lineStream >> token) )
        {
            // This line only consists of whitespace.
            continue;
        }
        if( token == "NAME" )
        {
            if( !(lineStream >> meta.name) )
                LogicError("Missing 'NAME' string");
            section = MPS_NAME;
            continue;
        }
        else if( token == "ROWS" )
        {
            section = MPS_ROWS;
            continue;
        }
        else if( token == "COLUMNS" )
        {
            section = MPS_COLUMNS;
            continue;
        }
        else if( token == "RHS" )
        {
            section = MPS_RHS;
            continue;
        }
        else if( token == "BOUNDS" )
        {
            section = MPS_BOUNDS;
            continue;
        }
        else if( token == "RANGES" )
        {
            section = MPS_RANGES;
            LogicError("MPS 'RANGES' section is not yet supported");
            continue;
        }
        else if( token == "ENDATA" )
        {
            section = MPS_END;
            continue;
        }
        else if( token == "MARKER" )
        {
            LogicError("MPS 'MARKER' section is not yet supported");
        }
        else if( token == "SOS" )
        {
            LogicError("MPS 'SOS' section is not yet supported");
        }

        // No section marker was found, so handle this data line.
        if( section == MPS_ROWS ) 
        {
            continue;
        }
        else if( section == MPS_COLUMNS )
        {
            // HERE
            variableName = token;
            auto variableIter = meta.variableDict.find( variableName );
            if( variableIter == meta.variableDict.end() )
            {
                MPSVariableData variableData;
                variableData.index = variableCounter++;    
                meta.variableDict[variableName] = variableData;
                variableIter = meta.variableDict.find( variableName );
            }
            MPSVariableData& variableData = variableIter->second;

            // There should be either one or two pairs of entries left to read
            // from this line. We will now read in the first such pair.
            if( !(lineStream >> rowName) )
                LogicError("Invalid 'COLUMNS' section");
            if( !(lineStream >> value) )
                LogicError("Invalid 'COLUMNS' section");
            ++variableData.numNonzeros;
            ++meta.numColumnEntries;

            // We will attempt to read in a second pair.
            if( !(lineStream >> rowName) )
            {
                // There was not a second pair.
                continue;
            }
            if( !(lineStream >> value) )
                LogicError("Invalid 'COLUMNS' section");
            ++variableData.numNonzeros;
            ++meta.numColumnEntries;
        }
        else if( section == MPS_RHS )
        {
            rhsNameCandidate = token;
            if( meta.numRHS == 0 )
            {
                // We should currently have that rhsName == "".
                meta.rhsName = rhsNameCandidate;
                meta.numRHS = 1;
            }
            else if( rhsNameCandidate != meta.rhsName )
                LogicError
                ("Only single problem instances are currently supported "
                 "(multiple right-hand side names were encountered)");
            continue;
        }
        else if( section == MPS_BOUNDS )
        {
            // We already have the first token of a bounding row, which should
            // be of the same general form as
            //
            //   FX BOUNDROW VARIABLENAME 1734.
            //
            // in the case of 'VARIABLENAME' being fixed ('FX') at the value
            // 1734 (with this problem's bound name being 'BOUNDROW').
            boundMark = token;
            if( !(lineStream >> boundNameCandidate) )
                LogicError("Invalid 'BOUNDS' section");
            if( meta.boundName == "" )
                meta.boundName = boundNameCandidate;
            else if( meta.boundName != boundNameCandidate )
                LogicError
                ("Only single problem instances are currently supported "
                 "(multiple bound names were encountered)");
            if( !(lineStream >> variableName) )
                LogicError("Invalid 'BOUNDS' section");
            auto variableIter = meta.variableDict.find( variableName );
            if( variableIter == meta.variableDict.end() )
                LogicError
                ("Invalid 'BOUNDS' section (variable name not found)");
            MPSVariableData& variableData = variableIter->second;
            if( boundMark == "UP" )
                variableData.upperBounded = true;
            else if( boundMark == "LO" )
                variableData.lowerBounded = true;
            else if( boundMark == "FX" )
                variableData.fixed = true;
            else if( boundMark == "FR" )
                variableData.free = true;
            else if( boundMark == "MI" )
                variableData.nonpositive = true;
            else if( boundMark == "PL" )
                variableData.nonnegative = true;
            else
                LogicError("Invalid 'BOUNDS' section (unknown bound mark)");
        }
        else if( section == MPS_RANGES )
        {
            LogicError("The 'RANGES' section is not yet supported");
        }
        else
        {
            LogicError("Invalid MPS file");
        }
    }
    if( meta.name == "" )
        LogicError("No nontrivial 'NAME' was found");
    if( meta.numRHS == 0 )
    {
        // Any unmentioned values are assumed to be zero.
        meta.numRHS = 1;
    }

    // Now iterate through the variable map and make use of the requested
    // conventions for counting the number of bounds of each type.
    // Also warn if there are possibly conflicting bound types.
    for( auto& entry : meta.variableDict )
    {
        // Handle explicit upper and lower bounds.
        if( upperBoundImplicitlyNonnegative )
        {
            if( entry.second.upperBounded && entry.second.lowerBounded )
            {
                ++meta.numUpperBounds;
                ++meta.numLowerBounds;
            }
            else if( entry.second.upperBounded )
            {
                ++meta.numUpperBounds;
                entry.second.nonnegative = true;
            }
            else if( entry.second.lowerBounded )
            {
                ++meta.numLowerBounds;
            }
        }
        else
        {
            if( entry.second.upperBounded )
            {
                ++meta.numUpperBounds;
            }
            if( entry.second.lowerBounded )
            {
                ++meta.numLowerBounds;
            }
        }
        
        // Handle fixed values.
        if( entry.second.fixed )
        {
            if( entry.second.upperBounded || entry.second.lowerBounded )
                LogicError("Invalid bound combination");
            ++meta.numFixedBounds;
        }

        // Handle free values.
        if( entry.second.free )
        {
            if( entry.second.upperBounded || entry.second.lowerBounded )
                LogicError("Invalid bound combination");
            ++meta.numFreeBounds;
        }

        // Handle non-positive values.
        if( entry.second.nonpositive )
        {
            if( entry.second.upperBounded )
                Output
                ("WARNING: Combined nonpositive constraint with upper bound");
            ++meta.numNonpositiveBounds;
        }

        // Handle non-negative values. 
        if( entry.second.nonnegative )
        {
            if( entry.second.lowerBounded )
                Output
                ("WARNING: Combined nonnegative constraint with lower bound");
            ++meta.numNonnegativeBounds;
        }
    }
}

template<typename Real>
void Helper
( AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const string& filename,
  bool compressed )
{
    EL_DEBUG_CSE
    if( compressed )
        LogicError("Compressed MPS is not yet supported");
    bool upperBoundImplicitlyNonnegative = false;

    // Perform an initial pass through the file to determine the metadata.
    MPSMeta meta;
    GetMetadata( filename, meta, upperBoundImplicitlyNonnegative );

    // We will formulate the LP as...

    // TODO(poulson): Test what we have so far on adlittle.

    Output("meta.name=",meta.name);
    Output("meta.costName=",meta.costName);
    Output("meta.numLesserRows=",meta.numLesserRows);
    Output("meta.numGreaterRows=",meta.numGreaterRows);
    Output("meta.numEqualityRows=",meta.numEqualityRows);
    Output("meta.numNonconstrainingRows=",meta.numNonconstrainingRows);
    Output("meta.numColumnEntries=",meta.numColumnEntries);
    Output("meta.boundName=",meta.boundName);
    Output("meta.numUpperBounds=",meta.numUpperBounds);
    Output("meta.numLowerBounds=",meta.numLowerBounds);
    Output("meta.numFixedBounds=",meta.numFixedBounds);
    Output("meta.numFreeBounds=",meta.numFreeBounds);
    Output("meta.numNonpositiveBounds=",meta.numNonpositiveBounds);
    Output("meta.numNonnegativeBounds=",meta.numNonnegativeBounds);
    Output("meta.rhsName=",meta.rhsName);
    Output("meta.numRHS=",meta.numRHS);

    FormProblem( problem, filename, meta );
}

template<typename Real>
void Helper
( AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
  const string& filename,
  bool compressed )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

template<typename Real>
void Helper
( AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const string& filename,
  bool compressed )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

template<typename Real>
void Helper
( AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
  const string& filename,
  bool compressed )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

} // namespace read_mps

template<class MatrixType,class VectorType>
void ReadMPS
( AffineLPProblem<MatrixType,VectorType>& problem,
  const string& filename,
  bool compressed )
{
    EL_DEBUG_CSE
    read_mps::Helper( problem, filename, compressed );
}

namespace write_mps {

template<typename Real>
void Helper
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const string& filename,
  bool compressed )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

template<typename Real>
void Helper
( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
  const string& filename,
  bool compressed )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

template<typename Real>
void Helper
( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const string& filename,
  bool compressed )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

template<typename Real>
void Helper
( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
  const string& filename,
  bool compressed )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

template<typename Real>
void Helper
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const string& filename,
  bool compressed )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

template<typename Real>
void Helper
( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
  const string& filename,
  bool compressed )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

template<typename Real>
void Helper
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const string& filename,
  bool compressed )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

template<typename Real>
void Helper
( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
  const string& filename,
  bool compressed )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

} // namespace write_mps

template<class MatrixType,class VectorType>
void WriteMPS
( const DirectLPProblem<MatrixType,VectorType>& problem,
  const string& filename,
  bool compressed )
{
    EL_DEBUG_CSE
    write_mps::Helper( problem, filename, compressed );
}

template<class MatrixType,class VectorType>
void WriteMPS
( const AffineLPProblem<MatrixType,VectorType>& problem,
  const string& filename,
  bool compressed )
{
    EL_DEBUG_CSE
    write_mps::Helper( problem, filename, compressed );
}

void CompressMPS
( const string& filename, const string& compressedFilename )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

void DecompressMPS
( const string& filename, const string& decompressedFilename )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

} // namespace El
