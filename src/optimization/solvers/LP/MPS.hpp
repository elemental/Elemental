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
  Int lowerBoundIndex=-1;

  bool upperBounded=false;
  Int upperBoundIndex=-1;

  bool fixed=false;
  Int fixedIndex=-1;

  bool free=false;
  Int freeIndex=-1;

  bool nonpositive=false;
  Int nonpositiveIndex=-1;

  bool nonnegative=false;
  Int nonnegativeIndex=-1;
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
  Int numEqualityEntries=0;
  Int numInequalityEntries=0;

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

  Int m=0, n=0, k=0;

  //
  //   | A0 | x = | b0 |
  //   | A1 |     | b1 |
  //
  Int equalityOffset=-1, fixedOffset=-1;

  //
  //   | G0 | x <= | h0 |
  //   | G1 |      | h1 |
  //   | G2 |      | h2 |
  //   | G3 |      | h3 |
  //   | G4 |      | h4 |
  //   | G5 |      | h5 |
  //
  Int lesserOffset=-1,
      greaterOffset=-1,
      upperBoundOffset=-1,
      lowerBoundOffset=-1,
      nonpositiveOffset=-1,
      nonnegativeOffset=-1;

  void PrintSummary() const
  {
      Output("MPSMeta summary:");
      Output("  name=",name);
      Output("  costName=",costName);
      Output("  numLesserRows=",numLesserRows);
      Output("  numGreaterRows=",numGreaterRows);
      Output("  numEqualityRows=",numEqualityRows);
      Output("  numNonconstrainingRows=",numNonconstrainingRows);
      Output("  numEqualityEntries=",numEqualityEntries);
      Output("  numInequalityEntries=",numInequalityEntries);
      Output("  boundName=",boundName);
      Output("  numUpperBounds=",numUpperBounds);
      Output("  numLowerBounds=",numLowerBounds);
      Output("  numFixedBounds=",numFixedBounds);
      Output("  numFreeBounds=",numFreeBounds);
      Output("  numNonpositiveBounds=",numNonpositiveBounds);
      Output("  numNonnegativeBounds=",numNonnegativeBounds);
      Output("  rhsName=",rhsName);
  }
};

enum AffineLPMatrixType {
  AFFINE_LP_COST_VECTOR, // The 'c' in 'c^T x'
  AFFINE_LP_EQUALITY_MATRIX, // The 'A' in 'A x = b'
  AFFINE_LP_EQUALITY_VECTOR, // The 'b' in 'A x = b'
  AFFINE_LP_INEQUALITY_MATRIX, // The 'G' in  'G x <= h'
  AFFINE_LP_INEQUALITY_VECTOR // The 'h' in 'G x <= h'
};

template<typename Real>
struct AffineLPEntry
{
  AffineLPMatrixType type;
  Int row;
  Int column;
  Real value;
};

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
// Furthermore, 'G x + s = h, s >= 0' takes the form
//
//   | G0 | x <= | h0 |,
//   | G1 |      | h1 |
//   | G2 |      | h2 |
//   | G3 |      | h3 |
//   | G4 |      | h4 |
//   | G5 |      | h5 |
//
// where 'G0 x <= h0' consists of the 'lesser' rows, 'G1 x <= h1' is the
// negation of the 'greater' rows, 'G2 x <= h2' consists of the upper
// bounds (each row of 'G0' is all zeros except for a single one),
// 'G3 x <= h3' is the negation of the lower bounds (each row of 'G3'
// is all zeros except for a single negative one), 'G4 x <= h4' is the
// set of nonpositive bounds, and 'G5 x <= h5' is the set of nonnegative
// bounds.
//
class MPSReader
{
public:
    MPSReader
    ( const string& filename,
      bool compressed=false,
      bool minimize=true,
      bool upperBoundImplicitlyNonnegative=false );

    // Attempt to enqueue another entry and return true if successful.
    bool QueuedEntry();

    // Only one call is allowed per call to 'QueuedEntry'.
    AffineLPEntry<double> GetEntry(); 

    const MPSMeta& Meta() const;

private:
    string filename_;
    std::ifstream file_;

    bool minimize_; 
    bool upperBoundImplicitlyNonnegative_;
    MPSMeta meta_;

    vector<AffineLPEntry<double>> queuedEntries_;
  
    MPSSection section_=MPS_NONE;

    // Temporaries for the metadata extraction process.
    string line_;
    string token_, rowType_, rowName_, variableName_, boundMark_;
    double value_;

    // For the final loop over the variable dictionary to extract the
    // nonpositive and nonnegative bounds.
    typename std::map<string,MPSVariableData>::const_iterator variableIter_;
};

MPSReader::MPSReader
( const string& filename,
  bool compressed,
  bool minimize,
  bool upperBoundImplicitlyNonnegative )
: filename_(filename),
  file_(filename.c_str()),
  minimize_(minimize),
  upperBoundImplicitlyNonnegative_(upperBoundImplicitlyNonnegative)
{
    EL_DEBUG_CSE
    if( compressed )
        LogicError("Compressed reads are not yet supported");
    if( !file_.is_open() )
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
    while( std::getline( file_, line ) )
    {
        std::stringstream lineStream( line );
        const char firstChar = lineStream.peek();
        const bool isDataLine = firstChar == ' ' || firstChar == '\t';

        // The first token on each line should be a string. We will check it for
        // equivalence with each section string.
        if( !(lineStream >> token) )
        {
            // This line only consists of whitespace.
            continue;
        }

        if( !isDataLine )
        {
            if( token == "NAME" )
            {
                if( meta_.name != "" )
                    Output("WARNING: Multiple 'NAME' sections");
                if( !(lineStream >> meta_.name) )
                    LogicError("Missing 'NAME' string");
                section = MPS_NAME;
            }
            else if( token == "OBJSENSE" )
            {
                LogicError("OBJSENSE is not yet handled");
            }
            else if( token == "ROWS" )
            {
                if( meta_.numLesserRows > 0 ||
                    meta_.numGreaterRows > 0 ||
                    meta_.numEqualityRows > 0 ||
                    meta_.numNonconstrainingRows > 0 )
                    Output("WARNING: Multiple ROWS sections");
                section = MPS_ROWS;
            }
            else if( token == "COLUMNS" )
            {
                if( meta_.numEqualityEntries > 0 ||
                    meta_.numInequalityEntries > 0 )
                    Output("WARNING: Multiple 'COLUMNS' sections");
                section = MPS_COLUMNS;
            }
            else if( token == "RHS" )
            {
                section = MPS_RHS;
            }
            else if( token == "BOUNDS" )
            {
                if( meta_.numUpperBounds > 0 ||
                    meta_.numLowerBounds > 0 ||
                    meta_.numFixedBounds > 0 ||
                    meta_.numFreeBounds > 0 ||
                    meta_.numNonpositiveBounds > 0 ||
                    meta_.numNonnegativeBounds > 0 )
                    Output("WARNING: Multiple 'BOUNDS' sections");
                section = MPS_BOUNDS;
            }
            else if( token == "RANGES" )
            {
                section = MPS_RANGES;
                LogicError("MPS 'RANGES' section is not yet supported");
            }
            else if( token == "ENDATA" )
            {
                section = MPS_END;
                break;
            }
            else if( token == "MARKER" )
            {
                LogicError("MPS 'MARKER' section is not yet supported");
            }
            else if( token == "SOS" )
            {
                LogicError("MPS 'SOS' section is not yet supported");
            }
            else
            {
                LogicError("Section token ",token," is not recognized");
            }
            continue;
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
                rowData.typeIndex = meta_.numLesserRows++;
                meta_.rowDict[rowName] = rowData;
            }
            else if( rowType == "G" )
            {
                rowData.type = MPS_GREATER_ROW;
                rowData.typeIndex = meta_.numGreaterRows++;
                meta_.rowDict[rowName] = rowData;
            }
            else if( rowType == "E" )
            {
                rowData.type = MPS_EQUALITY_ROW;
                rowData.typeIndex = meta_.numEqualityRows++;
                meta_.rowDict[rowName] = rowData;
            }
            else if( rowType == "N" )
            {
                rowData.type = MPS_NONCONSTRAINING_ROW;
                rowData.typeIndex = meta_.numNonconstrainingRows++;
                meta_.rowDict[rowName] = rowData;
                if( meta_.numNonconstrainingRows == 1 )
                    meta_.costName = rowName;
            }
            else
                LogicError("Invalid 'ROWS' section");
        }
        else if( section == MPS_COLUMNS )
        {
            variableName = token;
            auto variableIter = meta_.variableDict.find( variableName );
            if( variableIter == meta_.variableDict.end() )
            {
                MPSVariableData variableData;
                variableData.index = variableCounter++;
                meta_.variableDict[variableName] = variableData;
                variableIter = meta_.variableDict.find( variableName );
            }
            MPSVariableData& variableData = variableIter->second;

            // There should be either one or two pairs of entries left to read
            // from this line.
            for( Int pair=0; pair<2; ++pair )
            {
                if( !(lineStream >> rowName) )
                {
                    if( pair == 0 )
                        LogicError("Invalid 'COLUMNS' section");
                    else
                        break;
                }
                if( !(lineStream >> value) )
                    LogicError("Invalid 'COLUMNS' section");
                ++variableData.numNonzeros;


                auto rowIter = meta_.rowDict.find( rowName );
                if( rowIter == meta_.rowDict.end() )
                    LogicError("Could not find row ",rowName);
                const auto& rowData = rowIter->second;
                if( rowData.type == MPS_EQUALITY_ROW )
                {
                    // A(row,column) = value
                    ++meta_.numEqualityEntries;
                }
                else if( rowData.type == MPS_LESSER_ROW )
                {
                    // G(row,column) = value
                    ++meta_.numInequalityEntries;
                }
                else if( rowData.type == MPS_GREATER_ROW )
                {
                    // G(row,column) = -value
                    ++meta_.numInequalityEntries;
                }
                else if( rowData.type == MPS_NONCONSTRAINING_ROW )
                {
                    if( rowData.typeIndex == 0 )
                    {
                        // c(column) = value
                    }
                    // else?
                }
            }
        }
        else if( section == MPS_RHS )
        {
            rhsNameCandidate = token;
            if( meta_.numRHS == 0 )
            {
                // We should currently have that rhsName == "".
                meta_.rhsName = rhsNameCandidate;
                meta_.numRHS = 1;
            }
            else if( rhsNameCandidate != meta_.rhsName )
                LogicError
                ("Only single problem instances are currently supported "
                 "(multiple right-hand side names were encountered)");

            // There should be either one or two pairs of entries left to read
            // from this line.
            for( Int pair=0; pair<2; ++pair )
            {
                if( !(lineStream >> rowName) )
                {
                    if( pair == 0 )
                        LogicError("Invalid 'RHS' section");
                    else
                        break;
                }
                if( !(lineStream >> value) )
                    LogicError("Invalid 'RHS' section");
            }
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
            if( meta_.boundName == "" )
                meta_.boundName = boundNameCandidate;
            else if( meta_.boundName != boundNameCandidate )
                LogicError
                ("Only single problem instances are currently supported "
                 "(multiple bound names were encountered)");
            if( !(lineStream >> variableName) )
                LogicError("Invalid 'BOUNDS' section");
            auto variableIter = meta_.variableDict.find( variableName );
            if( variableIter == meta_.variableDict.end() )
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
    if( meta_.name == "" )
        LogicError("No nontrivial 'NAME' was found");
    if( meta_.numRHS == 0 )
    {
        // Any unmentioned values are assumed to be zero.
        meta_.numRHS = 1;
    }

    // Now iterate through the variable map and make use of the requested
    // conventions for counting the number of bounds of each type.
    // Also warn if there are possibly conflicting bound types.
    for( auto& entry : meta_.variableDict )
    {
        auto& data = entry.second;
        // Handle explicit upper and lower bounds.
        if( upperBoundImplicitlyNonnegative )
        {
            if( data.upperBounded && data.lowerBounded )
            {
                data.upperBoundIndex = meta_.numUpperBounds++;
                data.lowerBoundIndex = meta_.numLowerBounds++;
            }
            else if( data.upperBounded )
            {
                data.upperBoundIndex = meta_.numUpperBounds++;
                data.nonnegative = true;
            }
            else if( data.lowerBounded )
            {
                data.lowerBoundIndex = meta_.numLowerBounds++;
            }
        }
        else
        {
            if( data.upperBounded )
            {
                data.upperBoundIndex = meta_.numUpperBounds++;
            }
            if( data.lowerBounded )
            {
                data.lowerBoundIndex = meta_.numLowerBounds++;
            }
        }

        // Handle fixed values.
        if( data.fixed )
        {
            if( data.upperBounded || data.lowerBounded )
                LogicError("Invalid bound combination");
            data.fixedIndex = meta_.numFixedBounds++;
        }

        // Handle free values.
        if( data.free )
        {
            if( data.upperBounded || data.lowerBounded )
                LogicError("Invalid bound combination");
            data.freeIndex = meta_.numFreeBounds++;
        }

        // Handle non-positive values.
        if( data.nonpositive )
        {
            if( data.upperBounded )
                Output
                ("WARNING: Combined nonpositive constraint with upper bound");
            data.nonpositiveIndex = meta_.numNonpositiveBounds++;
        }

        // Default to nonnegative if there were not any markings.
        if( !data.upperBounded &&
            !data.lowerBounded &&
            !data.fixed &&
            !data.free &&
            !data.nonpositive )
        {
            data.nonnegative = true;
        }

        // Handle non-negative values.
        if( data.nonnegative )
        {
            if( data.lowerBounded )
                Output
                ("WARNING: Combined nonnegative constraint with lower bound");
            data.nonnegativeIndex = meta_.numNonnegativeBounds++;
        }
    }

    // Now that the initial pass over the file is done, we can set up for the
    // 'QueuedEntry'/'GetEntry' cycle.

    // Reset the file.
    file_.clear();
    file_.seekg( 0, std::ios::beg );

    // Extract the number of variables
    // (the matrix 'A' is 'm x n' and 'G' is 'k x n').
    meta_.n = meta_.variableDict.size();
    variableIter_ = meta_.variableDict.cbegin();

    //
    //   | A0 | x = | b0 |
    //   | A1 |     | b1 |
    //
    meta_.equalityOffset = 0;
    meta_.fixedOffset = meta_.numEqualityRows;
    meta_.m = meta_.fixedOffset + meta_.numFixedBounds;

    //
    //   | G0 | x <= | h0 |
    //   | G1 |      | h1 |
    //   | G2 |      | h2 |
    //   | G3 |      | h3 |
    //   | G4 |      | h4 |
    //   | G5 |      | h5 |
    //
    meta_.lesserOffset = 0;
    meta_.greaterOffset = meta_.numLesserRows;
    meta_.upperBoundOffset = meta_.greaterOffset + meta_.numGreaterRows;
    meta_.lowerBoundOffset = meta_.upperBoundOffset + meta_.numUpperBounds;
    meta_.nonpositiveOffset = meta_.lowerBoundOffset + meta_.numLowerBounds;
    meta_.nonnegativeOffset =
      meta_.nonpositiveOffset + meta_.numNonpositiveBounds;
    meta_.k = meta_.nonnegativeOffset + meta_.numNonnegativeBounds;
}

bool MPSReader::QueuedEntry()
{
    EL_DEBUG_CSE

    while( queuedEntries_.size() == 0 &&
           section_ != MPS_END &&
           std::getline( file_, line_ ) )
    {
        std::stringstream lineStream( line_ );
        const char firstChar = lineStream.peek();
        const bool isDataLine = firstChar == ' ' || firstChar == '\t';

        // The first token on each line should be a string. We will check it for
        // equivalence with each section string.
        if( !(lineStream >> token_) )
        {
            // This line only consists of whitespace.
            continue;
        }

        if( !isDataLine )
        {
            if( token_ == "NAME" )
            {
                section_ = MPS_NAME;
            }
            else if( token_ == "OBJSENSE" )
            {
                LogicError("OBJSENSE is not yet handled");
            }
            else if( token_ == "ROWS" )
            {
                section_ = MPS_ROWS;
            }
            else if( token_ == "COLUMNS" )
            {
                section_ = MPS_COLUMNS;
            }
            else if( token_ == "RHS" )
            {
                section_ = MPS_RHS;
            }
            else if( token_ == "BOUNDS" )
            {
                section_ = MPS_BOUNDS;
            }
            else if( token_ == "RANGES" )
            {
                section_ = MPS_RANGES;
            }
            else if( token_ == "ENDATA" )
            {
                section_ = MPS_END;
                break;
            }
            else
            {
                LogicError("Section token ",token_," is not recognized");
            }
            continue;
        }

        // No section marker was found, so handle this data line.
        double value;
        if( section_ == MPS_ROWS )
        {
            // We already have the row names.
        }
        else if( section_ == MPS_COLUMNS )
        {
            variableName_ = token_;
            auto variableIter = meta_.variableDict.find( variableName_ );
            if( variableIter == meta_.variableDict.end() )
                LogicError("Could not find variable ",variableName_);
            const MPSVariableData& variableData = variableIter->second;
            const Int column = variableData.index;

            // There should be either one or two pairs of entries left to read
            // from this line.
            for( Int pair=0; pair<2; ++pair )
            {
                if( !(lineStream >> rowName_) )
                {
                    if( pair == 0 )
                        LogicError("Invalid 'COLUMNS' section");
                    else
                        break;
                }
                if( !(lineStream >> value) )
                    LogicError("Invalid 'COLUMNS' section");
                auto rowIter = meta_.rowDict.find( rowName_ );
                if( rowIter == meta_.rowDict.end() )
                    LogicError("Could not find row ",rowName_);
                const auto& rowData = rowIter->second;
                AffineLPEntry<double> entry;
                if( rowData.type == MPS_EQUALITY_ROW )
                {
                    // A(row,column) = value
                    const Int row = meta_.equalityOffset + rowData.typeIndex;
                    entry.type = AFFINE_LP_EQUALITY_MATRIX;
                    entry.row = row;
                    entry.column = column;
                    entry.value = value;
                    queuedEntries_.push_back( entry );
                }
                else if( rowData.type == MPS_LESSER_ROW )
                {
                    // G(row,column) = value
                    const Int row = meta_.lesserOffset + rowData.typeIndex;
                    entry.type = AFFINE_LP_INEQUALITY_MATRIX;
                    entry.row = row;
                    entry.column = column;
                    entry.value = value;
                    queuedEntries_.push_back( entry );
                }
                else if( rowData.type == MPS_GREATER_ROW )
                {
                    // G(row,column) = -value
                    const Int row = meta_.greaterOffset + rowData.typeIndex;
                    entry.type = AFFINE_LP_INEQUALITY_MATRIX;
                    entry.row = row;
                    entry.column = column;
                    entry.value = -value;
                    queuedEntries_.push_back( entry );
                }
                else if( rowData.type == MPS_NONCONSTRAINING_ROW )
                {
                    if( rowData.typeIndex == 0 )
                    {
                        // c(column) = value
                        entry.type = AFFINE_LP_COST_VECTOR;
                        entry.row = column;
                        entry.column = 0;
                        entry.value = minimize_ ? value : -value;
                        queuedEntries_.push_back( entry );
                    }
                    // else?
                }
            }
        }
        else if( section_ == MPS_RHS )
        {
            // There should be either one or two pairs of entries left to read
            // from this line.
            for( Int pair=0; pair<2; ++pair )
            {
                if( !(lineStream >> rowName_) )
                {
                    if( pair == 0 )
                        LogicError("Invalid 'RHS' section");
                    else
                        break;
                }
                if( !(lineStream >> value) )
                    LogicError("Invalid 'RHS' section");
                auto rowIter = meta_.rowDict.find( rowName_ );
                if( rowIter == meta_.rowDict.end() )
                    LogicError("Could not find row ",rowName_);
                const auto& rowData = rowIter->second;
                AffineLPEntry<double> entry;
                if( rowData.type == MPS_EQUALITY_ROW )
                {
                    // b(row) = value
                    const Int row = meta_.equalityOffset + rowData.typeIndex;
                    entry.type = AFFINE_LP_EQUALITY_VECTOR;
                    entry.row = row;
                    entry.column = 0;
                    entry.value = value;
                    queuedEntries_.push_back( entry );
                }
                else if( rowData.type == MPS_LESSER_ROW )
                {
                    // h(row) = value
                    const Int row = meta_.lesserOffset + rowData.typeIndex;
                    entry.type = AFFINE_LP_INEQUALITY_VECTOR;
                    entry.row = row;
                    entry.column = 0;
                    entry.value = value;
                    queuedEntries_.push_back( entry );
                }
                else if( rowData.type == MPS_GREATER_ROW )
                {
                    // h(row) = -value
                    const Int row = meta_.greaterOffset + rowData.typeIndex;
                    entry.type = AFFINE_LP_INEQUALITY_VECTOR;
                    entry.row = row;
                    entry.column = 0;
                    entry.value = -value;
                    queuedEntries_.push_back( entry );
                }
                else if( rowData.type == MPS_NONCONSTRAINING_ROW )
                {
                    Output("WARNING: Nonsensical RHS for nonconstrained row");
                }
            }
        }
        else if( section_ == MPS_BOUNDS )
        {
            // We already have the first token of a bounding row, which should
            // be of the same general form as
            //
            //   FX BOUNDROW VARIABLENAME 1734.
            //
            // in the case of 'VARIABLENAME' being fixed ('FX') at the value
            // 1734 (with this problem's bound name being 'BOUNDROW').
            boundMark_ = token_;
            if( !(lineStream >> token_) )
                LogicError("Invalid 'BOUNDS' section");
            if( !(lineStream >> variableName_) )
                LogicError("Invalid 'BOUNDS' section");
            auto variableIter = meta_.variableDict.find( variableName_ );
            if( variableIter == meta_.variableDict.end() )
                LogicError
                ("Invalid 'BOUNDS' section (variable name not found)");
            const MPSVariableData& variableData = variableIter->second;
            const Int column = variableData.index;

            AffineLPEntry<double> entry;
            if( boundMark_ == "UP" )
            {
                if( !(lineStream >> value) )
                    LogicError("Invalid 'BOUNDS' section");
                const Int row =
                  meta_.upperBoundOffset + variableData.upperBoundIndex;

                // G(row,column) = 1
                entry.type = AFFINE_LP_INEQUALITY_MATRIX;
                entry.row = row;
                entry.column = column;
                entry.value = 1;
                queuedEntries_.push_back( entry );

                // h(row) = value
                entry.type = AFFINE_LP_INEQUALITY_VECTOR;
                entry.row = row;
                entry.column = 0;
                entry.value = value;
                queuedEntries_.push_back( entry );
            }
            else if( boundMark_ == "LO" )
            {
                if( !(lineStream >> value) )
                    LogicError("Invalid 'BOUNDS' section");
                const Int row =
                  meta_.lowerBoundOffset + variableData.lowerBoundIndex;

                // G(row,column) = -1
                entry.type = AFFINE_LP_INEQUALITY_MATRIX;
                entry.row = row;
                entry.column = column;
                entry.value = -1;
                queuedEntries_.push_back( entry );

                // h(row) = -value
                entry.type = AFFINE_LP_INEQUALITY_VECTOR;
                entry.row = row;
                entry.column = 0;
                entry.value = -value;
                queuedEntries_.push_back( entry );
            }
            else if( boundMark_ == "FX" )
            {
                if( !(lineStream >> value) )
                    LogicError("Invalid 'BOUNDS' section");
                const Int row = meta_.fixedOffset + variableData.fixedIndex;

                // A(row,column) = 1
                entry.type = AFFINE_LP_EQUALITY_MATRIX;
                entry.row = row;
                entry.column = column;
                entry.value = 1;
                queuedEntries_.push_back( entry );

                // h(row) = value
                entry.type = AFFINE_LP_INEQUALITY_VECTOR;
                entry.row = row;
                entry.column = 0;
                entry.value = value;
                queuedEntries_.push_back( entry );
            }
            else if( boundMark_ == "FR" )
            {
                continue;
            }
            else if( boundMark_ == "MI" )
            {
                // These will be handled later.
            }
            else if( boundMark_ == "PL" )
            {
                // These will be handled later.
            }
            else
                LogicError("Invalid 'BOUNDS' section (unknown bound mark)");
        }
        else if( section_ == MPS_RANGES )
        {
            LogicError("The 'RANGES' section is not yet supported");
        }
        else
        {
            LogicError("Invalid MPS file");
        }
    }
    if( queuedEntries_.size() > 0 )
        return true;

    // Now iterate through the variable map to handle the nonpositive and
    // nonnegative bounds (since many of the former were perhaps implicit).
    // We handle them both here for the sake of consistency.
    while( queuedEntries_.size() == 0 &&
           variableIter_ != meta_.variableDict.end() )
    {
        const auto& data = variableIter_->second;
        const Int column = data.index;

        // Handle non-positive values.
        AffineLPEntry<double> entry;
        if( data.nonpositive )
        {
            // G(row,column) = 1
            const Int row = meta_.nonpositiveOffset + data.nonpositiveIndex;
            entry.type = AFFINE_LP_INEQUALITY_MATRIX;
            entry.row = row;
            entry.column = column;
            entry.value = 1;
            queuedEntries_.push_back( entry );

            // There is no need to explicitly set h(row) to zero.
        }

        // Handle non-negative values.
        if( data.nonnegative )
        {
            // G(row,column) = -1
            const Int row = meta_.nonnegativeOffset + data.nonnegativeIndex;
            entry.type = AFFINE_LP_INEQUALITY_MATRIX;
            entry.row = row;
            entry.column = column;
            entry.value = -1;
            queuedEntries_.push_back( entry );

            // There is no need to explicitly set h(row) to zero.
        }

        ++variableIter_;
    }
    return queuedEntries_.size() > 0;
}

AffineLPEntry<double> MPSReader::GetEntry()
{
    EL_DEBUG_CSE
    if( queuedEntries_.size() == 0 )
        LogicError("No entries are currently enqueued");
    AffineLPEntry<double> entry = queuedEntries_.back();
    queuedEntries_.pop_back();
    return entry;
}

const MPSMeta& MPSReader::Meta() const
{
    EL_DEBUG_CSE
    return meta_;
}

namespace read_mps {

template<typename Real>
void Helper
( AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const string& filename,
  bool compressed,
  bool minimize )
{
    EL_DEBUG_CSE
    if( compressed )
        LogicError("Compressed MPS is not yet supported");
    const bool upperBoundImplicitlyNonnegative = false;
    const bool metadataSummary = true;

    MPSReader reader
      ( filename, compressed, minimize, upperBoundImplicitlyNonnegative );
    const MPSMeta& meta = reader.Meta();
    if( metadataSummary )
        meta.PrintSummary();

    Zeros( problem.c, meta.n, 1 );
    Zeros( problem.A, meta.m, meta.n );
    Zeros( problem.b, meta.m, 1 );
    Zeros( problem.G, meta.k, meta.n );
    Zeros( problem.h, meta.k, 1 );
    while( reader.QueuedEntry() )
    {
        const AffineLPEntry<double> entry = reader.GetEntry();
        if( entry.type == AFFINE_LP_COST_VECTOR )
            problem.c(entry.row) = entry.value;
        else if( entry.type == AFFINE_LP_EQUALITY_MATRIX )
            problem.A(entry.row,entry.column) = entry.value;
        else if( entry.type == AFFINE_LP_EQUALITY_VECTOR )
            problem.b(entry.row) = entry.value;
        else if( entry.type == AFFINE_LP_INEQUALITY_MATRIX )
            problem.G(entry.row,entry.column) = entry.value;
        else /* entry.type == AFFINE_LP_INEQUALITY_VECTOR */
            problem.h(entry.row) = entry.value;
    }
}

template<typename Real>
void Helper
( AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
  const string& filename,
  bool compressed,
  bool minimize )
{
    EL_DEBUG_CSE
   
    // TODO(poulson): Consider loading on a single process and then distributing
    // the results instead.

    if( compressed )
        LogicError("Compressed MPS is not yet supported");
    const bool upperBoundImplicitlyNonnegative = false;
    const bool metadataSummary = true;

    MPSReader reader
      ( filename, compressed, minimize, upperBoundImplicitlyNonnegative );
    const MPSMeta& meta = reader.Meta();
    if( metadataSummary && problem.A.Grid().Rank() == 0 )
        meta.PrintSummary();

    Zeros( problem.c, meta.n, 1 );
    Zeros( problem.A, meta.m, meta.n );
    Zeros( problem.b, meta.m, 1 );
    Zeros( problem.G, meta.k, meta.n );
    Zeros( problem.h, meta.k, 1 );

    while( reader.QueuedEntry() )
    {
        const AffineLPEntry<double> entry = reader.GetEntry();
        if( entry.type == AFFINE_LP_COST_VECTOR )
            problem.c.Set( entry.row, 0, entry.value );
        else if( entry.type == AFFINE_LP_EQUALITY_MATRIX )
            problem.A.Set( entry.row, entry.column, entry.value );
        else if( entry.type == AFFINE_LP_EQUALITY_VECTOR )
            problem.b.Set( entry.row, 0, entry.value );
        else if( entry.type == AFFINE_LP_INEQUALITY_MATRIX )
            problem.G.Set( entry.row, entry.column, entry.value );
        else /* entry.type == AFFINE_LP_INEQUALITY_VECTOR */
            problem.h.Set( entry.row, 0, entry.value );
    }
}

template<typename Real>
void Helper
( AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const string& filename,
  bool compressed,
  bool minimize )
{
    EL_DEBUG_CSE
    if( compressed )
        LogicError("Compressed MPS is not yet supported");
    const bool upperBoundImplicitlyNonnegative = false;
    const bool metadataSummary = true;

    MPSReader reader
      ( filename, compressed, minimize, upperBoundImplicitlyNonnegative );
    const MPSMeta& meta = reader.Meta();
    if( metadataSummary )
        meta.PrintSummary();

    Output("m=",meta.m,", n=",meta.n,", k=",meta.k);
    Zeros( problem.c, meta.n, 1 );
    Zeros( problem.A, meta.m, meta.n );
    Zeros( problem.b, meta.m, 1 );
    Zeros( problem.G, meta.k, meta.n );
    Zeros( problem.h, meta.k, 1 );

    problem.A.Reserve( meta.numEqualityEntries );
    problem.G.Reserve( meta.numInequalityEntries );
    while( reader.QueuedEntry() )
    {
        const AffineLPEntry<double> entry = reader.GetEntry();
        if( entry.type == AFFINE_LP_COST_VECTOR )
            problem.c.Set( entry.row, 0, entry.value );
        else if( entry.type == AFFINE_LP_EQUALITY_MATRIX )
            problem.A.QueueUpdate( entry.row, entry.column, entry.value );
        else if( entry.type == AFFINE_LP_EQUALITY_VECTOR )
            problem.b.Set( entry.row, 0, entry.value );
        else if( entry.type == AFFINE_LP_INEQUALITY_MATRIX )
            problem.G.QueueUpdate( entry.row, entry.column, entry.value );
        else /* entry.type == AFFINE_LP_INEQUALITY_VECTOR */
            problem.h.Set( entry.row, 0, entry.value );
    }
    problem.A.ProcessQueues();
    problem.G.ProcessQueues();
}

template<typename Real>
void Helper
( AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
  const string& filename,
  bool compressed,
  bool minimize )
{
    EL_DEBUG_CSE
    if( compressed )
        LogicError("Compressed MPS is not yet supported");
    const bool upperBoundImplicitlyNonnegative = false;
    const bool metadataSummary = false;

    MPSReader reader
      ( filename, compressed, minimize, upperBoundImplicitlyNonnegative );
    const MPSMeta& meta = reader.Meta();
    if( metadataSummary && problem.A.Grid().Rank() == 0 )
        meta.PrintSummary();

    Zeros( problem.c, meta.n, 1 );
    Zeros( problem.A, meta.m, meta.n );
    Zeros( problem.b, meta.m, 1 );
    Zeros( problem.G, meta.k, meta.n );
    Zeros( problem.h, meta.k, 1 );

    bool passive=true;
    problem.A.Reserve( meta.numEqualityEntries );
    problem.G.Reserve( meta.numInequalityEntries );
    while( reader.QueuedEntry() )
    {
        const AffineLPEntry<double> entry = reader.GetEntry();
        if( entry.type == AFFINE_LP_COST_VECTOR )
            problem.c.Set( entry.row, 0, entry.value );
        else if( entry.type == AFFINE_LP_EQUALITY_MATRIX )
            problem.A.QueueUpdate( entry.row, entry.column, entry.value, passive );
        else if( entry.type == AFFINE_LP_EQUALITY_VECTOR )
            problem.b.Set( entry.row, 0, entry.value );
        else if( entry.type == AFFINE_LP_INEQUALITY_MATRIX )
            problem.G.QueueUpdate( entry.row, entry.column, entry.value, passive );
        else /* entry.type == AFFINE_LP_INEQUALITY_VECTOR */
            problem.h.Set( entry.row, 0, entry.value );
    }
    problem.A.ProcessLocalQueues();
    problem.G.ProcessLocalQueues();
}

} // namespace read_mps

template<class MatrixType,class VectorType>
void ReadMPS
( AffineLPProblem<MatrixType,VectorType>& problem,
  const string& filename,
  bool compressed,
  bool minimize )
{
    EL_DEBUG_CSE
    read_mps::Helper( problem, filename, compressed, minimize );
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
