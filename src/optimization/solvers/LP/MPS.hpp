/*
   Copyright (c) 2009-2017, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

namespace {

const string invalidRowTokenString =
  "Invalid token in a line of the 'ROWS' section.";
const string tooFewRowTokensString =
  "Too few tokens in a line of the 'ROWS' section.";
const string tooManyRowTokensString =
  "Too many tokens in a line of the 'ROWS' section: please check that your "
  "variable names do not have spaces in them.";

const string invalidColumnTokenString =
  "Invalid token in a line of the 'COLUMNS' section.";
const string tooFewColumnTokensString =
  "Too few tokens in a line of the 'COLUMNS' section.";
const string tooManyColumnTokensString =
  "Too many tokens in a line of the 'COLUMNS' section: please check that your "
  "variable names do not have spaces in them.";

} // anonymous namespace

// We will follow the de facto convention of assuming that variables live within
// [0,\infty) if no bounds are specified and only modify the assumption as
// little as possible when explicit bounds are given. For example, if only a
// positive upper bound of 7 is specified, the interval will be changed to
// [0,7]. If a positive lower bound of 4 was instead specified, the interval
// would become [4,\infty).
//
// If the upper bound is negative, then the lower bound is reset to -\infty;
// if the upper bound is exactly zero, then there is disagreement in how to
// handle the lower bound. Many MPS readers reset the lower bound to -\infty,
// but ILOG famously does not. We provide an option to choose between the two
// behaviors.
//
// Handling 'RANGES' sections requires some extra care
// (Cf. http://lpsolve.sourceforge.net/5.5/mps-format.htm). In particular, the
// 'RANGES' section allowes 'GREATER', 'LESSER', and 'EQUALITY' rows to be
// transformed into two-sided bounds. That is to say, constraints of the form
//
//   g_i x >= d_i,
//   g_i x <= d_i, or
//   g_i x = d_i,
//
// become transformed into the form
//
//   d_i <= g_i x <= d_i + |r|, or
//   d_i - |r| <= g_i x <= d_i,
//
// depending upon the sign of the variable 'r' specified in the 'RANGES' section
// (in a format that mirrors the 'RHS' section). More specifically, the sign of
// 'r' is irrelevant for 'GREATER' and 'LESSER' rows, but a positive value for
// 'r' cases an 'EQUALITY' constraint to become the *lower* bound of an interval
// constraint of width 'r', whereas a negative value for 'r' transforms the
// 'EQUALITY' constraint into an *upper* bound on an interval of length '|r|'.
//
// It is of use to recognize that, if a 'RANGES' section modifies an 'EQUALITY'
// constraint, then we can henceforth treat said equality constraint as a
// 'GREATER' row if 'r' was positive, and as a 'LESSER' row if 'r' was negative.
// Thus, we are left with deciding how best to introduce additional inequality
// constraints.
//

enum MPSSection {
  MPS_NAME,
  MPS_ROWS,
  MPS_COLUMNS,
  MPS_RHS,
  MPS_BOUNDS,
  MPS_RANGES,
  MPS_OBJSENSE,
  MPS_MARKER,
  MPS_SOS,
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
  Int numNonzeros=0; // We will delete rows with no nonzeros (e.g., for tuff).

  // For handling range modifiers.
  bool hasRangeModifier=false;
  Int rangeIndex;
};

struct MPSRangeData
{
  string rowName;
  bool rhsIsLowerBound;
  double rangeSize;

  double lowerBound=0.;
  Int lowerBoundIndex=-1;

  double upperBound=0.;
  Int upperBoundIndex=-1;
};

struct MPSVariableData
{
  Int index;
  Int numNonzeros=0;

  bool lowerBounded=false;
  Int lowerBoundIndex=-1;
  double lowerBound;

  bool upperBounded=false;
  Int upperBoundIndex=-1;
  double upperBound;

  bool fixed=false;
  Int fixedIndex=-1;
  double fixedValue;

  bool free=false;
  Int freeIndex=-1;

  bool nonpositive=false;
  Int nonpositiveIndex=-1;

  bool nonnegative=false;
  Int nonnegativeIndex=-1;
};

struct MPSTrivialEquality
{
  string variableName;
  double singleNonzero;
};

void LPMPSMeta::PrintSummary() const
{
    Output("LPMPSMeta summary:");
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
    Output("  m=",m,", n=",n,", k=",k);
}

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
// where:
//
//   'G0 x <= h0' consists of the 'lesser' rows,
//
//   'G1 x <= h1' is the negation of the 'greater' rows,
//
//   'G2 x <= h2' consists of the upper bounds (each row of 'G0' is all zeros
//   except for a single one),
//
//   'G3 x <= h3' is the negation of the lower bounds (each row of 'G3' is all
//   zeros except for a single negative one),
//
//   'G4 x <= h4' is the set of nonpositive bounds, and
//
//   'G5 x <= h5' is the set of nonnegative bounds.
//
// The ordering of 'x' is: nontrivial variables, auxiliary variables
// (introduced by RANGE modifiers), and trivially fixed variables.
//
class MPSReader
{
public:
    MPSReader
    ( const string& filename,
      bool compressed=false,
      bool minimize=true,
      bool keepNonnegativeWithZeroUpperBound=true );
    // The PILOT netlib lp_data model appears to require
    // 'keepNonnegativeWithZeroUpperBound=true'.

    // Attempt to enqueue another entry and return true if successful.
    bool QueuedEntry();

    // Only one call is allowed per call to 'QueuedEntry'.
    AffineLPEntry<double> GetEntry();

    const LPMPSMeta& Meta() const;

    static bool IsDataLine( std::stringstream& lineStream );
    static MPSSection DecodeSection( const std::string& sectionToken );

private:
    string filename_;
    std::ifstream file_;

    bool minimize_;
    bool keepNonnegativeWithZeroUpperBound_;
    LPMPSMeta meta_;

    std::map<string,MPSRowData> rowDict_;
    std::map<string,MPSVariableData> variableDict_;
    vector<MPSRangeData> rangeList_;

    // A map from the original MPS_EQUALITY_ROW row name to the trivial equality
    // struct, which stores the variable name and the floating-point value
    // of the single nonzero in the row (eventually).
    std::map<string,MPSTrivialEquality> trivialEqualityDict_;

    std::set<string> emptyLesserRows_;
    std::set<string> emptyGreaterRows_;

    vector<AffineLPEntry<double>> queuedEntries_;

    MPSSection section_=MPS_NONE;

    // The RHS section typically has a name followed by either one or two pairs
    // per row, but some models (e.g., dfl001.mps) do not involve a name.
    bool initializedRHSSection_=false;
    bool rhsHasName_;

    // The RHS section typically has a name followed by either one or two pairs
    // per row, but some models do not involve a name.
    bool initializedRangesSection_=false;
    bool rangesHasName_;

    // The BOUNDS section typically has a bound type marker, followed by a
    // bound set name, followed by a variable name, and, if applicable,
    // a numeric value). But some models (e.g., dfl001.mps) do not involve a
    // bound set name.
    bool initializedBoundsSection_=false;
    bool boundsHasName_;

    // For the final loop over the variable dictionary to extract any bounds.
    typename std::map<string,MPSVariableData>::const_iterator variableIter_;

    // For the final loop over bounds on auxiliary variables.
    Int auxBoundIndex_=0;

    void ProcessRowDataLine( const string& line );
    void ProcessColumnDataLine( const string& line );
    void QueueColumnDataLine( const string& line );
    void ProcessRHSDataLine( const string& line );
    void QueueRHSDataLine( const string& line );
    void ProcessBoundsDataLine( const string& line );
    void ProcessRangesDataLine( const string& line );

    // If any remain, queue an entry associated with enforcing a variable bound.
    void QueueVariableBound();

    // After reading in the row data, this can be used to fast-forward to see if
    // apply any 'RANGES' modifiers before processing the column data.
    //
    // Any 'RANGES' induced constraint of the form
    //
    //   lowerBound <= g x <= upperBound
    //
    // can be expressed as
    //
    //   [g, -1] [x; zeta] = 0,
    //   lowerBound <= zeta <= upperBound,
    //
    // by introducing the auxiliary variable 'zeta'. And, if we are
    // to use this approach, the initial 'EQUALITY', 'GREATER', or
    // 'LESSER' row 'g' becomes an 'EQUALITY' row, and so we
    // temporarily fast-forward through the file to check for
    // 'RANGES' modifiers after processing the 'ROWS' section.
    //
    void PeekAtRanges( const string& sectionToken );

    // After reading in the column data, we can iterate through the row map to
    // delete any empty rows and convert any equality rows with a single nonzero
    // to a fixed state.
    void SimplifyRowTypes();

    // After having deleted empty rows, we can iterate over the MPSRange list
    // and assign the unique range index associated with each range modifier
    // within the MPSRowData.
    void SyncRangeData();

    // Once we are sure that the explict bound of each range has been filled in,
    // we can set the implicit bound to the explicit bound plus-or-minus the
    // interval size.
    void FillImplicitRangeBounds();

    // Make another pass through the COLUMNS section to store the variable name
    // associated with each trivial equality row.
    void StoreTrivialEqualityRowVariableNames
    ( const std::ifstream::pos_type& colSectionBeg );

    // Iterate through the variable map and make use of the requested
    // conventions for counting the number of bounds of each type.
    // Also warn if there are possibly conflicting bound types.
    void CountBoundTypes();
};

bool MPSReader::IsDataLine( std::stringstream& lineStream )
{
    const char firstChar = lineStream.peek();
    return firstChar == ' ' ||
           firstChar == '\t' ||
           firstChar == '*' ||
           firstChar == '#';
}

MPSSection MPSReader::DecodeSection( const string& sectionToken )
{
    EL_DEBUG_CSE
    MPSSection section;
    if( sectionToken == "NAME" )
    {
        section = MPS_NAME;
    }
    else if( sectionToken == "OBJSENSE" )
    {
        section = MPS_OBJSENSE;
    }
    else if( sectionToken == "ROWS" )
    {
        section = MPS_ROWS;
    }
    else if( sectionToken == "COLUMNS" )
    {
        section = MPS_COLUMNS;
    }
    else if( sectionToken == "RHS" )
    {
        section = MPS_RHS;
    }
    else if( sectionToken == "BOUNDS" )
    {
        section = MPS_BOUNDS;
    }
    else if( sectionToken == "RANGES" )
    {
        section = MPS_RANGES;
    }
    else if( sectionToken == "ENDATA" )
    {
        section = MPS_END;
    }
    else if( sectionToken == "MARKER" )
    {
        section = MPS_MARKER;
    }
    else if( sectionToken == "SOS" )
    {
        section = MPS_SOS;
    }
    else
    {
        LogicError("Section token ",sectionToken," is not recognized");
    }
    return section;
}

void MPSReader::ProcessRowDataLine( const string& line )
{
    EL_DEBUG_CSE
    std::stringstream stream( line );
    string rowType;
    if( !(stream >> rowType) )
        LogicError(tooFewRowTokensString);
    string rowName;
    if( !(stream >> rowName) )
        LogicError(tooFewRowTokensString);
    string token;
    if( stream >> token )
        LogicError(tooManyRowTokensString);

    MPSRowData rowData;
    // We set the 'typeIndex' fields later since it is not uncommon
    // (e.g., see tuff.mps) for rows to be empty.
    if( rowType == "L" )
    {
        rowData.type = MPS_LESSER_ROW;
        rowDict_[rowName] = rowData;
    }
    else if( rowType == "G" )
    {
        rowData.type = MPS_GREATER_ROW;
        rowDict_[rowName] = rowData;
    }
    else if( rowType == "E" )
    {
        rowData.type = MPS_EQUALITY_ROW;
        rowDict_[rowName] = rowData;
    }
    else if( rowType == "N" )
    {
        rowData.type = MPS_NONCONSTRAINING_ROW;
        rowDict_[rowName] = rowData;
        if( meta_.numNonconstrainingRows == 1 )
            meta_.costName = rowName;
    }
    else
        LogicError(invalidRowTokenString);
}

void MPSReader::ProcessColumnDataLine( const string& line )
{
    EL_DEBUG_CSE
    std::stringstream stream( line );
    string variableName;
    if( !(stream >> variableName) )
        LogicError(tooFewColumnTokensString);
    auto variableIter = variableDict_.find( variableName );
    if( variableIter == variableDict_.end() )
    {
        MPSVariableData variableData;
        variableDict_[variableName] = variableData;
        variableIter = variableDict_.find( variableName );
    }
    MPSVariableData& variableData = variableIter->second;

    // There should be either one or two pairs of entries left to read
    // from this line.
    for( Int pair=0; pair<2; ++pair )
    {
        string rowName;
        if( !(stream >> rowName) )
        {
            if( pair == 0 )
                LogicError(tooFewColumnTokensString);
            else
                break;
        }
        double value;
        if( !(stream >> value) )
            LogicError(tooManyColumnTokensString);
        ++variableData.numNonzeros;

        auto rowIter = rowDict_.find( rowName );
        if( rowIter == rowDict_.end() )
            LogicError("Could not find row ",rowName);
        auto& rowData = rowIter->second;
        if( rowData.type == MPS_EQUALITY_ROW )
        {
            // A(row,column) = value
            ++rowData.numNonzeros;
            ++meta_.numEqualityEntries;
        }
        else if( rowData.type == MPS_LESSER_ROW )
        {
            // G(row,column) = value
            ++rowData.numNonzeros;
            ++meta_.numInequalityEntries;
        }
        else if( rowData.type == MPS_GREATER_ROW )
        {
            // G(row,column) = -value
            ++rowData.numNonzeros;
            ++meta_.numInequalityEntries;
        }
        else if( rowData.type == MPS_NONCONSTRAINING_ROW )
        {
            // c(column) = value (assuming typeIndex==0)
            ++rowData.numNonzeros;
        }
        else
            LogicError("Unknown row data type");
    }
}

void MPSReader::QueueColumnDataLine( const string& line )
{
    EL_DEBUG_CSE
    double columnValue;
    std::stringstream stream( line );

    string variableName;
    if( !(stream >> variableName) )
        LogicError("Invalid 'COLUMNS' section");
    auto variableIter = variableDict_.find( variableName );
    if( variableIter == variableDict_.end() )
        LogicError("Could not find variable ",variableName);
    const MPSVariableData& variableData = variableIter->second;
    const Int column = variableData.index;

    // There should be either one or two pairs of entries left to read
    // from this line.
    for( Int pair=0; pair<2; ++pair )
    {
        string rowName;
        if( !(stream >> rowName) )
        {
            if( pair == 0 )
                LogicError("Invalid 'COLUMNS' section");
            else
                break;
        }
        if( !(stream >> columnValue) )
            LogicError("Invalid 'COLUMNS' section");

        AffineLPEntry<double> entry;
        auto rowIter = rowDict_.find( rowName );
        if( rowIter == rowDict_.end() )
        {
            auto trivialEqualityIter =
              trivialEqualityDict_.find( rowName );
            if( trivialEqualityIter == trivialEqualityDict_.end() )
            {
                LogicError("Could not find row ",rowName);
            }
            if( columnValue == 0. )
            {
                LogicError
                ("Trivial rows are not yet allowed to have their only "
                 "specified \"nonzero\" value be zero");
            }
            Output
            ("WARNING: Trivial equality row ",rowName,
             " only had a nonzero in column ",variableName);
            auto& trivialData = trivialEqualityIter->second;
            trivialData.variableName = variableName;
            trivialData.singleNonzero = columnValue;
            continue;
        }

        const auto& rowData = rowIter->second;
        if( rowData.type == MPS_EQUALITY_ROW )
        {
            // A(row,column) = columnValue
            const Int row = meta_.equalityOffset + rowData.typeIndex;
            entry.type = AFFINE_LP_EQUALITY_MATRIX;
            entry.row = row;
            entry.column = column;
            entry.value = columnValue;
            queuedEntries_.push_back( entry );
        }
        else if( rowData.type == MPS_LESSER_ROW )
        {
            // G(row,column) = columnValue
            const Int row = meta_.lesserOffset + rowData.typeIndex;
            entry.type = AFFINE_LP_INEQUALITY_MATRIX;
            entry.row = row;
            entry.column = column;
            entry.value = columnValue;
            queuedEntries_.push_back( entry );
        }
        else if( rowData.type == MPS_GREATER_ROW )
        {
            // G(row,column) = -columnValue
            const Int row = meta_.greaterOffset + rowData.typeIndex;
            entry.type = AFFINE_LP_INEQUALITY_MATRIX;
            entry.row = row;
            entry.column = column;
            entry.value = -columnValue;
            queuedEntries_.push_back( entry );
        }
        else if( rowData.type == MPS_NONCONSTRAINING_ROW )
        {
            if( rowData.typeIndex == 0 )
            {
                // c(column) = columnValue
                entry.type = AFFINE_LP_COST_VECTOR;
                entry.row = column;
                entry.column = 0;
                entry.value = minimize_ ? columnValue : -columnValue;
                queuedEntries_.push_back( entry );
            }
            // else?
        }
    }
}

void MPSReader::ProcessRHSDataLine( const string& line )
{
    EL_DEBUG_CSE
    if( !initializedRHSSection_ )
    {
        std::stringstream testStream( line );
        Int numTokens=0;
        string token;
        while( testStream >> token )
            ++numTokens;
        if( numTokens == 2 || numTokens == 4 )
        {
            // There were either one or two pairs with no name.
            rhsHasName_ = false;
            meta_.numRHS = 1;
        }
        else if( numTokens == 3 || numTokens == 5 )
        {
            // There were either one or two pairs with a name.
            rhsHasName_ = true;
        }
        else
        {
            LogicError("Invalid 'RHS' section (1)");
        }
        initializedRHSSection_ = true;
    }

    std::stringstream stream( line );
    if( rhsHasName_ )
    {
        string rhsNameCandidate;
        if( !(stream >> rhsNameCandidate) )
            LogicError("Invalid 'RHS' section (2)");
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
    }

    // There should be either one or two pairs of entries left to read
    // from this line.
    for( Int pair=0; pair<2; ++pair )
    {
        string rowName;
        if( !(stream >> rowName) )
        {
            if( pair == 0 )
                LogicError("Invalid 'RHS' section (3)");
            else
                break;
        }
        double value;
        if( !(stream >> value) )
            LogicError("Invalid 'RHS' section (4)");
        // TODO(poulson): Store fixed value?
    }
}

void MPSReader::QueueRHSDataLine( const string& line )
{
    double rhsValue;
    std::stringstream stream( line );

    if( rhsHasName_ )
    {
        string token;
        if( !(stream >> token) )
            LogicError("Invalid 'RHS' section");
    }

    // There should be either one or two pairs of entries left to read
    // from this line.
    for( Int pair=0; pair<2; ++pair )
    {
        string rowName;
        if( !(stream >> rowName) )
        {
            if( pair == 0 )
                LogicError("Invalid 'RHS' section");
            else
                break;
        }
        if( !(stream >> rhsValue) )
            LogicError("Invalid 'RHS' section");
        auto rowIter = rowDict_.find( rowName );
        if( rowIter == rowDict_.end() )
        {
            if( emptyLesserRows_.count(rowName) == 1 )
            {
                if( rhsValue < 0. )
                {
                    LogicError
                    ("Row ",rowName,
                     " has an invalid trivial upper bound of ",
                     rhsValue);
                }
                else
                {
                    Output
                    ("WARNING: Skipping trivial lesser row ",rowName,
                     " which has upper bound ",rhsValue);
                }
                continue;
            }
            else if( emptyGreaterRows_.count(rowName) == 1 )
            {
                if( rhsValue > 0. )
                {
                    LogicError
                    ("Row ",rowName,
                     " has an invalid trivial lower bound of ",
                     rhsValue);
                }
                else
                {
                    Output
                    ("WARNING: Skipping trivial greater row ",rowName,
                     " which has lower bound ",rhsValue);
                }
                continue;
            }

            auto trivialEqualityIter =
              trivialEqualityDict_.find( rowName );
            if( trivialEqualityIter == trivialEqualityDict_.end() )
            {
                LogicError("Could not find trivial row ",rowName);
            }
            const auto& trivialData = trivialEqualityIter->second;
            auto variableIter =
              variableDict_.find( trivialData.variableName );
            if( variableIter == variableDict_.end() )
            {
                LogicError
                ("Could not find variable ",trivialData.variableName,
                 " for trivial equality from row ",rowName);
            }
            variableIter->second.fixedValue =
              rhsValue / trivialData.singleNonzero;
            continue;
        }
        AffineLPEntry<double> entry;
        const auto& rowData = rowIter->second;
        if( rowData.hasRangeModifier )
        {
            // Fill in the explicit bound for the RANGE and queue the two-sided
            // bound later.
            auto& rangeData = rangeList_[rowData.rangeIndex];
            if( rangeData.rhsIsLowerBound )
            {
                rangeData.lowerBound = rhsValue;
            }
            else
            {
                rangeData.upperBound = rhsValue;
            }
            continue;
        }
        if( rowData.type == MPS_EQUALITY_ROW )
        {
            // b(row) = rhsValue
            const Int row = meta_.equalityOffset + rowData.typeIndex;
            entry.type = AFFINE_LP_EQUALITY_VECTOR;
            entry.row = row;
            entry.column = 0;
            entry.value = rhsValue;
            queuedEntries_.push_back( entry );
        }
        else if( rowData.type == MPS_LESSER_ROW )
        {
            // h(row) = rhsValue
            const Int row = meta_.lesserOffset + rowData.typeIndex;
            entry.type = AFFINE_LP_INEQUALITY_VECTOR;
            entry.row = row;
            entry.column = 0;
            entry.value = rhsValue;
            queuedEntries_.push_back( entry );
        }
        else if( rowData.type == MPS_GREATER_ROW )
        {
            // h(row) = -rhsValue
            const Int row = meta_.greaterOffset + rowData.typeIndex;
            entry.type = AFFINE_LP_INEQUALITY_VECTOR;
            entry.row = row;
            entry.column = 0;
            entry.value = -rhsValue;
            queuedEntries_.push_back( entry );
        }
        else if( rowData.type == MPS_NONCONSTRAINING_ROW )
        {
            Output("WARNING: Nonsensical RHS for nonconstrained row");
        }
    }
}

void MPSReader::ProcessBoundsDataLine( const string& line )
{
    EL_DEBUG_CSE
    if( !initializedBoundsSection_ )
    {
        std::stringstream testStream( line );
        string boundMark;
        if( !(testStream >> boundMark) )
            LogicError("Invalid 'BOUNDS' section");
        string token;
        if( !(testStream >> token) )
            LogicError("Invalid 'BOUNDS' section");

        Int numTokens=2;
        string rhsToken;
        while( testStream >> rhsToken )
            ++numTokens;
        if( boundMark == "LO" ||
            boundMark == "UP" ||
            boundMark == "FX" )
        {
            if( numTokens == 4 )
                boundsHasName_ = true;
            else if( numTokens == 3 )
                boundsHasName_ = false;
            else
                LogicError("Invalid ",boundMark," 'BOUNDS' line");
        }
        else if( boundMark == "FR" ||
                 boundMark == "MI" ||
                 boundMark == "PL" )
        {
            if( numTokens == 3 )
                boundsHasName_ = true;
            else if( numTokens == 2 )
                boundsHasName_ = false;
            else
                LogicError("Invalid ",boundMark," 'BOUNDS' line");
        }
        else
            LogicError("Unknown bound mark ",boundMark);

        // If there is a name, it occurred in the second position.
        if( boundsHasName_ )
            meta_.boundName = token;

        initializedBoundsSection_ = true;
    }

    std::stringstream stream( line );
    // We already have the first token of a bounding row, which should
    // be of the same general form as
    //
    //   FX BOUNDROW VARIABLENAME 1734.
    //
    // in the case of 'VARIABLENAME' being fixed ('FX') at the value
    // 1734 (with this problem's bound name being 'BOUNDROW').
    string boundMark;
    if( !(stream >> boundMark) )
        LogicError("Invalid 'BOUNDS' section");
    if( boundsHasName_ )
    {
        string boundNameCandidate;
        if( !(stream >> boundNameCandidate) )
            LogicError("Invalid 'BOUNDS' section");
        if( meta_.boundName != boundNameCandidate )
            LogicError
            ("Only single problem instances are currently supported "
             "(multiple bound names were encountered)");
    }
    string variableName;
    if( !(stream >> variableName) )
        LogicError("Invalid 'BOUNDS' section");
    auto variableIter = variableDict_.find( variableName );
    if( variableIter == variableDict_.end() )
        LogicError
        ("Invalid 'BOUNDS' section (name ",variableName," not found)");
    MPSVariableData& data = variableIter->second;
    if( boundMark == "UP" )
    {
        double value;
        if( !(stream >> value) )
            LogicError("Invalid 'BOUNDS' section");
        data.upperBounded = true;
        data.upperBound = value;
    }
    else if( boundMark == "LO" )
    {
        double value;
        if( !(stream >> value) )
            LogicError("Invalid 'BOUNDS' section");
        data.lowerBounded = true;
        data.lowerBound = value;
    }
    else if( boundMark == "FX" )
    {
        double value;
        if( !(stream >> value) )
            LogicError("Invalid 'BOUNDS' section");
        data.fixed = true;
        data.fixedValue = value;
    }
    else if( boundMark == "FR" )
        data.free = true;
    else if( boundMark == "MI" )
        data.nonpositive = true;
    else if( boundMark == "PL" )
        data.nonnegative = true;
    else
        LogicError("Invalid 'BOUNDS' section (unknown bound mark)");
}

void MPSReader::ProcessRangesDataLine( const string& line )
{
    EL_DEBUG_CSE
    if( !initializedRangesSection_ )
    {
        std::stringstream testStream( line );
        Int numTokens=0;
        string rangesToken;
        while( testStream >> rangesToken )
            ++numTokens;
        if( numTokens == 2 || numTokens == 4 )
        {
            // There were either one or two pairs with no name.
            rangesHasName_ = false;
            meta_.numRanges = 1;
        }
        else if( numTokens == 3 || numTokens == 5 )
        {
            // There were either one or two pairs with a name.
            rangesHasName_ = true;
        }
        else
        {
            LogicError("Invalid 'RANGES' section (1)");
        }
        initializedRangesSection_ = true;
    }

    std::stringstream stream( line );
    if( rangesHasName_ )
    {
        string rangesNameCandidate;
        if( !(stream >> rangesNameCandidate) )
            LogicError("Invalid 'RANGES' section (2)");
        if( meta_.numRanges == 0 )
        {
            // We should currently have that rangesName == "".
            meta_.rangesName = rangesNameCandidate;
            meta_.numRanges = 1;
        }
        else if( rangesNameCandidate != meta_.rangesName )
            LogicError
            ("Only single problem instances are currently supported "
             "(multiple range names were encountered)");
    }

    // There should be either one or two pairs of entries left to read
    // from this line.
    for( Int pair=0; pair<2; ++pair )
    {
        string rowName;
        if( !(stream >> rowName) )
        {
            if( pair == 0 )
                LogicError("Invalid 'RANGES' section (3)");
            else
                break;
        }
        double value;
        if( !(stream >> value) )
            LogicError("Invalid 'RANGES' section (4)");
        if( value == 0. )
            LogicError("Tried to force an empty RANGE");
        auto rowIter = rowDict_.find( rowName );
        if( rowIter == rowDict_.end() )
            LogicError("Could not find row ",rowName);
        MPSRangeData rangeData;
        rangeData.rowName = rowName;
        rangeData.rangeSize = Abs(value);

        auto& rowData = rowIter->second;
        if( rowData.type == MPS_EQUALITY_ROW )
        {
            rangeData.rhsIsLowerBound = value > 0.;
        }
        else if( rowData.type == MPS_GREATER_ROW )
        {
            rangeData.rhsIsLowerBound = true; 
        }
        else if( rowData.type == MPS_LESSER_ROW )
        {
            rangeData.rhsIsLowerBound = false;
        }
        else
        {
            LogicError
            ("Tried to apply a RANGE modifier to an invalid row type");
        }
        rangeList_.push_back( rangeData );

        rowData.type = MPS_EQUALITY_ROW;
        rowData.hasRangeModifier = true;
        // We will set the row index by iterating through the range data list
        // and either deleting entries that had deleted rows or assigning the
        // appropriate index when they exist.
    }
}

void MPSReader::QueueVariableBound()
{
    EL_DEBUG_CSE
    AffineLPEntry<double> entry;
    while( queuedEntries_.size() == 0 &&
           variableIter_ != variableDict_.end() )
    {
        const auto& data = variableIter_->second;
        const Int column = data.index;

        if( data.upperBounded )
        {
            const Int row = meta_.upperBoundOffset + data.upperBoundIndex;

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
            entry.value = data.upperBound;
            queuedEntries_.push_back( entry );
        }

        if( data.lowerBounded )
        {
            const Int row = meta_.lowerBoundOffset + data.lowerBoundIndex;

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
            entry.value = -data.lowerBound;
            queuedEntries_.push_back( entry );
        }

        if( data.fixed )
        {
            const Int row = meta_.fixedOffset + data.fixedIndex;

            // A(row,column) = 1
            entry.type = AFFINE_LP_EQUALITY_MATRIX;
            entry.row = row;
            entry.column = column;
            entry.value = 1;
            queuedEntries_.push_back( entry );

            // h(row) = value
            entry.type = AFFINE_LP_EQUALITY_VECTOR;
            entry.row = row;
            entry.column = 0;
            entry.value = data.fixedValue;
            queuedEntries_.push_back( entry );
        }

        // Handle non-positive values.
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
    if( queuedEntries_.size() > 0 )
        return;
    if( auxBoundIndex_ < rangeList_.size() )
    {
        const Int column =
          variableDict_.size() - meta_.numFixedBounds + auxBoundIndex_;
        const auto& data = rangeList_[auxBoundIndex_++];
        AffineLPEntry<double> entry;

        // Queue the auxiliary variable introductory entry of -1.
        {
            const auto rowIter = rowDict_.find( data.rowName );
            if( rowIter == rowDict_.end() )
                RuntimeError("Could not link range data to row dict");
            // G(row,column) = -1
            const Int row = meta_.equalityOffset + rowIter->second.typeIndex;
            entry.type = AFFINE_LP_EQUALITY_MATRIX;
            entry.row = row;
            entry.column = column;
            entry.value = -1;
            queuedEntries_.push_back( entry );
        }

        // Queue the lower bound
        {
            // G(row,column) = -1
            const Int row = meta_.lowerBoundOffset + data.lowerBoundIndex;
            entry.type = AFFINE_LP_INEQUALITY_MATRIX;
            entry.row = row;
            entry.column = column;
            entry.value = -1;
            queuedEntries_.push_back( entry );

            // h(row) = -lowerBound
            entry.type = AFFINE_LP_INEQUALITY_VECTOR;
            entry.row = row;
            entry.column = 0;
            entry.value = -data.lowerBound;
            queuedEntries_.push_back( entry );
        }

        // Queue the upper bound
        {
            // G(row,column) = 1
            const Int row = meta_.upperBoundOffset + data.upperBoundIndex;
            entry.type = AFFINE_LP_INEQUALITY_MATRIX;
            entry.row = row;
            entry.column = column;
            entry.value = 1;
            queuedEntries_.push_back( entry );

            // h(row) = upperBound
            entry.type = AFFINE_LP_INEQUALITY_VECTOR;
            entry.row = row;
            entry.column = 0;
            entry.value = data.upperBound;
            queuedEntries_.push_back( entry );
        }
    }
}

void MPSReader::PeekAtRanges( const string& sectionToken )
{
    EL_DEBUG_CSE
    std::ifstream::pos_type pos = file_.tellg();
    auto section = DecodeSection( sectionToken );

    std::string token, line, rowName;
    while( std::getline( file_, line ) )
    {
        std::stringstream stream( line );
        const bool isDataLine = IsDataLine(stream);

        // The first token on each line should be a string. We will check it for
        // equivalence with each section string.
        if( !(stream >> token) )
        {
            // This line only consists of whitespace.
            continue;
        }
        if( !isDataLine )
        {
            section = DecodeSection(token);
            if( section == MPS_END )
                break;
            continue;
        }
        if( section != MPS_RANGES )
            continue;

        ProcessRangesDataLine( line );
    }
    file_.seekg(pos);
}

void MPSReader::SimplifyRowTypes()
{
    EL_DEBUG_CSE
    for( auto iter=rowDict_.begin(); iter!=rowDict_.end(); )
    {
        if( iter->second.numNonzeros == 0 )
        {
            if( iter->second.type == MPS_NONCONSTRAINING_ROW )
            {
                Output("WARNING: Objective was entirely zero.");
                iter->second.typeIndex = 0;
                ++iter;
            }
            else
            {
                if( iter->second.type == MPS_EQUALITY_ROW )
                {
                    Output
                    ("WARNING: Deleting empty equality row ",
                     iter->first);
                }
                else if( iter->second.type == MPS_GREATER_ROW )
                {
                    emptyGreaterRows_.insert(iter->first);
                    Output
                    ("WARNING: Deleting empty greater row ",
                     iter->first);
                }
                else if( iter->second.type == MPS_LESSER_ROW )
                {
                    emptyLesserRows_.insert(iter->first);
                    Output
                    ("WARNING: Deleting empty lesser row ",
                     iter->first);
                }
                else
                {
                    LogicError("Unknown empty row type");
                }
                // Delete this entry.
                iter = rowDict_.erase(iter);
            }
        }
        else if( iter->second.numNonzeros == 1 &&
                 iter->second.type == MPS_EQUALITY_ROW &&
                 !iter->second.hasRangeModifier )
        {
            // Convert to a fixed value. We will divide the
            // corresponding RHS value by the single nonzero.
            MPSTrivialEquality trivialEquality;
            trivialEqualityDict_[iter->first] = trivialEquality;
            iter = rowDict_.erase( iter );
            // Subtract one from the equality entries
            // (which will be added back later when iterating over
            // all of the fixed bounds).
            --meta_.numEqualityEntries;

            // TODO(poulson): Exploit the fact that a trivial equality row that
            // involves a RANGE modifier would imply that the auxiliary variable
            // is proportional to the single nonzero column index's variable,
            // and so we could instead apply the (scaled) range requirements
            // on said variable. But perhaps this is rare enough to not
            // necessitate special handling.
        }
        else
        {
            if( iter->second.type == MPS_EQUALITY_ROW )
            {
                iter->second.typeIndex = meta_.numEqualityRows++;
            }
            else if( iter->second.type == MPS_GREATER_ROW )
            {
                iter->second.typeIndex = meta_.numGreaterRows++;
            }
            else if( iter->second.type == MPS_LESSER_ROW )
            {
                iter->second.typeIndex = meta_.numLesserRows++;
            }
            else if( iter->second.type == MPS_NONCONSTRAINING_ROW )
            {
                iter->second.typeIndex = 0;
            }
            else
                LogicError("Unknown row type");

            if( iter->second.hasRangeModifier )
            {
                // Account for a nonzero to introduce the auxiliary variable.
                ++meta_.numEqualityEntries;
                // Account for two nonzeros to bound the auxliary variable.
                meta_.numInequalityEntries += 2;
            }

            ++iter;
        }
    }
}

void MPSReader::SyncRangeData()
{
    EL_DEBUG_CSE
    Int rangeIndex = 0;
    for( auto iter=rangeList_.begin(); iter!=rangeList_.end(); )
    {
        // Attempt to find the row data associated with this range modifier.
        auto rowIter = rowDict_.find( iter->rowName );
        if( rowIter == rowDict_.end() )
        {
            // The row was deleted (due to not having any nonzeros). We will
            // therefore delete this member of the range list.
            iter = rangeList_.erase( iter ); 
        }
        else
        {
            rowIter->second.rangeIndex = rangeIndex++;
            ++iter;
        }
    }
}

void MPSReader::FillImplicitRangeBounds()
{
    EL_DEBUG_CSE
    for( auto& rangeData : rangeList_ )
    {
        if( rangeData.rhsIsLowerBound )
        {
            rangeData.upperBound = rangeData.lowerBound + rangeData.rangeSize;
        }
        else
        {
            rangeData.lowerBound = rangeData.upperBound - rangeData.rangeSize;
        }
    }
}

void MPSReader::StoreTrivialEqualityRowVariableNames
( const std::ifstream::pos_type& colSectionBeg )
{
    EL_DEBUG_CSE
    auto pos = file_.tellg();
    file_.seekg( colSectionBeg );
    string columnLine, columnToken;
    while( std::getline( file_, columnLine ) )
    {
        std::stringstream columnStream( columnLine );
        const bool isColumnLine = IsDataLine(columnStream);
        if( !isColumnLine )
            break;
        string variableName;
        if( !(columnStream >> variableName) )
            continue;
        auto variableIter = variableDict_.find( variableName );
        if( variableIter == variableDict_.end() )
            LogicError(invalidColumnTokenString);
        MPSVariableData& variableData = variableIter->second;
        for( Int pair=0; pair<2; ++pair )
        {
            string rowName;
            if( !(columnStream >> rowName) )
            {
                if( pair == 0 )
                    LogicError(tooFewColumnTokensString);
                else
                    break;
            }
            double value;
            if( !(columnStream >> value) )
                LogicError(tooManyColumnTokensString);
            auto trivialEqualityIter =
              trivialEqualityDict_.find( rowName );
            if( trivialEqualityIter ==
                trivialEqualityDict_.end() )
                continue;
            if( value == 0. )
                LogicError
                ("Trivial equality rows with a zero \"nonzero\""
                 " value are not yet supported");
            auto& trivialData = trivialEqualityIter->second;
            trivialData.variableName = variableName;
            trivialData.singleNonzero = value;
            Output
            ("WARNING: Storing nonzero of ",value,
             " for ",variableName," for trivial equality row ",rowName);
            // We initialize at zero and overwrite if there is
            // relevant RHS data.
            variableData.fixed = true;
            variableData.fixedValue = 0;
        }
    }
    file_.seekg( pos );
}

void MPSReader::CountBoundTypes()
{
    EL_DEBUG_CSE
    for( auto& entry : variableDict_ )
    {
        auto& data = entry.second;

        // Handle explicit upper and lower bounds.
        if( data.upperBounded )
        {
            if( data.fixed )
            {
                if( data.fixedValue > data.upperBound )
                {
                    LogicError
                    ("Incompatible fixed value of ",data.fixedValue,
                     " and upper bound of ",data.upperBound);
                }
                else
                {
                    data.upperBounded = false;
                    data.nonnegative = false;
                    data.nonpositive = false;
                }
            }
            else if( data.lowerBounded )
            {
                if( data.lowerBound < data.upperBound )
                {
                    // This should be the standard case.
                    data.upperBoundIndex = meta_.numUpperBounds++;
                    ++meta_.numInequalityEntries;
                }
                else if( data.lowerBound == data.upperBound )
                {
                    // We will instead fix the variable.
                    data.upperBounded = false;
                    data.lowerBounded = false;
                    data.fixed = true;
                    data.fixedValue = data.upperBound;
                    Output
                    ("WARNING: Fixing ",entry.first," since the lower and "
                     "upper bounds were both ",data.fixedValue);
                }
                else
                {
                    LogicError
                    ("Cannot enforce a lower bound of ",data.lowerBound,
                     " and an upper bound of ",data.upperBound," for ",
                     entry.first);
                }
            }
            else
            {
                // Handling the default non-negativity is somewhat subtle and
                // varies between different MPS readers.
                if( data.upperBound > 0. )
                {
                    // Preserve the default non-negativity assumption.
                    data.upperBoundIndex = meta_.numUpperBounds++;
                    ++meta_.numInequalityEntries;
                    data.nonnegative = true;
                }
                else if( data.upperBound == 0. && !data.nonnegative )
                {
                    // There is disagreement on how to handle the lower bounds
                    // in the case of a zero upper bound. We print a warning
                    // to signify this corner case.
                    if( keepNonnegativeWithZeroUpperBound_ )
                    {
                        // Preserve the default non-negativity and fix the
                        // variable at zero.
                        data.nonnegative = false;
                        data.upperBounded = false;
                        data.fixed = true;
                        data.fixedValue = 0.;
                        Output
                        ("WARNING: Fixing ",entry.first,
                         " at zero due to zero upper bound. If this is not "
                         "desired, please set "
                         "'keepNonnegativeWithZeroUpperBound=false'");
                    }
                    else
                    {
                        // Do not enforce non-negativity.
                        data.upperBoundIndex = meta_.numUpperBounds++;
                        ++meta_.numInequalityEntries;
                        data.nonnegative = false;
                        Output
                        ("WARNING: Removing default non-negativity of ",
                         entry.first," due to zero upper bound. If this is "
                         "not desired, please set "
                         "'keepNonnegativeWithZeroUpperBound=true'");
                    }
                }
                else
                {
                    data.upperBoundIndex = meta_.numUpperBounds++;
                    ++meta_.numInequalityEntries;
                    if( data.nonnegative )
                        LogicError
                        ("Cannot have non-negative variable with an upper bound"
                         " of ",data.upperBound);
                }
            }
        }
        if( data.lowerBounded )
        {
            // Any conflicts with an upper bound are already resolved.
            if( data.fixed )
            {
                if( data.fixedValue < data.lowerBound )
                {
                    LogicError
                    ("Incompatible fixed value of ",data.fixedValue,
                     " and lower bound of ",data.lowerBound);
                }
                else
                {
                    data.lowerBounded = false;
                    data.nonnegative = false;
                    data.nonpositive = false;
                }
            }
            else
            {
                data.lowerBoundIndex = meta_.numLowerBounds++;
                ++meta_.numInequalityEntries;
            }
        }

        // Handle fixed values.
        if( data.fixed )
        {
            if( data.upperBounded )
                LogicError("Tried to fix and upper bound ",entry.first);
            if( data.lowerBounded )
                LogicError("Tried to fix and lower bound ",entry.first);
            if( data.free )
            {
                Output
                ("WARNING: Tried to fix and free ",entry.first,
                 ", so we will only mark it as fixed");
                data.free = false;
            }
            if( data.nonpositive )
                LogicError("Tried to fix and nonpositive ",entry.first);
            if( data.nonnegative )
                LogicError("Tried to fix and nonnegative ",entry.first);
            data.fixedIndex = meta_.numFixedBounds++;
            ++meta_.numEqualityEntries;
        }

        // Handle free values.
        if( data.free )
        {
            if( data.upperBounded )
                LogicError("Tried to free and upper bound ",entry.first);
            if( data.lowerBounded )
                LogicError("Tried to free and lower bound ",entry.first);
            if( data.nonnegative )
                LogicError
                ("Tried to free and make nonnegative bound ",entry.first);
            if( data.nonpositive )
                LogicError
                ("Tried to free and make nonpositive bound ",entry.first);
            if( data.fixed )
                LogicError("Tried to free and fix bound ",entry.first);
            data.freeIndex = meta_.numFreeBounds++;
        }

        // Handle non-positive values.
        if( data.nonpositive )
        {
            if( data.upperBounded )
                Output
                ("WARNING: Combined nonpositive constraint with upper bound");
            data.nonpositiveIndex = meta_.numNonpositiveBounds++;
            ++meta_.numInequalityEntries;
        }

        // Default to nonnegative if there were not any markings.
        const bool hasAMark =
          data.upperBounded ||
          data.lowerBounded ||
          data.fixed ||
          data.free ||
          data.nonpositive ||
          data.nonnegative;
        if( !hasAMark )
            data.nonnegative = true;

        // Handle non-negative values.
        if( data.nonnegative )
        {
            if( data.lowerBounded )
                Output
                ("WARNING: Combined nonnegative constraint with lower bound");
            data.nonnegativeIndex = meta_.numNonnegativeBounds++;
            ++meta_.numInequalityEntries;
        }
    }
    // Now loop over the auxiliary variables introduced by the RANGE modifiers.
    for( auto& rangeData : rangeList_ )
    {
        rangeData.lowerBoundIndex = meta_.numLowerBounds++;
        rangeData.upperBoundIndex = meta_.numUpperBounds++;
    }
}

MPSReader::MPSReader
( const string& filename,
  bool compressed,
  bool minimize,
  bool keepNonnegativeWithZeroUpperBound )
: filename_(filename),
  file_(filename.c_str()),
  minimize_(minimize),
  keepNonnegativeWithZeroUpperBound_(keepNonnegativeWithZeroUpperBound)
{
    EL_DEBUG_CSE
    if( compressed )
        LogicError("Compressed reads are not yet supported");
    if( !file_.is_open() )
        RuntimeError("Could not open ",filename);

    // For quickly reparsing the column data section if there were any trivial
    // equality rows that required retrieving the associated single nonzero
    // column name.
    std::ifstream::pos_type colSectionBeg;

    // TODO(poulson): Convert each token to upper-case letters before each
    // comparison. While capital letters are used by convention, they are
    // not required.
    MPSSection section = MPS_NONE;
    string line, token;
    while( std::getline( file_, line ) )
    {
        // We first determine which section we are in
        // ------------------------------------------
        std::stringstream sectionStream( line );
        const bool isDataLine = IsDataLine(sectionStream);

        // The first token on each line should be a string. We will check it for
        // equivalence with each section string.
        if( !(sectionStream >> token) )
        {
            // This line only consists of whitespace.
            continue;
        }

        if( !isDataLine )
        {
            // Handle any special section exit procedures.
            if( section == MPS_ROWS )
            {
                PeekAtRanges( token );
            }
            else if( section == MPS_COLUMNS )
            {
                SimplifyRowTypes();
                SyncRangeData();
                if( trivialEqualityDict_.size() != 0 )
                {
                    StoreTrivialEqualityRowVariableNames( colSectionBeg );
                }
            }

            // Determine the new section
            section = DecodeSection( token );

            // Handle any special section entrance procedures.
            if( section == MPS_NAME )
            {
                if( meta_.name != "" )
                    LogicError("Multiple 'NAME' sections");
                if( !(sectionStream >> meta_.name) )
                    LogicError("Missing 'NAME' string");
            }
            else if( section == MPS_OBJSENSE )
            {
                LogicError("OBJSENSE is not yet handled");
            }
            else if( section == MPS_COLUMNS )
            {
                colSectionBeg = file_.tellg();
            }
            else if( section == MPS_END )
            {
                break;
            }
            else if( section == MPS_MARKER )
            {
                LogicError("MPS 'MARKER' section is not yet supported");
            }
            else if( section == MPS_SOS )
            {
                LogicError("MPS 'SOS' section is not yet supported");
            }

            continue;
        }

        // No section marker was found, so handle this data line.
        if( section == MPS_ROWS )
        {
            ProcessRowDataLine( line );
        }
        else if( section == MPS_COLUMNS )
        {
            ProcessColumnDataLine( line );
        }
        else if( section == MPS_RHS )
        {
            ProcessRHSDataLine( line );
        }
        else if( section == MPS_BOUNDS )
        {
            ProcessBoundsDataLine( line );
        }
        else if( section == MPS_RANGES )
        {
            // The RANGES section should have already been processed within
            // PeekAtRanges.
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

    CountBoundTypes();

    // Extract the number of variables
    // (the matrix 'A' is 'm x n' and 'G' is 'k x n').
    // Note that the variables from variableDict_ and rangeList_ are
    // potentially interwoven and are in the three sets of size:
    //
    //   variableDict_.size() - meta_.numFixedBounds,
    //   rangeList_.size(),
    //   meta_.numFixedBounds.
    //
    meta_.n = variableDict_.size() + rangeList_.size();

    // We now force the fixed variables to come last.
    Int numFixedSoFar = 0;
    Int numNotFixedSoFar = 0;
    for( auto& entry : variableDict_ )
    {
        auto& data = entry.second;
        if( data.fixed )
        {
            data.index = (meta_.n - meta_.numFixedBounds) + numFixedSoFar++;
        }
        else
        {
            data.index = numNotFixedSoFar++;
        }
    }

    // Now that the initial pass over the file is done, we can set up for the
    // 'QueuedEntry'/'GetEntry' cycle.

    // Reset the file.
    file_.seekg( 0, std::ios::beg );
    variableIter_ = variableDict_.cbegin();

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

    string line, token;
    while( queuedEntries_.size() == 0 &&
           section_ != MPS_END &&
           std::getline( file_, line ) )
    {
        std::stringstream sectionStream( line );
        const bool isDataLine = IsDataLine(sectionStream);

        // The first token on each line should be a string. We will check it for
        // equivalence with each section string.
        if( !(sectionStream >> token) )
        {
            // This line only consists of whitespace.
            continue;
        }

        if( !isDataLine )
        {
            // Handle any special section entrance procedures.
            if( section_ == MPS_RANGES )
            {
                FillImplicitRangeBounds();
            }

            // Determine the new section.
            section_ = DecodeSection( token );

            // Handle any special section entrance procedures.
            if( section_ == MPS_END )
                break;

            continue;
        }

        // No section marker was found, so handle this data line.
        if( section_ == MPS_ROWS )
        {
            // We already have the row names.
        }
        else if( section_ == MPS_COLUMNS )
        {
            QueueColumnDataLine( line );
        }
        else if( section_ == MPS_RHS )
        {
            QueueRHSDataLine( line );
        }
        else if( section_ == MPS_BOUNDS )
        {
            // All bounds are handled at the bottom of this routine.
        }
        else if( section_ == MPS_RANGES )
        {
            continue;
        }
        else
        {
            LogicError("Invalid MPS file");
        }
    }
    if( queuedEntries_.size() > 0 )
        return true;
    QueueVariableBound();
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

const LPMPSMeta& MPSReader::Meta() const
{
    EL_DEBUG_CSE
    return meta_;
}

namespace read_mps {

template<typename Real>
LPMPSMeta Helper
( AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const string& filename,
  bool compressed,
  bool minimize,
  bool keepNonnegativeWithZeroUpperBound )
{
    EL_DEBUG_CSE
    if( compressed )
        LogicError("Compressed MPS is not yet supported");

    MPSReader reader
      ( filename, compressed, minimize, keepNonnegativeWithZeroUpperBound );
    const LPMPSMeta& meta = reader.Meta();

    Zeros( problem.c, meta.n, 1 );
    Zeros( problem.A, meta.m, meta.n );
    Zeros( problem.b, meta.m, 1 );
    Zeros( problem.G, meta.k, meta.n );
    Zeros( problem.h, meta.k, 1 );
    while( reader.QueuedEntry() )
    {
        const AffineLPEntry<double> entry = reader.GetEntry();
        if( entry.type == AFFINE_LP_COST_VECTOR )
        {
            if( problem.c(entry.row) != Real(0) )
                LogicError
                ("c(",entry.row,") was already ",problem.c(entry.row));
            problem.c(entry.row) = entry.value;
        }
        else if( entry.type == AFFINE_LP_EQUALITY_MATRIX )
        {
            if( problem.A(entry.row,entry.column) != Real(0) )
                LogicError
                ("A(",entry.row,",",entry.column,") was already ",
                 problem.A(entry.row,entry.column));
            problem.A(entry.row,entry.column) = entry.value;
        }
        else if( entry.type == AFFINE_LP_EQUALITY_VECTOR )
        {
            if( problem.b(entry.row) != Real(0) )
                LogicError
                ("b(",entry.row,") was already ",problem.b(entry.row));
            problem.b(entry.row) = entry.value;
        }
        else if( entry.type == AFFINE_LP_INEQUALITY_MATRIX )
        {
            if( problem.G(entry.row,entry.column) != Real(0) )
                LogicError
                ("G(",entry.row,",",entry.column,") was already ",
                 problem.G(entry.row,entry.column));
            problem.G(entry.row,entry.column) = entry.value;
        }
        else /* entry.type == AFFINE_LP_INEQUALITY_VECTOR */
        {
            if( problem.h(entry.row) != Real(0) )
                LogicError
                ("h(",entry.row,") was already ",problem.h(entry.row));
            problem.h(entry.row) = entry.value;
        }
    }
    return meta;
}

template<typename Real>
LPMPSMeta Helper
( AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
  const string& filename,
  bool compressed,
  bool minimize,
  bool keepNonnegativeWithZeroUpperBound )
{
    EL_DEBUG_CSE

    // TODO(poulson): Consider loading on a single process and then distributing
    // the results instead.

    if( compressed )
        LogicError("Compressed MPS is not yet supported");

    MPSReader reader
      ( filename, compressed, minimize, keepNonnegativeWithZeroUpperBound );
    const LPMPSMeta& meta = reader.Meta();

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

    return meta;
}

template<typename Real>
LPMPSMeta Helper
( AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const string& filename,
  bool compressed,
  bool minimize,
  bool keepNonnegativeWithZeroUpperBound )
{
    EL_DEBUG_CSE
    if( compressed )
        LogicError("Compressed MPS is not yet supported");

    MPSReader reader
      ( filename, compressed, minimize, keepNonnegativeWithZeroUpperBound );
    const LPMPSMeta& meta = reader.Meta();

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

    return meta;
}

template<typename Real>
LPMPSMeta Helper
( AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
  const string& filename,
  bool compressed,
  bool minimize,
  bool keepNonnegativeWithZeroUpperBound )
{
    EL_DEBUG_CSE
    if( compressed )
        LogicError("Compressed MPS is not yet supported");

    MPSReader reader
      ( filename, compressed, minimize, keepNonnegativeWithZeroUpperBound );
    const LPMPSMeta& meta = reader.Meta();

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
            problem.A.QueueUpdate
            ( entry.row, entry.column, entry.value, passive );
        else if( entry.type == AFFINE_LP_EQUALITY_VECTOR )
            problem.b.Set( entry.row, 0, entry.value );
        else if( entry.type == AFFINE_LP_INEQUALITY_MATRIX )
            problem.G.QueueUpdate
            ( entry.row, entry.column, entry.value, passive );
        else /* entry.type == AFFINE_LP_INEQUALITY_VECTOR */
            problem.h.Set( entry.row, 0, entry.value );
    }
    problem.A.ProcessLocalQueues();
    problem.G.ProcessLocalQueues();

    return meta;
}

} // namespace read_mps

template<class MatrixType,class VectorType>
LPMPSMeta ReadMPS
( AffineLPProblem<MatrixType,VectorType>& problem,
  const string& filename,
  bool compressed,
  bool minimize,
  bool keepNonnegativeWithZeroUpperBound )
{
    EL_DEBUG_CSE
    return read_mps::Helper
    ( problem, filename, compressed,
      minimize, keepNonnegativeWithZeroUpperBound );
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
