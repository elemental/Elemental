/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

// Please see http://lpsolve.sourceforge.net/5.5/mps-format.htm for a very
// nuanced discussion of the MPS file format.
//
// An important, but seemingly not widely discussed, issue is that some of the
// lp_data LP examples (e.g., tuff.mps) are not well-formed. I found out the
// hard way that the fourth equality constraint of tuff.mps is empty.

// TODO(poulson): Allow the default lower and upper bounds to be configurable.

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
  Int numNonzeros=0; // We will delete rows with no nonzeros (e.g., for tuff).
};

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
  Int numFixedBounds=0; // The bottom rows of 'A' and 'x'
  Int numFreeBounds=0;
  Int numNonpositiveBounds=0;
  Int numNonnegativeBounds=0;

  // From the RANGES section
  // TODO(poulson)

  // From the RHS section
  string rhsName="";
  Int numRHS=0;

  // Pre-solve.

  // A map from the original MPS_EQUALITY_ROW row name to the trivial equality
  // struct, which stores the variable name and the floating-point value
  // of the single nonzero in the row (eventually).
  std::map<string,MPSTrivialEquality> trivialEqualityDict;

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
      Output("  m=",m,", n=",n,", k=",k);
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
      bool keepNonnegativeWithZeroUpperBound=true );
    // The PILOT netlib lp_data model appears to require
    // 'keepNonnegativeWithZeroUpperBound=true'.

    // Attempt to enqueue another entry and return true if successful.
    bool QueuedEntry();

    // Only one call is allowed per call to 'QueuedEntry'.
    AffineLPEntry<double> GetEntry();

    const MPSMeta& Meta() const;

private:
    string filename_;
    std::ifstream file_;

    bool minimize_;
    bool keepNonnegativeWithZeroUpperBound_;
    MPSMeta meta_;

    vector<AffineLPEntry<double>> queuedEntries_;

    MPSSection section_=MPS_NONE;

    // The RHS section typically has a name followed by either one or two pairs
    // per row, but some models (e.g., dfl001.mps) do not involve a name.
    bool initializedRHSSection_=false;
    bool rhsHasName_;

    // The BOUNDS section typically has a bound type marker, followed by a
    // bound set name, followed by a variable name, and, if applicable,
    // a numeric value). But some models (e.g., dfl001.mps) do not involve a
    // bound set name.
    bool initializedBoundsSection_=false;
    bool boundsHasName_;

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

    // Rather than assuming that std::map<string,Int>::size() is constant-time,
    // we can maintain counters for the sizes of the variable dictionary.
    Int variableCounter=0;

    // Temporaries for the metadata extraction process.
    string rhsNameCandidate, boundNameCandidate;

    // For quickly reparsing the column data section if there were any trivial
    // equality rows that required retrieving the associated single nonzero
    // column name.
    std::ifstream::pos_type colSectionBeg;

    // TODO(poulson): Convert each token to upper-case letters before each
    // comparison. While capital letters are used by convention, they are
    // not required.
    MPSSection section = MPS_NONE;
    while( std::getline( file_, line_ ) )
    {
        // We first determine which section we are in
        // ------------------------------------------
        std::stringstream sectionStream( line_ );
        const char firstChar = sectionStream.peek();
        const bool isDataLine =
          firstChar == ' ' ||
          firstChar == '\t' ||
          firstChar == '*' ||
          firstChar == '#';

        // The first token on each line should be a string. We will check it for
        // equivalence with each section string.
        if( !(sectionStream >> token_) )
        {
            // This line only consists of whitespace.
            continue;
        }

        if( !isDataLine )
        {
            if( section == MPS_COLUMNS )
            {
                // We are finished reading in the column data, so iterate
                // through the row map to delete any empty rows and convert any
                // equality rows with a single nonzero to a fixed state.
                for( auto iter=meta_.rowDict.begin();
                  iter!=meta_.rowDict.end(); )
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
                                Output
                                ("WARNING: Deleting empty equality row ",
                                 iter->first);
                            else if( iter->second.type == MPS_GREATER_ROW )
                                Output
                                ("WARNING: Deleting empty greater row ",
                                 iter->first);
                            else if( iter->second.type == MPS_LESSER_ROW )
                                Output
                                ("WARNING: Deleting empty lesser row ",
                                 iter->first);
                            else
                                LogicError("Unknown empty row type");
                            // Delete this entry.
                            iter = meta_.rowDict.erase(iter);
                        }
                    }
                    else if( iter->second.numNonzeros == 1 &&
                             iter->second.type == MPS_EQUALITY_ROW )
                    {
                        // Convert to a fixed value. We will divide the
                        // corresponding RHS value by the single nonzero.
                        MPSTrivialEquality trivialEquality;
                        meta_.trivialEqualityDict[iter->first] =
                          trivialEquality;
                        iter = meta_.rowDict.erase(iter);
                        // Subtract one from the equality entries
                        // (which will be added back later when iterating over
                        // all of the fixed bounds).
                        --meta_.numEqualityEntries;
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
                        ++iter;
                    }
                }

                if( meta_.trivialEqualityDict.size() != 0 )
                {
                    // Loop through the COLUMNS section again to store the
                    // variable name associated with each trivial equality row.
                    file_.seekg( colSectionBeg );
                    string columnLine, columnToken;
                    while( std::getline( file_, columnLine ) )
                    {
                        std::stringstream columnStream( columnLine );
                        const char firstColumnChar = columnStream.peek();
                        const bool isColumnLine =
                          firstColumnChar == ' ' ||
                          firstColumnChar == '\t' ||
                          firstColumnChar == '*' ||
                          firstColumnChar == '#';
                        if( !isColumnLine )
                            break;
                        if( !(columnStream >> variableName_) )
                            continue;
                        auto variableIter =
                            meta_.variableDict.find( variableName_ );
                        if( variableIter == meta_.variableDict.end() )
                            LogicError("Invalid 'COLUMNS' section");
                        MPSVariableData& variableData = variableIter->second;
                        for( Int pair=0; pair<2; ++pair )
                        {
                            if( !(columnStream >> rowName_) )
                            {
                                if( pair == 0 )
                                    LogicError("Invalid 'COLUMNS' section");
                                else
                                    break;
                            }
                            if( !(columnStream >> value_) )
                                LogicError("Invalid 'COLUMNS' section");
                            auto trivialEqualityIter =
                              meta_.trivialEqualityDict.find( rowName_ );
                            if( trivialEqualityIter ==
                                meta_.trivialEqualityDict.end() )
                                continue;
                            if( value_ == 0. )
                                LogicError
                                ("Trivial equality rows with a zero \"nonzero\""
                                 " value are not yet supported");
                            auto& trivialData = trivialEqualityIter->second;
                            trivialData.variableName = variableName_;
                            trivialData.singleNonzero = value_;
                            Output
                            ("WARNING: Storing nonzero of ",value_,
                             " for trivial equality row ",rowName_);
                            variableData.fixedIndex = meta_.numFixedBounds++;
                            // We initialize at zero and overwrite if there is
                            // relevant RHS data.
                            variableData.fixedValue = 0;
                        }
                    }
                }
            }

            if( token_ == "NAME" )
            {
                if( meta_.name != "" )
                    LogicError("Multiple 'NAME' sections");
                if( !(sectionStream >> meta_.name) )
                    LogicError("Missing 'NAME' string");
                section = MPS_NAME;
            }
            else if( token_ == "OBJSENSE" )
            {
                LogicError("OBJSENSE is not yet handled");
            }
            else if( token_ == "ROWS" )
            {
                section = MPS_ROWS;
            }
            else if( token_ == "COLUMNS" )
            {
                section = MPS_COLUMNS;
                colSectionBeg = file_.tellg();
            }
            else if( token_ == "RHS" )
            {
                section = MPS_RHS;
            }
            else if( token_ == "BOUNDS" )
            {
                section = MPS_BOUNDS;
            }
            else if( token_ == "RANGES" )
            {
                section = MPS_RANGES;
                LogicError("MPS 'RANGES' section is not yet supported");
            }
            else if( token_ == "ENDATA" )
            {
                section = MPS_END;
                break;
            }
            else if( token_ == "MARKER" )
            {
                LogicError("MPS 'MARKER' section is not yet supported");
            }
            else if( token_ == "SOS" )
            {
                LogicError("MPS 'SOS' section is not yet supported");
            }
            else
            {
                LogicError("Section token ",token_," is not recognized");
            }
            continue;
        }

        // No section marker was found, so handle this data line.
        if( section == MPS_ROWS )
        {
            std::stringstream rowStream( line_ );
            if( !(rowStream >> rowType_) )
                LogicError("Invalid 'ROWS' section");
            if( !(rowStream >> rowName_) )
                LogicError("Invalid 'ROWS' section");
            MPSRowData rowData;
            // We set the 'typeIndex' fields later since it is not uncommon
            // (e.g., see tuff.mps) for rows to be empty.
            if( rowType_ == "L" )
            {
                rowData.type = MPS_LESSER_ROW;
                meta_.rowDict[rowName_] = rowData;
            }
            else if( rowType_ == "G" )
            {
                rowData.type = MPS_GREATER_ROW;
                meta_.rowDict[rowName_] = rowData;
            }
            else if( rowType_ == "E" )
            {
                rowData.type = MPS_EQUALITY_ROW;
                meta_.rowDict[rowName_] = rowData;
            }
            else if( rowType_ == "N" )
            {
                rowData.type = MPS_NONCONSTRAINING_ROW;
                meta_.rowDict[rowName_] = rowData;
                if( meta_.numNonconstrainingRows == 1 )
                    meta_.costName = rowName_;
            }
            else
                LogicError("Invalid 'ROWS' section");
        }
        else if( section == MPS_COLUMNS )
        {
            std::stringstream columnStream( line_ );
            if( !(columnStream >> variableName_) )
                LogicError("Invalid 'COLUMNS' section");
            auto variableIter = meta_.variableDict.find( variableName_ );
            if( variableIter == meta_.variableDict.end() )
            {
                MPSVariableData variableData;
                variableData.index = variableCounter++;
                meta_.variableDict[variableName_] = variableData;
                variableIter = meta_.variableDict.find( variableName_ );
            }
            MPSVariableData& variableData = variableIter->second;

            // There should be either one or two pairs of entries left to read
            // from this line.
            for( Int pair=0; pair<2; ++pair )
            {
                if( !(columnStream >> rowName_) )
                {
                    if( pair == 0 )
                        LogicError("Invalid 'COLUMNS' section");
                    else
                        break;
                }
                if( !(columnStream >> value_) )
                    LogicError("Invalid 'COLUMNS' section");
                ++variableData.numNonzeros;

                auto rowIter = meta_.rowDict.find( rowName_ );
                if( rowIter == meta_.rowDict.end() )
                    LogicError("Could not find row ",rowName_);
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
        else if( section == MPS_RHS )
        {
            if( !initializedRHSSection_ )
            {
                std::stringstream rhsTestStream( line_ );
                Int numTokens=0;
                string rhsToken;
                while( rhsTestStream >> rhsToken )
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

            std::stringstream rhsStream( line_ );
            if( rhsHasName_ )
            {
                if( !(rhsStream >> rhsNameCandidate) )
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
                if( !(rhsStream >> rowName_) )
                {
                    if( pair == 0 )
                        LogicError("Invalid 'RHS' section (3)");
                    else
                        break;
                }
                if( !(rhsStream >> value_) )
                    LogicError("Invalid 'RHS' section (4)");
                // TODO(poulson): Store fixed value.
            }
        }
        else if( section == MPS_BOUNDS )
        {
            if( !initializedBoundsSection_ )
            {
                std::stringstream boundsTestStream( line_ );
                if( !(boundsTestStream >> boundMark_) )
                    LogicError("Invalid 'BOUNDS' section");
                if( !(boundsTestStream >> token_) )
                    LogicError("Invalid 'BOUNDS' section");

                Int numTokens=2;
                string rhsToken;
                while( boundsTestStream >> rhsToken )
                    ++numTokens;
                if( boundMark_ == "LO" ||
                    boundMark_ == "UP" ||
                    boundMark_ == "FX" )
                {
                    if( numTokens == 4 )
                        boundsHasName_ = true;
                    else if( numTokens == 3 )
                        boundsHasName_ = false;
                    else
                        LogicError("Invalid ",boundMark_," 'BOUNDS' line");
                }
                else if( boundMark_ == "FR" ||
                         boundMark_ == "MI" ||
                         boundMark_ == "PL" )
                {
                    if( numTokens == 3 )
                        boundsHasName_ = true;
                    else if( numTokens == 2 )
                        boundsHasName_ = false;
                    else
                        LogicError("Invalid ",boundMark_," 'BOUNDS' line");
                }
                else
                    LogicError("Unknown bound mark ",boundMark_);

                // If there is a name, it occurred in the second position.
                if( boundsHasName_ )
                    meta_.boundName = token_;

                initializedBoundsSection_ = true;
            }

            std::stringstream boundStream( line_ );
            // We already have the first token of a bounding row, which should
            // be of the same general form as
            //
            //   FX BOUNDROW VARIABLENAME 1734.
            //
            // in the case of 'VARIABLENAME' being fixed ('FX') at the value
            // 1734 (with this problem's bound name being 'BOUNDROW').
            if( !(boundStream >> boundMark_) )
                LogicError("Invalid 'BOUNDS' section");
            if( boundsHasName_ )
            {
                if( !(boundStream >> boundNameCandidate) )
                    LogicError("Invalid 'BOUNDS' section");
                if( meta_.boundName != boundNameCandidate )
                    LogicError
                    ("Only single problem instances are currently supported "
                     "(multiple bound names were encountered)");
            }
            if( !(boundStream >> variableName_) )
                LogicError("Invalid 'BOUNDS' section");
            auto variableIter = meta_.variableDict.find( variableName_ );
            if( variableIter == meta_.variableDict.end() )
                LogicError
                ("Invalid 'BOUNDS' section (name ",variableName_," not found)");
            MPSVariableData& data = variableIter->second;
            if( boundMark_ == "UP" )
            {
                if( !(boundStream >> value_) )
                    LogicError("Invalid 'BOUNDS' section");
                data.upperBounded = true;
                data.upperBound = value_;
            }
            else if( boundMark_ == "LO" )
            {
                if( !(boundStream >> value_) )
                    LogicError("Invalid 'BOUNDS' section");
                data.lowerBounded = true;
                data.lowerBound = value_;
            }
            else if( boundMark_ == "FX" )
            {
                if( !(boundStream >> value_) )
                    LogicError("Invalid 'BOUNDS' section");
                data.fixed = true;
                data.fixedValue = value_;
            }
            else if( boundMark_ == "FR" )
                data.free = true;
            else if( boundMark_ == "MI" )
                data.nonpositive = true;
            else if( boundMark_ == "PL" )
                data.nonnegative = true;
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
        if( data.upperBounded )
        {
            if( data.lowerBounded )
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
            data.lowerBoundIndex = meta_.numLowerBounds++;
            ++meta_.numInequalityEntries;
        }

        // Handle fixed values.
        if( data.fixed )
        {
            if( data.upperBounded || data.lowerBounded || data.free ||
                data.nonpositive || data.nonnegative )
                LogicError("Invalid bound combination");
            data.fixedIndex = meta_.numFixedBounds++;
            ++meta_.numEqualityEntries;
        }

        // Handle free values.
        if( data.free )
        {
            if( data.upperBounded || data.lowerBounded ||
                data.nonnegative || data.nonpositive ||
                data.fixed )
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
        std::stringstream sectionStream( line_ );
        const char firstChar = sectionStream.peek();
        const bool isDataLine =
          firstChar == ' ' ||
          firstChar == '\t' ||
          firstChar == '*' ||
          firstChar == '#';

        // The first token on each line should be a string. We will check it for
        // equivalence with each section string.
        if( !(sectionStream >> token_) )
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
        if( section_ == MPS_ROWS )
        {
            // We already have the row names.
        }
        else if( section_ == MPS_COLUMNS )
        {
           double columnValue;
           std::stringstream columnStream( line_ );

           if( !(columnStream >> variableName_) )
                LogicError("Invalid 'COLUMNS' section");
            auto variableIter = meta_.variableDict.find( variableName_ );
            if( variableIter == meta_.variableDict.end() )
                LogicError("Could not find variable ",variableName_);
            const MPSVariableData& variableData = variableIter->second;
            const Int column = variableData.index;

            // There should be either one or two pairs of entries left to read
            // from this line.
            for( Int pair=0; pair<2; ++pair )
            {
                if( !(columnStream >> rowName_) )
                {
                    if( pair == 0 )
                        LogicError("Invalid 'COLUMNS' section");
                    else
                        break;
                }
                if( !(columnStream >> columnValue) )
                    LogicError("Invalid 'COLUMNS' section");

                AffineLPEntry<double> entry;
                auto rowIter = meta_.rowDict.find( rowName_ );
                if( rowIter == meta_.rowDict.end() )
                {
                    auto trivialEqualityIter =
                      meta_.trivialEqualityDict.find( rowName_ );
                    if( trivialEqualityIter == meta_.trivialEqualityDict.end() )
                    {
                        LogicError("Could not find row ",rowName_);
                    }
                    if( columnValue == 0. )
                    {
                        LogicError
                        ("Trivial rows are not yet allowed to have their only "
                         "specified \"nonzero\" value be zero");
                    }
                    Output
                    ("WARNING: Trivial equality row ",rowName_,
                     " only had a nonzero in column ",variableName_);
                    auto& trivialData = trivialEqualityIter->second;
                    trivialData.variableName = variableName_; 
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
        else if( section_ == MPS_RHS )
        {
            double rhsValue;
            std::stringstream rhsStream( line_ );

            if( rhsHasName_ )
            {
                if( !(rhsStream >> token_) )
                    LogicError("Invalid 'RHS' section");
            }

            // There should be either one or two pairs of entries left to read
            // from this line.
            for( Int pair=0; pair<2; ++pair )
            {
                if( !(rhsStream >> rowName_) )
                {
                    if( pair == 0 )
                        LogicError("Invalid 'RHS' section");
                    else
                        break;
                }
                if( !(rhsStream >> rhsValue) )
                    LogicError("Invalid 'RHS' section");
                auto rowIter = meta_.rowDict.find( rowName_ );
                if( rowIter == meta_.rowDict.end() )
                {
                    auto trivialEqualityIter =
                      meta_.trivialEqualityDict.find( rowName_ );
                    if( trivialEqualityIter == meta_.trivialEqualityDict.end() )
                        LogicError("Could not find trivial row ",rowName_);
                    const auto& trivialData = trivialEqualityIter->second;
                    auto variableIter =
                      meta_.variableDict.find( trivialData.variableName );
                    if( variableIter == meta_.variableDict.end() )
                    {
                        LogicError
                        ("Could not find variable ",trivialData.variableName,
                         " for trivial equality from row ",rowName_);
                    }
                    variableIter->second.fixedValue =
                      rhsValue / trivialData.singleNonzero;
                    continue;
                }
                AffineLPEntry<double> entry;
                const auto& rowData = rowIter->second;
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
        else if( section_ == MPS_BOUNDS )
        {
            // All bounds are handled at the bottom of this routine.
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

    // Now iterate through the variable map to handle the bounds.
    AffineLPEntry<double> entry;
    while( queuedEntries_.size() == 0 &&
           variableIter_ != meta_.variableDict.end() )
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
  bool minimize,
  bool keepNonnegativeWithZeroUpperBound,
  bool metadataSummary )
{
    EL_DEBUG_CSE
    if( compressed )
        LogicError("Compressed MPS is not yet supported");

    MPSReader reader
      ( filename, compressed, minimize, keepNonnegativeWithZeroUpperBound );
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
}

template<typename Real>
void Helper
( AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
  const string& filename,
  bool compressed,
  bool minimize,
  bool keepNonnegativeWithZeroUpperBound,
  bool metadataSummary )
{
    EL_DEBUG_CSE

    // TODO(poulson): Consider loading on a single process and then distributing
    // the results instead.

    if( compressed )
        LogicError("Compressed MPS is not yet supported");

    MPSReader reader
      ( filename, compressed, minimize, keepNonnegativeWithZeroUpperBound );
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
  bool minimize,
  bool keepNonnegativeWithZeroUpperBound,
  bool metadataSummary )
{
    EL_DEBUG_CSE
    if( compressed )
        LogicError("Compressed MPS is not yet supported");

    MPSReader reader
      ( filename, compressed, minimize, keepNonnegativeWithZeroUpperBound );
    const MPSMeta& meta = reader.Meta();
    if( metadataSummary )
        meta.PrintSummary();

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
  bool minimize,
  bool keepNonnegativeWithZeroUpperBound,
  bool metadataSummary )
{
    EL_DEBUG_CSE
    if( compressed )
        LogicError("Compressed MPS is not yet supported");

    MPSReader reader
      ( filename, compressed, minimize, keepNonnegativeWithZeroUpperBound );
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
}

} // namespace read_mps

template<class MatrixType,class VectorType>
void ReadMPS
( AffineLPProblem<MatrixType,VectorType>& problem,
  const string& filename,
  bool compressed,
  bool minimize,
  bool keepNonnegativeWithZeroUpperBound,
  bool metadataSummary )
{
    EL_DEBUG_CSE
    read_mps::Helper
    ( problem, filename, compressed,
      minimize, keepNonnegativeWithZeroUpperBound, metadataSummary );
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
