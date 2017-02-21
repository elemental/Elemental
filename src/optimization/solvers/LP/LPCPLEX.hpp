/*
   Copyright (c) 2009-2017, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

// Please see http://lpsolve.sourceforge.net/5.5/CPLEX-format.htm for a
// reference.

// TODO(poulson): Add support for constants in the objective?

enum LPCPLEXSection {
  LPCPLEX_MAXIMIZE,
  LPCPLEX_MINIMIZE,
  LPCPLEX_SUBJECT_TO,
  LPCPLEX_BOUNDS,
  LPCPLEX_GENERAL,
  LPCPLEX_END,
  LPCPLEX_NONE
};

enum LPCPLEXRowType {
  LPCPLEX_LESSER_ROW,
  LPCPLEX_GREATER_ROW,
  LPCPLEX_EQUALITY_ROW
};

namespace lp_cplex {

static const string missingEndLabel = "There was no 'END' label.";

// Get the first character or return the failure condition of 'false'.
bool ExtractFirstCharacter( std::istream& stream, char& firstChar )
{
    EL_DEBUG_CSE
    while( true )
    {
        char firstCharCand = stream.peek();
        if( firstCharCand == std::char_traits<char>::eof() )
            return false;
        if( firstCharCand == '\\' )
        {
            char ignoredChar = ' ';
            while( ignoredChar != '\n' ) 
            {
                ignoredChar = stream.get();
                // TODO(poulson): Handle the last line being a comment and there
                // not being an 'END' label.
            }
            continue;
        }
        if( firstCharCand == ' ' ||
            firstCharCand == '\n' ||
            firstCharCand == '\t' )
        {
            stream.get();
        }
        else
        {
            firstChar = stream.get();
            break;
        }
    }
    return true;
}

bool IsNumberCharacter( const char& character )
{
    EL_DEBUG_CSE
    return (character >= '0' && character <= '9') ||
           character == 'e' || character == 'E' ||
           character == '.' ||
           character == '-' || character == '+';
}

bool IsInnerWordCharacter( const char& character )
{
    EL_DEBUG_CSE
    return (character >= 'a' && character <= 'z') ||
           (character >= 'A' && character <= 'Z') ||
           (character >= '0' && character <= '9') ||
           character == '!' ||
           character == '"' ||
           character == '#' ||
           character == '$' ||
           character == '%' ||
           character == '&' ||
           character == '(' ||
           character == ')' ||
           character == '/' ||
           character == '.' ||
           character == ';' ||
           character == '?' ||
           character == '@' ||
           character == '_' ||
           character == '`' ||
           character == '\'' ||
           character == '{' ||
           character == '}' ||
           character == '|' ||
           character == '~';
}

bool IsRelationToken( const string& token )
{
    EL_DEBUG_CSE
    return token == "<" ||
           token == "<=" ||
           token == "=" ||
           token == ">=" || 
           token == ">";
}

LPCPLEXRowType DecodeRowType( const string& token )
{
    EL_DEBUG_CSE
    if( !IsRelationToken( token ) )
        RuntimeError(token," is not a relation token");
    if( token == "<" || token == "<=" )
        return LPCPLEX_LESSER_ROW;
    else if( token == "=" )
        return LPCPLEX_EQUALITY_ROW;
    else
        return LPCPLEX_GREATER_ROW;
}

bool IsPositiveInfinity( const string& token )
{
    EL_DEBUG_CSE
    const string upperToken = ToUpper( token );
    return upperToken == "+INF" || upperToken == "+INFINITY";
}

bool IsNegativeInfinity( const string& token )
{
    EL_DEBUG_CSE
    const string upperToken = ToUpper( token );
    return upperToken == "-INF" || upperToken == "-INFINITY";
}

bool IsSignedInfinity( const string& token )
{
    EL_DEBUG_CSE
    return IsPositiveInfinity(token) || IsNegativeInfinity(token);
}

template<typename Real>
Real StringToReal( const string& token )
{
    EL_DEBUG_CSE
    Real number;
    std::stringstream stream( token );
    stream >> number;
    return number;
}

bool ExtractToken( std::istream& stream, string& token )
{
    EL_DEBUG_CSE
    token = "";

    char firstChar;
    bool gotFirstChar = ExtractFirstCharacter( stream, firstChar );
    if( !gotFirstChar )
        return false;
    token = firstChar;

    if( firstChar == 'e' || firstChar == 'E' )
    {
        auto pos = stream.tellg();
        const char char0 = stream.get();
        const char char1 = stream.get();
        if( ::toupper(char0) == 'N' && ::toupper(char1) == 'D' )
        {
            token = "END";
            return true;
        }
        else
            stream.seekg(pos);
    }

    if( firstChar == ':' )
    {
        while( stream.peek() == ':' )
            token += stream.get();
    }
    else if( firstChar == '>' )
    {
        if( stream.peek() == '=' )
            token += stream.get();
    }
    else if( firstChar == '<' )
    {
        if( stream.peek() == '=' )
            token += stream.get();
    }
    else if( IsNumberCharacter( firstChar ) )
    {
        if( firstChar == '-' || firstChar == '+' )
        {
            // See if we can extend to +-inf or +-infinty
            auto pos = stream.tellg();
            const char char0 = stream.get();
            const char char1 = stream.get();
            const char char2 = stream.get();
            if( ::toupper(char0) == 'I' &&
                ::toupper(char1) == 'N' &&
                ::toupper(char2) == 'F' )
            {
                token += "inf";
                // See if we can further extend to -infinity
                auto secondPos = stream.tellg();
                const char char3 = stream.get();
                const char char4 = stream.get();
                const char char5 = stream.get();
                const char char6 = stream.get();
                const char char7 = stream.get();
                if( ::toupper(char3) == 'I' &&
                    ::toupper(char4) == 'N' &&
                    ::toupper(char5) == 'I' &&
                    ::toupper(char6) == 'T' &&
                    ::toupper(char7) == 'Y' )
                {
                    token += "inity";
                }
                else
                {
                    stream.seekg( secondPos );
                }
            }
            else
            {
                stream.seekg( pos );
            }
        }
        else
        {
            while( IsNumberCharacter( stream.peek() ) )
                token += stream.get();
        }
    }
    else
    {
        while( IsInnerWordCharacter( stream.peek() ) )
            token += stream.get();
        string upperToken = ToUpper( token );
        if( upperToken == "SUBJECT" )
        {
            // See if we can extend to "SUBJECT TO".
            auto pos = stream.tellg();
            string nextToken;
            const bool gotNextToken = ExtractToken( stream, nextToken );
            if( gotNextToken && ToUpper(nextToken) == "TO" )
            {
                token += " " + nextToken;
            }
            else
            {
                stream.seekg( pos );
            }
        }
        else if( upperToken == "SUCH" )
        {
            // See if we can extend to "SUCH THAT".
            auto pos = stream.tellg();
            string nextToken;
            const bool gotNextToken = ExtractToken( stream, nextToken );
            if( gotNextToken && ToUpper(nextToken) == "THAT" )
            {
                token += " " + nextToken;
            }
            else
            {
                stream.seekg( pos );
            }
        }
        else if( upperToken == "ST" )
        {
            // See if we can extend to "ST.".
            if( stream.peek() == '.' )
                token += stream.get();
        }
        else if( upperToken == "S" )
        {
            // See if we can extend to s.t.
            auto pos = stream.tellg();
            const char char0 = stream.get();
            const char char1 = stream.get();
            const char char2 = stream.get();
            if( char0 == '.' && ::toupper(char1) == 'T' && char2 == '.' )
            {
                token += char0;
                token += char1;
                token += char2;
            }
            else
            {
                stream.seekg( pos );
            }
        }
    }

    return true;
}

LPCPLEXSection DecodeSection( const string& token )
{
    EL_DEBUG_CSE
    LPCPLEXSection section = LPCPLEX_NONE;
    const string upperToken = ToUpper( token );
    if( upperToken == "MAXIMIZE" )
        section = LPCPLEX_MAXIMIZE; 
    else if( upperToken == "MINIMIZE" )
        section = LPCPLEX_MINIMIZE;
    else if( upperToken == "SUBJECT TO" ||
             upperToken == "SUCH THAT" ||
             upperToken == "ST" ||
             upperToken == "ST." ||
             upperToken == "S.T." )
        section = LPCPLEX_SUBJECT_TO;
    else if( upperToken == "BOUND" )
        section = LPCPLEX_BOUNDS;
    else if( upperToken == "GENERAL" || upperToken == "GENERALS" )
        section = LPCPLEX_GENERAL;
    else if( upperToken == "END" )
        section = LPCPLEX_END;
    return section;
}

} // namespace lp_cplex

struct LPCPLEXRowData
{
  string name="";
  LPCPLEXRowType type;
  Int typeIndex;
  Int numNonzeros=0; // We will delete rows with no nonzeros (e.g., for tuff).
};

template<typename Real>
struct LPCPLEXVariableData
{
  Int index;
  Int numNonzeros=0;

  bool lowerBounded=false;
  Int lowerBoundIndex=-1;
  Real lowerBound;

  bool upperBounded=false;
  Int upperBoundIndex=-1;
  Real upperBound;

  bool fixed=false;
  Int fixedIndex=-1;
  Real fixedValue;

  bool free=false;
  Int freeIndex=-1;

  bool nonpositive=false;
  Int nonpositiveIndex=-1;

  bool nonnegative=true;
  Int nonnegativeIndex=-1;
};

void LPCPLEXMeta::PrintSummary() const
{
    Output("LPCPLEXMeta summary:");
    Output("  objectiveName=",objectiveName);
    Output("  numLesserRows=",numLesserRows);
    Output("  numGreaterRows=",numGreaterRows);
    Output("  numEqualityRows=",numEqualityRows);
    Output("  numEqualityEntries=",numEqualityEntries);
    Output("  numInequalityEntries=",numInequalityEntries);
    Output("  numUpperBounds=",numUpperBounds);
    Output("  numLowerBounds=",numLowerBounds);
    Output("  numFixedBounds=",numFixedBounds);
    Output("  numFreeBounds=",numFreeBounds);
    Output("  numNonpositiveBounds=",numNonpositiveBounds);
    Output("  numNonnegativeBounds=",numNonnegativeBounds);
    Output("  m=",m,", n=",n,", k=",k);
}

// We form the primal problem
//
//   arginf_{x,s} { c^T x | A x = b, G x + s = h, s >= 0 },
//
// where 'A x = b' consists of the 'equality' LP CPLEX rows and the 'fixed'
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
// The ordering of 'x' is: nontrivial variables, and trivially fixed variables.
//
template<typename Real>
class LPCPLEXReader
{
public:
    LPCPLEXReader( const string& filename );

    // Attempt to enqueue another entry and return true if successful.
    bool QueuedEntry();

    // Only one call is allowed per call to 'QueuedEntry'.
    AffineLPEntry<Real> GetEntry();

    const LPCPLEXMeta& Meta() const;

private:
    string filename_;
    std::ifstream file_;

    bool minimize_;
    LPCPLEXMeta meta_;

    // TODO(poulson): Expose this somehow.
    bool keepNonnegativeWithZeroUpperBound_=true;

    // TODO(poulson): Expose this somehow.
    bool progress_=false;

    vector<LPCPLEXRowData> rowData_;
    std::map<string,LPCPLEXVariableData<Real>> variableDict_;

    vector<AffineLPEntry<Real>> queuedEntries_;

    LPCPLEXSection section_=LPCPLEX_NONE;

    // For the final loop over the variable dictionary to extract any bounds.
    typename std::map<string,LPCPLEXVariableData<Real>>::const_iterator
      variableIter_;

    Int rowCounter_=0;

    // If any remain, queue an entry associated with enforcing a variable bound.
    void QueueVariableBound();

    // Iterate through the variable map and make use of the requested
    // conventions for counting the number of bounds of each type.
    // Also warn if there are possibly conflicting bound types.
    void CountBoundTypes();

    void SimplifyRowTypes();
};

template<typename Real>
void LPCPLEXReader<Real>::QueueVariableBound()
{
    EL_DEBUG_CSE
    using namespace El::lp_cplex;
    AffineLPEntry<Real> entry;
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
}

template<typename Real>
void LPCPLEXReader<Real>::SimplifyRowTypes()
{
    EL_DEBUG_CSE
    using namespace El::lp_cplex;
    // TODO(poulson): Add support for trivial equalities...
    for( auto& rowData : rowData_ )
    {
        if( rowData.type == LPCPLEX_EQUALITY_ROW )
        {
            rowData.typeIndex = meta_.numEqualityRows++;
            meta_.numEqualityEntries += rowData.numNonzeros;
        }
        else if( rowData.type == LPCPLEX_GREATER_ROW )
        {
            rowData.typeIndex = meta_.numGreaterRows++;
            meta_.numInequalityEntries += rowData.numNonzeros;
        }
        else if( rowData.type == LPCPLEX_LESSER_ROW )
        {
            rowData.typeIndex = meta_.numLesserRows++;
            meta_.numInequalityEntries += rowData.numNonzeros;
        }
        else
            LogicError("Unknown row type");
    }
}

template<typename Real>
void LPCPLEXReader<Real>::CountBoundTypes()
{
    EL_DEBUG_CSE
    using namespace El::lp_cplex;
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
                // varies between different readers.
                if( data.upperBound > Real(0) )
                {
                    // Preserve the default non-negativity assumption.
                    data.upperBoundIndex = meta_.numUpperBounds++;
                    ++meta_.numInequalityEntries;
                    data.nonnegative = true;
                }
                else if( data.upperBound == Real(0) && !data.nonnegative )
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
                        data.fixedValue = 0;
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

        // Handle non-positive values.
        if( data.nonpositive )
        {
            if( data.upperBounded )
                Output
                ("WARNING: Combined nonpositive constraint with upper bound");
            data.nonpositiveIndex = meta_.numNonpositiveBounds++;
            ++meta_.numInequalityEntries;
        }

        // Handle non-negative values.
        if( data.nonnegative )
        {
            if( data.lowerBounded )
                Output
                ("WARNING: Combined nonnegative constraint with lower bound");
            data.nonnegativeIndex = meta_.numNonnegativeBounds++;
            ++meta_.numInequalityEntries;
        }

        // Default to free if no bounds are marked and the 'nonnegative' mark
        // has been removed.
        const bool hasAMark =
          data.upperBounded ||
          data.lowerBounded ||
          data.fixed ||
          data.free ||
          data.nonpositive ||
          data.nonnegative;
        if( !hasAMark )
            data.free = true;

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
    }
}

template<typename Real>
LPCPLEXReader<Real>::LPCPLEXReader( const string& filename )
: filename_(filename), file_(filename.c_str())
{
    EL_DEBUG_CSE
    using namespace El::lp_cplex;
    if( !file_.is_open() )
        RuntimeError("Could not open ",filename);

    LPCPLEXSection section = LPCPLEX_NONE;
    string token;

    // TODO(poulson): Expose this parameter.
    const bool progress_ = false;

    // Scan the objective information
    // TODO(poulson): Support constants?
    {
        string directionToken;
        if( !ExtractToken( file_, directionToken ) )
            RuntimeError("Could not extract first token");
        section = DecodeSection( directionToken );
        if( section == LPCPLEX_MAXIMIZE )
            minimize_ = false;
        else if( section == LPCPLEX_MINIMIZE )
            minimize_ = true;
        else
            RuntimeError("First token (",directionToken,") was invalid");

        string secondToken;
        if( !ExtractToken( file_, secondToken ) )
            RuntimeError("There was not second token!");
        auto secondTokenSection = DecodeSection( secondToken );
        if( secondTokenSection != LPCPLEX_NONE )
        {
            if( progress_ )
                Output("Jumping from OBJECTIVE to ",secondToken);
            section = secondTokenSection;
        }
        else
        {
            if( file_.peek() == ':' )
            {
                meta_.objectiveName = secondToken;
                file_.get();
            }

            // Loop until the token stores the next section label.
            while( ExtractToken( file_, token ) )
            {
                auto tokenSection = DecodeSection( token );
                if( tokenSection != LPCPLEX_NONE )
                {
                    if( progress_ )
                        Output("Jumping from OBJECTIVE to ",token);
                    section = tokenSection;
                    break;
                }
            }
        }
    }

    if( section == LPCPLEX_SUBJECT_TO )
    {
        if( progress_ )
            Output("Entered Subject To section");
        // Each of the following outer loops will handle an (in)equality.
        Int unnamedRowCounter = 0;
        while( true )
        {
            // Get the first token.
            if( !ExtractToken( file_, token ) )
                RuntimeError(missingEndLabel);
            auto tokenSection = DecodeSection( token );
            if( tokenSection != LPCPLEX_NONE )
            {
                if( progress_ )
                    Output("Jumping from SUBJECT TO to ",token);
                section = tokenSection;
                break;
            }

            // Handle another (in)equality.
            LPCPLEXRowData rowData;
            if( file_.peek() == ':' )
            {
                rowData.name = token;
                file_.get();
                if( !ExtractToken( file_, token ) )
                    RuntimeError(missingEndLabel);
            }
            else
            {
                rowData.name = "R" + std::to_string(unnamedRowCounter++);
            }

            bool pastRelation = false;
            while( true )
            {
                if( !pastRelation )
                {
                    if( IsRelationToken( token ) )
                    {
                        rowData.type = DecodeRowType( token );
                        pastRelation = true;
                    }
                    else if( IsNumberCharacter(token[0]) )
                    {
                        // We would accumulate here.
                    }
                    else
                    {
                        ++rowData.numNonzeros;

                        // Look up the token in the variable dictionary.
                        auto varIter = variableDict_.find( token );
                        if( varIter == variableDict_.end() )
                        {
                            LPCPLEXVariableData<Real> varData;
                            varData.numNonzeros = 1;
                            variableDict_[token] = varData;
                        }
                        else
                        {
                            auto& varData = varIter->second;
                            ++varData.numNonzeros;
                        }

                        // We would reset the accumulation here...
                    }
                }
                else
                {
                    // We would store the right-hand side value here...
                    break;
                }
                if( !ExtractToken( file_, token ) )
                    RuntimeError(missingEndLabel);
            }
            rowData_.push_back( rowData );
        }
    }

    if( section == LPCPLEX_BOUNDS )
    {
        if( progress_ )
            Output("Entered Bounds section");
        string line;

        // The initial call to getline should not involve any nontrivial tokens
        // and should finish off the whitespace on the 'BOUNDS' line.
        if( !std::getline( file_, line ) ) 
            RuntimeError(missingEndLabel);

        // Loop over the BOUNDS lines. We will end with the line that contains
        // the next section label.
        while( true )
        {
            if( !std::getline( file_, line ) ) 
                RuntimeError(missingEndLabel);
            vector<string> lineTokens;
            {
                std::stringstream lineStream( line );
                while( ExtractToken( lineStream, token ) )
                    lineTokens.push_back( token );
            }
            auto tokenSection = DecodeSection( lineTokens[0] );
            if( tokenSection != LPCPLEX_NONE )
            {
                if( progress_ )
                    Output("Jumping from BOUNDS to ",lineTokens[0]);
                section = tokenSection;
                break;
            }
            if( lineTokens.size() == 2 )
            {
                // We should have a line of the form 'x free'.
                auto varIter = variableDict_.find( lineTokens[0] );
                if( varIter == variableDict_.end() )
                    RuntimeError("Could not find variable ",lineTokens[0]);
                if( ToUpper(lineTokens[1]) != "FREE" )
                    RuntimeError
                    ("Unrecognized bound line of (",lineTokens[0],",",
                     lineTokens[1],")");
                auto& varData = varIter->second;
                varData.nonnegative = false;
                varData.free = true;
            }
            else if( lineTokens.size() == 3 )
            {
                // We should have a three-token line of the form
                // 'lowerBound <= x' or 'x <= upperBound'.
                if( !IsRelationToken(lineTokens[1]) )
                    RuntimeError
                    ("Unrecognized bound line of (",lineTokens[0],",", 
                     lineTokens[1],",",lineTokens[2],")");
                if( IsNumberCharacter(lineTokens[0][0]) )
                {
                    if( IsSignedInfinity(lineTokens[0]) )
                        continue;
                    const Real number = StringToReal<Real>( lineTokens[0] );

                    auto varIter = variableDict_.find( lineTokens[2] );
                    if( varIter == variableDict_.end() )
                        RuntimeError("Could not find variable ",lineTokens[2]);
                    auto& varData = varIter->second;

                    if( lineTokens[1] == "<" || lineTokens[1] == "<=" )
                    {
                        // number <= x.
                        varData.nonnegative = false;
                        varData.lowerBounded = true;
                        varData.lowerBound = number;
                    }
                    else if( lineTokens[1] == ">" || lineTokens[1] == ">=" )
                    {
                        // number >= x.
                        varData.upperBounded = true;
                        varData.upperBound = number;
                    }
                    else
                    {
                        // number = x.
                        varData.nonnegative = false;
                        varData.fixed = true;
                        varData.fixedValue = number;
                    }
                }
                else
                {
                    const bool boundIsInfinite =
                      IsSignedInfinity( lineTokens[2] );

                    auto varIter = variableDict_.find( lineTokens[0] );
                    if( varIter == variableDict_.end() )
                        RuntimeError("Could not find variable ",lineTokens[0]);
                    auto& varData = varIter->second;

                    if( lineTokens[1] == "<" || lineTokens[1] == "<=" )
                    {
                        // x <= number.
                        if( !boundIsInfinite )
                        {
                            const Real number =
                              StringToReal<Real>( lineTokens[2] );
                            varData.upperBounded = true;
                            varData.upperBound = number;
                        }
                    }
                    else if( lineTokens[1] == ">" || lineTokens[1] == ">=" )
                    {
                        // x >= number.
                        varData.nonnegative = false;
                        if( !boundIsInfinite )
                        {
                            const Real number =
                              StringToReal<Real>( lineTokens[2] );
                            varData.lowerBounded = true;
                            varData.lowerBound = number;
                        }
                        // else: variable is implicitly free unless another
                        // bound is applied.
                    }
                    else
                    {
                        // x = number.
                        varData.nonnegative = false;
                        if( boundIsInfinite )
                            RuntimeError("Tried to fix variable to infinity");
                        const Real number =
                          StringToReal<Real>( lineTokens[2] );
                        varData.fixed = true;
                        varData.fixedValue = number; 
                    }
                }
            }
            else if( lineTokens.size() == 5 )
            {
                // We should have a five-token line of the form
                // 'lowerBound <= x <= upperBound'.
                if( !IsRelationToken(lineTokens[1]) ||
                    !IsRelationToken(lineTokens[3]) )
                    RuntimeError
                    ("Unrecognized bound line of (",lineTokens[0],",", 
                     lineTokens[1],",",lineTokens[2],",",lineTokens[3],",",
                     lineTokens[4],")");

                auto varIter = variableDict_.find( lineTokens[2] );
                if( varIter == variableDict_.end() )
                    RuntimeError("Could not find variable ",lineTokens[2]);
                auto& varData = varIter->second;

                if( lineTokens[1] != "<" &&
                    lineTokens[1] != "<=" )
                    RuntimeError("Expected ",lineTokens[1]," to be < or <=");
                if( lineTokens[3] != "<" &&
                    lineTokens[3] != "<=" )
                    RuntimeError("Expected ",lineTokens[3]," to be < or <=");

                const bool lowerBoundIsNegativeInfinity =
                  IsNegativeInfinity( lineTokens[0] );
                const bool upperBoundIsPositiveInfinity =
                  IsPositiveInfinity( lineTokens[4] );

                varData.nonnegative = false;
                if( lowerBoundIsNegativeInfinity &&
                    upperBoundIsPositiveInfinity )
                {
                    varData.free = true;
                }
                else if( lowerBoundIsNegativeInfinity )
                {
                    const Real upperBound = StringToReal<Real>( lineTokens[4] );
                    varData.upperBounded = true;
                    varData.upperBound = upperBound;
                }
                else if( upperBoundIsPositiveInfinity )
                {
                    const Real lowerBound = StringToReal<Real>( lineTokens[0] );
                    varData.lowerBounded = true;
                    varData.lowerBound = lowerBound;
                }
                else
                {
                    const Real lowerBound = StringToReal<Real>( lineTokens[0] );
                    const Real upperBound = StringToReal<Real>( lineTokens[4] );
                    varData.lowerBounded = true;
                    varData.upperBounded = true;
                    varData.upperBound = upperBound;
                    varData.lowerBound = lowerBound;
                }
            }
        }
    }

    if( section == LPCPLEX_GENERAL )
    {
        // We do not support MIPs yet.
        Output("WARNING: GENERAL section is being skipped.");
        while( true )
        {
            // Get the first token.
            if( !ExtractToken( file_, token ) )
                RuntimeError(missingEndLabel);
            auto tokenSection = DecodeSection( token );
            if( tokenSection != LPCPLEX_NONE )
            {
                if( progress_ )
                    Output("Jumping from GENERAL to ",token);
                section = tokenSection;
                break;
            }
        }
    }

    if( section != LPCPLEX_END )
        RuntimeError("Improperly handled section");

    SimplifyRowTypes();
    CountBoundTypes();

    // Extract the number of variables
    // (the matrix 'A' is 'm x n' and 'G' is 'k x n').
    //
    // The entries will fall into two sets of sizes:
    //   variableDict_.size() - meta_.numFixedBounds, and
    //   meta_.numFixedBounds.
    //
    meta_.n = variableDict_.size();

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
    section_ = LPCPLEX_NONE;

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

template<typename Real>
bool LPCPLEXReader<Real>::QueuedEntry()
{
    EL_DEBUG_CSE
    using namespace El::lp_cplex;
    string token;
    if( queuedEntries_.size() > 0 )
        return true;

    // Scan the objective information
    // TODO(poulson): Support constants?
    if( section_ == LPCPLEX_NONE )
    {
        if( !ExtractToken( file_, token ) )
            RuntimeError("Could not extract first token");
        section_ = DecodeSection( token );

        if( !ExtractToken( file_, token ) )
            RuntimeError("There was not second token!");
        auto secondTokenSection = DecodeSection( token );
        if( secondTokenSection != LPCPLEX_NONE )
        {
            if( progress_ )
                Output("Jumping from OBJECTIVE to ",token);
            section_ = secondTokenSection;
        }
        else
        {
            if( file_.peek() == ':' )
            {
                meta_.objectiveName = token;
                file_.get();
            }

            // Loop until the token stores the next section label.
            Real multiple = 1;
            while( ExtractToken( file_, token ) )
            {
                auto tokenSection = DecodeSection( token );
                if( tokenSection != LPCPLEX_NONE )
                {
                    if( progress_ )
                        Output("Jumping from OBJECTIVE to ",token);
                    section_ = tokenSection;
                    break;
                }
                if( IsNumberCharacter(token[0]) )
                {
                    if( token == "+" )
                    {
                        multiple = 1;
                    }
                    else if( token == "-" )
                    {
                        multiple = -1;
                    }
                    else
                    {
                        const Real number = StringToReal<Real>( token );
                        multiple *= number;
                    }
                }
                else
                {
                    auto varIter = variableDict_.find( token );
                    if( varIter == variableDict_.end() )
                        RuntimeError("Couldn't find variable ",token);
                    auto& varData = varIter->second;
                    AffineLPEntry<Real> entry;
                    entry.type = AFFINE_LP_COST_VECTOR;
                    entry.row = varData.index;
                    entry.column = 0;
                    entry.value = minimize_ ? multiple : -multiple;
                    queuedEntries_.push_back( entry );

                    multiple = 1;
                }
            }
        }
    }
    if( queuedEntries_.size() > 0 )
        return true;

    if( section_ == LPCPLEX_SUBJECT_TO )
    {
        if( progress_ )
            Output("Entered Subject To section");
        // Each of the following outer loops will handle an (in)equality.
        while( true )
        {
            // Get the first token.
            if( !ExtractToken( file_, token ) )
                RuntimeError(missingEndLabel);
            auto tokenSection = DecodeSection( token );
            if( tokenSection != LPCPLEX_NONE )
            {
                if( progress_ )
                    Output("Jumping from SUBJECT TO to ",token);
                section_ = tokenSection;
                break;
            }

            const auto& rowData = rowData_[rowCounter_++];
            if( file_.peek() == ':' )
            {
                file_.get();
                if( !ExtractToken( file_, token ) )
                    RuntimeError(missingEndLabel);
            }

            bool pastRelation = false;
            Real multiple = 1;
            while( true )
            {
                if( !pastRelation )
                {
                    if( IsRelationToken( token ) )
                    {
                        pastRelation = true;
                    }
                    else if( IsNumberCharacter(token[0]) )
                    {
                        // We would accumulate here.
                        if( token == "+" )
                        {
                            multiple = 1;
                        }
                        else if( token == "-" )
                        {
                            multiple = -1;
                        }
                        else
                        {
                            const Real number = StringToReal<Real>( token );
                            multiple *= number;
                        }
                    }
                    else
                    {
                        // Look up the token in the variable dictionary.
                        auto varIter = variableDict_.find( token );
                        if( varIter == variableDict_.end() )
                            RuntimeError("Could not find variable ",token);
                        const auto& varData = varIter->second;
                        const Int column = varData.index;
                        
                        if( rowData.type == LPCPLEX_GREATER_ROW )
                        {
                            // G(row,col) = -multiple
                            const Int row =
                              meta_.greaterOffset + rowData.typeIndex;
                            AffineLPEntry<Real> entry;
                            entry.type = AFFINE_LP_INEQUALITY_MATRIX;
                            entry.row = row;
                            entry.column = column;
                            entry.value = -multiple;
                            queuedEntries_.push_back( entry );
                        }
                        else if( rowData.type == LPCPLEX_LESSER_ROW )
                        {
                            // G(row,col) = multiple
                            const Int row =
                              meta_.lesserOffset + rowData.typeIndex;
                            AffineLPEntry<Real> entry;
                            entry.type = AFFINE_LP_INEQUALITY_MATRIX;
                            entry.row = row;
                            entry.column = column;
                            entry.value = multiple;
                            queuedEntries_.push_back( entry );
                        }
                        else if( rowData.type == LPCPLEX_EQUALITY_ROW )
                        {
                            // A(row,col) = multiple
                            const Int row =
                              meta_.equalityOffset + rowData.typeIndex;
                            AffineLPEntry<Real> entry;
                            entry.type = AFFINE_LP_EQUALITY_MATRIX;
                            entry.row = row;
                            entry.column = column;
                            entry.value = multiple;
                            queuedEntries_.push_back( entry );
                        }
                        else
                            RuntimeError("Unknown row type");

                        multiple = 1;
                    }
                }
                else
                {
                    const Real number = StringToReal<Real>( token );
                    if( rowData.type == LPCPLEX_GREATER_ROW )
                    {
                        // h(row) = -number
                        const Int row =
                          meta_.greaterOffset + rowData.typeIndex;
                        AffineLPEntry<Real> entry;
                        entry.type = AFFINE_LP_INEQUALITY_VECTOR;
                        entry.row = row;
                        entry.column = 0;
                        entry.value = -number;
                        queuedEntries_.push_back( entry );
                    }
                    else if( rowData.type == LPCPLEX_LESSER_ROW )
                    {
                        // h(row) = number
                        const Int row =
                          meta_.lesserOffset + rowData.typeIndex;
                        AffineLPEntry<Real> entry;
                        entry.type = AFFINE_LP_INEQUALITY_VECTOR;
                        entry.row = row;
                        entry.column = 0;
                        entry.value = number;
                        queuedEntries_.push_back( entry );
                    }
                    else
                    {
                        // b(row) = number
                        const Int row =
                          meta_.equalityOffset + rowData.typeIndex;
                        AffineLPEntry<Real> entry;
                        entry.type = AFFINE_LP_EQUALITY_VECTOR;
                        entry.row = row;
                        entry.column = 0;
                        entry.value = number;
                        queuedEntries_.push_back( entry );
                    }
                    // We would store the right-hand side value here...
                    break;
                }
                if( !ExtractToken( file_, token ) )
                    RuntimeError(missingEndLabel);
            }
            break;
        }
    }
    if( queuedEntries_.size() > 0 )
        return true;

    QueueVariableBound();

    return queuedEntries_.size() > 0;
}

template<typename Real>
AffineLPEntry<Real> LPCPLEXReader<Real>::GetEntry()
{
    EL_DEBUG_CSE
    using namespace El::lp_cplex;
    if( queuedEntries_.size() == 0 )
        LogicError("No entries are currently enqueued");
    AffineLPEntry<Real> entry = queuedEntries_.back();
    queuedEntries_.pop_back();
    return entry;
}

template<typename Real>
const LPCPLEXMeta& LPCPLEXReader<Real>::Meta() const
{
    EL_DEBUG_CSE
    return meta_;
}

namespace read_lp_cplex {

template<typename Real>
LPCPLEXMeta Helper
( AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const string& filename )
{
    EL_DEBUG_CSE
    using namespace El::lp_cplex;

    LPCPLEXReader<Real> reader( filename );
    const auto& meta = reader.Meta();

    Zeros( problem.c, meta.n, 1 );
    Zeros( problem.A, meta.m, meta.n );
    Zeros( problem.b, meta.m, 1 );
    Zeros( problem.G, meta.k, meta.n );
    Zeros( problem.h, meta.k, 1 );
    while( reader.QueuedEntry() )
    {
        const auto entry = reader.GetEntry();
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
LPCPLEXMeta Helper
( AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
  const string& filename )
{
    EL_DEBUG_CSE
    using namespace El::lp_cplex;

    // TODO(poulson): Consider loading on a single process and then distributing
    // the results instead.

    LPCPLEXReader<Real> reader( filename );
    const auto& meta = reader.Meta();

    Zeros( problem.c, meta.n, 1 );
    Zeros( problem.A, meta.m, meta.n );
    Zeros( problem.b, meta.m, 1 );
    Zeros( problem.G, meta.k, meta.n );
    Zeros( problem.h, meta.k, 1 );

    while( reader.QueuedEntry() )
    {
        const auto entry = reader.GetEntry();
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
LPCPLEXMeta Helper
( AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const string& filename )
{
    EL_DEBUG_CSE
    using namespace El::lp_cplex;

    LPCPLEXReader<Real> reader( filename );
    const auto& meta = reader.Meta();

    Zeros( problem.c, meta.n, 1 );
    Zeros( problem.A, meta.m, meta.n );
    Zeros( problem.b, meta.m, 1 );
    Zeros( problem.G, meta.k, meta.n );
    Zeros( problem.h, meta.k, 1 );

    problem.A.Reserve( meta.numEqualityEntries );
    problem.G.Reserve( meta.numInequalityEntries );
    while( reader.QueuedEntry() )
    {
        const auto entry = reader.GetEntry();
        if( entry.type == AFFINE_LP_COST_VECTOR )
        {
            problem.c.Set( entry.row, 0, entry.value );
        }
        else if( entry.type == AFFINE_LP_EQUALITY_MATRIX )
        {
            problem.A.QueueUpdate( entry.row, entry.column, entry.value );
        }
        else if( entry.type == AFFINE_LP_EQUALITY_VECTOR )
        {
            problem.b.Set( entry.row, 0, entry.value );
        }
        else if( entry.type == AFFINE_LP_INEQUALITY_MATRIX )
        {
            problem.G.QueueUpdate( entry.row, entry.column, entry.value );
        }
        else /* entry.type == AFFINE_LP_INEQUALITY_VECTOR */
        {
            problem.h.Set( entry.row, 0, entry.value );
        }
    }
    problem.A.ProcessQueues();
    problem.G.ProcessQueues();

    return meta;
}

template<typename Real>
LPCPLEXMeta Helper
( AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
  const string& filename )
{
    EL_DEBUG_CSE
    using namespace El::lp_cplex;

    LPCPLEXReader<Real> reader( filename );
    const auto& meta = reader.Meta();

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
        const auto entry = reader.GetEntry();
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

} // namespace read_lp_cplex

template<class MatrixType,class VectorType>
LPCPLEXMeta ReadLPCPLEX
( AffineLPProblem<MatrixType,VectorType>& problem,
  const string& filename )
{
    EL_DEBUG_CSE
    return read_lp_cplex::Helper( problem, filename );
}

namespace write_lp_cplex {

template<typename Real>
void Helper
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const string& filename )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

template<typename Real>
void Helper
( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
  const string& filename )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

template<typename Real>
void Helper
( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const string& filename )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

template<typename Real>
void Helper
( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
  const string& filename )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

template<typename Real>
void Helper
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const string& filename )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

template<typename Real>
void Helper
( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
  const string& filename )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

template<typename Real>
void Helper
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const string& filename )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

template<typename Real>
void Helper
( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
  const string& filename )
{
    EL_DEBUG_CSE
    LogicError("This routine is not yet written");
}

} // namespace write_lp_cplex

template<class MatrixType,class VectorType>
void WriteLPCPLEX
( const DirectLPProblem<MatrixType,VectorType>& problem,
  const string& filename )
{
    EL_DEBUG_CSE
    write_lp_cplex::Helper( problem, filename );
}

template<class MatrixType,class VectorType>
void WriteLPCPLEX
( const AffineLPProblem<MatrixType,VectorType>& problem,
  const string& filename )
{
    EL_DEBUG_CSE
    write_lp_cplex::Helper( problem, filename );
}

} // namespace El
