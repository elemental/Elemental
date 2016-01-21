/*
   Copyright (c) 2013, Jack Poulson
   All rights reserved.

   This file is a modification of Choice, a simple command-line option library.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_IMPORTS_MPICHOICE_HPP
#define EL_IMPORTS_MPICHOICE_HPP

#include <algorithm>

namespace El {
namespace choice {

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;
using std::string;
using std::vector;

class MpiArgs
{
public:
    MpiArgs
    ( int argc, char** argv, 
      mpi::Comm comm=mpi::COMM_WORLD, ostream& error=cerr );
    virtual ~MpiArgs() { }

    template<typename T>
    T Input( string name, string desc );
    template<typename T>
    T Input( string name, string desc, T defaultVal );

    void Process( ostream& os=cout ) const;
    void PrintReport( ostream& os=cout ) const;

protected:
    int argc_;
    char** argv_;
    vector<bool> usedArgs_;
    ostream& error_;
    mpi::Comm comm_;

    virtual void HandleVersion( ostream& os=cout ) const { }
    virtual void HandleBuild( ostream& os=cout ) const { }

    struct RequiredArg
    { 
        string name, desc, typeInfo, usedVal; 
        bool found;

        RequiredArg( string n, string d, string t, string uv, bool f ) 
        : name(n), desc(d), typeInfo(t), usedVal(uv), found(f) { };
    };

    struct OptionalArg
    { 
        string name, desc, typeInfo, defaultVal, usedVal; 
        bool found;

        OptionalArg
        ( string n, string d, string t, string dv, string uv, bool f )
        : name(n), desc(d), typeInfo(t), 
          defaultVal(dv), usedVal(uv), found(f) { } 
    };

    vector<RequiredArg> requiredArgs_;
    vector<OptionalArg> optionalArgs_;
};

inline
MpiArgs::MpiArgs( int argc, char** argv, mpi::Comm comm, ostream& error )
: argc_(argc), argv_(argv), usedArgs_(argc,false), error_(error), comm_(comm)
{ }

template<typename T>
inline T
MpiArgs::Input( string name, string desc )
{
    const int commRank = mpi::Rank( comm_ );

    char** arg = std::find( argv_, argv_+argc_, name );
    const bool found = ( arg != argv_+argc_ );
    const bool invalidFound = ( arg == argv_+argc_-1 );
    if( invalidFound )
    {
        if( commRank == 0 )
            error_ << "Missing value for last command-line argument" << endl;
        throw ArgException();
    }

    const string typeInfo = TypeName<T>();
    string usedVal = ( found ? arg[1] : "N/A" );
    requiredArgs_.push_back( RequiredArg(name,desc,typeInfo,usedVal,found) );

    // Before returning, store the used indices and check for duplication
    if( commRank == 0 && found )
    {
        const int offset = arg - argv_;
        if( usedArgs_[offset] || usedArgs_[offset+1] )
        {
            error_ << "WARNING: conflict with " << name << " detected at ";
            if( usedArgs_[offset] && usedArgs_[offset+1] )
                error_ << "arguments " << offset << " and " << offset+1
                       << endl;
            else if( usedArgs_[offset] )
                error_ << "argument " << offset << endl;
            else
                error_ << "argument " << offset+1 << endl;
            error_ << "Please ensure that you did request argument " 
                   << name << " multiple times" << endl;
        }
        usedArgs_[offset+0] = true;
        usedArgs_[offset+1] = true;

        arg = std::find( arg+1, argv_+argc_, name );
        if( arg != argv_+argc_ )
            error_ << "WARNING: " << name << " was specified twice and only "
                   << "the first instance is used" << endl;
    }

    return Cast<T>( usedVal );
}

template<typename T>
inline T
MpiArgs::Input( string name, string desc, T defaultVal )
{
    const int commRank = mpi::Rank( comm_ );

    char** arg = std::find( argv_, argv_+argc_, name );
    const bool found = ( arg != argv_+argc_ );
    const bool invalidFound = ( arg == argv_+argc_-1 );
    if( invalidFound )
    {
        if( commRank == 0 )
            error_ << "Missing value for last command-line argument" 
                   << endl;
        throw ArgException();
    }

    const string typeInfo = TypeName<T>();
    string defValString = Cast<string>( defaultVal );
    string usedVal = ( found ? arg[1] : defValString );

    optionalArgs_.push_back
    ( OptionalArg(name,desc,typeInfo,defValString,usedVal,found) );

    // Before returning, store the used indices and check for duplication
    if( commRank == 0 && found )
    {
        const int offset = arg - argv_;
        if( usedArgs_[offset] || usedArgs_[offset+1] )
        {
            error_ << "WARNING: conflict with " << name << " detected at ";
            if( usedArgs_[offset] && usedArgs_[offset+1] )
                error_ << "arguments " << offset << " and " << offset+1
                       << endl;
            else if( usedArgs_[offset] )
                error_ << "argument " << offset << endl;
            else
                error_ << "argument " << offset+1 << endl;
            error_ << "Please ensure that you did request argument " 
                   << name << " multiple times" << endl;
        }
        usedArgs_[offset+0] = true;
        usedArgs_[offset+1] = true;

        arg = std::find( arg+1, argv_+argc_, name );
        if( arg != argv_+argc_ )
            error_ << "WARNING: " << name << " was specified twice and only "
                   << "the first instance is used" << endl;
    }

    if( found )
        return Cast<T>( usedVal );
    else
        return defaultVal; // avoid the double-cast
}

inline void 
MpiArgs::Process( ostream& os ) const
{
    HandleVersion( os );
    HandleBuild( os );

    string help = "--help";
    char** arg = std::find( argv_, argv_+argc_, help );
    const bool foundHelp = ( arg != argv_+argc_ );
    
    int numFailed = 0;
    const int numRequired = requiredArgs_.size();
    for( int i=0; i<numRequired; ++i )
        if( !requiredArgs_[i].found )
            ++numFailed;
    if( numFailed > 0 || foundHelp )
    {
        PrintReport( os );
        throw ArgException(); 
    }
}

inline void 
MpiArgs::PrintReport( ostream& os ) const
{
    const int commRank = mpi::Rank( comm_ );
    if( commRank != 0 )
        return;

    const int numRequired = requiredArgs_.size();
    const int numOptional = optionalArgs_.size();

    if( numRequired > 0 )
        os << "Required arguments:\n";
    int numReqFailed = 0;
    for( int i=0; i<numRequired; ++i )
    {
        const RequiredArg& reqArg = requiredArgs_[i];
        if( !reqArg.found )
            ++numReqFailed;
        string foundString = ( reqArg.found ? "found" : "NOT found" );
        os << "  " << reqArg.name
           << " [" << reqArg.typeInfo << "," << reqArg.usedVal << ","
           << foundString << "]\n"
           << "    " << reqArg.desc << "\n\n";
    }

    if( numOptional > 0 )
        os << "Optional arguments:\n";
    int numOptFailed = 0;
    for( int i=0; i<numOptional; ++i )
    {
        const OptionalArg& optArg = optionalArgs_[i];
        if( !optArg.found )
            ++numOptFailed;
        string foundString = ( optArg.found ? "found" : "NOT found" );
        os << "  " << optArg.name
           << " [" << optArg.typeInfo
           << "," << optArg.defaultVal << "," << optArg.usedVal << ","
           << foundString << "]\n"
           << "    " << optArg.desc << "\n\n";
    }

    os << "Out of " << numRequired << " required arguments, " 
       << numReqFailed << " were not specified." << endl;

    os << "Out of " << numOptional << " optional arguments, "
       << numOptFailed << " were not specified.\n" << endl;
}

} // namespace choice
} // namespace El

#endif // ifndef EL_IMPORTS_MPICHOICE_HPP
