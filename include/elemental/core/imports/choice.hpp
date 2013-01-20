/*
   Copyright (c) 2012-2013, Jack Poulson
   All rights reserved.

   This file is a modification of Choice, a simple command-line option library.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CORE_CHOICE_HPP
#define CORE_CHOICE_HPP

#include <typeinfo>

namespace elem {

template<typename TOut,typename TIn>
inline TOut Cast( const TIn& input )
{
    std::stringstream stream;
    TOut output;

    stream << input;
    stream >> output;

    return output;
}

template<>
inline bool Cast( const std::string& input )
{
    std::string trueString("true");
    std::string falseString("false");
    if( input.compare(trueString) == 0 )
        return true;
    else if( input.compare(falseString) == 0 )
        return false;
    else
    {
        bool output;
        std::stringstream stream;
        stream << input;
        stream >> output;
        return output;
    }
}

class ArgException : public std::logic_error
{
public:
    ArgException( const char* msg="Argument exception" )
    : std::logic_error( msg ) { }
};

class Args
{
public:
    Args( int argc, char** argv, std::ostream& error=std::cerr );

    template<typename T>
    T Input( std::string name, std::string desc );
    template<typename T>
    T Input( std::string name, std::string desc, T defaultVal );

    void Process( std::ostream& output=std::cout ) const;
    void PrintReport( std::ostream& output=std::cout ) const;

private:
    int argc_;
    char** argv_;
    std::vector<bool> usedArgs_;
    std::ostream& error_;

    struct RequiredArg
    { 
        std::string name, desc, typeInfo, usedVal; 
        bool found;

        RequiredArg
        ( std::string n, std::string d, std::string t, std::string uv, bool f ) 
        : name(n), desc(d), typeInfo(t), usedVal(uv), found(f) { };
    };

    struct OptionalArg
    { 
        std::string name, desc, typeInfo, defaultVal, usedVal; 
        bool found;

        OptionalArg
        ( std::string n, std::string d, std::string t, 
          std::string dv, std::string uv, bool f )
        : name(n), desc(d), typeInfo(t), 
          defaultVal(dv), usedVal(uv), found(f) { } 
    };

    std::vector<RequiredArg> requiredArgs_;
    std::vector<OptionalArg> optionalArgs_;
};

inline
Args::Args( int argc, char** argv, std::ostream& error )
: argc_(argc), argv_(argv), usedArgs_(argc,false), error_(error)
{ }

template<typename T>
inline T
Args::Input( std::string name, std::string desc )
{
    char** arg = std::find( argv_, argv_+argc_, name );
    const bool found = ( arg != argv_+argc_ );
    const bool invalidFound = ( arg == argv_+argc_-1 );
    if( invalidFound )
    {
        error_ << "Missing value for last command-line argument" << std::endl;
        throw ArgException();
    }

    std::string typeInfo( typeid(T).name() );
    std::string usedVal = ( found ? arg[1] : "N/A" );
    requiredArgs_.push_back( RequiredArg(name,desc,typeInfo,usedVal,found) );

    // Before returning, store the used indices and check for duplication
    if( found )
    {
        const int offset = arg - argv_;
        if( usedArgs_[offset] || usedArgs_[offset+1] )
        {
            error_ << "WARNING: conflict with " << name << " detected at ";
            if( usedArgs_[offset] && usedArgs_[offset+1] )
                error_ << "arguments " << offset << " and " << offset+1
                       << std::endl;
            else if( usedArgs_[offset] )
                error_ << "argument " << offset << std::endl;
            else
                error_ << "argument " << offset+1 << std::endl;
            error_ << "Please ensure that you did request argument "
                   << name << " multiple times" << std::endl;
        }
        usedArgs_[offset+0] = true;
        usedArgs_[offset+1] = true;

        arg = std::find( arg+1, argv_+argc_, name );
        if( arg != argv_+argc_ )
            error_ << "WARNING: " << name << " was specified twice and only "
                   << "the first instance is used" << std::endl;
    }

    return Cast<T>( usedVal );
}

template<typename T>
inline T
Args::Input( std::string name, std::string desc, T defaultVal )
{
    char** arg = std::find( argv_, argv_+argc_, name );
    const bool found = ( arg != argv_+argc_ );
    const bool invalidFound = ( arg == argv_+argc_-1 );
    if( invalidFound )
    {
        error_ << "Missing value for last command-line argument" << std::endl;
        throw ArgException();
    }

    std::string typeInfo( typeid(T).name() );

    std::string defValString = Cast<std::string>( defaultVal );
    std::string usedVal = ( found ? arg[1] : defValString );

    optionalArgs_.push_back
    ( OptionalArg(name,desc,typeInfo,defValString,usedVal,found) );

    // Before returning, store the used indices and check for duplication
    if( found )
    {
        const int offset = arg - argv_;
        if( usedArgs_[offset] || usedArgs_[offset+1] )
        {
            error_ << "WARNING: conflict with " << name << " detected at ";
            if( usedArgs_[offset] && usedArgs_[offset+1] )
                error_ << "arguments " << offset << " and " << offset+1
                       << std::endl;
            else if( usedArgs_[offset] )
                error_ << "argument " << offset << std::endl;
            else
                error_ << "argument " << offset+1 << std::endl;
            error_ << "Please ensure that you did request argument "
                   << name << " multiple times" << std::endl;
        }
        usedArgs_[offset+0] = true;
        usedArgs_[offset+1] = true;

        arg = std::find( arg+1, argv_+argc_, name );
        if( arg != argv_+argc_ )
            error_ << "WARNING: " << name << " was specified twice and only "
                   << "the first instance is used" << std::endl;
    }

    if( found )
        return Cast<T>( usedVal );
    else
        return defaultVal; // avoid the double-cast
}

inline void
Args::Process( std::ostream& output ) const
{
    std::string help = "--help";
    char** arg = std::find( argv_, argv_+argc_, help );
    const bool foundHelp = ( arg != argv_+argc_ );

    int numFailed = 0;
    const int numRequired = requiredArgs_.size();
    for( int i=0; i<numRequired; ++i )
        if( !requiredArgs_[i].found )
            ++numFailed;
    if( numFailed > 0 || foundHelp )
    {
        PrintReport( output );
        throw ArgException();
    }
}

inline void
Args::PrintReport( std::ostream& output ) const
{
    const int numRequired = requiredArgs_.size();
    const int numOptional = optionalArgs_.size();

    if( numRequired > 0 )
        output << "Required arguments:\n";
    int numReqFailed = 0;
    for( int i=0; i<numRequired; ++i )
    {
        const RequiredArg& reqArg = requiredArgs_[i];
        if( !reqArg.found )
            ++numReqFailed;
        std::string foundString = ( reqArg.found ? "found" : "NOT found" );
        output << "  " << reqArg.name
               << " [" << reqArg.typeInfo << "," << reqArg.usedVal << ","
               << foundString << "]\n"
               << "    " << reqArg.desc << "\n\n";
    }

    if( numOptional > 0 )
        output << "Optional arguments:\n";
    int numOptFailed = 0;
    for( int i=0; i<numOptional; ++i )
    {
        const OptionalArg& optArg = optionalArgs_[i];
        if( !optArg.found )
            ++numOptFailed;
        std::string foundString = ( optArg.found ? "found" : "NOT found" );
        output << "  " << optArg.name
               << " [" << optArg.typeInfo
               << "," << optArg.defaultVal << "," << optArg.usedVal << ","
               << foundString << "]\n"
               << "    " << optArg.desc << "\n\n";
    }

    output << "Out of " << numRequired << " required arguments, " 
           << numReqFailed << " were not specified." << std::endl;

    output << "Out of " << numOptional << " optional arguments, "
           << numOptFailed << " were not specified.\n" << std::endl;
}

} // namespace elem

#endif // ifndef CORE_CHOICE_HPP
