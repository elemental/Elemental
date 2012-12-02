/*
   Copyright (c) 2012, Jack Poulson
   All rights reserved.

   This file is a modification of Choice, a simple command-line option library.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef MPI_CHOICE_H
#define MPI_CHOICE_H 1

namespace elem {

class MpiArgs
{
public:
    MpiArgs
    ( int argc, char** argv, 
      mpi::Comm comm=mpi::COMM_WORLD, std::ostream& error=std::cerr );

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
    mpi::Comm comm_;

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

    template<typename TOut,typename TIn>
    TOut Cast( const TIn& input )
    {
        std::stringstream stream;
        TOut output;

        stream << input;
        stream >> output;

        return output;
    }

    std::vector<RequiredArg> requiredArgs_;
    std::vector<OptionalArg> optionalArgs_;
};

inline
MpiArgs::MpiArgs( int argc, char** argv, mpi::Comm comm, std::ostream& error )
: argc_(argc), argv_(argv), usedArgs_(argc,false), error_(error), comm_(comm)
{ }

template<typename T>
inline T
MpiArgs::Input( std::string name, std::string desc )
{
    const int commRank = mpi::CommRank( comm_ );

    char** arg = std::find( argv_, argv_+argc_, name );
    const bool found = ( arg != argv_+argc_ );
    const bool invalidFound = ( arg == argv_+argc_-1 );
    if( invalidFound )
    {
        if( commRank == 0 )
            error_ << "Missing value for last command-line argument" 
                   << std::endl;
        throw ArgException();
    }

    std::string typeInfo( typeid(T).name() );
    std::string usedVal = ( found ? arg[1] : "N/A" );
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
MpiArgs::Input( std::string name, std::string desc, T defaultVal )
{
    const int commRank = mpi::CommRank( comm_ );

    char** arg = std::find( argv_, argv_+argc_, name );
    const bool found = ( arg != argv_+argc_ );
    const bool invalidFound = ( arg == argv_+argc_-1 );
    if( invalidFound )
    {
        if( commRank == 0 )
            error_ << "Missing value for last command-line argument" 
                   << std::endl;
        throw ArgException();
    }

    std::string typeInfo( typeid(T).name() );

    std::string defValString = Cast<std::string>( defaultVal );
    std::string usedVal = ( found ? arg[1] : defValString );

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
MpiArgs::Process( std::ostream& output ) const
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
MpiArgs::PrintReport( std::ostream& output ) const
{
    const int commRank = mpi::CommRank( comm_ );
    if( commRank != 0 )
        return;

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
        output << "  " << reqArg.name << "\n"
               << "    description: " << reqArg.desc << "\n"
               << "    type string: " << reqArg.typeInfo << "\n"
               << "    used value:  " << reqArg.usedVal << "\n";
        if( reqArg.found )
            output << "    found\n\n";
        else
            output << "    NOT found\n\n";
    }

    if( numOptional > 0 )
        output << "Optional arguments:\n";
    int numOptFailed = 0;
    for( int i=0; i<numOptional; ++i )
    {
        const OptionalArg& optArg = optionalArgs_[i];
        if( !optArg.found )
            ++numOptFailed;
        output << "  " << optArg.name << "\n"
               << "    description:   " << optArg.desc << "\n"
               << "    type string:   " << optArg.typeInfo << "\n"
               << "    default value: " << optArg.defaultVal << "\n"
               << "    used value:    " << optArg.usedVal << "\n";
        if( optArg.found )
            output << "    found\n\n";
        else
            output << "    NOT found\n\n";
    }

    output << "Out of " << numRequired << " required arguments, " 
           << numReqFailed << " were not specified." << std::endl;

    output << "Out of " << numOptional << " optional arguments, "
           << numOptFailed << " were not specified.\n" << std::endl;
}

} // namespace elem

#endif // ifndef MPI_CHOICE_H
