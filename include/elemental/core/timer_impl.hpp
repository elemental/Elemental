/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

inline 
Timer::Timer()
: running_(false), time_(0), name_("[blank]")
{ }

inline 
Timer::Timer( const std::string name )
: running_(false), time_(0), name_(name)
{ }

inline void 
Timer::Start()
{
#ifndef RELEASE
    PushCallStack("Timer::Start");
    if( running_ )
        throw std::logic_error("Forgot to stop timer before restarting");
#endif
    lastStartTime_ = mpi::Time();
    running_ = true;
    running_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void 
Timer::Stop()
{
#ifndef RELEASE
    PushCallStack("Timer::Stop");
    if( !running_ )
        throw std::logic_error("Tried to stop a timer before starting it");
#endif
    time_ += mpi::Time()-lastStartTime_;
    running_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void 
Timer::Reset()
{ time_ = 0; }

inline const std::string 
Timer::Name() const
{ return name_; }

inline double 
Timer::Time() const
{
#ifndef RELEASE
    PushCallStack("Timer::Time");
    if( running_ )
        throw std::logic_error("Asked for time while still timing");
    PopCallStack();
#endif
    return time_;
}

} // namespace elem
