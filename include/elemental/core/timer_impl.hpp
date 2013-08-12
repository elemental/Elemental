/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_TIMER_IMPL_HPP
#define ELEM_CORE_TIMER_IMPL_HPP

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
    CallStackEntry entry("Timer::Start");
    if( running_ )
        LogicError("Forgot to stop timer before restarting");
#endif
    lastStartTime_ = mpi::Time();
    running_ = true;
}

inline double
Timer::Stop()
{
#ifndef RELEASE
    CallStackEntry entry("Timer::Stop");
    if( !running_ )
        LogicError("Tried to stop a timer before starting it");
#endif
    const double partial = Partial();  
    time_ += partial;
    running_ = false;
    return partial;
}

inline void 
Timer::Reset()
{ time_ = 0; }

inline const std::string 
Timer::Name() const
{ return name_; }

inline double 
Timer::Partial() const
{
#ifndef RELEASE
    CallStackEntry entry("Timer::Partial");
#endif
    return mpi::Time()-lastStartTime_;
}

inline double 
Timer::Total() const
{
#ifndef RELEASE
    CallStackEntry entry("Timer::Total");
#endif
    return time_;
}

} // namespace elem

#endif // ifndef ELEM_CORE_TIMER_IMPL_HPP
