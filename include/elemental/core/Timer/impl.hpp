/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TIMER_IMPL_HPP
#define ELEM_TIMER_IMPL_HPP

namespace elem {

inline Timer::Timer( const std::string& name )
: name_(name)
{ }

inline void
Timer::Start()
{
    DEBUG_ONLY(
        if( running_ )
            LogicError("Forgot to stop timer before restarting.");
    )
#ifdef ELEM_HAVE_STEADYCLOCK
    lastTime_ = steady_clock::now();
#else
    lastTime_ = high_resolution_clock::now();
#endif
    running_ = true;
}

inline double
Timer::Stop()
{
    DEBUG_ONLY(
        if( !running_ )
            LogicError("Tried to stop a timer before starting it.");
    )
    lastPartialTime_ = Partial();
    running_ = false;
    totalTime_ += lastPartialTime_;
    return lastPartialTime_;
}

inline void
Timer::Reset( const std::string& name )
{ 
    name_ = name;
    running_ = false;
    totalTime_ = 0; 
    lastPartialTime_ = 0;
}

inline const std::string&
Timer::Name() const
{ return name_; }

inline double
Timer::Partial() const
{ 
    if( running_ )
    {
#ifdef ELEM_HAVE_STEADYCLOCK
        auto now = steady_clock::now();
#else
        auto now = high_resolution_clock::now();
#endif
        auto timeSpan = duration_cast<duration<double>>(now-lastTime_);
        return timeSpan.count();
    }
    else
        return lastPartialTime_; 
}

inline double
Timer::Total() const
{
    if( running_ )
        return totalTime_ + Partial();
    else
        return totalTime_;
}

} // namespace elem

#endif // ifndef ELEM_TIMER_IMPL_HPP
