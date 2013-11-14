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

inline Timer::Timer( const std::string& name )
: name_(name)
{ }

inline void
Timer::Start()
{
#ifndef RELEASE
    if( running_ )
        throw std::logic_error("Forgot to stop timer before restarting.");
#endif
    lastTime_ = steady_clock::now();
    running_ = true;
}

inline double
Timer::Stop()
{
#ifndef RELEASE
    if( !running_ )
        throw std::logic_error("Tried to stop a timer before starting it.");
#endif
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
        auto now = steady_clock::now();
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

#endif // ifndef ELEM_CORE_TIMER_IMPL_HPP
